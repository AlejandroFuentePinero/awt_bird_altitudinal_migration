# Packages
library(tidyverse)
library(jagsUI)
library(AHMbook)
library(bayestestR)

# Read data
df <- read.csv("01 - Data/pop/data_migration_expanded.csv")

# Encode season for the model
df <- df %>% mutate(season_num = ifelse(season == "summer", 0.5, -0.5))

# Constants
nsites <- length(unique(df$site_name))
ntime <- length(unique(df$season_step))
nreps <- length(unique(df$rep))
nmountain <- length(unique(df$mountain))

# Select species to model
spp_list <- unique(df$species)

#spp_list <- spp_list[41:49]

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
#(spp_sel <- spp_list[[1]])
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

for(s in 1:length(spp_list)){
  
  (spp_sel <- spp_list[[s]])
  
  
  df_spp <- df %>% filter(species == spp_sel) %>% dplyr::select(-c(1,2))
  
  # Detection covatiates
  
  # Wind
  wind <- df_spp %>% dplyr::select(site_name, season_step, rep, wind)
  wind$clump <- paste0("X", sprintf("%02d", wind$season_step), ".", wind$rep)
  wind <- wind[,-c(2,3)]
  wind_spread <- wind %>% spread(clump, wind)
  wind_spread <- wind_spread[, c(1, order(names(wind_spread)[-1]) + 1)]
  wind_matrix <- as.matrix(wind_spread[, -1])
  wind_matrix <- standardize(wind_matrix)
  wind_array <- array(wind_matrix, dim = c(nsites, nreps, ntime))
  wind_array[is.na(wind_array)] <- 0 # inputation
  
  # Rain
  rain <- df_spp %>% dplyr::select(site_name, season_step, rep, rain)
  rain$clump <- paste0("X", sprintf("%02d", rain$season_step), ".", rain$rep)
  rain <- rain[,-c(2,3)]
  rain_spread <- rain %>% spread(clump, rain)
  rain_spread <- rain_spread[, c(1, order(names(rain_spread)[-1]) + 1)]
  rain_matrix <- as.matrix(rain_spread[, -1])
  rain_matrix <- standardize(rain_matrix)
  rain_array <- array(rain_matrix, dim = c(nsites, nreps, ntime))
  rain_array[is.na(rain_array)] <- 0
  
  # Date
  date <- df_spp %>% dplyr::select(site_name, season_step, rep, julian)
  date$clump <- paste0("X", sprintf("%02d", date$season_step), ".", date$rep)
  date <- date[,-c(2,3)]
  date_spread <- date %>% spread(clump, julian)
  date_spread <- date_spread[, c(1, order(names(date_spread)[-1]) + 1)]
  date_matrix <- as.matrix(date_spread[, -1])
  date_matrix <- standardize(date_matrix)
  date_array <- array(date_matrix, dim = c(nsites, nreps, ntime))
  date_array[is.na(date_array)] <- 0
  
  
  # Count array
  df_count <- df_spp %>% dplyr::select(site_name, season_step, rep, count)
  df_count$clump <- paste0("X", sprintf("%02d", df_count$season_step), ".", df_count$rep)
  df_count <- df_count[,-c(2,3)]
  df_spread <- df_count %>% spread(clump, count)
  df_spread <- df_spread[, c(1, order(names(df_spread)[-1]) + 1)]
  counts_hl <- as.matrix(df_spread[, -1])
  C <- array(counts_hl, dim = c(nsites, nreps, ntime))
  
  # Count check
  C[6,,12]
  C[7,,20]
  C[1,,3]
  
  
  # Other variables
  mountain <- tibble(site_name = df_spread$site_name,
                     mountain = substr(df_spread$site_name, 1, 2),
                     mountain_id = as.integer(as.factor(mountain)))
  
  elev <- df %>% 
    dplyr::select(site_name, elev) %>% 
    distinct() %>% 
    mutate(elev_std = standardize(elev))
  
  
  season <- df %>%
    dplyr::select(bio_year, season, season_step, season_num) %>% 
    distinct()
  
  elev_mean <- mean(elev$elev)
  elev_sd <- sd(elev$elev)
  
  elev_pred <- seq(min(df$elev), max(df$elev), by = 50)
  nelev_pred <- length(elev_pred)
  
  elev_pred_scaled <- (elev_pred - elev_mean) / elev_sd
  
  
  range(df$elev)
  range(elev_pred)
  # Compile data
  str(bdata <- list(C = C,
                    wind = wind_array,
                    rain = rain_array,
                    date = date_array,
                    mountain = mountain$mountain_id,
                    season = season$season_num,
                    bio_year = standardize(season$bio_year),
                    elev = elev$elev_std,
                    elev_pred = elev_pred_scaled,
                    nelev_pred = nelev_pred,
                    nsites = dim(C)[1],
                    nreps = dim(C)[2],
                    ntimes = dim(C)[3],
                    nmountain = nmountain
  ))
  
  # Model
  cat(file = "altitudinal_migration_awt.txt","
  
  model {

  #############
  # 1. Priors #
  #############
  
  ########################
  #### Random effects ####
  ########################
  
  
  #=== Detectability
  
    #=== Random season intercept

  for(t in 1:ntimes){
  
    p0[t] ~ dunif(0,1)
    
    alpha0[t] <- log(p0[t]/(1-p0[t])) # intercept detection
  
  } 
  
    tau.lp <- 1 / (sigma.p * sigma.p)
    sigma.p ~ dunif(0, 100)
  
  
  #=== Abundance
  
    #=== Site random effect

    for(i in 1:nsites){
  
      eps[i] ~ dnorm(0, tau.site) # overdispersion
      
    }
    
      tau.site <- 1 / (sigma.site * sigma.site)
      sigma.site ~ dunif(0, 100)
      
      
    #=== Mountain-segregated intercept
    
    for(m in 1:nmountain){
    
      beta0_mountain[m] ~ dnorm(0,0.01)
      
    }


  ######################
  #### Slopes priors ###
  ######################
  
  
  #=== Abundance covs
  
  beta.elev ~ dnorm(0,0.1) 
  beta.elev2 ~ dnorm(0,0.1) 
  beta.season ~ dnorm(0,0.1) 
  beta.trend ~ dnorm(0,0.1) 
  beta.trend2 ~ dnorm(0,0.1) 
  beta.elev_season ~ dnorm(0, 0.1)
  beta.elev_season2 ~ dnorm(0, 0.1)
  
  #=== Detectability covs

  alpha.wind ~ dnorm(0,0.1) 
  alpha.rain ~ dnorm(0,0.1)                                                    
  

  ##############
  # Likelihood #
  ##############
  
  #########################################
  ### Ecological model for true abundance #
  #########################################
  
  for (i in 1:nsites){
  
      for(t in 1:ntimes){
    
    
      N[i,t] ~ dpois(lambda[i,t])
      
      log(lambda[i,t]) <- beta0_mountain[mountain[i]] +                         # mountain-segregated interecept
                          beta.elev * elev[i] +                                 # Spatial elevation effect
                          beta.elev2 * pow(elev[i],2) +                         # Spatial elevation quadratic effect
                          beta.season * season[t] +                             # Temporal season effect
                          beta.elev_season * elev[i] * season[t] +              # Interaction effect (altitudinal migration spotter)
                          beta.elev_season2 * pow(elev[i],2) * season[t] +      # Interaction quadratic effect (altitudinal migration spotter)
                          beta.trend * bio_year[t] +                            # Abundance trend effect
                          beta.trend2 * pow(bio_year[t],2) +                    # Abundance trend effect (quadratic)
                          eps[i]                                                # Site random effect

      
  #############################################
  ### Observation model for replicated counts #
  #############################################
      
       for (j in 1:nreps){
       
        C[i,j,t] ~ dbin(p[i,j,t], N[i,t])
        p[i,j,t] <- exp(lp[i,j,t])/(1+exp(lp[i,j,t]))
        lp[i,j,t] ~ dnorm(mu.lp[i,j,t], tau.lp)
        mu.lp[i,j,t] <- alpha0[t] + alpha.wind * wind[i,j,t] +                  # Wind effect on detection
                                    alpha.rain * rain[i,j,t]                    # Rain effect on detection

      }
    }
  }
  
  ###################
  ### Predictions ###
  ###################
  
  beta_0_avg <- mean(beta0_mountain[])
  
  for (e in 1:nelev_pred) {
  for (t in 1:ntimes) {
    
    log(lambda_pred[e,t]) <- beta_0_avg +        
                             beta.elev * elev_pred[e] +
                             beta.elev2 * pow(elev_pred[e],2) +
                             beta.season * season[t] +
                             beta.elev_season * elev_pred[e] * season[t] +
                             beta.elev_season2 * pow(elev_pred[e],2) * season[t] +
                             beta.trend * bio_year[t] +
                             beta.trend2 * pow(bio_year[t], 2)

    N_pred[e,t] ~ dpois(lambda_pred[e,t])
  }
}

  ########################
  ### Bayesian p-value ###
  ########################


   for(t in 1:ntimes){
    
     for(i in 1:nsites){
     
      for(j in 1:nreps){
      
        
      
     #=== Actual data
     
     eval[i,j,t] <-N[i,t]*p[i,j,t] # Expected value
     sd.resi[i,j,t]<-sqrt(eval[i,j,t]*(1-p[i,j,t])) +0.5
     E[i,j,t]<-(C[i,j,t]-eval[i,j,t])/ sd.resi[i,j,t]
     E2[i,j,t] <- pow(E[i,j,t],2) 
     
     
     
     #=== Replicate data sets
     
     C.new[i,j,t]~dbin(p[i,j,t],N[i,t])
     E.new[i,j,t]<-(C.new[i,j,t]-eval[i,j,t])/sd.resi[i,j,t]
     E2.new[i,j,t] <- pow(E.new[i,j,t], 2)
     
  
      }
    }
  }    
    fit <- sum(E2[,,])                                                          # Sum up squared resi for actual data set
    fit.new <- sum(E2.new[,,])                                                  # Sum up for replicate data sets 
    
    
  ##########################
  ### Derived parameters ###
  ##########################
  
  #=== Backtransformed intercept
  
  for(m in 1:nmountain){
  
    beta0_exp[m] <- exp(beta0_mountain[m])
    
  }
  
  # Average deteection probability
  
  avg_p <- mean(p0[])
  
  # Decriphing direction
  
  slope_winter <- beta.elev - 0.5 * beta.elev_season
  slope_summer <- beta.elev + 0.5 * beta.elev_season

  
}
")
  
  
  
  # Initial values
  Nst <- apply(C, c(1,3), max, na.rm = TRUE)+1 # Inits for latent N
  Nst[Nst == '-Inf'] <- 1
  inits <- function() list(N = Nst,
                           ###
                           beta0_mountain = runif(bdata$nmountain),
                           ###
                           beta.elev = runif(1),
                           beta.elev2 = runif(1),
                           beta.season = runif(1),
                           beta.elev_season = runif(1),
                           beta.elev_season2 = runif(1),
                           beta.trend = runif(1),
                           beta.trend2 = runif(1),
                           ###
                           alpha.wind = runif(1),
                           alpha.rain = runif(1)
  )
  
  # Parameters monitored
  params <- c("N_pred","fit", "fit.new")
  
  #na <- 1000 ; ni <- 200000 ; nt <- 100 ; nb <- 100000 ; nc <- 3
  
  na <- 1000 ; ni <- 100000 ; nt <- 100 ; nb <- 50000 ; nc <- 3
  
  # Run model
  
  set.seed(123 + s)
  out3 <- jags(bdata, 
               inits, 
               params, 
               "altitudinal_migration_awt.txt", 
               n.adapt = na,
               n.chains = nc, 
               n.thin = nt, 
               n.iter = ni, 
               n.burnin = nb, 
               parallel = TRUE)
  
  pp.check(out3, observed = "fit", simulated = "fit.new")
  
  options(max.print = 2000)
  
  
  print(out3, 3)
  
  save.image(paste0("03 - Output/R_output_pred_pop/",spp_sel,"pred_pop.RData"))
  
}

