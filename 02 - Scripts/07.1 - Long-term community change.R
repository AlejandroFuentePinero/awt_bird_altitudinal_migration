library(tidyverse)
library(jagsUI)
library(bayestestR)
library(Hmisc)
library(vegan)
library(betapart)


# Compile posterior predictions -------------------------------------------

df <- read.csv("01 - Data/pop/data_migration_expanded.csv")

covs <- read.csv("01 - Data/covs/covs.csv")

bioclim <- read.csv("01 - Data/covs/bioclim_extract.csv")

covs2 <- left_join(covs,bioclim)

covs2 %>% ggplot(aes(elevation, Annual.Mean.Temperature))+
  geom_point()


# Create elevational lookup table
elev_pred <- seq(min(df$elev), max(df$elev), by = 50)

elev_lookup <- data.frame(
  elev_band = 1:length(elev_pred),
  elevation = elev_pred
)


# Select species to model
spp_list <- unique(df$species)
spp_list

# Create empty list to hold results
all_N_list <- list()

# Loop through each species
for (spp in spp_list) {
  message(paste("🔄 Processing", spp))
  
  model_path <- file.path("03 - Output/R_output_pred_pop", paste0(spp, "pred_pop.RData"))
  
  if (!file.exists(model_path)) {
    message(paste("⚠️ Skipping", spp, "- model file not found"))
    next
  }
  
  # Load into isolated environment to avoid variable overwrites
  tmp_env <- new.env()
  load(model_path, envir = tmp_env)
  
  if (!"out3" %in% names(tmp_env)) {
    message(paste("❌ Skipping", spp, "- 'out3' not found in file"))
    next
  }
  
  if (!"N_pred" %in% names(tmp_env$out3$sims.list)) {
    message(paste("❌ Skipping", spp, "- 'N'_pred not found in out3$sims.list"))
    next
  }
  
  N_array <- tmp_env$out3$sims.list$N_pred
  
  if (length(dim(N_array)) != 3) {
    message(paste("❌ Skipping", spp, "- 'N'_pred is not a 3D array (dim:", paste(dim(N_array), collapse = ", "), ")"))
    next
  }
  
  dimnames(N_array) <- list(
    sample = 1:dim(N_array)[1],
    elev_band = 1:dim(N_array)[2],
    time   = 1:dim(N_array)[3]
  )
  
  # Convert to tidy data frame
  N_df <- as.data.frame.table(N_array, responseName = "population") %>%
    mutate(
      sample    = as.integer(as.character(sample)),
      elev_band = as.integer(as.character(elev_band)),
      time      = as.integer(as.character(time)),
      species   = spp
    ) %>%
    left_join(elev_lookup, by = "elev_band")
  
  if (any(is.na(N_df$site_name))) {
    message(paste("⚠️ Warning for", spp, "- unmatched site IDs in site_lookup"))
  }
  
  all_N_list[[spp]] <- N_df
  message(paste("✅ Done:", spp, "with", nrow(N_df), "rows"))
}

# Combine all into one dataframe
all_N <- bind_rows(all_N_list)
head(all_N)
nrow(all_N) == (1500*49*34*23)

# Join with metadata
meta <- df[, c("species", "season_step", "bio_year", "season")]
names(meta)[names(meta) == "season_step"] <- "time"
meta <- distinct(meta)

all_N <- left_join(all_N, meta)
unique(all_N$species)


all_N %>% 
  group_by(elevation, species, season) %>% 
  summarise(mean = mean(population)) %>% 
  ggplot(aes(elevation,mean, col = season))+
  geom_point()+
  geom_line()+
  facet_wrap(~species, scales = "free")

all_N %>% filter(species == "ABT") %>% 
  group_by(elevation, species, season) %>% 
  summarise(mean = mean(population)) %>% 
  ggplot(aes(elevation,mean, col = season))+
  geom_point()+
  geom_line()+
  facet_wrap(~species, scales = "free")

#all_N <- all_N %>% filter(!species %in% c("ABT"))

# Checks

unique(all_N$species)
length(unique(all_N$species))
sum(is.na(all_N)) # should be zero
nrow(all_N) == (1500*49*34*23)


# Community ---------------------------------------------------------------

ci_level <- 0.89

# Normalise within species × sample × bio_year (distribute each species' yearly abundance across elevations)
norm_all_N_year <- all_N %>%
  group_by(sample, species, bio_year) %>%                      # <- was season
  mutate(
    total_pop = sum(population, na.rm = TRUE),
    norm_pop  = if_else(total_pop == 0, 0, population / total_pop)
  ) %>%
  ungroup() %>%
  select(-total_pop)

sum(is.na(norm_all_N_year))

# Sum across species at each elevation × sample × year (community, species-weighted)
community_norm_year <- norm_all_N_year %>%
  group_by(sample, elevation, bio_year) %>%
  summarise(total_norm = sum(norm_pop), .groups = "drop")


# Mean & CI across posterior samples for each elevation × year
summary_norm_year <- community_norm_year %>%
  group_by(elevation, bio_year) %>%
  summarise(
    mean = mean(total_norm, na.rm = TRUE),
    CI   = list(ci(total_norm, ci = 0.89)),
    .groups = "drop"
  ) %>%
  unnest_wider(CI)


# Plot: community normalised abundance by elevation for each year (can be busy; facet/year optional)
(p_norm_total_abund_year <- ggplot(summary_norm_year,
                                  aes(x = elevation, y = mean, group = factor(bio_year), col = factor(bio_year))) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high),
                  fatten = 4, shape = 21, linewidth = 0.4) +
  labs(x = "Elevation (m)",
       y = "Normalised total abundance (species-weighted)",
       col = "Year") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10)))



# Wide by year, then compute delta to next year within each elevation × sample
delta_norm_year <- community_norm_year %>%
  arrange(sample, elevation, bio_year) %>%
  group_by(sample, elevation) %>%
  mutate(delta = total_norm - lag(total_norm)) %>%           # year_t − year_(t−1)
  ungroup() %>%
  filter(!is.na(delta))

# Summarise mean delta (across samples) per elevation (overall) or per elevation×year (optional)
delta_norm_year_summary <- delta_norm_year %>%
  group_by(elevation) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    CI = list(ci(delta, ci = 0.89)),
    .groups = "drop"
  ) %>%
  unnest_wider(CI)


(p_delta_community_norm_year <- ggplot(delta_norm_year_summary, aes(elevation, mean_delta)) +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high),
                  fatten = 4, shape = 21, fill = "white", linewidth = 0.5) +
  labs(x = "Elevation (m)", y = "Year-to-year change in normalised abundance") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1),
        axis.title = element_text(size = 14)))

#========================================================
# 3) BETA DIVERSITY (BASELGA) OVER TIME:
#    CONSECUTIVE YEAR COMPARISONS WITHIN ELEVATION BAND
#========================================================

ci_level <- 0.89
# all_N must have: species, sample, bio_year, elev_band, elevation, population

#========================================================
# 1) COMMUNITY TOTALS BY ELEVATION × YEAR (RAW ABUND)
#========================================================

# Sum raw abundance across species at each elevation × sample × year
community_raw_year <- all_N %>%
  group_by(sample, elevation, bio_year) %>%
  summarise(total_raw = sum(population, na.rm = TRUE), .groups = "drop")

# Summarise across posterior samples to mean & CI per elevation × year
summary_raw_year <- community_raw_year %>%
  group_by(elevation, bio_year) %>%
  summarise(
    mean = mean(total_raw, na.rm = TRUE),
    CI   = list(ci(total_raw, ci = ci_level)),
    .groups = "drop"
  ) %>%
  unnest_wider(CI)

# Plot: raw community abundance by elevation for each year (facet for readability)
(p_raw_total_abund_year_facet <- ggplot(summary_raw_year,
                                       aes(x = bio_year, y = mean, col = elevation, group = elevation)) +
  # geom_pointrange(aes(ymin = CI_low, ymax = CI_high),
  #                 fatten = 3, shape = 21, linewidth = 0.35) +
    geom_line()+
  #facet_wrap(~ bio_year, scales = "free_y") +
  labs(x = "Elevation (m)", y = "Total abundance (raw)") +
  theme_classic() +
  theme(strip.text = element_text(size = 9, face = "bold"),
        axis.text = element_text(colour = "black", size = 9),
        axis.title = element_text(size = 12)))

#========================================================
# 2) YEAR-TO-YEAR Δ IN RAW COMMUNITY ABUNDANCE
#========================================================

delta_raw_year <- community_raw_year %>%
  arrange(sample, elevation, bio_year) %>%
  group_by(sample, elevation) %>%
  mutate(delta = total_raw - lag(total_raw)) %>%  # change to next year
  ungroup() %>%
  filter(!is.na(delta))

# Optional: Δ per transition year (to_year label)
delta_raw_year_by_transition <- delta_raw_year %>%
  mutate(to_year = bio_year) %>%
  group_by(elevation, to_year) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    CI = list(ci(delta, ci = ci_level)),
    .groups = "drop"
  ) %>% unnest_wider(CI)

(p <- delta_raw_year_by_transition %>% ggplot(aes(to_year, mean_delta, group = elevation, fill = elevation))+
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.035)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "black", size = 1)+
  geom_line(aes(, col = elevation), size = 1)+
  scale_fill_viridis_c()+
  scale_colour_viridis_c()+
    labs(x = 'Year to year comparison\n(t vs t-1)',
         y = "Mean abundance difference",
         fill = "Elevation",
         colour = "Elevation")+
  theme_classic() +
    theme(axis.text = element_text(colour = "black", size = 12),
                         axis.title = element_text(size = 14),
                         legend.text = element_text(size =12), 
                         legend.title = element_text(size = 14)))

ggsave(plot = p, filename = "03 - Output/Figures/net_abundance_change_long_term.tiff", width = 9, height = 10, dpi = 300)


#========================================================
# 3) BASELGA β-DIVERSITY YEAR→YEAR (RAW ABUND, BRAY)
#========================================================
# Build species-by-year matrices per elev_band × sample (use raw or normalised abundance)
abund_matrix_year <- all_N %>%
  group_by(sample, elev_band, bio_year, species) %>%
  summarise(abund = sum(population, na.rm = TRUE), .groups = "drop") %>%  # or norm_pop if you prefer
  arrange(sample, elev_band, bio_year) %>%
  tidyr::pivot_wider(names_from = species, values_from = abund, values_fill = 0)

beta_year <- abund_matrix_year %>%
  group_by(elev_band, sample) %>%
  arrange(bio_year, .by_group = TRUE) %>%
  nest() %>%
  mutate(
    beta_parts = map(data, ~{
      df <- .x
      yrs <- df$bio_year
      comm <- df %>% select(-bio_year)
      
      if (nrow(comm) < 2) return(tibble(from_year = NA_integer_, to_year = NA_integer_,
                                        turnover = NA_real_, nestedness = NA_real_, total = NA_real_))
      
      map2_dfr(seq_len(nrow(comm)-1), seq_len(nrow(comm)-1)+1, \(i,j){
        pair <- comm[c(i,j), , drop = FALSE]
        b <- beta.pair.abund(as.matrix(pair), index.family = "bray")
        tibble(
          from_year  = yrs[i],
          to_year    = yrs[j],     # ← use this on the x-axis
          turnover   = as.numeric(b$beta.bray.bal[1]),
          nestedness = as.numeric(b$beta.bray.gra[1]),
          total      = as.numeric(b$beta.bray[1])
        )
      })
    })
  ) %>%
  unnest(beta_parts) %>%
  select(-data) %>%
  filter(!is.na(from_year))

beta_year_by_transition <- beta_year %>%
  group_by(elev_band, to_year) %>%
  summarise(
    mean_turnover   = mean(turnover, na.rm = TRUE),
    ci_turnover     = list(ci(turnover,   ci = ci_level)),
    mean_nestedness = mean(nestedness,   na.rm = TRUE),
    ci_nestedness   = list(ci(nestedness, ci = ci_level)),
    mean_total      = mean(total,        na.rm = TRUE),
    ci_total        = list(ci(total,      ci = ci_level)),
    .groups = "drop"
  ) %>%
  unnest_wider(ci_turnover, names_sep = "_") %>%
  unnest_wider(ci_nestedness, names_sep = "_") %>%
  unnest_wider(ci_total, names_sep = "_") %>%
  rename(
    turn_low  = ci_turnover_CI_low,  turn_high  = ci_turnover_CI_high,
    nest_low  = ci_nestedness_CI_low,nest_high  = ci_nestedness_CI_high,
    total_low = ci_total_CI_low,     total_high = ci_total_CI_high
  )

beta_year_by_transition <- beta_year_by_transition %>%
  left_join(elev_lookup, by = "elev_band")

(p_beta_total_lines <- ggplot(beta_year_by_transition,
                             aes(x = to_year, y = mean_total, group = elev_band)) +
  geom_ribbon(aes(ymin = total_low, ymax = total_high), alpha = 0.15) +
  geom_line() +
  facet_wrap(~ elev_band) +
  labs(x = "Year (t)", y = "Bray β-diversity: total (t vs t-1)") +
  theme_classic())

(p_beta_total_heat <- ggplot(beta_year_by_transition,
                            aes(x = to_year, y = elevation %||% elev_band, fill = mean_total)) +
  geom_tile() +
  scale_fill_viridis_c(name = "β_total") +
  labs(x = "Year (t)", y = "Elevation", title = "Year-to-year β-diversity (t vs t-1)") +
  theme_classic())

