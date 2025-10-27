library(tidyverse)
library(jagsUI)


# Read data
df <- read.csv("01 - Data/pop/data_migration_expanded.csv")

# Select species to model
spp_list <- unique(df$species)



# Create an empty list to store the dataframes
posterior_list <- list()

# Loop over the species list
for (spp in spp_list) {
  # Construct the file path
  file_path <- file.path("03 - Output/R_output", paste0(spp, ".RData"))
  
  # Load the RData file (this will load 'out3' into the environment)
  load(file_path)
  
  # Extract the posterior samples for beta.elev_season
  posterior <- out3$sims.list$beta.elev_season
  
  # Create a dataframe with the samples and species code
  df <- data.frame(species = spp, posterior = posterior)
  
  # Append to the list
  posterior_list[[spp]] <- df
}

# Combine all into one dataframe
posterior_all <- do.call(rbind, posterior_list)

posterior_all %>% group_by(species) %>% summarise(mean = mean(posterior))

#write.csv(posterior_all, "03 - Output/Posterior traits data/migration_posterior.csv")
