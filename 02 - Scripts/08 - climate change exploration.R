library(tidyverse)

#=== Create posterior compilation for all species

folder_path <- "03 - Output/Posterior climate change/"
file_list <- list.files(path = folder_path, full.names = TRUE)
combined_df <- file_list %>%
  lapply(read.csv) %>%
  bind_rows()
head(combined_df)

#=== Filter species with significant migration
migrating_spp <- read.csv("03 - Output/Posterior summary/combined_posteriors.csv") %>% 
  filter(significance == "yes" & parameter == "Elev * Season")

combined_df$parameter <- factor(combined_df$parameter,
                                levels = c("Early migration (1999-2004)", "Mid migration (2005-2010)", "Late migration (2011-2015)"))

spp_list <- unique(combined_df$species)

combined_df %>% filter(species %in% spp_list[1:10]) %>% 
  ggplot(aes(parameter, average))+
  geom_hline(yintercept = 0)+
  geom_pointrange(aes(ymin = lower_cri, ymax = upper_cri))+
  facet_wrap(~species, scales = "free")


#=== Full posteriors

posterior <- read.csv("03 - Output/Posterior traits data/migration_posterior_cc.csv")

pos_summary <- posterior[,-c(1)] %>%
  filter(species %in% migrating_spp$species) %>% 
  dplyr::select(-species) %>% 
  gather("par", "val", 1:3) %>% 
  group_by(par) %>% 
  summarise(mean = mean(val),
            low = bayestestR::ci(val, ci = 0.89)[[2]],
            hi = bayestestR::ci(val, ci = 0.89)[[3]]) 

pos_summary$par <- factor(pos_summary$par, 
                          levels = c("posterior_early", "posterior_mid", "posterior_late"), 
                          labels = c("Early", "Mid", "Late"))

pos_summary %>% ggplot(aes(par, mean))+
  geom_pointrange(aes(ymin = low, ymax = hi))

