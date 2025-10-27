library(tidyverse)

#=== Create posterior compilation for all species

folder_path <- "03 - Output/Posterior summary/"
file_list <- list.files(path = folder_path, full.names = TRUE)
combined_df <- file_list %>%
  lapply(read.csv) %>%
  bind_rows()
head(combined_df)

parameter_list <- unique(combined_df$parameter)

#=== Remove species

combined_df <- combined_df %>% filter(!species %in% c("MSTAR"))

length(unique(combined_df$species))
#=== Check p-value consistency

(pvalue <- combined_df %>% 
  filter(parameter == "Bayesian p-value") %>% 
  ggplot(aes(species, average))+
  geom_hline(yintercept = 0.35)+
  geom_hline(yintercept = 0.65)+
  geom_point(shape = 21, stroke = 1)+
  labs(y = "Bayesian p-value", x = "Species code")+
  ylim(0,1)+
  coord_flip()+
  theme_classic()+
  theme(axis.text = element_text(colour = "black")))

ggsave(plot = pvalue, filename = "03 - Output/Figures/bayesian_pvalue.tiff", width = 9, height = 10, dpi = 300)

#=== Create significance column

combined_df <- combined_df %>% mutate(significance = case_when(
  average > 0 & lower_cri >0 ~ "yes",
  average < 0 & upper_cri< 0 ~ "yes",
  TRUE ~ "no"
),
sign = case_when(
  average > 0 ~ "positive",
  average < 0 ~ "negative"
))

#=== Migration slope
(migration_slope_plot <- combined_df %>% filter(parameter == parameter_list[[4]]) %>% 
  mutate(migration_type = case_when(
    significance == "no" ~ "No migration",
    average > 0 ~ "Uphill summer migration",
    average < 0 ~ "Uphill winter migration"
  )) %>% group_by(mountain) %>%
  mutate(species = fct_reorder(species, average)) %>% 
  ggplot(aes(species, average, col = migration_type))+
  geom_hline(yintercept = 0)+
  geom_pointrange(aes(ymin = lower_cri, ymax = upper_cri))+
  scale_color_manual(values = c("grey80", "steelblue1", "tan2"))+
  labs(col = "Migration pattern:", y = "Effect size", x = "Species code")+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom",
        panel.grid.major.y = element_line(color = "grey90"),
        axis.text = element_text(colour = "black")))

ggsave(plot = migration_slope_plot, filename = "03 - Output/Figures/migration_slope.tiff", width = 9, height = 10, dpi = 300)



#=== Identify species with at least one significant Elev * Season somewhere
species_with_migration <- combined_df %>%
  filter(parameter == "Elev * Season", significance == "yes") %>%
  pull(species) %>%
  unique()

migration_df <- combined_df %>% 
  filter(species %in% species_with_migration) 

migration_df %>% 
  filter(parameter == "Elev * Season" & significance == "yes") %>% 
  ggplot(aes(species, average))+
  geom_hline(yintercept = 0)+
  geom_pointrange(aes(ymin = lower_cri, ymax = upper_cri))+
  coord_flip()+
  theme_bw()

#=== Elevational effect

(elev_effect <- combined_df %>% filter(parameter == parameter_list[[2]])%>% 
  ggplot(aes(reorder(species, average), average))+
  geom_hline(yintercept = 0)+
  geom_pointrange(aes(ymin = lower_cri, ymax = upper_cri))+
  labs(y = "Elevation effect", x = "Species code")+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90"),
        axis.text = element_text(colour = "black"))+
  annotate(geom = "text", label = "Lowland preference", x = 10, y = -3.1, size = 4)+
  annotate(geom = "text", label = "Upland preference", x = 50, y = 3, size = 4)+
  annotate(geom = "text", label = "Midland preference", x = 28, y = -0.85, size = 4, angle = 70))

ggsave(plot = elev_effect, filename = "03 - Output/Figures/elev_effect.tiff", width = 9, height = 10, dpi = 300)


#=== Season effect

(season_effect <- combined_df %>% filter(parameter == parameter_list[[3]])%>% 
  ggplot(aes(reorder(species, average), average))+
  geom_hline(yintercept = 0)+
  geom_pointrange(aes(ymin = lower_cri, ymax = upper_cri))+
  labs(y = "Season effect", x = "Species code")+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90"),
        axis.text = element_text(colour = "black"))+
  annotate(geom = "text", label = "Abundance increases\nin winter", x = 10, y = -2, size = 4)+
  annotate(geom = "text", label = "Abundance increases\nin summer", x = 40, y = 1.4, size = 4))

ggsave(plot = season_effect, filename = "03 - Output/Figures/season_effect.tiff", width = 9, height = 10, dpi = 300)

#=== Stack plot migration

migration_summary_table <- combined_df %>% 
  filter(parameter == parameter_list[[4]]) %>% 
  mutate(migration_type = case_when(
    significance == "no" ~ "No migration",
    average > 0 ~ "Uphill summer migration",
    average < 0 ~ "Uphill winter migration"
  )) %>% 
  group_by(migration_type) %>%
  summarise(n_species = n(), .groups = "drop") %>% 
  mutate(migration_type = factor(migration_type, levels = c(
    "Uphill summer migration",
    "No migration",
    "Uphill winter migration"
  ))) 



migration_summary_table %>% 
  ggplot(aes(x = NA, y = n_species, fill = migration_type)) +
  geom_col(position = "stack") +
  theme_minimal() +
  scale_fill_manual(
    values = c(
      "Uphill summer migration" = "steelblue1",
      "Uphill winter migration" = "tan2",
      "No migration" = "grey80"
    )
  ) +
  labs(
    x = "",
    y = "Number of species",
    fill = "Migration pattern:"
  )+
  annotate(geom = "text", label = "34", x = 1, y = 20, size = 5)+
  annotate(geom = "text", label = "3", x = 1, y = 1.5, size = 5)+
  annotate(geom = "text", label = "22", x = 1, y = 47, size = 5)+
  theme_classic()+
  theme(legend.position = "bottom",
        axis.text = element_text(colour = "black"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

#=== Save dataset

#write.csv(combined_df, "03 - Output/Posterior summary/combined_posteriors.csv")


#=== Marginal slopes


combined_df %>% 
  filter(species %in% species_with_migration) %>% 
  filter(parameter %in% c(parameter_list[[4]], parameter_list[[5]], parameter_list[[6]])) %>% 
  select(species, parameter, average, significance) %>% 
  ggplot(aes(species, average, fill = parameter))+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("grey", "tan2", "steelblue1"))

