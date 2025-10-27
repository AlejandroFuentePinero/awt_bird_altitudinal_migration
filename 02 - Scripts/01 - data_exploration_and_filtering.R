library(tidyverse)


df <- read.csv("01 - Data/pop/data_expanded.csv")

covs <- read.csv("01 - Data/covs/covs.csv")

elev <- tibble(site_name = covs$site,
               elev = covs$elev_round,
               mountain = substr(covs$site, 1, 2),
               elev_cat = case_when(
                 elev %in% c(100:400) ~ "low",
                 elev %in% c(600:900) ~ "mid",
                 elev %in% c(1000:1300) ~ "high"
               ))

df <- left_join(df, elev)

df <- df[,-1] %>% filter(!is.na(count)) %>% 
  mutate(season = cut(date,
                      breaks = c(0, 90, 182, 273, 366),
                      labels = c("DJF", "MAM", "JJA", "SON"),
                      include.lowest = TRUE))



survey_coverage <- df %>% filter(species == "ASW") %>% 
  group_by(site_name, elev, season, mountain) %>% 
  summarise(n_survey = n(), .groups = "drop")



survey_coverage %>% filter(!is.na(season)) %>% ggplot(aes(season, elev, fill = n_survey))+
  geom_tile()+
  facet_wrap(~mountain)+
  scale_fill_viridis_c()+
  theme_bw()

# NOT ENOUGH SPATIO-TEMPORAL COVERAGE FOR WU, DROP IT FOR ANALYSES

df <- df %>% 
  filter(mountain %in% c("AU", "CU")) %>% 
  filter(!is.na(date)) # also filter out surveys without date



species_summary <- df %>% 
  group_by(species) %>% 
  filter(count > 0) %>% 
  summarise(
    total_detections = n(),
    n_elev_bands = n_distinct(elev),
    n_seasons = n_distinct(season),
    .groups = "drop"
  )


species_keep <- species_summary %>%
  filter(total_detections >= 30, n_elev_bands >= 3, n_seasons >= 3)

species_keep <- unique(species_keep$species)

(species_rmv <- setdiff(unique(df$species), species_keep))



df <- df %>% filter(!species %in% species_rmv)  


# Remove birds for which rainforest is not core habitat 

traits <- read.csv("/Users/alejandrofp/Library/CloudStorage/OneDrive-JamesCookUniversity/Postdoc/Projects/10 - Bird elevational migration/awt_birds_altitudinal_migration/01 - Data/traits/bird_awt_traits.csv")

traits_rainforest <- traits %>% filter(rainforest_specialization..1.occasionally...6.obligate. >=4) %>% select(spp, rainforest_specialization..1.occasionally...6.obligate.)

rainforest_spp <- unique(sort(traits_rainforest$spp))

(not_rainforest <- setdiff(species_keep,sort(rainforest_spp)))
  
df <- df %>% filter(species %in% rainforest_spp)

unique(df$species)
unique(df$mountain)

#write.csv(df,"01 - Data/filtered_df.csv")


# expand data -------------------------------------------------------------

df_prep <- df %>% dplyr::select(1,2,3,5:11)


str(df_prep)

unique(df_prep$mountain)

df_season <- df_prep %>% 
  mutate(
    julian = date,
    bio_year = ifelse(julian < 121, year - 1, year),
    season = case_when(
      julian >= 121 & julian <= 304 ~ "winter",
      TRUE ~ "summer"),
    year_season = paste0(bio_year,"_", season)) %>% 
  group_by(species, site_name, year_season) %>% 
  arrange(date, .by_group = TRUE) %>% 
  mutate(rep = row_number()) %>% 
  ungroup() %>% 
  mutate(season_id = as.integer(factor(year_season, levels = sort(unique(year_season))))
  )


df_fill <- df_season %>% dplyr::select(1,2,10,9,12,13,11,8,5:7,15)



max_rep <- df_fill %>%
  count(species, site_name, bio_year, season) %>%
  summarise(max_rep = max(n)) %>%
  pull(max_rep)

full_grid <- df_fill %>%
  distinct(species, site_name) %>%
  crossing(bio_year = 1999:2015,
           season = c("winter", "summer"),
           rep = 1:max_rep)


# Add back season-level metadata
meta_cols <- elev[,-4]

# Join metadata
full_grid <- left_join(full_grid, meta_cols)

# Join observations
complete_df <- left_join(full_grid, df_fill)


season_id <- complete_df %>%
  mutate(season_order = ifelse(season == "summer", 1, 2)) %>%
  distinct(bio_year, season, season_order) %>%
  arrange(bio_year, season_order) %>%
  mutate(season_step = row_number()) %>% 
  dplyr::select(-season_order)

df_final <- left_join(complete_df, season_id)

nrow(df_final %>% filter(!is.na(count))) == nrow(df) # CHECK POINT - should be TRUE

unique(df_final$species)
unique(df_final$site_name)
unique(df_final$bio_year)
unique(df_final$mountain)
#write.csv(df_final,"01 - Data/pop/data_migration_expanded.csv")




# -------------------------------------------------------------------------


# -------------------------------------------------------------------------


# -------------------------------------------------------------------------



species_sel <- "RBBE"

df_plot <- df_final %>%
  filter(species == species_sel) %>%
  group_by(mountain, season, elev) %>%
  summarise(mean_count = mean(count, na.rm = TRUE),
            n = n()) %>%
  ungroup()

ggplot(df_plot, aes(x = elev, y = mean_count, color = season)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ mountain) +
  labs(x = "Elevation (m)", y = "Mean Count", color = "Season",
       title = "Seasonal changes in abundance across elevation") +
  theme_bw()


library(dplyr)
library(tidyr)
library(broom)
library(purrr)

# STEP 1: Prepare the data
df_model <- df_final %>%
  filter(!is.na(count)) %>%
  mutate(
    season_bin = ifelse(season == "summer", 0.5, -0.5),  # effect-coded for interpretability
    species = as.factor(species),
    mountain = as.factor(mountain)
  )

# STEP 2: Fit linear model for each species
# Model: count ~ season_bin * elev + (1 | mountain) — but using fixed effect mountain for simplicity
model_results <- df_model %>%
  group_by(species) %>%
  group_split() %>%
  map_df(~ {
    data <- .
    if (length(unique(data$season_bin)) > 1 && length(unique(data$elev)) > 1) {
      mod <- lm(count ~ season_bin * elev + mountain, data = data)
      tidy(mod) %>%
        filter(term == "season_bin:elev") %>%
        select(term, estimate, std.error, p.value) %>%
        mutate(species = unique(data$species))
    } else {
      tibble(term = "season_bin:elev", estimate = NA, std.error = NA, p.value = NA, species = unique(data$species))
    }
  })

# STEP 3: Arrange species by interaction strength
migration_summary <- model_results %>%
  arrange(desc(abs(estimate))) %>%
  mutate(
    direction = case_when(
      estimate > 0 ~ "uphill in summer",
      estimate < 0 ~ "downhill in summer",
      TRUE ~ "no signal"
    )
  )

# View top results
migration_summary

