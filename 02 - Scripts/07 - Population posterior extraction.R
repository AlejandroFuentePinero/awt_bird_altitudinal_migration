library(tidyverse)
library(jagsUI)
library(bayestestR)
library(Hmisc)
library(vegan)
library(betapart)
library(patchwork)


# Compile posterior predictions -------------------------------------------

df <- read.csv("01 - Data/pop/data_migration_expanded.csv")

df %>% ggplot(aes(elev, reorder(site_name, elev), col = mountain))+
  geom_point(shape = 21, stroke = 1)+
  theme_classic()+
  labs(x = "Elevation (m)", y = "Site name", color = "Mountain range")+
  scale_color_manual(values = c("black", "orange"))+
  theme()+
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size =12), 
        legend.title = element_text(size = 14),
        legend.position = c(0.2,0.8))


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

all_N <- all_N %>% filter(!species %in% c("ABT", "FTCUC", "SBCUC"))

# Checks

unique(all_N$species)
length(unique(all_N$species))
sum(is.na(all_N)) # should be zero
nrow(all_N) == (1500*46*34*23)

# Species level analyses --------------------------------------------------

#=== Compute centroid and range
season_stats <- all_N %>%
  group_by(species, sample, season) %>%
  filter(sum(population, na.rm = TRUE) > 0) %>%
  summarise(
    centroid = weighted.mean(elevation, population),
    q10 = Hmisc::wtd.quantile(elevation, weights = population, probs = 0.10, na.rm = TRUE),
    q90 = Hmisc::wtd.quantile(elevation, weights = population, probs = 0.9, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(range_width = q90 - q10)

head(season_stats)

ci_level <- 0.89

#=== Compute overall elevation range
overall_stats <- season_stats %>%
  group_by(species, sample) %>%
  summarise(
    centroid = mean(centroid),
    q10 = mean(q10),
    q90 = mean(q90),
    .groups = "drop"
  ) %>%
  group_by(species) %>%
  summarise(
    centroid_mean = mean(centroid),
    centroid_ci = list(ci(centroid, ci = ci_level)),
    q10_mean = mean(q10),
    q10_ci = list(ci(q10, ci = ci_level)),
    q90_mean = mean(q90),
    q90_ci = list(ci(q90, ci = ci_level)),
    .groups = "drop"
  )

(centroid_range <- overall_stats %>% ggplot(aes(reorder(species, centroid_mean), centroid_mean))+
  geom_pointrange(aes(ymin = q10_mean, ymax = q90_mean), fatten = 4, shape = 21, fill = "white", linewidth = 0.5)+
  labs(x = "Species code", y = "Elevation (m)")+
  coord_flip()+
  theme_classic()+
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1),
        axis.title = element_text(size = 14)))

#ggsave(plot = centroid_range, filename = "03 - Output/Figures/mean_centroid_range.tiff", width = 9, height = 10, dpi = 300)


#=== Seasonal range

seasonal_stats <- season_stats %>%
  group_by(species, season) %>%
  summarise(
    centroid_mean = mean(centroid),
    centroid_ci = list(ci(centroid, ci = ci_level)),
    q10_mean = mean(q10),
    q10_ci = list(ci(q10, ci = ci_level)),
    q90_mean = mean(q90),
    q90_ci = list(ci(q90, ci = ci_level)),
    .groups = "drop") 

season_stats %>% filter(species == "BCS") %>% summarise(mean = mean(range_width))

# range overlap

# Wide format: one row per sample per species
range_overlap <- season_stats %>%
  group_by(species, sample) %>%
  filter(n_distinct(season) == 2) %>%
  ungroup() %>%
  select(species, sample, season, q10, q90) %>%
  pivot_wider(names_from = season, values_from = c(q10, q90)) %>%
  mutate(
    lower_overlap = pmax(q10_summer, q10_winter, na.rm = TRUE),
    upper_overlap = pmin(q90_summer, q90_winter, na.rm = TRUE),
    overlap_width = upper_overlap - lower_overlap,
    union_width = pmax(q90_summer, q90_winter, na.rm = TRUE) - pmin(q10_summer, q10_winter, na.rm = TRUE),
    prop_overlap = ifelse(overlap_width > 0, overlap_width / union_width, 0)
  )

overlap_summary <- range_overlap %>%
  group_by(species) %>%
  summarise(
    mean_overlap = mean(prop_overlap),
    CI = list(ci(prop_overlap, ci = 0.89)),  # your existing style
    .groups = "drop"
  ) %>%
  unnest_wider(CI) 

seasonal_stats_labeled <- seasonal_stats %>%
  left_join(overlap_summary, by = "species")

label_data <- seasonal_stats_labeled %>%
  filter(season == "summer") %>%
  mutate(overlap_label = paste0(round(mean_overlap * 100), "%"))


# 1) Rectangles for the 10–90% range (as bars)
ranges <- seasonal_stats_labeled %>%
  mutate(y_mid = (q10_mean + q90_mean)/2,   # bar centre
         y_h   = (q90_mean - q10_mean))     # bar height

# 2) One row per species with winter & summer centroids for the arrow
arrow_df <- seasonal_stats_labeled %>%
  select(species, season, centroid_mean) %>%
  pivot_wider(names_from = season, values_from = centroid_mean) %>%
  # keep only species with both seasons present
  filter(!is.na(winter), !is.na(summer))

(seasonal_centroid_range <-
    ggplot(ranges, aes(x = reorder(species, centroid_mean))) +
    # translucent bars for ranges
    geom_tile(aes(y = y_mid, height = y_h, width = 0.42, fill = season, col = season),
              alpha = 0.25) +
    # bigger centroid points
    geom_point(aes(y = centroid_mean, colour = season), size = 1.5, shape = 21, stroke = 1.4) +
    # arrow from winter -> summer centroid with smaller head + neutral colour
    # geom_segment(data = arrow_df,
    #              aes(x = reorder(species, winter),
    #                  xend = reorder(species, winter),
    #                  y = winter, yend = summer),
    #              arrow = arrow(length = unit(1, "mm"), type = "closed"), # smaller head
    #              linewidth = 0.5, colour = "black") +  # neutral colour
    coord_flip(clip = "off") +
    labs(x = "Species code", y = "Elevation (m)",
         colour = "Season", fill = "Season") +
    scale_colour_manual(values = c(winter = "steelblue", summer = "darkorange3")) +
    scale_fill_manual(values   = c(winter = "steelblue", summer = "darkorange3")) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1),
          axis.title  = element_text(size = 14),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.position = c(0.2, 0.9),
          legend.background = element_blank()))


#ggsave(plot = seasonal_centroid_range, filename = "03 - Output/Figures/seasonal_centroid_range.tiff", width = 9, height = 10, dpi = 300)


#=== Changes in centroid and range

delta_stats <- season_stats %>%
  select(species, sample, season, centroid, range_width) %>%
  pivot_wider(
    names_from = season,
    values_from = c(centroid, range_width),
    names_sep = "."
  ) %>%
  filter(!is.na(centroid.summer) & !is.na(centroid.winter)) %>%
  mutate(
    delta_centroid = centroid.summer - centroid.winter,
    delta_range = range_width.summer - range_width.winter
  ) %>%
  group_by(species) %>%
  summarise(
    delta_centroid_mean = mean(delta_centroid),
    delta_centroid_ci = list(ci(delta_centroid, ci = ci_level)),
    delta_range_mean = mean(delta_range),
    delta_range_ci = list(ci(delta_range, ci = ci_level)),
    .groups = "drop"
  )

# global average
season_stats %>%
  select(species, sample, season, centroid, range_width) %>%
  pivot_wider(
    names_from = season,
    values_from = c(centroid, range_width),
    names_sep = "."
  ) %>%
  filter(!is.na(centroid.summer) & !is.na(centroid.winter)) %>%
  mutate(
    delta_centroid = centroid.summer - centroid.winter,
    delta_range = range_width.summer - range_width.winter
  ) %>% summarise(mean_centroid = mean(delta_centroid),
                  centroid_low = ci(delta_centroid, ci = 0.89)[[2]],
                  centroid_high = ci(delta_centroid, ci = 0.89)[[3]],
                  mean_range = mean(delta_range),
                  range_low = ci(delta_range, ci = 0.89)[[2]],
                  range_high = ci(delta_range, ci = 0.89)[[3]])


delta_plot_data <- delta_stats %>%
  unnest_wider(delta_centroid_ci, names_sep = "_") %>%
  rename(
    delta_centroid_lower = delta_centroid_ci_CI_low,
    delta_centroid_upper = delta_centroid_ci_CI_high
  ) %>%
  unnest_wider(delta_range_ci, names_sep = "_") %>%
  rename(
    delta_range_lower = delta_range_ci_CI_low,
    delta_range_upper = delta_range_ci_CI_high
  )

delta_plot_data <- delta_plot_data %>%
  mutate(centroid_cat = case_when(
    delta_centroid_mean > 50   ~ "bright_blue",
    delta_centroid_mean > 25   ~ "faint_blue",
    delta_centroid_mean >= -25 ~ "grey",
    delta_centroid_mean >= -50 ~ "faint_red",
    TRUE                       ~ "bright_red"
  ))

delta_plot_data <- left_join(delta_plot_data, seasonal_stats_labeled)

# Plot
(delta_centroid <- delta_plot_data %>%
  ggplot(aes(reorder(species, centroid_mean),
             delta_centroid_mean,
             color = centroid_cat)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey60", linewidth = 0.4) +
  geom_pointrange(aes(ymin = delta_centroid_lower,
                      ymax = delta_centroid_upper),
                  fatten = 4, shape = 21, fill = "white", linewidth = 0.8, stroke = 1.5) +
  scale_color_manual(values = c(
    "bright_blue" = "darkorange3",  # bright blue
    "faint_blue"  = "darkorange3",  # faint blue
    "grey"        = "grey80",   # neutral grey
    "faint_red"   = "steelblue",  # faint red
    "bright_red"  = "steelblue"   # bright red
  )) +
  labs(x = "Species code",
       y = "Centroid shift (m)\nSummer - Winter") +
    scale_y_continuous(breaks = c(-400,-200,0,200,400),
                       limits = c(-450,450))+
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank()))

#ggsave(plot = delta_centroid, filename = "03 - Output/Figures/delta_centroid.tiff", width = 9, height = 10, dpi = 300)


delta_plot_data <- delta_plot_data %>%
  mutate(range_cat = case_when(
    delta_range_mean > 100   ~ "bright_blue",
    delta_range_mean > 50   ~ "faint_blue",
    delta_range_mean >= -50 ~ "grey",
    delta_range_mean >= -100 ~ "faint_red",
    TRUE                    ~ "bright_red"
  ))

# Plot
(delta_range <- delta_plot_data %>%
  ggplot(aes(reorder(species, centroid_mean),
             delta_range_mean,
             color = range_cat)) +
    geom_hline(yintercept = 0, linetype = "dashed", col = "grey60", linewidth = 0.4) +
    geom_pointrange(aes(ymin = delta_range_lower,
                      ymax = delta_range_upper),
                  fatten = 4, shape = 21, fill = "white", linewidth = 0.8, stroke = 1.5) +
  scale_color_manual(values = c(
    "bright_blue" = "darkorange3",  # bright blue
    "faint_blue"  = "darkorange3",  # faint blue
    "grey"        = "grey80",   # neutral grey
    "faint_red"   = "steelblue",  # faint red
    "bright_red"  = "steelblue"   # bright red
  )) +
  labs(x = "Species code",
       y = "Range width shift (m)\nSummer - Winter") +
    scale_y_continuous(breaks = c(-400,-200,0,200,400),
                       limits = c(-450,450))+
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank()))

#ggsave(plot = delta_range, filename = "03 - Output/Figures/delta_range.tiff", width = 9, height = 10, dpi = 300)




library(cowplot)

(wrapped_plot <- plot_grid(seasonal_centroid_range, delta_centroid, delta_range, 
                           labels = c("A)", "B)", "C)"),
                           label_size = 14,
                           hjust = -0.6,
                           rel_widths = c(1.3, 1, 1),
          ncol = 3))

ggsave(plot = wrapped_plot, filename = "03 - Output/Figures/species_level_result.tiff", width = 12, height = 10, dpi = 300)


# Traits ------------------------------------------------------------------

traits <- read.csv("01 - Data/traits/traits_curated.csv")
traits_raw <- read.csv("01 - Data/traits/bird_awt_traits.csv")

colnames(traits)[colnames(traits) == "spp"] <- "species"
colnames(traits_raw)[colnames(traits_raw) == "spp"] <- "species"

#=== Prepare and join posterior summary to traits
df <- delta_stats
spp_list <- unique(df$species)

traits <- traits %>% filter(species %in% spp_list)
traits_raw <- traits_raw %>% filter(species %in% spp_list)

df <- left_join(df, traits)
df <- left_join(df, traits_raw)


(centroid_diet_plot <- df %>% ggplot(aes(x = diet, y = delta_centroid_mean))+
  geom_point(shape = 21, stroke = 2, size = 3)+
  labs(x = "Diet", y = "Δ centroid\n(Summer - Winter ; m)")+
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 14)))


(centroid_mass_plot <- df %>% ggplot(aes(x = log(mass), y = delta_centroid_mean))+
    geom_point(shape = 21, stroke = 2, size = 3)+
    labs(x = "Mass (log scale)", y = "Δ centroid\n(Summer - Winter ; m)")+
    theme_classic() +
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(size = 14)))

(centroid_foraging_plot <- df %>% ggplot(aes(x = forage, y = delta_centroid_mean))+
    geom_point(shape = 21, stroke = 2, size = 3)+
    labs(x = "Foraging behaviour", y = "Δ centroid\n(Summer - Winter ; m)")+
    theme_classic() +
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(size = 14)))

(centroid_habitat_plot <- df %>% ggplot(aes(x = habitat_special, y = delta_centroid_mean))+
    geom_point(shape = 21, stroke = 2, size = 3)+
    labs(x = "Habitat specialisation", y = "Δ centroid\n(Summer - Winter ; m)")+
    theme_classic() +
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(size = 14)))


(centroid_plot <- plot_grid(centroid_diet_plot, centroid_mass_plot, centroid_foraging_plot, centroid_habitat_plot,
                           labels = c("A)", "B)", "C)", "D)"),
                           label_size = 14,
                           hjust = -0.6,
                           ncol = 2))


summary(lm(delta_range_mean ~ diet, df))

ggsave(plot = centroid_plot, filename = "03 - Output/Figures/centroid_traits.tiff", width = 12, height = 10, dpi = 300)


(range_diet_plot <- df %>% ggplot(aes(x = diet, y = delta_range_mean))+
    geom_point(shape = 21, stroke = 2, size = 3)+
    labs(x = "Diet", y = "Δ range\n(Summer - Winter ; m)")+
    theme_classic() +
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(size = 14)))


(range_mass_plot <- df %>% ggplot(aes(x = log(mass), y = delta_range_mean))+
    geom_point(shape = 21, stroke = 2, size = 3)+
    labs(x = "Mass (log scale)", y = "Δ range\n(Summer - Winter ; m)")+
    theme_classic() +
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(size = 14)))

(range_foraging_plot <- df %>% ggplot(aes(x = forage, y = delta_range_mean))+
    geom_point(shape = 21, stroke = 2, size = 3)+
    labs(x = "Foraging behaviour", y = "Δ range\n(Summer - Winter ; m)")+
    theme_classic() +
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(size = 14)))

(range_habitat_plot <- df %>% ggplot(aes(x = habitat_special, y = delta_range_mean))+
    geom_point(shape = 21, stroke = 2, size = 3)+
    labs(x = "Habitat specialisation", y = "Δ range\n(Summer - Winter ; m)")+
    theme_classic() +
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(size = 14)))


(range_plot <- plot_grid(range_diet_plot, range_mass_plot, range_foraging_plot, range_habitat_plot,
                            labels = c("A)", "B)", "C)", "D)"),
                            label_size = 14,
                            hjust = -0.6,
                            ncol = 2))

ggsave(plot = range_plot, filename = "03 - Output/Figures/range_traits.tiff", width = 12, height = 10, dpi = 300)



(trend_range <- df %>% ggplot(aes(x = delta_range_mean, y = pop_trend))+
    geom_point(shape = 21, stroke = 2, size = 3)+
    labs(x = "Δ range", y = "Long-term population trend")+
    #geom_smooth(method = "lm", col = "red")+
    theme_classic() +
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(size = 14)))

summary(lm(pop_trend ~ delta_range_mean, data = df))

(trend_centroid <- df %>% ggplot(aes(x = delta_centroid_mean, y = pop_trend))+
    geom_point(shape = 21, stroke = 2, size = 3)+
    labs(x = "Δ centroid", y = "Long-term population trend")+
    #geom_smooth(method = "lm", col = "red")+
    theme_classic() +
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(size = 14)))

(trend_plot <- plot_grid(trend_range, trend_centroid, 
                         labels = c("A)", "B)",
                         label_size = 14,
                         hjust = -0.6,
                         ncol = 2)))

ggsave(plot = trend_plot, filename = "03 - Output/Figures/migration_trend.tiff", width = 14, height = 10, dpi = 300)


# Community level analyses ------------------------------------------------

#=== Total abundance
# head(all_N)
# 
# abund_summary <- all_N %>%
#   group_by(sample, elevation, season) %>%
#   summarise(total_abund = sum(population), .groups = "drop")
# 
# summary_stats <- abund_summary %>%
#   group_by(elevation, season) %>%
#   summarise(
#     mean = mean(total_abund),
#     CI = list(ci(total_abund, ci = 0.89)),
#     .groups = "drop"
#   ) %>%
#   unnest_wider(CI)
# 
# summary_wide <- summary_stats %>%
#   pivot_wider(names_from = season,
#               values_from = c(mean, CI_low, CI_high),
#               names_sep = "_")
# 
# calc_overlap_pct <- function(low1, high1, low2, high2) {
#   overlap = max(0, min(high1, high2) - max(low1, low2))
#   union_range = max(high1, high2) - min(low1, low2)
#   if (union_range == 0) return(0)
#   return(round(100 * overlap / union_range, 0))
# }
# 
# 
# summary_wide <- summary_wide %>%
#   mutate(overlap_pct = mapply(
#     calc_overlap_pct,
#     CI_low_summer, CI_high_summer,
#     CI_low_winter, CI_high_winter
#   ))
# 
# delta_abund <- abund_summary %>%
#   pivot_wider(names_from = season, values_from = total_abund) %>%
#   mutate(delta = summer - winter)
# 
# delta_summary <- delta_abund %>%
#   group_by(elevation) %>%
#   summarise(
#     mean_delta = mean(delta, na.rm = TRUE),
#     CI = list(ci(delta, ci = 0.89)),  # Using bayestatr CI as requested
#     .groups = "drop"
#   ) %>%
#   unnest_wider(CI)
# 
# (total_abundance_elev <- ggplot(summary_stats, aes(x = elevation, y = mean, col = season, fill = season)) +
#   geom_pointrange(aes(ymin = CI_low, ymax = CI_high),
#                   fatten = 4, shape = 21, linewidth = 0.5) +
#   # geom_text(data = summary_wide,
#   #           aes(x = elevation,
#   #               y = pmax(mean_summer, mean_winter) + 1000,
#   #               label = paste0(overlap_pct, "%")),
#   #           inherit.aes = FALSE, size = 3) +
#   scale_colour_manual(values = c("winter" = "steelblue1", "summer" = "tan2")) +
#   scale_fill_manual(values = c("winter" = "steelblue1", "summer" = "tan2")) +
#   labs(x = "Elevation (m)", y = "Total bird abundance",
#        col = "Season", fill = "Season") +
#   theme_classic() +
#   theme(axis.text = element_text(colour = "black")))
# 
# ggsave(plot = total_abundance_elev, filename = "03 - Output/Figures/total_abundance_elev.tiff", width = 9, height = 10, dpi = 300)
# 
# #=== Delta abundance
# (delta_community_abundance <- delta_summary %>% ggplot(aes(elevation, mean_delta))+
#   geom_hline(yintercept = 0, colour = "red", linetype = "dashed")+
#   geom_pointrange(aes(ymin = CI_low, ymax = CI_high),
#                   fatten = 4, shape = 21, fill = "white", linewidth = 0.5) +
#   labs(x = "Elevation (m)", y = "Seasonal change in abundance\n(Summer-Winter)") +
#   theme_classic() +
#   theme(axis.text = element_text(colour = "black")))
# 
# ggsave(plot = delta_community_abundance, filename = "03 - Output/Figures/delta_community_abundance.tiff", width = 9, height = 10, dpi = 300)


#=== Normalised abundance

norm_all_N <- all_N %>%
  group_by(sample, species, season) %>%
  mutate(
    total_pop = sum(population, na.rm = TRUE),
    norm_pop = if_else(total_pop == 0, 0, population / total_pop)
  ) %>%
  ungroup() %>%
  select(-total_pop)

sum(is.na(norm_all_N))

# Sum normalized population across species for each elevation × sample × season
community_norm <- norm_all_N %>%
  group_by(sample, elevation, season) %>%
  summarise(total_norm = sum(norm_pop), .groups = "drop")

# Then calculate mean and credible intervals
summary_norm <- community_norm %>%
  group_by(elevation, season) %>%
  summarise(
    mean = mean(total_norm, na.rm = T),
    CI = list(ci(total_norm, ci = 0.89)),
    .groups = "drop"
  ) %>% unnest_wider(CI)

delta_norm_abund <- community_norm %>%
  pivot_wider(names_from = season, values_from = total_norm) %>%
  mutate(delta = summer - winter)

delta_norm_summary <- delta_norm_abund %>%
  group_by(elevation) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    CI = list(ci(delta, ci = 0.89)),
    .groups = "drop"
  ) %>%
  unnest_wider(CI)

# compute tile centre/height and a sensible tile width from your elevation bands
norm_tiles <- summary_norm %>%
  mutate(y_mid = (CI_low + CI_high)/2,
         y_h   = (CI_high - CI_low)) %>%
  arrange(elevation)

band_width <- 0.4 * median(diff(sort(unique(norm_tiles$elevation))))  # ~80% of band gap

(normalised_total_abund <-
  ggplot(norm_tiles, aes(x = elevation, colour = season)) +
  # translucent tiles for the CI (overlap shows clearly)
  geom_tile(aes(y = y_mid, height = y_h, width = band_width, fill = season),
            alpha = 0.25, colour = NA) +
  # mean point on top (outlined so it pops on both tiles)
  geom_point(aes(y = mean, colour = season), size = 2.5, shape = 21, stroke = 1.6) +
  geom_vline(xintercept = 550, linetype = "dashed", col = "grey60", linewidth = 0.4) +
  scale_colour_manual(values = c(winter = "steelblue", summer = "darkorange")) +
  scale_fill_manual(values   = c(winter = "steelblue", summer = "darkorange")) +
  labs(x = "Elevation (m)", y = "Total abundance\n(Normalised)",
       colour = "Season", fill = "Season") +
  theme_classic() +
  theme(axis.text  = element_text(colour = "black", size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 14),
        legend.text  = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.position = c(0.80, 0.15),
        legend.background = element_blank(),
        legend.direction  = "horizontal",
        legend.box        = "horizontal"))


# Crossover in panel A
# per-sample deltas already in delta_norm_abund (one row per sample×elev)
find_zero <- function(df) {
  d <- df[order(df$elevation), c("elevation","delta")]
  idx <- which(diff(sign(d$delta)) != 0)
  if (!length(idx)) return(NA_real_)
  i <- idx[1]
  x1<-d$elevation[i]; x2<-d$elevation[i+1]
  y1<-d$delta[i];     y2<-d$delta[i+1]
  x1 - y1*(x2-x1)/(y2-y1)  # linear interpolation
}
xovers <- delta_norm_abund %>% group_by(sample) %>% summarise(xover = find_zero(cur_data()))
c(mean = mean(xovers$xover, na.rm=TRUE),
  ci = bayestestR::ci(xovers$xover, ci=.89))



#ggsave(plot = normalised_total_abund, filename = "03 - Output/Figures/norm_total_abundance_elev.tiff", width = 9, height = 10, dpi = 300)


(delta_community_norm_plot <- ggplot(delta_norm_summary, aes(elevation, mean_delta)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey60", linewidth = 0.4) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high),
                  fatten = 4, shape = 21, fill = "white", linewidth = 0.8, stroke = 1.5) +
  labs(x = "Elevation (m)", y = "Abundance change\n(Summer - Winter)") +
  theme_classic()  +
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 14),
          legend.text = element_text(size =12), 
          legend.title = element_text(size = 14)))


delta_norm_summary2 <- delta_norm_summary %>%
  mutate(dir = ifelse(mean_delta >= 0, "pos", "neg"))

(delta_community_norm_plot <-
  ggplot(delta_norm_summary2, aes(elevation, mean_delta)) +
  # reference line
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +
  # # CI + mean
  # geom_linerange(aes(ymin = CI_low, ymax = CI_high, colour = dir),
  #                linewidth = 0.9) +
  # geom_point(aes(fill = dir),
  #            shape = 21, size = 3.2, stroke = 0.6) +
    
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high, colour = dir),
                  fatten = 4, shape = 21, fill = "white", linewidth = 0.8, stroke = 1.5)+
  # symmetric axis, small headroom
  scale_y_continuous(limits = c(-0.6, 0.6),
                     breaks = seq(-0.5, 0.5, 0.25),
                     expand = expansion(mult = c(0.02, 0.05))) +
  # harmonised colours (match rest of fig)
  scale_colour_manual(values = c(pos = "black",      # summer > winter
                                 neg = "black")) +
  # scale_fill_manual(values   = c(pos = "darkorange",
  #                                neg = "steelblue")) +
  labs(x = "Elevation (m)",
       y = "Abundance change") +
  theme_classic() +
  theme(axis.text  = element_text(colour = "black", size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 14),
        legend.position = "none"))

lm(mean_delta ~ elevation, data = delta_norm_summary) %>% broom::tidy()

#ggsave(plot = delta_community_norm_plot, filename = "03 - Output/Figures/norm_delta_community_abundance.tiff", width = 9, height = 10, dpi = 300)

#=== betapart

abund_matrix <- norm_all_N %>%  # from your normalized abundance pipeline
  group_by(sample, elev_band, season, species) %>%
  summarise(abund = sum(norm_pop), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = abund, values_fill = 0)

# abund_matrix <- all_N %>% # for non-normalised version
#   group_by(sample, elev_band, season, species) %>%
#   summarise(abund = sum(population), .groups = "drop")%>%
#   pivot_wider(names_from = species, values_from = abund, values_fill = 0)

beta_abund <- abund_matrix %>%
  group_by(elev_band, sample) %>%
  nest() %>%
  mutate(
    beta_parts = map(data, ~ {
      df <- .x %>%
        arrange(season) %>%
        select(where(is.numeric))
      
      if (nrow(df) < 2) return(tibble(turnover = NA, nestedness = NA, total = NA))
      
      beta <- beta.pair.abund(df, index.family = "bray")
      tibble(
        turnover = beta$beta.bray.bal[1],
        nestedness = beta$beta.bray.gra[1],
        total = beta$beta.bray[1]
      )
    })
  ) %>%
  unnest(beta_parts) %>%
  select(-data)

beta_abund %>% ungroup() %>%  summarise(mean_total = mean(total),
                         total_low = ci(total, ci = 0.89)[[2]],
                         total_high = ci(total, ci = 0.89)[[3]],
                         mean_turn = mean(turnover),
                         turn_low = ci(turnover, ci = 0.89)[[2]],
                         turn_high = ci(turnover, ci = 0.89)[[3]],
                         mean_nest = mean(nestedness),
                         nest_low = ci(nestedness, ci = 0.89)[[2]],
                         nest_high = ci(nestedness, ci = 0.89)[[3]])

beta_summary <- beta_abund %>%
  group_by(elev_band) %>%
  summarise(
    mean_turnover = mean(turnover, na.rm = TRUE),
    ci_turnover = list(ci(turnover, ci = 0.89)),
    mean_nestedness = mean(nestedness, na.rm = TRUE),
    ci_nestedness = list(ci(nestedness, ci = 0.89)),
    mean_total = mean(total, na.rm = TRUE),
    ci_total = list(ci(total, ci = 0.89)),
    .groups = "drop"
  ) %>%
  mutate(elevation = elev_pred) %>%
  unnest_wider(ci_turnover, names_sep = "_") %>%
  unnest_wider(ci_nestedness, names_sep = "_") %>%
  unnest_wider(ci_total, names_sep = "_")


means <- beta_summary %>%
  select(elevation, starts_with("mean_")) %>%
  pivot_longer(
    cols = -elevation,
    names_to = "metric",
    names_prefix = "mean_",
    values_to = "value"
  )


ci_lows <- beta_summary %>%
  select(elevation, starts_with("ci_"), ends_with("_CI_low")) %>%
  pivot_longer(
    cols = -elevation,
    names_to = "metric",
    names_pattern = "ci_(.*)_CI_low",
    values_to = "low"
  )

ci_highs <- beta_summary %>%
  select(elevation, starts_with("ci_"), ends_with("_CI_high")) %>%
  pivot_longer(
    cols = -elevation,
    names_to = "metric",
    names_pattern = "ci_(.*)_CI_high",
    values_to = "high"
  )


beta_long <- reduce(
  list(means, ci_lows, ci_highs),
  left_join,
  by = c("elevation", "metric")
) %>%
  mutate(
    metric = recode(metric,
                    turnover = "Turnover",
                    nestedness = "Nestedness",
                    total = "Total"),
    metric = factor(metric, levels = c("Turnover", "Nestedness", "Total"))
  )



(beta_facet_plot <- ggplot(beta_long, aes(x = elevation, y = value, col = metric)) +
  geom_pointrange(aes(ymin = low, ymax = high),
                  fatten = 4, shape = 21, fill = "white", linewidth = 0.8, stroke = 1.5) +
    #geom_line()+
  labs(x = "Elevation (m)",
       y = expression(atop(beta*"-diversity", "")),
       fill = "Metric",
       col = "Metric") +
    scale_color_manual(values = c("#1B9E77", "#7570B3", "black"))+
    scale_fill_manual(values = c("#1B9E77", "#7570B3", "black"))+
  theme_classic() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        panel.spacing = unit(1, "lines"),
        axis.text = element_text(colour = "black", size = 12),
        axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size =12), 
        legend.title = element_blank(),
        legend.position = c(0.5,0.8),
        legend.direction = "horizontal",
        legend.box = "horizontal"))

#ggsave(plot = beta_facet_plot, filename = "03 - Output/Figures/beta_diversity_baselga.tiff", width = 9, height = 10, dpi = 300)



(wrapped_plot_comm <- plot_grid(normalised_total_abund, delta_community_norm_plot, beta_facet_plot, 
                                labels = c("A)", "B)", "C)"),
                                label_size = 14,
                                ncol = 1))



# 1) Give them all a consistent left margin (optional but helps)
padL <- margin(5.5, 5.5, 5.5, 20)  # t, r, b, l
normalised_total_abund    <- normalised_total_abund    + theme(plot.margin = padL)
delta_community_norm_plot <- delta_community_norm_plot + theme(plot.margin = padL)
beta_facet_plot           <- beta_facet_plot           + theme(plot.margin = padL)


# 2) Align by the left axis, then stack
al <- align_plots(
  normalised_total_abund,
  delta_community_norm_plot,
  beta_facet_plot,
  align = "v",
  axis  = "l"
)

(wrapped_plot_comm <- plot_grid(
  al[[1]], al[[2]], al[[3]],
  labels = c("A)", "B)", "C)"),
  label_size = 14,
  ncol = 1
))



ggsave(plot = wrapped_plot_comm, filename = "03 - Output/Figures/community_level_result.tiff", width = 12, height = 10, dpi = 300)




