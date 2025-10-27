library(tidyverse)

#=== Read data

pos <- read.csv("03 - Output/Posterior summary/combined_posteriors.csv")
traits <- read.csv("01 - Data/traits/traits_curated.csv")
traits_raw <- read.csv("01 - Data/traits/bird_awt_traits.csv")
posterior_migration <- read.csv("03 - Output/Posterior traits data/migration_posterior.csv")

#=== Standardise col names for joining
colnames(traits)[colnames(traits) == "spp"] <- "species"
colnames(traits_raw)[colnames(traits_raw) == "spp"] <- "species"

#=== Prepare and join posterior summary to traits
df <- posterior_migration[,-1]
spp_list <- unique(df$species)

traits <- traits %>% filter(species %in% spp_list)
traits_raw <- traits_raw %>% filter(species %in% spp_list)

df <- left_join(df, traits)
df <- left_join(df, traits_raw)

df <- df %>% mutate(forage_ground = forage_G + forage_V,
                    forage_canopy = forage_B + forage_C,
                    diet_protein = diet_invert + diet_insect,
                    diet_veg = diet_frug + diet_flower + diet_seed,
                    rainforest_spec = rainforest_specialization..1.occasionally...6.obligate.)

df_ecology <- df[,c(4,1,2,5,6,8,10,17,23,27,48:50)]

df_ecology <- df_ecology %>% drop_na()


cor(df_ecology[,-c(1:3)])

df_ecology %>% ggplot(aes(diet_protein, posterior))+
  geom_point()+
  geom_smooth(se = F, method = "lm")+
  facet_wrap(~elev_cat, scales = "free")

df_ecology %>% ggplot(aes(diet_veg, posterior))+
  geom_point()+
  geom_smooth(se = F, method = "lm")+
  facet_wrap(~elev_cat)


# BAYESIAN META-ANALYSIS
