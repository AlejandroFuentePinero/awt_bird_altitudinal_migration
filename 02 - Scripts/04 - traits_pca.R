library(tidyverse)
library(factoextra)
library(FactoMineR)


traits <- read.csv("01 - Data/traits/bird_awt_traits.csv")

traits <- traits %>% filter(!spp %in% c("LKOOK"))

# diet

diet <- traits[,c(32:36)]

rownames(diet) <- traits$spp

pca_diet <- PCA(diet)

summary(pca_diet)

diet_pca <- pca_diet$ind$coord[,1]

fviz_pca_ind(pca_diet)


# forage

forage <- traits[,c(27:31)]

rownames(forage) <- traits$spp

pca_forage <- PCA(forage)

summary(pca_forage)

forage_pca <- pca_forage$ind$coord[,1]

fviz_pca_ind(pca_forage)

# habitat

habitat <- traits[,c(15,17,18)]

rownames(habitat) <- traits$spp


pca_habitat <- PCA(habitat)


habitat_pca <- pca_habitat$ind$coord[,1]

fviz_pca_ind(pca_habitat)


# Curated traits data

traits_curated <- tibble(family = traits$FAMILY,
                         spp = traits$spp,
                         pop_trend = traits$population_trend,
                         elev_cat = traits$elevation.pref.1...low...5.up.,
                         repro_seasonality = traits$rep_seasonality..estimated.number.of.births.per.year.,
                         pot_dispersal = traits$pot_dispersal,
                         habitat_special = habitat_pca,
                         forage = forage_pca,
                         diet = diet_pca)

write.csv(traits_curated, "01 - Data/traits_curated.csv")

