# --- Packages ---
library(tidyverse)
library(sf)
library(terra)
library(ggspatial)
library(png)
library(grid)
library(viridis)

# --- Output device ---
ragg::agg_tiff(
  "/Users/alejandrofp/Library/CloudStorage/OneDrive-JamesCookUniversity/Postdoc/Projects/10 - Bird elevational migration/awt_birds_altitudinal_migration/03 - Output/map.tiff",
  width = 168, height = 130, units = "mm", res = 350
)

# --- Paths (unchanged) ---
dem_path <- "/Users/alejandrofp/Library/CloudStorage/OneDrive-JamesCookUniversity/PhD OneDrive folder/PhD - Finished Projects/community_reshuffling/simulation_setup_and_first_model_with_all_spp/spatial_data/QLD_DEM_3SEC/awt_3s_DEM.tif"
r_path   <- "/Users/alejandrofp/Library/CloudStorage/OneDrive-JamesCookUniversity/PhD OneDrive folder/PhD - Finished Projects/abundance_suitability_hypothesis/spatial_data/predictors_1km/predictors/asc.tif"
awt_path <- "/Users/alejandrofp/Library/CloudStorage/OneDrive-JamesCookUniversity/PhD OneDrive folder/PhD - Spatial data/subregion_shapefile/disolved.shp"
pts_path <- "/Users/alejandrofp/Library/CloudStorage/OneDrive-JamesCookUniversity/PhD active folder/PhD - projects/climate_change_and_cyclones/BOM_climate/data/sites_geolocation.csv"
aus_png  <- "/Users/alejandrofp/Library/CloudStorage/OneDrive-JamesCookUniversity/PhD OneDrive folder/PhD - Finished Projects/community_reshuffling/Manuscript/delaFuente_et_al_2022_Diversity_and_Distributions/DDI_final_versions/aus_awt2.png"

# --- Read data ---
dem <- terra::rast(dem_path)
ref <- terra::rast(r_path)                         # for CRS
awt <- sf::st_read(awt_path, quiet = TRUE)

site_sel <- c("AU10A","AU10B","AU1A","AU2A","AU4A","AU6A","AU8A",
              "CU10A","CU10B","CU12A","CU1A","CU4A","CU6A","CU8A")

points_raw <- read.csv(pts_path) |>
  dplyr::filter(site %in% site_sel)
pts_wgs <- sf::st_as_sf(points_raw, coords = c("long","lat"), crs = 4326)

# --- Project to target CRS (this is very likely lon/lat degrees) ---
target_crs <- terra::crs(ref, proj = TRUE)
dem <- terra::project(dem, target_crs, method = "bilinear")
awt <- sf::st_transform(awt, target_crs)
pts <- sf::st_transform(pts_wgs, target_crs)   # retains 'site' col

# --- Zoom/crop DEM around sites ---
buffer_km <- 25
buf <- sf::st_buffer(sf::st_union(pts), dist = buffer_km * 1000)
dem_crop <- terra::crop(dem, terra::vect(buf))
dem_mask <- terra::mask(dem_crop, terra::vect(awt))
dem_fast <- terra::aggregate(dem_mask, fact = 3, fun = mean, na.rm = TRUE)

dem_df <- as.data.frame(dem_fast, xy = TRUE, na.rm = TRUE)
names(dem_df)[3] <- "elev"

pts_df <- pts |>
  sf::st_coordinates() |>
  as.data.frame() |>
  dplyr::bind_cols(sf::st_drop_geometry(pts))

# --- Label positions (no polygons) ---
# Tag groups: AU* = Atherton (south); CU* = Carbine (north)
pts_grouped <- pts |>
  dplyr::mutate(
    region = ifelse(grepl("^AU", site), "Atherton Uplands", "Carbine Uplands")
  )

# A point inside each cluster
label_pts <- pts_grouped |>
  dplyr::group_by(region) |>
  summarise(geometry = sf::st_union(geometry), .groups = "drop") |>
  sf::st_point_on_surface()

# Convert to coords
label_df <- label_pts |>
  sf::st_coordinates() |>
  as.data.frame() |>
  dplyr::bind_cols(region = label_pts$region) |>
  dplyr::rename(X = X, Y = Y)

# --- Compute km->degree offsets and nudge labels inside the frame ---
# Map bounds (degrees)
xmin <- min(dem_df$x); xmax <- max(dem_df$x)
ymin <- min(dem_df$y); ymax <- max(dem_df$y)
pad  <- 0.01  # degree padding from edges (~1.1 km N-S)

# km we want to nudge
dx_km <- 3   # ~3 km horizontally
dy_km <- 14   # ~2 km vertically

# Convert km to degrees at each label's latitude
km_to_deg_lat <- function(km) km / 111.0
km_to_deg_lon <- function(km, lat_deg) km / (111.0 * cos(lat_deg * pi/180))

label_df <- label_df |>
  mutate(
    dx = km_to_deg_lon(dx_km, Y),
    dy = km_to_deg_lat(dy_km),
    X = case_when(
      region == "Atherton Uplands" ~ X + dx,  # east
      region == "Carbine Uplands"   ~ X - dx,  # west
      TRUE ~ X
    ),
    Y = case_when(
      region == "Atherton Uplands" ~ Y - dy,  # south
      region == "Carbine Uplands"   ~ Y + dy,  # north
      TRUE ~ Y
    ),
    # keep inside plotting window
    X = pmin(pmax(X, xmin + pad), xmax - pad),
    Y = pmin(pmax(Y, ymin + pad), ymax - pad)
  )

# --- Theme ---
base_theme <- theme_bw() +
  theme(
    axis.text  = element_text(colour = "black", size = 8),
    axis.title = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.key.height = unit(3.0, "mm"),
    legend.key.width  = unit(3.5, "mm"),
    plot.margin = margin(5.5, 12, 5.5, 5.5)
  )

# --- Map ---
elev_plot <-
  ggplot() +
  geom_raster(data = dem_df, aes(x = x, y = y, fill = elev)) +
  geom_sf(data = awt, fill = NA, colour = "grey25", linewidth = 0.3) +
  geom_point(data = pts_df, aes(x = X, y = Y),
             shape = 21, size = 1.6, stroke = 0.6, fill = "white") +
  geom_label(data = label_df, aes(x = X, y = Y, label = region),
             size = 3.2, fontface = "bold",
             label.padding = unit(0.15, "lines"),
             label.size = 0.2, fill = scales::alpha("white", 0.85)) +
  labs(x = "Longitude", y = "Latitude", fill = "Elevation\n(m a.s.l.)") +
  coord_sf(xlim = range(dem_df$x), ylim = range(dem_df$y), expand = FALSE) +
  annotation_scale() +
  annotation_north_arrow(
    height = unit(0.8, "cm"),
    width  = unit(0.8, "cm"),
    location = "bl",
    pad_x = unit(0.15, "in"),   # move arrow left
    pad_y = unit(0.4, "in")
  ) +
  # Brighter = lowlands, darker = uplands
  scale_fill_viridis_c(option = "rocket", direction = -1) +
  guides(fill = guide_colorbar(title.position = "top",
                               barheight = unit(28, "mm"))) +
  base_theme

# --- Draw plot, then overlay AUS inset PNG ---
print(elev_plot)

aus_img <- png::readPNG(aus_png)
grid::grid.raster(aus_img, x = 0.82, y = 0.82, width = 0.22)

dev.off()
