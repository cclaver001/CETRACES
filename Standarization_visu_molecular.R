####################################
# Standardize sighting & eDNA data #
####################################

library(readxl)
library(ggspatial)
library(ggplot2)
library(sf)
library(ggpubr)

# 1. Create the spatial grid 
# Define grid size in degrees
lat_grid_size <- 0.7  #0.8, 0.4
lon_grid_size <- 0.6 #1, 0.5

# Create a grid using `sf` package
xmin <- -7.5
xmax <- -1.25
ymin <- 43.2
ymax <- 48

# Number of grid cells in x and y direction
n_lon <- ceiling((xmax - xmin) / lon_grid_size)
n_lat <- ceiling((ymax - ymin) / lat_grid_size)

# Create grid cells and set the same CRS as data
grid <- st_make_grid(st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)),cellsize = c(lon_grid_size, lat_grid_size), what = "polygons")

# Create an `sf` object and set CRS
grid_sf <- st_as_sf(grid) %>% st_set_crs(4326)

# Ensure unique IDs: Label each cell with row and column indices
grid_df <- grid_sf %>% mutate(
  row = rep(1:ceiling((ymax - ymin) / lat_grid_size), each = ceiling((xmax - xmin) / lon_grid_size)),  # row index
  column = rep(1:ceiling((xmax - xmin) / lon_grid_size), times = ceiling((ymax - ymin) / lat_grid_size)),  # column index
  grid_id = paste0("R", row, "C", column)  # unique identifier
)

# View the resulting data frame to ensure proper grid_id assignment
grid_df$geometry <- st_geometry(grid_df)
grid_centroids <- grid_df %>%   st_transform(4326) %>%  mutate(centroid = st_centroid(geometry)) #CALCULATE CENTROIDS TO NAME EACH CELL

ggplot() +
  annotation_map_tile(zoom = 8, type = "osm") + 
  geom_sf(data = grid_df, fill = "transparent", color = "grey30") +  # Setting color to NA removes cell outlines
  ggspatial::geom_spatial_point(data = obs_data[obs_data$code_esp=="water",],  aes(x = Longitude, y = Latitude), color = "grey30", size = 0.05, shape = 19) +
  ggspatial::geom_spatial_point(data = obs_data[obs_data$code_esp!="water"& obs_data$type =="target",],  aes(x = Longitude, y = Latitude), color = "green", size = 2, shape = 19) +
  geom_spatial_point(data = effort_mol, aes(x = Longitud, y = Latitud), size=2, color="red", alpha=0.5)+ theme_minimal()+
  geom_sf_text(data = grid_centroids, aes(geometry = centroid, label = grid_id)) +
  theme_minimal()

# Join data with grid to get the grid_id for each data point

# OBSERVATION DATA
obs_data_sf <- st_as_sf(obs_data, coords = c("Longitude", "Latitude"), crs = 4326)
obs_data_sf <- st_transform(obs_data_sf, st_crs(grid_df))
obs_data_sf <- st_join(grid_df, obs_data_sf)
obs_data_df <- obs_data_sf %>% st_drop_geometry() |> select(grid_id, year, code_esp) |>  (\(x) filter(x, complete.cases(x)))()
obs_data_df$Class = ifelse(obs_data_df$code_esp =="water", "effort", "sighting")
obs_data_count <- obs_data_df |>  group_by(grid_id, year, Class) |>
  summarise(n = n(), .groups = "drop") |> tidyr::pivot_wider(names_from = Class, values_from = n, values_fill = 0) |> mutate( ratio = sighting*100 / effort)
obs_data_count = merge(grid_df, obs_data_count, by= "grid_id")
obs_data_count_tot <- obs_data_df |>  group_by(grid_id, Class) |>  summarise(n = n(), .groups = "drop") |> tidyr::pivot_wider(names_from = Class, values_from = n, values_fill = 0) |> mutate( ratio = sighting*100 / effort)
obs_data_count_tot = merge(grid_df, obs_data_count_tot, by= "grid_id")
# plot visual effort
ggplot() +
  annotation_map_tile(zoom = 8, type = "cartolight") +  # OpenStreetMap as the base map
  geom_sf(data = obs_data_count_tot, aes(fill = ratio), color = "transparent", alpha = 0.7) +  # Plot the grid
  scale_fill_viridis_c(option = "F", na.value = "transparent", begin = 1, end = 0) +  # Viridis palette
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", title ="Encounter rate (sightings per unit effort) fo the 8-year period (2018-2025")
ggplot() +
  annotation_map_tile(zoom = 8, type = "cartolight") +  # OpenStreetMap as the base map
  geom_sf(data = obs_data_count, aes(fill = ratio), color = "transparent", alpha = 0.7) +  # Plot the grid
  scale_fill_viridis_c(option = "F", na.value = "transparent", begin = 1, end = 0) +  # Viridis palette
  theme_bw() +
  facet_wrap(~year, ncol=2)+
  labs(x = "Longitude", y = "Latitude")

# MOLECULAR DATA
effort_mol <- read_excel("Species specific assays/dPCR analysis/Muestras JUVENA18-25/dPCR_Cetacea_details.xlsx")
effort_mol = effort_mol |> select(Campaña, Latitud, Longitud, Estacion)
effort_mol = unique(effort_mol)
effort_mol <- effort_mol %>%  mutate(an = gsub("Juvena_", "", Campaña))
effort_mol$code_esp = paste0("WATER_",c(1:dim(effort_mol)[1]))
effort_mol_sf <- st_as_sf(effort_mol, coords = c("Longitud", "Latitud"), crs = 4326)
effort_mol_sf <- st_transform(effort_mol_sf, st_crs(grid_df))
effort_mol_sf <- st_join(grid_df, effort_mol_sf)
effort_mol_sf$Class = "molecular_effort"
effort_mol_df <- effort_mol_sf %>% st_drop_geometry() |> select(grid_id, an, code_esp, Class) |>  (\(x) filter(x, complete.cases(x)))()

molecular <- read_excel("Species specific assays/dPCR analysis/dPCR_cetacean_results.xlsx")
molecular = merge(molecular, effort_mol[,2:4], by= "Estacion")
molecular <- molecular %>%  mutate(an = gsub("Juvena_", "", Campaña))
molecular_sf <- st_as_sf(molecular, coords = c("Longitud", "Latitud"), crs = 4326)
molecular_sf <- st_transform(molecular_sf, st_crs(grid_df))
molecular_sf <- st_join(grid_df, molecular_sf)
molecular_sf$Class = "molecular_data"
molecular_df <- molecular_sf %>% st_drop_geometry() |> select(grid_id, an, Class) |> unique() |>  (\(x) filter(x, complete.cases(x)))()
molecular_df$code_esp = "any detection" #solo contamos una detecion por muestra para mantener el porcentaje <100%

molecular_data = rbind(effort_mol_df, molecular_df)
molecular_count <- molecular_data |>
  group_by(grid_id, an, Class) |>
  summarise(n = n(), .groups = "drop") |> 
  tidyr::pivot_wider(names_from = Class, values_from = n, values_fill = 0) |> 
  mutate(
    ratio = molecular_data*100 / molecular_effort
  )

molecular_count = merge(grid_df, molecular_count, by= "grid_id")
ggplot() +
  annotation_map_tile(zoom = 8, type = "stamenterrain") +  # "cartolight", "stamenbw"
  geom_sf(data = molecular_count, aes(fill = ratio), color = "transparent", alpha = 0.7) +  # Plot the grid
  scale_fill_viridis_c(option = "F", na.value = "transparent", begin = 1, end = 0.2) +  # Viridis palette
  theme_minimal() +
  facet_grid(~an)+
  labs(x = "Longitude", y = "Latitude", title = "Molecular detections (observation/effort)")


