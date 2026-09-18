library(ggplot2)
library(dplyr)
library(ggspatial)
library(ggOceanMaps)
library(readxl)
library(viridis)
library(tidyr)
library(scales)
library(sf)

metaprobes <- read_excel("~/OneDrive - AZTI/PROYECTOS/BIOcean5D/BayofBiscay_eDNA_Cetaceos/Tables.xlsx", sheet = "Table S2 aux")
basemap(data = dt, bathymetry = T, bathy.style = "poly_grays") + 
  ggspatial::geom_spatial_point(data = metaprobes, aes(x = Longitude, y = Latitude, color = Material), size = 10, shape = 20, alpha = 0.7) + theme(legend.position = "right")

##############################
# Process cetacean sightings #
##############################

obs_data <- read.csv("C:/Use/OneDrive - AZTI/ESCRITORIO/Cristina/Datos/BIOcean5D/Marine_mammals/Sigthing_data/All_sighting_data_AZTI/obs_effort.csv", header = F)
names(obs_data) = obs_data[1,]
obs_data=obs_data[-1,]
obs_data$year.y = as.numeric(obs_data$year.y)
obs_data <- obs_data %>%  mutate(number = as.numeric(as.character(number)), Longitude = as.numeric(as.character(Longitude)), Latitude = as.numeric(as.character(Latitude)))
obs_data = obs_data |> filter(Survey=="JUVENA") |> filter(year.y > 2017)
unique(obs_data$code_esp)
obs_data <- obs_data |>  mutate(code_esp = ifelse(is.na(code_esp), "water", code_esp))
obs_data <- obs_data |>  mutate(number = ifelse(is.na(number), "0", number))
obs_data <- obs_data %>%   mutate(number = as.numeric(number))
names(obs_data)[16] <- "year"
names(obs_data)[15] <- "code_leg"

target = c("GLOMEL", "TURTRU", "ZIPCAV", "DELDEL", "STECOE", "BALPHY", "water")
othertarget = c("DELSPP", "BALSPP", "ZIPSPP", "SMADEL", "BALACU")

obs_data <- obs_data |>  filter(code_esp %in% c(target, othertarget)) |>
  mutate(
    type = case_when(
      code_esp %in% target      ~ "target",
      code_esp %in% othertarget ~ "other"
    )) |>  select(Survey, year, code_leg, Latitude, Longitude, code_esp, number, type, leg_heure)

# because dPCR assays fails to differentiate these two dolphin species
obs_data$code_esp2 <- ifelse(as.character(obs_data$code_esp) %in% c("DELDEL", "STECOE"), "DELDEL&STECOE",as.character(obs_data$code_esp))

# export clean dataset 
write.csv(obs_data, file = "Visu_data_JUVENA_18-25.csv", row.names = F)

# Figure 1a. Effort and sampling area
basemap(data = dt, bathymetry = F ) + 
  ggspatial::geom_spatial_point(data = obs_data[obs_data$code_esp =="water", ], aes(x = Longitude, y = Latitude, color = factor(year)), size = 1, shape = 20, alpha = 0.3) + 
  scale_color_viridis_d(name = "Year", option = "turbo") + theme(legend.position = "right")

# Figure S1. Plot visual effort (black) and sightings (red)
basemap(data = dt, bathymetry = T) + #bathy.style = "poly_grays"
  ggspatial::geom_spatial_point(data = obs_data[obs_data$code_esp =="water", ], aes(x = Longitude, y = Latitude), color = "black", size = 0.05, shape = 20, alpha = 0.5) + 
  ggspatial::geom_spatial_point(data = obs_data[obs_data$type =="target" & obs_data$code_esp!= "water", ], aes(x = Longitude, y = Latitude), color = "red", size = 2, shape = 8, alpha = 0.5) + 
  facet_wrap(~year, ncol = 2) + theme(legend.position = "none")

# Figure 3. Plot sightings per target species
basemap(data = dt, bathymetry = T, bathy.style = "contour_gray") + 
  ggspatial::geom_spatial_point(data = obs_data[obs_data$type =="target" & obs_data$code_esp!= "water", ], aes(x = Longitude, y = Latitude), color = "red", size = 2, shape = 8, alpha = 0.5) + 
  facet_grid(year ~ code_esp) +
  theme_classic()+
  theme(legend.position = "none")
# Figure 3. Plot sightings per target species WITH effort
# Datos de esfuerzo (sin columna code_esp)
obs_data$code_esp <- factor(obs_data$code_esp,  levels = c("BALPHY", "GLOMEL", "DELDEL", "STECOE", "TURTRU", "ZIPCAV", "water"))
obs_data$code_esp2 <- factor(obs_data$code_esp2,  levels = c("BALPHY", "GLOMEL", "DELDEL&STECOE", "TURTRU", "ZIPCAV", "water"))
obs_water <- obs_data[obs_data$code_esp2 == "water", ]
obs_water$code_esp2 <- NULL
basemap(data = dt, bathymetry = TRUE, bathy.style = "contour_gray") +
  ggspatial::geom_spatial_point(
    data = obs_water,
    aes(x = Longitude, y = Latitude),
    color = "gray10", size = 0.1, shape = 20, alpha = 0.3,
    inherit.aes = FALSE
  ) +
  ggspatial::geom_spatial_point(
    data = obs_data[obs_data$type == "target" & obs_data$code_esp2 != "water", ],
    aes(x = Longitude, y = Latitude),
    color = "red", size = 2, shape = 8, alpha = 0.5
  ) +
  facet_grid(year ~ code_esp2) +
  theme_classic() +
  theme(legend.position = "none")

# metrics about sightings
obs_data |> filter(type == "target") |> filter(number > 0) %>% group_by(code_esp) %>%
  summarise(n_records = n(), sum_number = sum(number, na.rm = TRUE) ) %>% arrange(desc(sum_number))   

rows_per_year <- obs_data %>%  filter(type == "target") %>%  count(year, name = "n_rows")
ggplot(rows_per_year, aes(x = factor(year), y = n_rows, group = 1)) +   geom_line(linewidth = 1, color = "black") +  geom_point(size = 2.5, color = "black") +  labs(x = "Year", y = "Effort (leg)") +  theme_minimal(base_size = 12) + theme(panel.grid.major.x = element_blank())

by_year_species <- obs_data %>%   filter(type == "target") %>%  group_by(year, code_esp) %>%  summarise(    n_records = n(),    sum_number = sum(number, na.rm = TRUE),    .groups = "drop")
ggplot( by_year_species %>% filter(code_esp != "water"),  aes(x = factor(year), y = n_records, fill = code_esp)) +  geom_col(color = "grey30", width = 0.8) +
  labs(x = "Year", y = "Number of sightings", fill = "Species") +
  scale_fill_viridis_d( name = "Species", option = "turbo", na.translate = FALSE) +
  theme_minimal(base_size = 12) +  theme(panel.grid.major.x = element_blank())

species_levels <- obs_data %>%  filter(type == "target", code_esp != "water") %>%  distinct(code_esp) %>% pull(code_esp) %>% sort()
years_levels <- sort(unique(obs_data$year))
by_year_species_complete <- by_year_species %>%  filter(code_esp != "water") %>%  complete(    year = years_levels,  code_esp = species_levels, fill = list(n_records = 0, sum_number = 0))
ggplot(by_year_species_complete, aes(x = factor(year), y = n_records, fill = code_esp)) + geom_col(width = 0.75, color = "grey30", show.legend = FALSE) +  facet_wrap(~ code_esp, scales = "free_y", ncol=1) +
  scale_y_continuous(  breaks = function(x) {
      m <- suppressWarnings(max(x, na.rm = TRUE))
      if (!is.finite(m) || m <= 0) {
        return(0)  # especie sin datos > 0
      }
      b <- c(0, round(m/3), round(2*m/3), m)
      unique(b)
    }, labels = scales::number_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) + 
  scale_fill_viridis_d( name = "Species", option = "turbo", na.translate = FALSE) +
  labs( x = "Year", y = "Number of sightings") + theme_minimal(base_size = 12) + theme( panel.grid.major.x = element_blank(),  axis.text.x = element_text(angle = 45, hjust = 1))

# Plot sightings per species (other)
basemap(data = dt, bathymetry = T, bathy.style = "contour_gray") + 
  ggspatial::geom_spatial_point(data = obs_data[obs_data$type =="other", ], aes(x = Longitude, y = Latitude), color = "red", size = 2, shape = 8, alpha = 0.5) + 
  facet_grid(year ~ code_esp) + theme(legend.position = "none")

# Plot rorquals
metab_data = metabarcoding_eDNA |> filter(Species == "Balaenoptera_acutorostrata") |> select(Sample, Latitud, Longitud, Campaña)
metab_data$code_esp = "BALACU"
metab_data$method = "molecular"  
names(metab_data) = c("code_leg", "Latitude",  "Longitude", "Campaña", "code_esp", "method")
metab_data$year = c(2019,2023)
metab_data$Campaña <-NULL

dPCR_data = molecular_data |> mutate(year = as.numeric(unlist(an)))  |> filter(code_esp == "BALPHY", an %in% c(2019, 2023)) |> select(Estacion, Latitude, Longitude, an, code_esp)
names(dPCR_data) = c("code_leg", "Latitude",  "Longitude", "year", "code_esp")
dPCR_data$method = "molecular"
  
plot_data <- obs_data |>  dplyr::filter(code_esp %in% c("BALACU", "BALPHY"), year %in% c(2019, 2023) ) |>
  dplyr::mutate(code_esp = factor(code_esp, levels = c("BALACU", "BALPHY")),  year= factor(year,     levels = c(2019, 2023))) |> select(code_leg, code_esp ,Longitude, Latitude, year)
plot_data$method = "visu"
plot_data <- bind_rows(mutate(metab_data, year = as.numeric(year)),  mutate(dPCR_data,  year = as.numeric(year)),plot_data %>% mutate(year = as.numeric(as.character(year))) )
metab_data <- NULL
dPCR_data <- NULL

plot_data <- plot_data %>%  mutate( year = factor(as.character(year), levels = c("2019", "2023")),code_esp = factor(code_esp, levels = c("BALACU", "BALPHY")) )
ggOceanMaps::basemap( data= dt, bathymetry  = TRUE, bathy.style = "contour_gray") +
  ggspatial::geom_spatial_point(data = plot_data,  aes(x = Longitude, y = Latitude, color = method, shape = method, alpha = method, size = method)) +
  scale_color_manual(values = c(molecular = "dodgerblue", visu = "#d62728")) + 
  scale_shape_manual(values = c(molecular = 16,  visu = 8)) + 
  scale_size_manual(values  = c(molecular = 3.5, visu = 2)) +
  scale_alpha_manual(values = c(molecular = 0.8, visu = 0.5)) +
  facet_grid(rows = vars(year), cols = vars(code_esp), drop = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold")
  ) + labs( x = NULL, y = NULL, shape = "Método")

# Plot pilot whales
dPCR_data =  molecular_data  |> mutate(year = as.numeric(unlist(an))) |> filter(Species == "Pilot whale", an %in% c(2019, 2020, 2021, 2022, 2024)) |> select(Estacion, Latitud, Longitud, an, Species)
names(dPCR_data) = c("code_leg", "Latitude",  "Longitude", "year", "code_esp")
dPCR_data$method = "molecular"
plot_data <- obs_data |>  dplyr::filter(code_esp %in% c("GLOMEL"), year %in% c(2019, 2020, 2021, 2022, 2024) ) |>
  dplyr::mutate(code_esp = factor(code_esp, levels = c("GLOMEL")),  year= factor(year,     levels = c(2019,  2020, 2021, 2022, 2024))) |> select(code_leg, code_esp ,Longitude, Latitude, year)
plot_data$method = "visu"
plot_data <- bind_rows(mutate(dPCR_data,  year = as.numeric(year)), plot_data %>% mutate(year = as.numeric(as.character(year))) )
dPCR_data <- NULL

plot_data <- plot_data %>%  mutate( year = factor(as.character(year), levels = c("2019", "2020", "2021", "2022", "2024")),code_esp = factor(code_esp, levels = c("GLOMEL")) )
plot_data$code_esp = "GLOMEL"
ggOceanMaps::basemap( data= dt, bathymetry  = TRUE, bathy.style = "contour_gray") +
  ggspatial::geom_spatial_point(data = plot_data,  aes(x = Longitude, y = Latitude, color = method, shape = method, alpha = method, size = method)) +
  scale_color_manual(values = c(molecular = "dodgerblue", visu = "#d62728")) + 
  scale_shape_manual(values = c(molecular = 16,  visu = 8)) + 
  scale_size_manual(values  = c(molecular = 3.5, visu = 2)) +
  scale_alpha_manual(values = c(molecular = 0.8, visu = 0.5)) +
  facet_grid(rows = vars(year), cols = vars(code_esp), drop = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold")
  ) + labs( x = NULL, y = NULL, shape = "Método")

# Plot bottlenose dolphin
dPCR_data =  molecular_data  |> mutate(year = as.numeric(unlist(an))) |> filter(Species == "Bottlenose dolphin", an %in% c(2019)) |> select(Estacion, Latitud, Longitud, an, Species)
names(dPCR_data) = c("code_leg", "Latitude",  "Longitude", "year", "code_esp")
dPCR_data$method = "molecular"
plot_data <- obs_data |>  dplyr::filter(code_esp %in% c("TURTRU"), year %in% c(2019) ) |>
  dplyr::mutate(code_esp = factor(code_esp, levels = c("TURTRU")),  year= factor(year,     levels = c(2019))) |> select(code_leg, code_esp ,Longitude, Latitude, year)
plot_data$method = "visu"
plot_data <- bind_rows(mutate(dPCR_data,  year = as.numeric(year)), plot_data %>% mutate(year = as.numeric(as.character(year))) )
dPCR_data <- NULL

plot_data <- plot_data %>%  mutate( year = factor(as.character(year), levels = c("2019")),code_esp = factor(code_esp, levels = c("TURTRU")) )
plot_data$code_esp = "TURTRU"
ggOceanMaps::basemap( data= dt, bathymetry  = TRUE, bathy.style = "contour_gray") +
  ggspatial::geom_spatial_point(data = plot_data,  aes(x = Longitude, y = Latitude, color = method, shape = method, alpha = method, size = method)) +
  scale_color_manual(values = c(molecular = "dodgerblue", visu = "#d62728")) + 
  scale_shape_manual(values = c(molecular = 16,  visu = 8)) + 
  scale_size_manual(values  = c(molecular = 3.5, visu = 2)) +
  scale_alpha_manual(values = c(molecular = 0.8, visu = 0.5)) +
  facet_grid(rows = vars(year), cols = vars(code_esp), drop = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold")
  ) + labs( x = NULL, y = NULL, shape = "Método")

#################################
# Process molecular data (dPCR) #
#################################

eDNA_detections <- read_excel("Species specific assays/dPCR analysis/dPCR_cetacean_results.xlsx")
eDNA_detections = eDNA_detections |> select(Estacion, Campaña, Species, Concentration)
names(eDNA_detections) = c("Estacion", "Survey", "code_esp", "Concentration")
eDNA_effort <- read_excel("Species specific assays/dPCR analysis/Muestras JUVENA18-25/dPCR_Cetacea_details.xlsx")
eDNA_effort = eDNA_effort |> select(Campaña, Estacion, Latitud, Longitud) |> unique()
names(eDNA_effort) = c( "Survey","Estacion", "Latitude", "Longitude")
molecular_data0 = merge(eDNA_effort, eDNA_detections, by=c("Estacion", "Survey"), all=T)
molecular_data0$code_esp = ifelse(is.na(molecular_data0$code_esp), "water", molecular_data0$code_esp)
molecular_data0$Concentration = ifelse(is.na(molecular_data0$Concentration), 0 , molecular_data0$Concentration)
molecular_data0 <- molecular_data0 %>%  mutate(an = sub(".*_", "", Survey))
molecular_data0$code_esp = ifelse(molecular_data0$code_esp =="Balaenoptera", "BALPHY", molecular_data0$code_es)
molecular_data0$code_esp = ifelse(molecular_data0$code_esp =="Pilot whale", "GLOMEL", molecular_data0$code_es)
molecular_data0$code_esp = ifelse(molecular_data0$code_esp =="Beaked whale", "ZIPCAV", molecular_data0$code_es)
molecular_data0$code_esp = ifelse(molecular_data0$code_esp =="Bottlenose dolphin", "TURTRU", molecular_data0$code_es)
molecular_data0$code_esp = ifelse(molecular_data0$code_esp =="Dolphin (Common and stripped)", "DELDEL&STECOE", molecular_data0$code_es)
molecular_data0$code_esp = factor(molecular_data0$code_esp, levels=c("BALPHY", "ZIPCAV","TURTRU","GLOMEL", "DELDEL&STECOE", "water"))
names(molecular_data0)
molecular_data0$code_esp2 = molecular_data0$code_esp
names(molecular_data0)[7] <- "year"

# Figure 1b. Effort and sampling area
basemap(data = dt, bathymetry = F ) + 
  ggspatial::geom_spatial_point(data = molecular_data0, aes(x = Longitude, y = Latitude, color = factor(year)), size = 4, shape = 20, alpha = 0.7) + 
  scale_color_viridis_d(name = "Year", option = "turbo") + theme(legend.position = "right")

# all eDNA
basemap(data = dt, bathymetry = T) +
  ggspatial::geom_spatial_point(data = molecular_data0[molecular_data0$code_esp!="water",],  aes(x = Longitude, y = Latitude), color = "red", size = 2, alpha = 0.9, shape = 19) + 
  facet_grid( ~ code_esp) + theme(legend.position = "none")

basemap(data = dt, bathymetry = T) + 
  ggspatial::geom_spatial_point(data = molecular_data0[molecular_data0$code_esp!="water",],  aes(x = Longitude, y = Latitude), color = "red", size = 2, shape = 19) + 
  facet_grid(code_esp ~ year) +  theme(legend.position = "none")

# both data
basemap(data = dt, bathymetry = T, bathy.style = "contour_gray") + 
  ggspatial::geom_spatial_point(data = obs_data[obs_data$code_esp2!="water"&obs_data$type=="target",],  aes(x = Longitude, y = Latitude), color = "darkred", size = 1, alpha =0.9, shape = 18) + 
  ggspatial::geom_spatial_point(data = molecular_data0[molecular_data0$code_esp2!="water",],  aes(x = Longitude, y = Latitude), color = "dodgerblue", size = 2, alpha=.9, shape = 19) + 
  facet_grid(code_esp2 ~ year) +  theme(legend.position = "none")

basemap(data = dt, bathymetry = T, bathy.style = "poly_grays") + 
  ggspatial::geom_spatial_point(data = obs_data[obs_data$code_esp2!="water",],  aes(x = Longitude, y = Latitude), color = "darkred", size = 1.5, alpha =0.9, shape = 18) + 
  ggspatial::geom_spatial_point(data = molecular_data0[molecular_data0$code_esp2!="water",],  aes(x = Longitude, y = Latitude), color = "red", size = 2, alpha=.9, shape = 19) + 
  facet_grid( ~ code_esp2) +  theme(legend.position = "none")

#############################
# visu vs molecular effort #
#############################

basemap(data = dt, bathymetry = F) + 
  ggspatial::geom_spatial_point(data = obs_data[obs_data$code_esp2=="water",],  aes(x = Longitude, y = Latitude), color = "grey35", size = 0.5, alpha =0.9, shape = 18) + 
  ggspatial::geom_spatial_point(data = molecular_data0,  aes(x = Longitude, y = Latitude), color = "dodgerblue", size = 2, alpha=.7, shape = 19) + 
  facet_wrap( ~ year, ncol= 2) +  theme(legend.position = "none")

#get info about environmental conditions of sightings
env_data <- read_excel("C:/Use/OneDrive - AZTI/ESCRITORIO/Cristina/Datos/BIOcean5D/Marine_mammals/Sigthing_data/All_sighting_data_AZTI/environment_data.xlsx")
unique(env_data$beaufort) #viento
unique(env_data$houle_hauteur) #ola
env_data$houle_hauteur = as.numeric(env_data$houle_hauteur)
env_data <- env_data %>%  mutate(
    douglas = case_when(
      is.na(houle_hauteur)            ~ NA_integer_,
      houle_hauteur < 0.10            ~ 0L,  # Calma (glassy/rippled)
      houle_hauteur < 0.50            ~ 1L,  # Rizada / ripples
      houle_hauteur < 1.25            ~ 2L,  # Marejadilla
      houle_hauteur < 2.50            ~ 3L,  # Marejada
      houle_hauteur < 4.00            ~ 4L,  # Fuerte marejada
      houle_hauteur < 6.00            ~ 5L,  # Gruesa
      houle_hauteur < 9.00            ~ 6L,  # Muy gruesa
      houle_hauteur < 14.0            ~ 7L,  # Arbolada
      houle_hauteur < 20.0            ~ 8L,  # Montañosa
      TRUE                            ~ 9L   # Enorme 
    ))
env_data = env_data |> filter(year > 2016) |> filter(month!=5) |> select("code_leg", "day","month", "year", "beaufort", "douglas", "cond_generale_babord_MM", "cond_generale_tribord_MM")
env_data %>%  filter(cond_generale_babord_MM != cond_generale_tribord_MM)

proof = merge(env_data, obs_data, by= c("code_leg", "year"), all.y= T)
proof$cond_generale_babord_MM  = ifelse(proof$cond_generale_babord_MM=="NA", "Bonne", proof$cond_generale_babord_MM)
proof$cond_generale_tribord_MM = ifelse(proof$cond_generale_tribord_MM =="NA","Bonne", proof$cond_generale_tribord_MM)
unique(proof$cond_generale_babord_MM)
unique(proof$cond_generale_tribord_MM)

pal_cond <- c("Excellente" = "#00A300", "Bonne"= "#00AA55", "Moyenne"= "#FFA500", "Mauvaise"= "#CC0000")
ord <- c("Excellente", "Bonne", "Moyenne", "Mauvaise")

proof <- proof %>% mutate(cond_bab = factor(cond_generale_babord_MM,  levels = ord, ordered = TRUE), cond_tri = factor(cond_generale_tribord_MM, levels = ord, ordered = TRUE), cond_final = factor(pmax(cond_bab, cond_tri), levels = ord, ordered = TRUE))

basemap(data = dt, bathymetry = FALSE) +
  ggspatial::geom_spatial_point(data = proof[proof$code_esp2 == "water", ], aes(x = Longitude, y = Latitude, color = cond_final), size = 1, alpha = 0.3, shape = 16) +
  #ggspatial::geom_spatial_point(data = molecular_data0[molecular_data0$code_esp=="water",],  aes(x = Longitude, y = Latitude), color = "black", size = 2, alpha=.7, shape = 8) + 
  ggspatial::geom_spatial_point(data = obs_data[obs_data$code_esp!="water",],  aes(x = Longitude, y = Latitude), color = "black", size = 2, alpha=.1, shape = 8) + 
  ggspatial::geom_spatial_point(data = molecular_data[molecular_data$code_esp!="water",],  aes(x = Longitude, y = Latitude), color = "dodgerblue", size = 2, alpha=.7, shape = 17) + 
  facet_wrap(~ year, ncol = 2) +  scale_color_manual(values = pal_cond, drop = FALSE) +  theme(legend.position = "right") +  labs(color = "Environmental / conditions")

# Associate hydrological stations to weather conditions
no_match_ids <- c(paste0("S", sprintf("%03d", 46:51)), paste0("S", sprintf("%03d", 57:64))) #no visual effort available
edna_ok <- molecular_data %>%  filter(!(Estacion %in% no_match_ids))
edna_excluded <- molecular_data %>%  filter(Estacion %in% no_match_ids) %>% mutate(dist_m = NA_real_)
proof_sf <- st_as_sf(proof, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>% st_transform(3857)
names(edna_ok)[4] = "Latitude"
names(edna_ok)[5] = "Longitude"
edna_sf  <- st_as_sf(edna_ok, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%  st_transform(3857)

proof_by_year <- split(proof_sf, proof_sf$year) #Para cada muestra eDNA, nearest point en proof dentro del MISMO año

nearest_idx <- vapply(seq_len(nrow(edna_sf)), function(i) {
  yy <- edna_sf$an[i]
  pyy <- proof_by_year[[as.character(yy)]]
  if (is.null(pyy) || nrow(pyy) == 0) return(NA_integer_)
  st_nearest_feature(edna_sf[i, ], pyy)[1]
}, integer(1))

proof_match <- bind_rows(lapply(seq_len(nrow(edna_sf)), function(i) {
  yy <- edna_sf$an[i]
  pyy <- proof_by_year[[as.character(yy)]]
  if (is.na(nearest_idx[i])) return(NULL)
  pyy[nearest_idx[i], ] %>%
    st_drop_geometry() %>%
    mutate(.row_edna = i)
}))

edna_joined <- edna_sf %>%
  mutate(.row_edna = seq_len(n()),
         .nearest_in_year = nearest_idx) %>%
  left_join(proof_match, by = ".row_edna", suffix = c("_edna", "_proof"))

# 5) Distancia (m) entre eDNA y el punto proof asignado (si existe)
#    Calculamos la geometría del punto matched en el mismo orden de edna_joined
geom_matched <- st_sfc(lapply(seq_len(nrow(edna_sf)), function(i) {
  yy <- edna_sf$an[i]
  pyy <- proof_by_year[[as.character(yy)]]
  if (is.na(nearest_idx[i])) return(st_point(c(NA_real_, NA_real_)))
  st_geometry(pyy[nearest_idx[i], ])[1][[1]]
}), crs = st_crs(edna_sf))

edna_joined <- edna_joined %>%
  mutate(dist_m = as.numeric(st_distance(st_geometry(edna_sf), geom_matched, by_element = TRUE))) %>%
  st_drop_geometry() %>%
  select(-.row_edna)  

# Cambia aquí si tu columna se llama distinto:
cond_col <- "cond_final"
cond_col <- "beaufort"
cond_col <- "douglas"

library(scales)
df_plot = edna_joined |> select(Estacion, Species, cond_final) |> unique() |> filter(Species !="water")
ggplot(df_plot, aes(x = cond_final)) +
  geom_bar(aes(y = after_stat(count / sum(count) * 100)),
           fill = "darkgreen") +  labs(
             x = "Environmental condition",
             y = "Molecular detections (%)"
           ) +   ylim(0, 100) +  theme_bw(base_size = 12)

df_plot <- edna_joined %>%  mutate(detected = Species != "water") %>%  filter(!is.na(.data[[cond_col]])) %>%
  group_by(.data[[cond_col]]) %>%
  summarise(
    n = n(),
    n_detect = sum(detected, na.rm = TRUE),
    prop_detect = n_detect / n,
    .groups = "drop"
  )

ggplot(df_plot, aes(x = .data[[cond_col]], y = prop_detect)) +
  geom_col(fill = "#2C7FB8") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Environmental conditions", y = "Water samples with positive cetacean DNA detections (%)") +
  theme_minimal(base_size = 12)

df_plot <- edna_joined %>%
  mutate(status = ifelse(Species == "water", "water", "no water")) %>%
  filter(!is.na(.data[[cond_col]])) %>%
  count(.data[[cond_col]], status, name = "n") %>%
  group_by(.data[[cond_col]]) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(df_plot, aes(x = .data[[cond_col]], y = prop, fill = status)) +
  geom_col(width = 0.8, color = "grey20") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = c("water" = "#1f78b4", "no water" = "#33a02c"),
    labels = c("water" = "Non detection", "no water" = "Detection")
  ) + 
  labs(x = "Environmental conditions", y = "Water samples (%)", fill = "") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#lo mismo para avistamientos 
visu = proof |> select(leg_heure, code_esp, cond_final) |> unique() |> filter(code_esp !="water")
ggplot(visu, aes(x = cond_final)) +
  geom_bar(aes(y = after_stat(count / sum(count) * 100)),
           fill = "darkgreen") +  labs(
    x = "Environmental condition",
    y = "Sightings (%)"
  ) +   ylim(0, 100) +  theme_bw(base_size = 12)

visu = proof |> select(code_leg, code_esp, cond_final) |> unique()
library(scales)
df_clean <- visu %>%  group_by(code_leg) %>%
  filter(
    # mantener filas no-water
    code_esp != "water" |
      # o mantener water solo si no hay ninguna detección
      all(code_esp == "water")
  ) %>% ungroup()
cond_col <- "cond_final"
df_plot <- df_clean %>%
  mutate(
    status = ifelse(code_esp == "water", "water", "detection")
  ) %>%
  filter(!is.na(.data[[cond_col]])) %>%
  distinct(code_leg, .data[[cond_col]], status) %>%  # 1 fila por code_leg
  count(.data[[cond_col]], status, name = "n") %>%
  group_by(.data[[cond_col]]) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(df_plot, aes(x = .data[[cond_col]], y = prop, fill = status)) +
  geom_col(width = 0.8, color = "grey20") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = c("water" = "#1f78b4", "detection" = "#33a02c"),
    labels = c("water" = "No detection", "detection" = "Detection")
  ) +
  labs(
    x = "Environmental conditions",
    y = "Legs (%)",
    fill = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


