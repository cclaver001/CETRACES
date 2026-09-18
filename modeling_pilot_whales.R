
###########################################
########## Bayesian models ###########
###########################################

pilot_visu <- effort_grid_year  %>% left_join(sightings_grid_year[sightings_grid_year$code_esp2== "GLOMEL", ], by = c("grid_id", "year")) %>%  mutate(ratio_sightings_km   = n_sightings   / effort_km, ratio_individuals_km = n_individuals / effort_km)
pilot_visu <- pilot_visu %>% mutate(code_esp2 = ifelse(is.na(code_esp2), "NA", code_esp2),  n_sightings   = ifelse(is.na(n_sightings),   0L, n_sightings),  n_individuals = ifelse(is.na(n_individuals), 0,  n_individuals), ratio_sightings_km   = ifelse(is.na(ratio_sightings_km),   0, ratio_sightings_km), ratio_individuals_km = ifelse(is.na(ratio_individuals_km), 0, ratio_individuals_km))
pilot_visu = pilot_visu %>% st_drop_geometry() %>% left_join(grid_with_effort_visu, by = c("grid_id")) %>% st_as_sf() 
ggplot() + geom_sf(data = pilot_visu, aes(fill = ratio_individuals_km), color = NA) +
  scale_fill_viridis_c(option = "A", direction = 1, begin = 1, end = 0.3, name = "Individual encounter rate
  (individuals·km⁻¹)") +
  { if (exists("land_ea")) geom_sf(data = land_ea, fill = "grey85", color = "grey70", linewidth = 0.2) } +
  geom_sf(data = trs_ea, shape = 21, fill = "grey50", color = "grey50", size = 0.05) +
  coord_sf(crs = st_crs(grid_ea)) +  theme_bw() +  facet_wrap (~year, ncol= 2) +   labs(title = "pilot whales per km", x = NULL, y = NULL)

### PILOT WHALE ###
pilot_molecular = pts_ea_join 
pilot_molecular$code_esp = ifelse(pilot_molecular$code_esp!="GLOMEL", "water", "GLOMEL")
pilot_molecular$Concentration = ifelse(pilot_molecular$code_esp=="water", 0, pilot_molecular$Concentration)

# Prepare data for the joint model
all_cells <- sort(unique(c(globi_visu$grid_id, pilot_molecular$grid_id)))
cell_id <- setNames(seq_along(all_cells), all_cells)
# Sightings
z <- log1p(globi_visu$ratio_individuals_km)
site_z <- cell_id[as.character(globi_visu$grid_id)]
M <- length(z)
# eDNA
y <- log1p(pilot_molecular$Concentration)
site_y <- cell_id[as.character(pilot_molecular$grid_id)]
N <- length(y)
pilot_data <- list(
  S = length(all_cells),
  N = N,
  y = y,
  site_y = site_y,
  M = M,
  z = z,
  site_z = site_z
)

pilot_visu_fit <- stan(file = "sightings_only.stan", data = pilot_data, chains  = 4, iter    = 6000,  warmup  = 3000,  seed    = 123,  control = list(adapt_delta = 0.9))
plot(pilot_visu_fit, par = c("X"))
pairs(pilot_visu_fit, pars = c("kappa", "omega", "tau"))

post <- rstan::extract(pilot_visu_fit)
ratio_mean <- apply(post$ratio_cell, 2, mean)
ratio_sd   <- apply(post$ratio_cell, 2, sd)
X_mean <- apply(post$X, 2, mean)
X_sd   <- apply(post$X, 2, sd)
X_z <- as.numeric(scale(X_mean))
cell_lookup_visu <- data.frame(
  grid_id        = all_cells,
  ratio_mean = ratio_mean,
  ratio_sd   = ratio_sd,
  X_mean         = X_mean,
  X_sd           = X_sd,
  X_z            = X_z)
pilot_cells_sf_visu <- grid_ea |> right_join(cell_lookup_visu, by = "grid_id")
ggplot(pilot_cells_sf_visu) +  geom_sf(aes(fill = ifelse(X_mean <= 0, NA_real_, X_mean)), color = NA) +
  scale_fill_gradient(
    low  = "#deebf7",
    high = "#08519c",
    limits = c(
      0,
      max(pilot_cells_sf_visu$X_mean, na.rm = TRUE)
    ),
    na.value = "#deebf7",
    name = "X posterior mean"
  ) +
  { if (exists("land_ea"))
    geom_sf(
      data = land_ea,
      fill = "grey85",
      color = "grey70",
      linewidth = 0.2
    ) } +  geom_sf(data = iso200_sf, color = "grey40",  linewidth = 0.3) +  theme_classic() +
  labs(title    = "Estimated abundance of pilots (visu model)")

###########################
## JOINT MODEL pilot ####
###########################

pilot_joint_fit <- stan(file = "joint_eDNA_sightings_new.stan", data = pilot_data, chains  = 4, iter    = 6000,  warmup  = 3000,  seed    = 123,  control = list(adapt_delta = 0.95, max_treedepth =12))
plot(pilot_joint_fit, par = c("X"))
post <- rstan::extract(pilot_joint_fit)
ratio_mean <- apply(post$ratio_cell, 2, mean) # escala observacional de la estima
ratio_sd   <- apply(post$ratio_cell, 2, sd) # escala observacional de la estima
X_mean <- apply(post$X, 2, mean)
X_sd   <- apply(post$X, 2, sd)
X_z <- as.numeric(scale(X_mean))
cell_lookup_joint <- data.frame(
  grid_id        = all_cells,
  ratio_mean = ratio_mean,
  ratio_sd   = ratio_sd,
  X_mean         = X_mean,
  X_sd           = X_sd,
  X_z            = X_z)

pilot_cells_sf_joint <- grid_ea |> right_join(cell_lookup_joint, by = "grid_id")
ggplot(pilot_cells_sf_joint) +  geom_sf(aes(fill = ifelse(X_mean <= 0, NA_real_, X_mean)),
                                          color = NA) +
  scale_fill_gradient(
    low  = "#deebf7",
    high = "#08519c",
    limits = c(
      0,
      max(pilot_cells_sf_joint$X_mean, na.rm = TRUE)
    ),
    na.value = "#deebf7",
    name = "X posterior mean"
  ) +
  { if (exists("land_ea"))
    geom_sf(
      data = land_ea,
      fill = "grey85",
      color = "grey70",
      linewidth = 0.2
    )
  } +
  geom_sf(data = iso200_sf, color = "grey40",  linewidth = 0.3) +  theme_classic() +  labs(title    = "Estimated abundance of pilots (joint model)")

###########################
### eDNA MODEL pilot ####
###########################
pilot_dna_fit <- stan(file = "eDNA_only_2.stan", data = pilot_data, chains  = 4, iter    = 6000,  warmup  = 3000,  seed    = 123,  control = list(adapt_delta = 0.9))

#calculate beta to set the priors for the joint model
post_dna <- rstan::extract(pilot_dna_fit)
beta_dna <- post_dna$beta
beta_hat <- mean(beta_dna)
beta_sd  <- sd(beta_dna)

plot(pilot_dna_fit, par = c("X"))
pairs(pilot_dna_fit, pars = c("alpha", "beta","sigma", "X[1]"))
post <- rstan::extract(pilot_dna_fit)
X_mean <- apply(post$X, 2, mean)
X_sd   <- apply(post$X, 2, sd)
X_z <- as.numeric(scale(X_mean))
cell_lookup_dna <- data.frame(
  grid_id        = all_cells,
  X_mean         = X_mean,
  X_sd           = X_sd,
  X_z            = X_z)

pilot_cells_sf_dna <- grid_ea |> right_join(cell_lookup_dna, by = "grid_id")
cells_with_edna <-  unique(all_cells[pilot_data$site_y])
pilot_cells_sf_dna <- pilot_cells_sf_dna |> filter(grid_id %in% cells_with_edna)
ggplot(pilot_cells_sf_dna) +
  geom_sf(
    aes(fill = ifelse(X_mean <= 0, NA_real_, X_mean)),
    color = NA
  ) +
  scale_fill_gradient(
    low  = "#deebf7",
    high = "#08519c",
    limits = c(
      0,
      max(pilot_cells_sf_dna$X_mean, na.rm = TRUE)
    ),
    na.value = "#deebf7",
    name = "X posterior mean"
  ) +
  { if (exists("land_ea"))
    geom_sf(
      data = land_ea,
      fill = "grey85",
      color = "grey70",
      linewidth = 0.2
    )
  } +
  geom_sf(
    data = iso200_sf,
    color = "grey40",
    linewidth = 0.3
  ) +
  theme_classic() + labs( title = "Estimated abundance of pilots (eDNA model)")

#######################################
df_joint <- pilot_cells_sf_joint |>  sf::st_drop_geometry() |>  select(grid_id, X_mean, X_sd) |>  mutate(Model = "Joint")
df_dna <- pilot_cells_sf_dna |>  sf::st_drop_geometry() |>  select(grid_id, X_mean, X_sd) |>  mutate(Model = "DNA")
df_visu <- pilot_cells_sf_visu |>  sf::st_drop_geometry() |>  select(grid_id, X_mean, X_sd) |>  mutate(Model = "Visual")
df_X <- bind_rows(df_joint, df_dna, df_visu)
cells_3models <- Reduce(intersect,  list(df_joint$grid_id, df_dna$grid_id,  df_visu$grid_id))
df_X_3models <- df_X |>  filter(grid_id %in% cells_3models)
ordered_ids <- df_X_3models |>filter(Model == "Joint") |>  arrange(X_mean) |>  pull(grid_id)
df_X_3models$grid_id <- factor( df_X_3models$grid_id,  levels = ordered_ids)
ggplot(df_X_3models, aes(x = grid_id, y = X_mean, color = Model)) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +  geom_errorbar(
    aes(ymin = X_mean - X_sd,  ymax = X_mean + X_sd ),
    position = position_dodge(width = 0.6),    width = 0.3, alpha = 0.7) + coord_flip() +
  scale_color_manual(values = c("Joint"  = "darkorange",  "DNA"    = "steelblue3", "Visual" = "olivedrab3")) +  theme_classic() +
  labs( x = "Grid cell", y = "Latent abundance (X)", color = "Model", subtitle = "Pilot whales")
