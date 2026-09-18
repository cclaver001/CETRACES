################################################
###################### BLAST ###################
################################################

library(dplyr)
library(tidyr)  
library(rlang)
library(ggplot2)
library(forcats)
library(scales)
library(tidyverse)
library(scatterpie)
library(sf)
library(rnaturalearth)

# Load data
# Mock samples (true proportions)
mock <- read_excel("New Metabarcoding Primers - CetAZTI/BLAST y SWARM/BLAST_results.xlsx", sheet = "mock_true")
# Unique sequences
sequs <- read_excel("New Metabarcoding Primers - CetAZTI/BLAST y SWARM/BLAST_results.xlsx", sheet = "count_table")
# Tax assignment
blast <- read_excel("New Metabarcoding Primers - CetAZTI/BLAST y SWARM/BLAST_results.xlsx", sheet = "BLAST_1_hit")

# Merge tax assignations with unique sequences
 All_seqs = merge(blast[,1:5], sequs, by= "Representative_Sequence", all = T)
 
 All_seqs$Match_rfDB[is.na(All_seqs$Match_rfDB)] <- "unclassified"
 
 All_seqs = All_seqs |> select(!c("Representative_Sequence", "perc_iden", "length", "mismatch"))
 
 # 4) Agrupar por el taxón (Match_rfDB / match_rf) y sumar las columnas numéricas
 All_seqs <- All_seqs %>% group_by(Match_rfDB) %>%
   summarise(
     across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
     .groups = "drop"
   ) %>%
   arrange(Match_rfDB)
 
All_seqs$Match_rfDB <- sub("^[^.]*\\..{2}", "", All_seqs$Match_rfDB)

# Prepare color pallette
pal <- c(
  # Balaenopteridae
  Balaenoptera_musculus        = "#1B9E77",
  Balaenoptera_physalus        = "#33A02C",
  Balaenoptera_borealis        = "#66C266",
  Balaenoptera_acutorostrata   = "#A6D854",
  # Physeteridae + Kogiidae
  Physeter_catodon             = "#6A3D9A",
  Kogia_breviceps              = "#9E63C9",
  # Ziphiidae
  Hyperoodon_ampullatus        = "#1F78B4",
  Ziphius_cavirostris          = "#3182BD",
  Mesoplodon_bidens            = "#6BAED6",
  Mesoplodon_densirostris      = "#9ECAE1",
  Mesoplodon_europaeus         = "#C6DBEF",
  Mesoplodon_minus             = "#DEEBF7",
  # Delphinidae 
  Orcinus_orca                 = "darkred",
  Pseudorca_crassidens         = "red",
  Globicephala_melas           = "#E6550D",
  Globicephala_macrorhynchus   = "#FC8D59",
  Grampus_griseus              = "#FEE0D2",
  Tursiops_truncatus           = "#FB6A4A",
  Stenella_coeruleoalba        = "#FCAE91",
  Delphinus_delphis            = "#F16975",
  # Phocoenidae – amarillo
  Phocoena_phocoena            = "#FFD92F",
  # Otros
  Homo_sapiens                 = "#F0F0F0",
  unclassified                 = "#BDBDBD"
)

# plotear las mock (pre sequencing)
mock_long <- mock %>% pivot_longer(cols = where(is.numeric), names_to = "muestra", values_to = "proporcion")
names(mock_long)[1] <- "Species"
ggplot(mock_long, aes(x = muestra, y = proporcion, fill = Species)) +
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  #scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = pal, breaks = names(pal) ) +
  labs(
    x = "Mock sample",
    y = "Proportion",
    fill = "Species",
    title = "Mock samples(true proportions)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# plotear las mock (post sequencing)
mock_post = All_seqs |> select(!total)
mock_post = mock_post[,1:11]
mock_post_long <- mock_post %>% pivot_longer(cols = where(is.numeric), names_to = "muestra", values_to = "proporcion")
names(mock_post_long)[1] <- "Species"

ggplot(mock_post_long, aes(x = muestra, y = proporcion, fill = Species)) +
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  #scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = pal, breaks = names(pal) ) +
  labs(
    x = "Mock sample",
    y = "Proportion",
    fill = "Species",
    title = "Mock samples (post sequencing)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# plot eDNA samples 
Count_table_JUVENA = All_seqs |> select(!total) |>  select(!starts_with("25CETMock"))
colSums(Count_table_JUVENA[ , -1, drop = FALSE], na.rm = TRUE)
Count_table_JUVENA = Count_table_JUVENA |> filter(Match_rfDB != "unclassified")
colSums(Count_table_JUVENA[ , -1, drop = FALSE], na.rm = TRUE)
Count_table_JUVENA = Count_table_JUVENA |> filter(Match_rfDB != "Homo_sapiens")
colSums(Count_table_JUVENA[ , -1, drop = FALSE], na.rm = TRUE)

# only phylotypes with more than 10 reads per sample will be considered
Count_table_JUVENA[, c(TRUE, colSums(Count_table_JUVENA[, -1, drop = FALSE], na.rm = TRUE) > 10)]

metabarcoding_eDNA <- Count_table_JUVENA %>%
  select(1, where(~ sum(as.numeric(.x), na.rm = TRUE) > 10)) %>%
  pivot_longer(-1, names_to = "Sample", values_to = "Reads") %>%
  rename(Species = 1) %>%
  filter(Reads > 0) %>%
  mutate(
    Estacion = gsub("25SMixJUV", "S0", Sample)
  ) %>%
  merge(eDNA_effort, by = "Estacion") %>%
  group_by(Sample) %>%
  mutate(RelAbund = Reads / sum(Reads)) %>%
  ungroup()

pie_df <- metabarcoding_eDNA %>%  select(Latitud, Longitud, Species, RelAbund) %>%  pivot_wider(names_from = Species, values_from = RelAbund, values_fill = 0)
cols <- setNames(c("#1f77b4", "#ff7f0e", "#2ca02c"),  names(pie_df)[-(1:2)])

# Set coastline (fast version)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
lon_rng <- range(pie_df$Longitud, na.rm = TRUE) + c(-3, 1)
lat_rng <- range(pie_df$Latitud,  na.rm = TRUE) + c(-0.6, 1.2)
ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70", linewidth = 0.2) +
  coord_sf(xlim = lon_rng, ylim = lat_rng, expand = FALSE) +
  geom_scatterpie(
    data = pie_df,
    aes(x = Longitud, y = Latitud, r = 0.3), 
    cols = names(cols),
    color = "black"
  ) +
  scale_fill_manual(values = cols, name = "Species") +
  theme_bw() +   labs(x = "Longitud", y = "Latitud")

