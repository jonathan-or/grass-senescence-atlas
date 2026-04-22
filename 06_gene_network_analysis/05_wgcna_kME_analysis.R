############################################################
## kME Structure Preservation Analysis
## - Computes intramodular connectivity (kME) across species
## - Correlates Annual (Maize/Sorghum) vs Perennial hub-ness
## - Generates preservation grids and integrated hub heatmaps
############################################################

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(ggpubr)
})

# --- Reproducibility ---
set.seed(42)

# =========================
# 0) Paths & Functions
# =========================
data_dir <- file.path(getwd(), "06_gene_network_analysis", "data")
out_dir  <- file.path(getwd(), "06_gene_network_analysis", "results")

# =========================
# 1) Load kME Master Table
# =========================
kme_master_path <- file.path(data_dir, "Supplemental_Data_kME_Master_Table_v2.csv")
if (file.exists(kme_master_path)) {
  cat("Loading pre-computed kME master table...\n")
  kME_final <- read_csv(kme_master_path, show_col_types = FALSE)
} else {
  stop("kME master table not found. Please ensure Supplemental_Data_kME_Master_Table_v2.csv is in the data directory.")
}

# =========================
# 2) Preservation Grid Plots
# =========================
cat("Generating kME preservation grids...\n")

plot_kme_grid <- function(df, plot_title) {
  ggplot(df, aes(x = kME_per, y = kME_ann, fill = subnetwork)) +
    geom_point(size = 1, alpha = 0.6, color = "grey14", shape = 21) +
    geom_smooth(method = "lm", color = "red", linetype = "dashed", se = FALSE, linewidth = 0.5) +
    stat_cor(method = "pearson", label.x = -0.4, label.y = 1.0, size = 3) +
    facet_grid(Annual ~ Perennial) +
    scale_fill_identity() +
    labs(x = "Perennial connectivity (kME)", y = "Annual connectivity (kME)", title = plot_title) +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          strip.background = element_rect(fill = "grey95", color = NA),
          strip.text = element_text(face = "italic"))
}

# Format for grid: Annuals (Zm, Sb) x Perennials (Ag, Sn, Td)
df_grid <- kME_final %>%
  mutate(species_label = case_when(
    Species == "Zea_mays" ~ "Zea mays",
    Species == "Sorghum_bicolor" ~ "S. bicolor",
    Species == "Andropogon_gerardii" ~ "A. gerardi",
    Species == "Sorghastrum_nutans" ~ "S. nutans",
    Species == "Tripsacum_dactyloides" ~ "T. dactyloides",
    TRUE ~ Species
  ))

# Separate Annuals and Perennials
perennials <- c("A. gerardi", "S. nutans", "T. dactyloides")
annuals    <- c("Zea mays", "S. bicolor")

df_grid_paired <- df_grid %>%
  filter(species_label %in% c(perennials, annuals)) %>%
  pivot_wider(id_cols = c(Orthogroup, Network, subnetwork), names_from = species_label, values_from = kME) %>%
  pivot_longer(cols = any_of(perennials), names_to = "Perennial", values_to = "kME_per") %>%
  pivot_longer(cols = any_of(annuals), names_to = "Annual", values_to = "kME_ann") %>%
  filter(!is.na(kME_per), !is.na(kME_ann))

# Leaf Blue Grid
p_leaf_blue <- plot_kme_grid(df_grid_paired %>% filter(Network == "leaf", subnetwork == "blue"), "Leaf m2_blue Preservation")
ggsave(file.path(out_dir, "15a_grid_leaf_blue_kME.png"), p_leaf_blue, width = 5, height = 3.5, dpi = 300)

# Root Blue Grid
p_root_blue <- plot_kme_grid(df_grid_paired %>% filter(Network == "underground", subnetwork == "blue"), "Root um2_blue Preservation")
ggsave(file.path(out_dir, "15b_grid_root_blue_kME.png"), p_root_blue, width = 5, height = 3.5, dpi = 300)

# =========================
# 3) Integrated Hub Heatmap
# =========================
cat("Generating integrated hub heatmap...\n")

og_hubs <- read_tsv(file.path(data_dir, "OG_symbols_hubs.txt"), col_types = "cccc") %>%
  rename(OG_symbol_hub = OG_symbol) # Avoid conflict if any

bwr_palette <- rev(c("#ca1875", "#ff1744", "#f06848", "#fceba5", "#e4f2f8", "#0dd2e8"))

plot_df_hubs <- kME_final %>%
  inner_join(og_hubs, by = c("Orthogroup" = "OG", "Network" = "network")) %>%
  mutate(species_label = factor(case_when(
    Species == "Andropogon_gerardii" ~ "A. gerardi",
    Species == "Sorghastrum_nutans" ~ "S. nutans",
    Species == "Tripsacum_dactyloides" ~ "T. dactyloides",
    Species == "Sorghum_bicolor" ~ "S. bicolor",
    Species == "Zea_mays" ~ "Zea mays",
    TRUE ~ Species
  ), levels = c("A. gerardi", "S. nutans", "T. dactyloides", "S. bicolor", "Zea mays"))) %>%
  mutate(subnetwork_f = factor(subnetwork.y, levels = c("m2_blue", "um1_turquoise", "um2_blue")))

p_hubs <- ggplot(plot_df_hubs, aes(x = OG_symbol_hub, y = species_label, fill = kME)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(colors = bwr_palette, limits = c(0, 1), oob = scales::squish, name = "kME") +
  facet_grid(. ~ subnetwork_f, scales = "free_x", space = "free_x") +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.y = element_text(face = "italic"), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, face = "italic"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "15c_heatmap_kME_hubs_integrated.png"), p_hubs, width = 5, height = 3, dpi = 300)

cat("\n==== kME preservation analysis complete ====\n")
