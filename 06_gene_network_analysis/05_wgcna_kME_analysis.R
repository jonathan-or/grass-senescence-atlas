############################################################
## kME Structure Preservation Analysis
## - Computes intramodular connectivity (kME) across species
## - Correlates Annual (Maize/Sorghum) vs Perennial hub-ness
## - Generates preservation grids and integrated hub heatmaps
##
## kME master table (Supplementary Data SD13): per-accession kME for
## leaf m2_blue (leaf VST, as 01_wgcna_leaf.R) and underground um1_turquoise /
## um2_blue (ROOT-ONLY VST: root libraries of all accessions). Loaded if present,
## otherwise rebuilt here.
############################################################

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(ggpubr)
  library(DESeq2)
})

# --- Reproducibility ---
set.seed(42)

source(file.path(getwd(), "03_cross-species_transcriptome_analysis", "Functions.R"))

# Set to TRUE to rebuild the kME master table even if it exists
force_recompute <- FALSE

# =========================
# 0) Paths & Functions
# =========================
ps_dir    <- file.path(getwd(), "04_photosynthetic_index", "data")
annot_dir <- file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")
data_dir  <- file.path(getwd(), "06_gene_network_analysis", "data")
out_dir   <- file.path(getwd(), "06_gene_network_analysis", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# kME within modules for ONE accession: module eigengenes are computed from this
# accession's samples only; kME = Pearson cor(gene, own-module eigengene).
# Genes with zero variance within the accession (e.g. no reads in any of its
# libraries) are left out of the eigengene and get kME = NA; if kept, they shift
# the eigengene in a WGCNA-version-dependent way.
get_kME_within_modules_per_accession <- function(meta, full_vst, ref_colors_named, accession, species,
                                                 drop_controls = TRUE, min_genes_per_module = 5) {
  full_vst <- as.matrix(full_vst); mode(full_vst) <- "numeric"

  samp_ids <- meta %>%
    filter(species == !!species, accession == !!accession) %>%
    pull(sample_id) %>%
    intersect(colnames(full_vst))
  if (!length(samp_ids)) {
    warning(sprintf("No samples for accession '%s' (species '%s').", accession, species))
    return(tibble())
  }

  keep_genes <- intersect(names(ref_colors_named), rownames(full_vst))
  Expr <- t(full_vst[keep_genes, samp_ids, drop = FALSE])  # samples x genes
  colors_here <- ref_colors_named[keep_genes]

  if (drop_controls) {
    keep <- !(colors_here %in% c("grey", "gold"))
    Expr <- Expr[, keep, drop = FALSE]
    colors_here <- colors_here[keep]
  }

  var_ok <- setNames(matrixStats::colSds(Expr) > 0, colnames(Expr))
  ME_w <- moduleEigengenes(Expr[, var_ok, drop = FALSE], colors = colors_here[var_ok])$eigengenes
  ME_w <- orderMEs(ME_w)

  imap_dfr(split(names(colors_here), colors_here), function(genes_in_mod, mod_col) {
    me_name <- paste0("ME", mod_col)
    if (length(genes_in_mod) < min_genes_per_module || !me_name %in% colnames(ME_w)) return(NULL)
    kme_vals <- suppressWarnings(cor(Expr[, genes_in_mod, drop = FALSE], ME_w[, me_name], use = "pairwise"))
    kme_vals[!var_ok[genes_in_mod]] <- NA
    tibble(gene_id = genes_in_mod, module_color = mod_col, kME = as.numeric(kme_vals),
           accession = accession, species = species)
  })
}

# kME for the focal accessions in one network, restricted to the exported modules
kme_accessions <- c("Ag_Niagara", "Sn_WI", "Td_Pete", "Sb_BTx623", "Zm_lh287Xlh244")

get_network_kME <- function(meta, full_vst, net, DatExpr, network, modules) {
  ref_colors <- setNames(labels2colors(net$colors), colnames(DatExpr))
  map_dfr(kme_accessions, function(acc) {
    sp <- unique(meta$species[meta$accession == acc])
    get_kME_within_modules_per_accession(meta, full_vst, ref_colors, acc, sp)
  }) %>%
    filter(module_color %in% modules) %>%
    mutate(network = network)
}

# =========================
# 1) Load or build kME Master Table
# =========================
kme_master_path <- file.path(data_dir, "Supplemental_Data_kME_Master_Table_v2.csv")
if (file.exists(kme_master_path) && !force_recompute) {
  cat("Loading pre-computed kME master table...\n")
  kME_final <- read_csv(kme_master_path, show_col_types = FALSE)
} else {
  cat("Building kME master table...\n")

  # Leaf: same VST as 01_wgcna_leaf.R (all leaf libraries)
  leaf_metadata_ps_index <- readRDS(file.path(ps_dir, "leaf_metadata_ps_index.rds")) %>%
    mutate(sample_id = as.character(sample_id))
  vst_leaf <- plot_pca_ggplot2(
    count_data = readRDS(file.path(ps_dir, "leaf_og_counts_dps.rds")),
    metadata = leaf_metadata_ps_index,
    color_by = "species",
    species_palette = species_palette,
    accession_palette = accession_palette
  )$vst_matrix
  kME_leaf <- get_network_kME(leaf_metadata_ps_index, vst_leaf,
                              readRDS(file.path(data_dir, "net_leaf.rds")),
                              readRDS(file.path(data_dir, "DatExpr_leaf.rds")),
                              "leaf", "blue")

  # Underground: ROOT-ONLY VST (all root libraries). Note: not the root + stem
  # VST the underground network was built on (02_wgcna_underground.R)
  root_metadata_ps_index <- readRDS(file.path(ps_dir, "root_metadata_ps_index.rds")) %>%
    mutate(sample_id = as.character(sample_id))
  vst_root <- plot_pca_ggplot2(
    count_data = as.matrix(readRDS(file.path(ps_dir, "root_og_counts_dps.rds"))),
    metadata = root_metadata_ps_index,
    color_by = "species",
    species_palette = species_palette,
    accession_palette = accession_palette
  )$vst_matrix
  kME_root <- get_network_kME(root_metadata_ps_index, vst_root,
                              readRDS(file.path(data_dir, "net_und.rds")),
                              readRDS(file.path(data_dir, "DatExpr_und.rds")),
                              "underground", c("blue", "turquoise"))

  # Maize gene members per OG
  zm_id_map <- readRDS(file.path(annot_dir, "gene_og_mapping.rds")) %>%
    filter(Species == "Zea_mays") %>%
    group_by(Orthogroup) %>%
    summarise(Maize_GeneID = paste(unique(GeneID), collapse = "; "), .groups = "drop")

  kME_final <- bind_rows(kME_leaf, kME_root) %>%
    left_join(zm_id_map, by = c("gene_id" = "Orthogroup")) %>%
    select(Orthogroup = gene_id, Network = network, subnetwork = module_color,
           Species = species, Accession = accession, kME, Maize_GeneID) %>%
    arrange(Network, subnetwork, Species, desc(kME))

  write_csv(kME_final, kme_master_path)
  cat("  kME master table saved:", nrow(kME_final), "rows\n")
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

# Leaf Blue Grid (Figure 6B)
p_leaf_blue <- plot_kme_grid(df_grid_paired %>% filter(Network == "leaf", subnetwork == "blue"), "Leaf m2_blue Preservation")
ggsave(file.path(out_dir, "15a_grid_leaf_blue_kME.png"), p_leaf_blue, width = 5, height = 3.5, dpi = 300)

# Root Turquoise Grid (Figure 6E)
p_root_turq <- plot_kme_grid(df_grid_paired %>% filter(Network == "underground", subnetwork == "turquoise"), "Root um1_turquoise Preservation")
ggsave(file.path(out_dir, "15d_grid_root_turquoise_kME.png"), p_root_turq, width = 5, height = 3.5, dpi = 300)

# Root Blue Grid (Figure 6F)
p_root_blue <- plot_kme_grid(df_grid_paired %>% filter(Network == "underground", subnetwork == "blue"), "Root um2_blue Preservation")
ggsave(file.path(out_dir, "15b_grid_root_blue_kME.png"), p_root_blue, width = 5, height = 3.5, dpi = 300)

# =========================
# 3) Integrated Hub Heatmap
# =========================
cat("Generating integrated hub heatmap...\n")

og_hubs <- read_tsv(file.path(data_dir, "OG_symbols_hubs.txt"), col_types = "cccc") %>%
  dplyr::rename(OG_symbol_hub = OG_symbol) # Avoid conflict if any; S4Vectors (DESeq2) masks rename()

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
