############################################################
## Module Preservation Plotting
## - Generates Zsummary vs module size plots for annuals (Maize/Sorghum)
## - Faceted by species and tissue type
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

# =========================
# 0) Paths
# =========================
res_dir <- file.path(getwd(), "06_gene_network_analysis", "results")
out_dir <- file.path(getwd(), "06_gene_network_analysis", "results")

# Annual accessions of interest
annuals_map <- c(
  "Sb_BTx623"      = "Sorghum",
  "Zm_lh287Xlh244" = "Maize"
)

# =========================
# 1) Plotting Function
# =========================
plot_preservation_z <- function(pres_df, me_ps_df, tissue_name, label_prefix = "m") {
  
  # Wide table for annuals
  df_wide <- pres_df %>%
    filter(accession %in% names(annuals_map)) %>%
    mutate(species = annuals_map[accession]) %>%
    filter(module_color != "grey") %>%
    group_by(module_num, module_color, module_size, species) %>%
    summarise(Zsummary = mean(Zsummary.pres, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = species, values_from = Zsummary)
  
  # Join with ps-index correlation for ordering
  df_plot <- df_wide %>%
    left_join(me_ps_df %>% select(module_num, module_label, rho_ps), by = "module_num") %>%
    mutate(module_label = factor(module_label, levels = unique(module_label[order(rho_ps, decreasing = TRUE)]))) %>%
    pivot_longer(cols = any_of(c("Maize", "Sorghum")), names_to = "species", values_to = "Zsummary") %>%
    filter(!is.na(Zsummary))
  
  # Threshold lines: 2 (weak), 10 (strong)
  ggplot(df_plot, aes(x = module_label, y = Zsummary)) +
    geom_hline(yintercept = c(2, 10), linetype = "dashed", linewidth = 0.5, color = "coral") +
    geom_point(aes(size = module_size, color = module_color), alpha = 0.85) +
    facet_wrap(~ species, ncol = 2) +
    scale_color_identity() +
    scale_size_continuous(range = c(1, 6)) +
    scale_y_continuous(breaks = c(0, 2, 10, 20, 30, 40, 60)) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    ) +
    labs(x = paste(tissue_name, "modules"), y = expression(Z[summary]~"Preservation"), size = "module size")
}

# =========================
# 2) Leaf Plot
# =========================
cat("Plotting leaf preservation...\n")
pres_leaf <- read_tsv(file.path(res_dir, "pres_summary_accessions_leaf.txt"), show_col_types = FALSE)
me_ps_leaf <- read_csv(file.path(res_dir, "04_ME_vs_ps_index_leaf.csv"), show_col_types = FALSE)

p_leaf <- plot_preservation_z(pres_leaf, me_ps_leaf, "leaf", "m")
ggsave(file.path(out_dir, "14a_module_preservation_leaf.pdf"), p_leaf, width = 5, height = 4)

# =========================
# 3) Root Plot
# =========================
cat("Plotting root preservation...\n")
pres_root <- read_tsv(file.path(res_dir, "pres_summary_accessions_root.txt"), show_col_types = FALSE)
me_ps_und  <- read_csv(file.path(res_dir, "04_ME_vs_ps_index_und.csv"), show_col_types = FALSE)

p_root <- plot_preservation_z(pres_root, me_ps_und, "root", "um")
ggsave(file.path(out_dir, "14b_module_preservation_root.pdf"), p_root, width = 5, height = 4)

# =========================
# 4) Stem Plot
# =========================
cat("Plotting stem preservation...\n")
pres_stem <- read_tsv(file.path(res_dir, "pres_summary_accessions_stem.txt"), show_col_types = FALSE)

p_stem <- plot_preservation_z(pres_stem, me_ps_und, "stem", "um")
ggsave(file.path(out_dir, "14c_module_preservation_stem.pdf"), p_stem, width = 3.5, height = 4)

cat("\n==== Preservation plotting complete ====\n")
