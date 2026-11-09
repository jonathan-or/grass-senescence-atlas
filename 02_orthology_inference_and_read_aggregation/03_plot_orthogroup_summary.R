# 03_plot_orthogroup_summary.R
# Supplementary Figure 3: definition of the core set of orthogroups.
#   A  genes per species: in core orthogroups, in any orthogroup, and annotated in total (bars;
#      written as published and with the analysis core set, see below)
#   B  number of species represented per orthogroup (core set: all 14 species)
#   C  number of genes per orthogroup, all orthogroups (dashed line: 75-gene limit); inset: core set
#   D  share of each library's reads assigned to core orthogroups (2,685 libraries)
# The species tree and schematic in panel A were drawn separately.
#
# Core set (as in 02_aggregate_counts_to_orthogroups.R): orthogroups with at least one gene in
# each of the 14 species and at most 75 genes in total (16,778 orthogroups).
#
# Usage (from the repository root):
#   Rscript 02_orthology_inference_and_read_aggregation/03_plot_orthogroup_summary.R [counts_rds] [out_dir]
#   counts_rds: per-species gene-level counts (list of genes x libraries data frames), as written by
#               01_import_counts_and_metadata.R (default: count_data_frames.rds in this directory).
#               Only needed for panel D, which is skipped if the file is missing.

args <- commandArgs(trailingOnly = TRUE)
mod_dir    <- file.path(getwd(), "02_orthology_inference_and_read_aggregation")
counts_rds <- if (length(args) > 0) args[1] else file.path(mod_dir, "count_data_frames.rds")
out_dir    <- if (length(args) > 1) args[2] else file.path(mod_dir, "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})
source(file.path(getwd(), "03_cross-species_transcriptome_analysis", "Functions.R"))  # species_palette

# ==========================================
# Orthogroup membership (genes per species per orthogroup)
# ==========================================
OGs <- read.table(gzfile(file.path(mod_dir, "Orthogroups.tsv.gz")), sep = "\t", header = TRUE,
                  quote = "", stringsAsFactors = FALSE) %>%
  dplyr::rename(Anatherum_virginicum = Andropogon_virginicus)
gene_num <- sapply(OGs[, -1], function(x) ifelse(x == "", 0L, lengths(strsplit(x, ", "))))
rownames(gene_num) <- OGs$Orthogroup

n_species <- rowSums(gene_num > 0)
n_genes   <- rowSums(gene_num)
is_core   <- n_species == ncol(gene_num) & n_genes <= 75
cat("Orthogroups:", nrow(gene_num), "| in all", ncol(gene_num), "species:", sum(n_species == ncol(gene_num)),
    "| core set:", sum(is_core), "\n")

# ==========================================
# A) Genes per species: core / assigned to an orthogroup / annotated
# ==========================================
# Two versions of the core segment are drawn:
#   OG_plot.tiff       as published: the core segment counts genes in the first 16,778 of the 17,329
#                      orthogroups present in all 14 species, in OrthoFinder order (largest first).
#                      This set includes 551 orthogroups with more than 75 genes and leaves out 551
#                      core orthogroups with 18-26 genes.
#   OG_plot_core.tiff  the core set used in the analysis (all 14 species, <= 75 genes).
total_genes <- read.delim(file.path(mod_dir, "total_gene_counts_named.tsv"))
species_order <- c("Panicum_virgatum", "Elionorus_tripsacoides", "Tripsacum_dactyloides", "Zea_mays",
                   "Zea_diploperennis", "Chrysopogon_zizanioides", "Sorghastrum_nutans",
                   "Miscanthus_giganteus", "Sorghum_bicolor", "Themeda_triandra",
                   "Bothriochloa_laguroides", "Andropogon_gerardii", "Anatherum_virginicum",
                   "Schizachyrium_scoparium")  # top to bottom, as in the species tree

in_all_species <- which(n_species == ncol(gene_num))
as_published   <- seq_len(nrow(gene_num)) %in% in_all_species[seq_len(sum(is_core))]

gene_sets <- tibble(Species = colnames(gene_num),
                    core = colSums(gene_num[is_core, ]),
                    core_as_published = colSums(gene_num[as_published, ]),
                    assigned = colSums(gene_num)) %>%
  inner_join(total_genes, by = "Species") %>%
  mutate(Species = factor(Species, levels = rev(species_order)))
write.table(gene_sets, file.path(out_dir, "orthogroup_gene_sets_per_species.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(gene_sets))

plot_gene_sets <- function(core_col) {
  bars <- gene_sets %>%
    transmute(Species, Filtered = .data[[core_col]], Assigned = assigned, Total) %>%
    pivot_longer(-Species, names_to = "Category", values_to = "Count")
  labels <- gene_sets %>% transmute(Species, x_pos = .data[[core_col]] * 0.05, Label = .data[[core_col]])
  ggplot(bars, aes(y = Species, x = Count, fill = Species, alpha = Category)) +
    geom_bar(stat = "identity", position = "identity") +
    geom_text(data = labels, aes(x = x_pos, y = Species, label = Label), hjust = 0, size = 2,
              color = "white", inherit.aes = FALSE) +
    scale_fill_manual(values = species_palette) +
    scale_alpha_manual(values = c(Total = 0.2, Assigned = 0.5, Filtered = 1)) +
    scale_x_continuous(breaks = seq(0, 225000, 25000), labels = label_number(scale_cut = cut_short_scale())) +
    theme_classic(base_size = 8) +
    labs(x = "gene count", y = NULL) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), legend.position = "none")
}
for (v in list(c("core_as_published", "OG_plot.tiff"), c("core", "OG_plot_core.tiff"))) {
  tiff(file.path(out_dir, v[2]), width = 2.5, height = 2.5, units = "in", res = 600)
  print(plot_gene_sets(v[1]))
  dev.off()
}

# ==========================================
# B) Number of species per orthogroup
# ==========================================
species_hist <- tibble(n_species) %>%
  dplyr::count(n_species, name = "n") %>%
  mutate(all_species = n_species == ncol(gene_num))
print(as.data.frame(species_hist))

p_species <- ggplot(species_hist, aes(factor(n_species), n, fill = all_species)) +
  geom_col(width = 1, color = "black", linewidth = 0.6) +
  scale_fill_manual(values = c(`FALSE` = "#A4C3D5", `TRUE` = "#0E6190"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "number of species in orthogroup", y = "number of orthogroups") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "OG_hist.tiff"), p_species, width = 5, height = 3.4, units = "in", dpi = 300,
       compression = "lzw")

# ==========================================
# C) Genes per orthogroup: all orthogroups (75-gene limit) and the core set (inset)
# ==========================================
p_size_all <- ggplot(tibble(n_genes), aes(n_genes)) +
  # Bin width chosen to match the published panel (the original setting was not recorded)
  geom_histogram(binwidth = 19.5, fill = "#9ECFE0", color = "black", linewidth = 0.2) +
  geom_vline(xintercept = 75, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = seq(0, 1000, 250)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "number of genes in orthogroup", y = "frequency") +
  theme_classic(base_size = 16)
ggsave(file.path(out_dir, "OG_genehist.tiff"), p_size_all, width = 5, height = 3.4, units = "in",
       dpi = 300, compression = "lzw")

p_size_core <- ggplot(tibble(n_genes = n_genes[is_core]), aes(n_genes)) +
  # bins = 56 reproduces the published inset (bins slightly wider than one gene)
  geom_histogram(bins = 56, fill = "#0E6190", color = "black", linewidth = 0.2) +
  scale_x_continuous(breaks = seq(15, 75, 15)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 18)
ggsave(file.path(out_dir, "filtOGhist.tiff"), p_size_core, width = 3, height = 2, units = "in", dpi = 300,
       compression = "lzw")

# ==========================================
# D) Reads assigned to core orthogroups per library
# ==========================================
if (file.exists(counts_rds)) {
  counts    <- readRDS(counts_rds)
  og_counts <- readRDS(file.path(getwd(), "03_cross-species_transcriptome_analysis", "data", "og_counts.rds"))
  assigned_reads <- colSums(og_counts)

  per_library <- bind_rows(lapply(names(counts), function(sp) {
    total <- colSums(counts[[sp]], na.rm = TRUE)
    libs  <- intersect(names(total), names(assigned_reads))
    tibble(Library = libs, Species = sp, Total = total[libs], Assigned = assigned_reads[libs])
  })) %>%
    mutate(Assigned_Pct = Assigned / Total, Unassigned_Pct = 1 - Assigned_Pct)
  stopifnot(setequal(per_library$Library, colnames(og_counts)))
  mean_assigned_pct <- mean(per_library$Assigned_Pct, na.rm = TRUE)
  cat(sprintf("Panel D: %d libraries; mean share of reads in core orthogroups %.1f%%\n",
              nrow(per_library), 100 * mean_assigned_pct))

  plot_data <- per_library %>%
    arrange(desc(Assigned_Pct)) %>%
    mutate(Library = factor(Library, levels = unique(Library))) %>%
    pivot_longer(c(Assigned_Pct, Unassigned_Pct), names_to = "Category", values_to = "Percent") %>%
    mutate(Category = factor(recode(Category, Assigned_Pct = "Assigned", Unassigned_Pct = "Unassigned"),
                             levels = c("Unassigned", "Assigned")))

  p_assigned <- ggplot(plot_data, aes(Library, Percent, fill = Category)) +
    geom_col(width = 0.85, color = NA) +
    geom_hline(yintercept = mean_assigned_pct, linetype = "dashed", color = "red") +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = c(Unassigned = "grey80", Assigned = "steelblue")) +
    labs(x = paste0("libraries (n = ", nlevels(plot_data$Library), ")"), y = "% of reads assigned to core OGs") +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_blank(), legend.position = "none", panel.grid.major.x = element_blank()) +
    coord_cartesian(clip = "off")
  tiff(file.path(out_dir, "assigned_percentage.tiff"), width = 2, height = 2, units = "in", res = 600)
  print(p_assigned)
  dev.off()
} else {
  message("Panel D skipped: ", counts_rds, " not found (run 01_import_counts_and_metadata.R first).")
}

cat("Supplementary Figure 3 panels saved to", out_dir, "\n")
