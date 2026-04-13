# 05_N_quant/02_amino_acid_analysis.R
# Evaluates Free Amino Acids vs Photosynthetic Index for Tripsacum dactyloides

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) > 0) args[1] else file.path(getwd(), "05_N_quant", "data")
plot_dir <- if (length(args) > 1) args[2] else file.path(getwd(), "05_N_quant", "plots")
meta_dir <- file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

source(file.path(getwd(), "03_cross-species_transcriptome_analysis", "Functions.R"))

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(ggpubr)
library(tidyr)
library(DESeq2)

# --- 1. SETUP DATA & MAPPING ---
# 8/4/23 is Week 2, 10/27/23 is Week 14
date_map <- data.frame(
  date = c("8/4/23", "8/11/23", "8/18/23", "8/25/23", "9/1/23", "9/8/23", 
           "9/15/23", "9/22/23", "9/29/23", "10/6/23", "10/13/23", "10/20/23", "10/27/23"),
  week = 2:14
)

aa_raw <- read_tsv(file.path(in_dir, "aa_data_ng_mg.txt"), show_col_types = FALSE)
leaf_metadata_ps_index <- readRDS(file.path(meta_dir, "leaf_metadata_ps_index.rds"))

td_aa_data <- aa_raw %>%
  filter(species == "Tripsacum_dactyloides", tissue == "leaves", total_aa > 0) %>%
  left_join(date_map, by = "date") %>%
  filter(!is.na(week))

supplementary_aa_tripsacum <- td_aa_data %>%
  select(species, accession, date, week, total_aa, everything()) %>%
  arrange(week, accession)

write_csv(supplementary_aa_tripsacum, file.path(in_dir, "Supplementary_Dataset_Tripsacum_AA.csv"))

td_ps_data <- leaf_metadata_ps_index %>%
  filter(species == "Tripsacum_dactyloides") %>%
  mutate(week_num = as.numeric(as.character(week)))

# --- 2. PLOT: TIMELINE (PHOTOSYNTHESIS VS. NITROGEN BURST) ---
total_aa_scale_ug <- 5

timeline_plot <- ggplot() +
  geom_jitter(data = td_ps_data, aes(x = week_num, y = ps_index), 
              color = "green4", alpha = 0.3, width = 0.25, size = 0.8, stroke = 0.3) +
  geom_jitter(data = td_aa_data, aes(x = week, y = (total_aa / 1000) / total_aa_scale_ug), 
              color = "deeppink3", alpha = 0.3, width = 0.25, size = 0.8, stroke = 0.3) +
  stat_summary(data = td_ps_data, aes(x = week_num, y = ps_index, group = 1), 
               fun = mean, geom = "line", color = "green4", linewidth = 1) +
  stat_summary(data = td_aa_data, aes(x = week, y = (total_aa / 1000) / total_aa_scale_ug, group = 1), 
               fun = mean, geom = "line", color = "deeppink3", linewidth = 1) +
  stat_summary(data = td_ps_data, aes(x = week_num, y = ps_index), 
               fun = mean, geom = "point", shape = 21, color = "black", fill = "green3", size = 2.5, stroke = 0.8) +
  stat_summary(data = td_aa_data, aes(x = week, y = (total_aa / 1000) / total_aa_scale_ug), 
               fun = mean, geom = "point", shape = 21, color = "black", fill = "deeppink3", size = 2.5, stroke = 0.8) +
  scale_y_continuous(
    name = "photosynthetic index", limits = c(0, 1.1), breaks = seq(0, 1, 0.2),
    sec.axis = sec_axis(~ . * total_aa_scale_ug, name = expression(paste("free amino acids (", mu, "g/mg DW)")), breaks = seq(0, 5, 0.5))
  ) +
  scale_x_continuous(limits = c(1, 14), breaks = seq(1, 14, 1)) +
  theme_classic(base_size = 10) +
  labs(x = "week")

tiff(file.path(plot_dir, "total_aas_Td_ug_v2_filtered.tiff"), width = 3.5, height = 2.5, units = "in", res = 600)
print(timeline_plot)
dev.off()

# --- 3. AMINO ACID VS. PS_INDEX SCATTERPLOT ---
td_plot_data <- td_aa_data %>%
  inner_join(td_ps_data %>% group_by(accession, week) %>% summarize(weekly_ps_index = mean(ps_index, na.rm = TRUE), .groups = "drop"), by = c("accession", "week"))

aa_ps_corr <- ggplot(td_plot_data, aes(x = weekly_ps_index, y = total_aa / 1000)) +
  geom_point(shape = 21, size = 2, color = "black", fill = "deeppink3", alpha = 0.7, stroke = 0.5) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE, linewidth = 0.8) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 1)) +
  scale_x_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +
  theme_classic(base_size = 10) +
  theme(axis.title = element_text(face = "bold")) +
  labs(x = "photosynthetic index", y = expression(paste("free amino acids (", mu, "g/mg DW)")))

tiff(file.path(plot_dir, "aa_ps_index_correlation_filtered.tiff"), width = 2.5, height = 2.5, units = "in", res = 600)
print(aa_ps_corr)
dev.off()

# --- 4. INDIVIDUAL AMINO ACIDS ---
aa_list <- c("Alanine", "Arginine", "Asparagine", "Aspartate", "Glutamine", 
             "Glutamate", "Glycine", "Histidine", "Isoleucine", "Leucine", 
             "Lysine", "Methionine", "Phenylalanine", "Proline", "Serine", 
             "Threonine", "Tryptophan", "Tyrosine", "Valine", "Cysteine")

aa_long_data <- td_plot_data %>%
  pivot_longer(cols = all_of(aa_list), names_to = "amino_acid", values_to = "concentration_ng_mg") %>%
  mutate(log2_conc = log2((concentration_ng_mg / 1000) + 1), detected = concentration_ng_mg > 0)

aa_faceted_corr <- ggplot(aa_long_data, aes(x = weekly_ps_index, y = log2_conc)) +
  geom_point(shape = 21, size = 1.2, color = "black", fill = "deeppink3", alpha = 0.6, stroke = 0.3) +
  geom_smooth(data = filter(aa_long_data, detected == TRUE),
              method = "lm", color = "black", linewidth = 0.6, se = FALSE, linetype = "dashed") +
  stat_cor(data = filter(aa_long_data, detected == TRUE),
           method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 2.5) +
  facet_wrap(~amino_acid, ncol = 5, scales = "fixed") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
        strip.text = element_text(face = "bold", size = 8),
        legend.position = "none",
        panel.spacing = unit(0.5, "lines")) +
  labs(x = "photosynthetic index", y = expression(bold(paste(log[2], "(", mu, "g / mg DW + 1)"))))

tiff(file.path(plot_dir, "aa_correlations_final_clean.tiff"), width = 7, height = 5, units = "in", res = 600)
print(aa_faceted_corr)
dev.off()

# --- 5. EXPRESSION (SAG12 / SAG15) VS AMINO ACIDS ---
cat("Preparing marker gene VST models...\n")

leaf_og_counts <- readRDS(file.path(meta_dir, "og_counts.rds"))
metadata <- readRDS(file.path(meta_dir, "metadata.rds"))
leaf_metadata <- metadata %>% filter(tissue == "leaf")
leaf_og_counts <- leaf_og_counts[, colnames(leaf_og_counts) %in% leaf_metadata$sample_id]

dds <- DESeqDataSetFromMatrix(countData = leaf_og_counts, colData = leaf_metadata, design = ~ 1)
dds <- estimateSizeFactors(dds)
vst_mat <- assay(vst(dds, blind = TRUE))

target_genes <- c("SAG12" = "OG0008818", "SAG15" = "OG0012961")

for (marker_name in names(target_genes)) {
  marker_og <- target_genes[[marker_name]]
  expr_vector <- as.numeric(vst_mat[marker_og, ])
  
  scaled_marker <- data.frame(marker_zscore = as.vector(scale(expr_vector)), sample_id = colnames(vst_mat))
  
  td_marker_data <- scaled_marker %>%
    inner_join(leaf_metadata %>% select(sample_id, species, week), by = "sample_id") %>%
    filter(species == "Tripsacum_dactyloides") %>%
    mutate(week_num = as.numeric(as.character(week)))
  
  aa_shift <- 1
  aa_scale_factor <- 1.0 
  
  aligned_plot <- ggplot() +
    geom_jitter(data = td_marker_data, aes(x = week_num, y = marker_zscore), 
                color = "purple4", alpha = 0.2, width = 0.2, size = 0.6) +
    geom_jitter(data = td_aa_data, aes(x = week, y = (total_aa / 1000) - aa_shift), 
                color = "deeppink3", alpha = 0.2, width = 0.2, size = 0.6) +
    stat_summary(data = td_marker_data, aes(x = week_num, y = marker_zscore, group = 1), 
                 fun = mean, geom = "line", color = "purple4", linewidth = 0.8) +
    stat_summary(data = td_aa_data, aes(x = week, y = (total_aa / 1000) - aa_shift, group = 1), 
                 fun = mean, geom = "line", color = "deeppink3", linewidth = 0.8) +
    stat_summary(data = td_marker_data, aes(x = week_num, y = marker_zscore), 
                 fun = mean, geom = "point", shape = 21, color = "black", fill = "purple", size = 2, stroke = 0.5) +
    stat_summary(data = td_aa_data, aes(x = week, y = (total_aa / 1000) - aa_shift), 
                 fun = mean, geom = "point", shape = 21, color = "black", fill = "deeppink3", size = 2, stroke = 0.5) +
    scale_y_continuous(
      name = paste(marker_name, "expression (z-score)"),
      limits = c(-2.1, 3.1),
      breaks = seq(-2, 3, 1),
      sec.axis = sec_axis(~ (. + aa_shift) * aa_scale_factor, name = expression(paste("free amino acids (", mu, "g/mg DW)")), breaks = seq(0, 5, 1))
    ) +
    scale_x_continuous(limits = c(1, 14), breaks = seq(1, 14, 1)) +
    theme_classic(base_size = 10) +
    theme(axis.title.y.left = element_text(color = "purple4", face = "bold"),
          axis.title.y.right = element_text(color = "deeppink3", face = "bold")) +
    labs(x = "week")
  
  tiff(file.path(plot_dir, paste0("Td_", marker_name, "_vs_AA_aligned.tiff")), width = 3.5, height = 2.5, units = "in", res = 600)
  print(aligned_plot)
  dev.off()
}

cat("Amino acid analysis module completed.\n")
