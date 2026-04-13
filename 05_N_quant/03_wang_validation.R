# 05_N_quant/03_wang_validation.R
# Validates the photosynthetic index against the Wang et al. 2014 maize leaf gradient dataset

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) > 0) args[1] else file.path(getwd(), "05_N_quant", "data")
plot_dir <- if (length(args) > 1) args[2] else file.path(getwd(), "05_N_quant", "plots")

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

source(file.path(getwd(), "03_cross-species_transcriptome_analysis", "Functions.R"))

library(ggplot2)
library(dplyr)
library(DESeq2)
library(tidyr)

# --- 1. Load data ---
metadata <- read.csv(file.path(in_dir, "metadata.txt"), header = TRUE, sep = "\t")
count_data <- read.delim(file.path(in_dir, "Zea_mays_counts.txt"), row.names = 1, stringsAsFactors = FALSE)

# Clean transcript IDs
rownames(count_data) <- sub("_T0\\d{2}$", "", rownames(count_data))
count_data <- aggregate(count_data, by = list(Gene = rownames(count_data)), FUN = sum)
rownames(count_data) <- count_data$Gene
count_data$Gene <- NULL

# Filter for Wang dataset and limit stages
wang_metadata <- metadata %>% filter(batch == "wang", dev_stage != 16)
wang_counts <- count_data[, colnames(count_data) %in% wang_metadata$sample_id, drop = FALSE]

photo_zm_genes <- read.table(file.path(in_dir, "photo_genes_Zm.txt"), header = FALSE, stringsAsFactors = FALSE)$V1

# --- 2. Functional Calculation of PS_Index ---
calculate_local_photo_proportion <- function(counts, photo_genes, meta) {
  photo_genes <- intersect(photo_genes, rownames(counts))
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ 1)
  dds <- estimateSizeFactors(dds)
  normalized_counts <- counts(dds, normalized = TRUE)
  
  total_counts <- colSums(normalized_counts)
  photo_sums <- colSums(normalized_counts[photo_genes, , drop = FALSE])
  photo_percent_counts <- (photo_sums / total_counts) * 100
  
  lower_bound <- quantile(photo_percent_counts, 0.05, na.rm = TRUE)
  upper_bound <- quantile(photo_percent_counts, 0.95, na.rm = TRUE)
  
  data.frame(
    sample_id = colnames(counts),
    photo_percent_counts = photo_percent_counts
  ) %>%
    mutate(ps_index = (photo_percent_counts - lower_bound) / (upper_bound - lower_bound)) %>%
    mutate(ps_index = pmin(pmax(ps_index, 0), 1))
}

wang_proportions <- calculate_local_photo_proportion(wang_counts, photo_zm_genes, wang_metadata)

wang_metadata_with_photo <- wang_metadata %>% inner_join(wang_proportions, by = "sample_id")

# --- 3. Panel A: Gradient ps_index plot ---
wp <- ggplot(wang_metadata_with_photo, aes(x = dev_stage, y = ps_index)) +
  geom_point(shape = 21, size = 1.2, fill = "green4", color = "black", alpha = 0.6, stroke = 0.3) +
  scale_x_continuous(breaks = 0:16) +   
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = NULL, x = "leaf section", y = "photosynthetic index") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text = element_text(color = "black"))

tiff(file.path(plot_dir, "wang_ps_index_solid_green.tiff"), width = 2.5, height = 2.5, units = "in", res = 600)
print(wp)
dev.off()

# --- 4. Metabolites Panels B-G ---
metabolites_mean <- read.table(file.path(in_dir, "metabolites_wang.txt"), header = TRUE, sep = "\t", check.names = FALSE)
metabolites_sd   <- read.table(file.path(in_dir, "metabolites_wang_sd.txt"), header = TRUE, sep = "\t", check.names = FALSE)

metabolites_long <- left_join(
  metabolites_mean %>% pivot_longer(cols = -metabolite, names_to = "dev_stage", values_to = "mean_value") %>% mutate(dev_stage = as.integer(dev_stage)),
  metabolites_sd %>% pivot_longer(cols = -metabolite, names_to = "dev_stage", values_to = "sd_value") %>% mutate(dev_stage = as.integer(dev_stage)),
  by = c("metabolite", "dev_stage")
) %>% distinct()

plot_metabolite_vs_photo <- function(metabolite_name, sec_axis_label, sec_axis_breaks, sec_axis_max, meta_color, output_file) {
  meta_df <- metabolites_long %>%
    mutate(metabolite = trimws(metabolite)) %>%
    filter(metabolite == metabolite_name) %>%
    mutate(n = 6, se_value = sd_value / sqrt(n))
  
  photo_df <- wang_metadata_with_photo %>% select(dev_stage, ps_index) %>% mutate(metric = "ps_index")
  
  photo_min <- -0.1
  photo_max <- 1.1
  meta_min <- 0
  meta_max <- sec_axis_max
  scale_factor <- (photo_max - photo_min) / (meta_max - meta_min)
  offset <- photo_min - scale_factor * meta_min
  
  meta_plot <- meta_df %>%
    mutate(y = offset + scale_factor * mean_value, ymin = offset + scale_factor * (mean_value - se_value), ymax = offset + scale_factor * (mean_value + se_value))
  
  p <- ggplot() +
    geom_point(data = photo_df, aes(x = dev_stage, y = ps_index), color = "green4", size = 0.8, alpha = 0.8) +
    geom_ribbon(data = meta_plot, aes(x = dev_stage, ymin = ymin, ymax = ymax), fill = meta_color, alpha = 0.1) +
    geom_line(data = meta_plot, aes(x = dev_stage, y = y), color = meta_color, linewidth = 0.5) +
    geom_point(data = meta_plot, aes(x = dev_stage, y = y), color = meta_color, size = 0.8, alpha = 1) +
    scale_x_continuous(breaks = 0:16) +
    scale_y_continuous(name = "ps_index", limits = c(-0.1, 1.1), breaks = seq(0, 1, 0.25),
                       sec.axis = sec_axis(trans = ~ (. - offset) / scale_factor, name = sec_axis_label, breaks = sec_axis_breaks)) +
    theme_classic(base_size = 8) +
    labs(x = "leaf section", y = NULL) +
    theme(legend.position = "none")
  
  tiff(file.path(plot_dir, output_file), width = 2.5, height = 2.5, units = "in", res = 600)
  print(p)
  dev.off()
}

plot_metabolite_vs_photo("Amino acids", "free amino acids (µmol/g FW)", seq(0, 45, by = 5), 48, "deeppink", "aa_photo_wang.tiff")
plot_metabolite_vs_photo("Protein", "protein (mg/g FW)", seq(0, 11, by = 1), 12, "deepskyblue", "protein_photo_wang.tiff")
plot_metabolite_vs_photo("Chlorophyll a", "Chlorophyll a (mg/g FW)", seq(0, 1.8, by = 0.2), 1.8, "darkolivegreen4", "chlorophyllA_photo_wang.tiff")
plot_metabolite_vs_photo("Rubisco (maximal)", "Rubisco activity (nmol/min/g FW)", seq(0, 10000, by = 1000), 10000, "green", "rubisco_photo_wang.tiff")
plot_metabolite_vs_photo("Glucose", "Glucose (µmol/g FW)", seq(0, 25, by = 5), 25, "orange", "glucose_photo_wang.tiff")
plot_metabolite_vs_photo("Starch", "Starch (µmol/g FW)", seq(0, 14, by = 2), 13, "coral", "starch_photo_wang.tiff")

# --- 5. Pearson Correlations ---
correlate_metabolite_with_photo <- function(metabolite_name, dev_stage_min = 3, dev_stage_max = 15) {
  photo_df <- wang_metadata_with_photo %>%
    filter(dev_stage >= dev_stage_min, dev_stage <= dev_stage_max) %>%
    group_by(dev_stage) %>%
    summarise(mean_photo_percent = mean(photo_percent_counts, na.rm = TRUE), .groups = "drop")
  
  metabolite_df <- metabolites_long %>%
    mutate(metabolite = trimws(metabolite)) %>%
    filter(metabolite == metabolite_name, dev_stage >= dev_stage_min, dev_stage <= dev_stage_max) %>%
    group_by(dev_stage) %>%
    summarise(mean_metabolite = mean(mean_value, na.rm = TRUE), .groups = "drop")
  
  cor_df <- inner_join(photo_df, metabolite_df, by = "dev_stage")
  cor_result <- cor.test(cor_df$mean_photo_percent, cor_df$mean_metabolite, method = "pearson")
  
  cat(sprintf("\n📊 Correlation between ps_index and %s: r = %.3f (p = %.3e)\n", 
              metabolite_name, cor_result$estimate, cor_result$p.value))
}

cat("\n--- Pearson Correlations (Dataset Verification) ---\n")
correlate_metabolite_with_photo("Chlorophyll a")
correlate_metabolite_with_photo("Rubisco (maximal)")
correlate_metabolite_with_photo("Starch")
correlate_metabolite_with_photo("Glucose")
correlate_metabolite_with_photo("Amino acids")
correlate_metabolite_with_photo("Protein")

cat("\nWang validation analysis completed successfully.\n")
