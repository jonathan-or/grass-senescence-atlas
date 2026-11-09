# 04_n_allocation_all_accessions.R
# Supplementary Figure 10:
#   A  underground / aboveground N allocation ratio, weeks 7 and 14, for all 8 accessions
#      (7 perennials + Zm_lh287Xlh244), with Tukey HSD (week 7 vs week 14 per accession)
#   B  Zea mays %N concentration ratios: underground/aboveground, kernel/aboveground,
#      kernel/underground, with one-way ANOVA and Tukey HSD
# Inputs: 05_N_quant/data/{sep,oct}_n_id.txt, 05_N_quant/data/Zm_kernel_N_quant.txt,
#         03_cross-species_transcriptome_analysis/data/metadata.rds
#
# Usage (from the repository root):
#   Rscript 05_N_quant/04_n_allocation_all_accessions.R [plot_dir]

library(tidyverse)
library(lubridate)
library(ggpubr)
library(rstatix)
library(glue)
library(scales)

cat("1. Loading N combustion data, metadata, and Maize kernel N quantification...\n")
in_dir <- file.path(getwd(), "05_N_quant", "data")
args <- commandArgs(trailingOnly = TRUE)
plot_dir <- if (length(args) > 0) args[1] else file.path(getwd(), "05_N_quant", "plots")
meta_dir <- file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

source(file.path(getwd(), "03_cross-species_transcriptome_analysis", "Functions.R"))

sep_n <- read_tsv(file.path(in_dir, "sep_n_id.txt"), show_col_types = FALSE) %>% mutate(timepoint = "Sep")
oct_n <- read_tsv(file.path(in_dir, "oct_n_id.txt"), show_col_types = FALSE) %>% mutate(timepoint = "Oct")

n_data <- bind_rows(sep_n, oct_n)

metadata <- readRDS(file.path(meta_dir, "metadata.rds"))
metadata_clean <- metadata %>%
  mutate(accession = as.character(accession), date = as.Date(date)) %>%
  distinct(accession, date, week, life_history, species)

n_with_week <- n_data %>%
  mutate(accession = as.character(accession), date = as.Date(date)) %>%
  left_join(metadata_clean, by = c("accession", "date")) %>%
  mutate(
    week = as.factor(week),
    biomass_ratio = udg_g / abg_g,
    n_total_ratio = total_n_udg / total_n_abg
  )

# ==============================================================================
# Task 1: 8-Accession Plot (7.2 in x 2.2 in, y-axis limit = 5, base_size = 8)
# ==============================================================================
cat("\n2. Generating 8-Accession n_remob_ratio plots (7.2 x 2.2 in, y-limit = 5)...\n")

p_8acc_base <- ggplot(n_with_week, aes(x = week, y = n_total_ratio)) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", width = 0.4, linewidth = 0.4) +
  geom_jitter(aes(color = species), width = 0.2, alpha = 0.65, size = 1.1) +
  facet_wrap(~accession, nrow = 1) +
  scale_color_manual(values = species_palette, name = "Species") +
  labs(x = "week", y = "underground / aboveground N allocation") +
  scale_y_continuous(limits = c(0, 5.0), breaks = seq(0, 5, 1), oob = oob_keep) +
  theme_classic(base_size = 8) +
  theme(
    strip.text = element_text(size = 6.5, face = "plain"),
    axis.title = element_text(size = 8, face = "plain"),
    axis.text.x = element_text(size = 6, color = "black"),
    axis.text.y = element_text(size = 6.5, color = "black"),
    plot.margin = margin(3, 5, 3, 5, "pt")
  )

# A. Plot WITH Legend
p_8acc_legend <- p_8acc_base +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 7, face = "plain"),
    legend.text = element_text(size = 6),
    legend.key.size = unit(9, "pt")
  )

ggsave(file.path(plot_dir, "n_remob_ratio_8acc_legend.png"), p_8acc_legend, width = 8.7, height = 2.2, units = "in", dpi = 300)
ggsave(file.path(plot_dir, "n_remob_ratio_8acc_legend.tiff"), p_8acc_legend, width = 8.7, height = 2.2, units = "in", dpi = 600)

# B. Plot WITHOUT Legend (Exact 7.2 in x 2.2 in)
p_8acc_nolegend <- p_8acc_base +
  theme(legend.position = "none")

ggsave(file.path(plot_dir, "n_remob_ratio_8acc_nolegend.png"), p_8acc_nolegend, width = 7.2, height = 2.2, units = "in", dpi = 300)
ggsave(file.path(plot_dir, "n_remob_ratio_8acc_nolegend.tiff"), p_8acc_nolegend, width = 7.2, height = 2.2, units = "in", dpi = 600)

# ==============================================================================
# Task 2: Tukey HSD Test for All 8 Accessions (Week 7 vs Week 14)
# ==============================================================================
cat("\n3. Performing Tukey HSD test on all 8 accessions (Week 7 vs Week 14)...\n")

aov_all <- aov(n_total_ratio ~ week * accession, data = n_with_week)
tukey_all <- TukeyHSD(aov_all)

tukey_8acc_df <- as.data.frame(tukey_all$`week:accession`) %>%
  rownames_to_column(var = "comparison") %>%
  filter(grepl("14:", comparison) & grepl("7:", comparison)) %>%
  mutate(
    acc1 = gsub("14:", "", sub("-.*", "", comparison)),
    acc2 = gsub("7:", "", sub(".*-", "", comparison))
  ) %>%
  filter(acc1 == acc2) %>%
  rename(accession = acc1, mean_diff_14_vs_7 = diff, tukey_p_adj = `p adj`) %>%
  mutate(
    mean_diff_14_vs_7 = round(mean_diff_14_vs_7, 3),
    tukey_p_adj = round(tukey_p_adj, 5),
    significance = case_when(
      tukey_p_adj < 0.001 ~ "***",
      tukey_p_adj < 0.01  ~ "**",
      tukey_p_adj < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

summary_8acc <- n_with_week %>%
  group_by(accession, species, life_history, week) %>%
  summarise(
    mean_val = round(mean(n_total_ratio), 3),
    sd_val = round(sd(n_total_ratio), 3),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = week,
    values_from = c(mean_val, sd_val, n),
    names_glue = "week_{week}_{.value}"
  )

final_8acc_df <- summary_8acc %>%
  left_join(tukey_8acc_df %>% select(accession, mean_diff_14_vs_7, tukey_p_adj, significance), by = "accession") %>%
  arrange(tukey_p_adj)

write_csv(final_8acc_df, file.path(plot_dir, "n_allocation_8acc_tukey_stats.csv"))

# ==============================================================================
# Task 3: Zea mays %N Concentration Ratio Figure (2.0 in x 2.0 in) & Stats
# ==============================================================================
cat("\n4. Processing Maize kernel N quantification data (%N concentration ratios)...\n")

zm_kernel_raw <- read_tsv(file.path(in_dir, "Zm_kernel_N_quant.txt"), show_col_types = FALSE) %>%
  rename(plot = `Sample ID`, kernel_n_pct = `N%`) %>%
  mutate(plot = as.numeric(plot))

# Merge kernel N data with Oct Maize samples and calculate %N concentration ratios
zm_oct_merged <- oct_n %>%
  filter(accession == "Zm_lh287Xlh244") %>%
  left_join(zm_kernel_raw %>% select(plot, kernel_n_pct), by = "plot") %>%
  mutate(
    udg_to_abg_n    = udg_n / abg_n,
    kernel_to_abg_n = kernel_n_pct / abg_n,
    kernel_to_udg_n = kernel_n_pct / udg_n
  )

# Reshape long for 3-ratio Maize plot
zm_ratios_long <- zm_oct_merged %>%
  select(plot, udg_to_abg_n, kernel_to_abg_n, kernel_to_udg_n) %>%
  pivot_longer(
    cols = c(udg_to_abg_n, kernel_to_abg_n, kernel_to_udg_n),
    names_to = "Ratio_Type_Raw",
    values_to = "Ratio_Value"
  ) %>%
  mutate(
    Ratio_Type = case_when(
      Ratio_Type_Raw == "udg_to_abg_n"    ~ "underground /\naboveground",
      Ratio_Type_Raw == "kernel_to_abg_n" ~ "kernel /\naboveground",
      Ratio_Type_Raw == "kernel_to_udg_n" ~ "kernel /\nunderground"
    ),
    Ratio_Type = factor(Ratio_Type, levels = c("underground /\naboveground",
                                               "kernel /\naboveground",
                                               "kernel /\nunderground")),
    dot_color = ifelse(Ratio_Type_Raw == "udg_to_abg_n", "darkslategray", "#FFD700")
  )

cat("\n5. Generating Zea mays 2.0 x 2.0 in %N concentration ratio plot...\n")

p_zm_ratios <- ggplot(zm_ratios_long, aes(x = Ratio_Type, y = Ratio_Value)) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", width = 0.4, linewidth = 0.4) +
  geom_jitter(aes(fill = dot_color), shape = 21, color = "black", size = 1.5, width = 0.12, alpha = 0.9) +
  scale_fill_identity() +
  labs(
    x = "concentration ratio",
    y = "N concentration ratio"
  ) +
  scale_y_continuous(limits = c(0, 4.2), breaks = seq(0, 4, 1)) +
  theme_classic(base_size = 8) +
  theme(
    axis.title = element_text(size = 8, face = "plain"),
    axis.text.x = element_text(size = 5.0, color = "black"),
    axis.text.y = element_text(size = 6.0, color = "black"),
    plot.margin = margin(2, 2, 2, 2, "pt")
  )

ggsave(file.path(plot_dir, "zea_mays_kernel_allocation_ratios_2.0in.png"), p_zm_ratios, width = 2.0, height = 2.0, units = "in", dpi = 300)
ggsave(file.path(plot_dir, "zea_mays_kernel_allocation_ratios_2.0in.tiff"), p_zm_ratios, width = 2.0, height = 2.0, units = "in", dpi = 600)

# Maize One-Way ANOVA & Tukey HSD
aov_zm <- aov(Ratio_Value ~ Ratio_Type, data = zm_ratios_long)
tukey_zm <- TukeyHSD(aov_zm)

tukey_zm_df <- as.data.frame(tukey_zm$Ratio_Type) %>%
  rownames_to_column(var = "comparison") %>%
  mutate(
    mean_difference = round(diff, 3),
    ci_lower = round(lwr, 3),
    ci_upper = round(upr, 3),
    tukey_p_adj = round(`p adj`, 5),
    significance = case_when(
      `p adj` < 0.001 ~ "***",
      `p adj` < 0.01  ~ "**",
      `p adj` < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  select(comparison, mean_difference, ci_lower, ci_upper, tukey_p_adj, significance)

write_csv(tukey_zm_df, file.path(plot_dir, "zea_mays_kernel_ratios_tukey_stats.csv"))

cat("Revision N quantification analysis complete!\n")
