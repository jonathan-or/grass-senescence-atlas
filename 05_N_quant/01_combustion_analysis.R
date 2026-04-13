# 05_N_quant/01_combustion_analysis.R
# Evaluates total nitrogen recycling and allocation between aboveground and underground structures

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) > 0) args[1] else file.path(getwd(), "05_N_quant", "data")
plot_dir <- if (length(args) > 1) args[2] else file.path(getwd(), "05_N_quant", "plots")
meta_dir <- file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

source(file.path(getwd(), "03_cross-species_transcriptome_analysis", "Functions.R"))

library(tidyverse)
library(lubridate)
library(ggpubr)
library(rstatix)
library(glue)

# === Step 1: Read in N data ===
sep_n <- read_tsv(file.path(in_dir, "sep_n_id.txt"), show_col_types = FALSE)
oct_n <- read_tsv(file.path(in_dir, "oct_n_id.txt"), show_col_types = FALSE)

n_data <- bind_rows(
  sep_n %>% mutate(timepoint = "Sep"),
  oct_n %>% mutate(timepoint = "Oct")
)

# === Step 2: Prepare metadata ===
metadata <- readRDS(file.path(meta_dir, "metadata.rds"))
metadata_clean <- metadata %>%
  mutate(
    accession = as.character(accession),
    date = as.Date(date)
  ) %>%
  distinct(accession, date, week, life_history, species)

# === Step 3: Merge week info into N data ===
n_data_clean <- n_data %>%
  mutate(
    accession = as.character(accession),
    date = as.Date(date)
  )

n_with_week <- left_join(n_data_clean, metadata_clean, by = c("accession", "date")) %>%
  mutate(week = as.factor(week))

# === Step 4: Create derived ratios ===
n_with_week <- n_with_week %>%
  mutate(
    biomass_ratio = udg_g / abg_g,
    n_total_ratio = total_n_udg / total_n_abg,
    n_alloc_ratio = udg_n / abg_n
  )

# === Step 5: Plotting `n_remob_ratio.tiff` ===
n_remob_ratio <- ggplot(n_with_week, aes(x = week, y = n_total_ratio)) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  geom_jitter(aes(color = species), width = 0.2, alpha = 0.6) +
  facet_wrap(~accession, nrow = 2) +
  scale_color_manual(values = species_palette, guide = "legend") +
  labs(x = "week", y = "undeground / aboveground N allocation") +
  scale_y_continuous(limits = c(0,5), breaks = seq(0, 4, 1))+
  theme_classic(base_size = 10)

tiff(file.path(plot_dir, "n_remob_ratio.tiff"), width = 4.1, height = 3.5, units = "in", res = 600)
print(n_remob_ratio)
dev.off()

# === Step 6: Link with ps_index trends for 2024 season ===
leaf_metadata_ps_index <- readRDS(file.path(meta_dir, "leaf_metadata_ps_index.rds"))
ps_2024 <- leaf_metadata_ps_index %>% filter(season == "2024")

ps_n_all <- ps_2024 %>%
  mutate(accession = as.character(accession)) %>%
  select(sample_id, accession, species, ps_index, week, photo_percent_counts) %>%
  inner_join(n_with_week %>% select(sample_id, n_total_ratio), by = "sample_id")

ps_n_all_exp <- ps_n_all %>% filter(n_total_ratio > 0)

# Exponential model
exp_fit <- nls(n_total_ratio ~ a * exp(b * ps_index), data = ps_n_all_exp,
               start = list(a = 0.1, b = 1))

a_exp <- coef(exp_fit)["a"]
b_exp <- coef(exp_fit)["b"]

# Add R-squared-like pseudo metric for nonlinear model
y <- ps_n_all_exp$n_total_ratio
y_hat <- predict(exp_fit)
ss_res <- sum((y - y_hat)^2)
ss_tot <- sum((y - mean(y))^2)
pseudo_r2 <- 1 - ss_res / ss_tot

fit_label <- glue("y = {round(a_exp, 2)} * exp({round(b_exp, 2)} * x)
pseudo R² = {round(pseudo_r2, 3)}")

ps_n_plot_exp <- ggplot(ps_n_all_exp, aes(x = ps_index, y = n_total_ratio)) +
  geom_point(aes(color = species), size = 2, alpha = 0.8) +
  stat_function(fun = function(x) a_exp * exp(b_exp * x), color = "black", linewidth = 1.2) +
  scale_color_manual(values = species_palette, guide = "none") +
  labs(x = "ps_index", y = "N allocation ratio (udg/abg)") +
  theme_classic(base_size = 10) +
  annotate("text", x = 0.95, y = max(y), hjust = 1, vjust = 1,
           label = fit_label, size = 2.5)

tiff(file.path(plot_dir, "ps_n_plot.tiff"), width = 2.5, height = 2.5, units = "in", res = 600)
print(ps_n_plot_exp)
dev.off()

cat("Combustion N quantification analysis complete.\n")
