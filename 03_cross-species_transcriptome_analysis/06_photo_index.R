# 06_photo_index.R
# Computes the photosythesis index (ps_index) metric and evaluates environmental correlation.

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) > 0) args[1] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")
plot_dir <- if (length(args) > 1) args[2] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "plots")

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

source(file.path(getwd(), "03_cross-species_transcriptome_analysis", "Functions.R"))

library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(lubridate)
library(DESeq2)

# --- 1. Identify 115 Photosynthetic Associated OGs --- #
cat("Identifying the core photosynthetic orthogroup set...\n")
photo_go_terms <- c(
  "GO:0009535", "GO:0015979", "GO:0009523", "GO:0009522", "GO:0009765",
  "GO:0055035", "GO:0009941", "GO:0010287", "GO:0016168",
  "GO:0009416", "GO:0009654", "GO:0018298", "GO:0009532", "GO:0019904",
  "GO:0019253", "GO:0022900", "GO:0010207", "GO:0009773", "GO:0010598",
  "GO:0015995", "GO:0009570", "GO:0031977", "GO:0006094", "GO:0030388",
  "GO:0010206", "GO:0042170", "GO:0045156", "GO:0010319",
  "GO:0009579", "GO:0009644", "GO:0010196", "GO:0016620", "GO:0009534",
  "GO:0050661", "GO:0031969", "GO:0009658", "GO:0009642", "GO:0009543",
  "GO:0051537", "GO:0051536", "GO:0042651", "GO:0009055","GO:0015977"
)

# Load the negative enrichment results from pseudotime modeling
neg_go_res <- read.table(file.path(plot_dir, "leaf_neg_OG_GO_enrichment.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter for the target photosynthesis GO IDs and extract matching genes
photo_go_results <- neg_go_res %>%
  filter(ID %in% photo_go_terms) %>%
  pull(geneID)

photo_og_vector <- unique(unlist(strsplit(photo_go_results, "/")))
cat("Number of unique photosynthesis-associated orthogroups identified:", length(photo_og_vector), "\n")

# --- 2. Calculate ps_index for leaf libraries --- #
cat("Loading matrices and calculating ps_index scaling...\n")
metadata <- readRDS(file.path(in_dir, "metadata.rds"))
og_counts <- readRDS(file.path(in_dir, "og_counts.rds"))

leaf_metadata <- metadata %>% filter(tissue == "leaf")
leaf_og_counts <- og_counts[, colnames(og_counts) %in% leaf_metadata$sample_id]

# Calculate Photo OGs percentage arrays and scale index bounding via Functions.R
photo_og_proportions <- calculate_photo_hog_proportion(
  hog_counts = leaf_og_counts,
  photo_hogs = photo_og_vector,
  metadata = leaf_metadata
)

# Join results and load pseudotime model results to get timeline mappings
leaf_pseudotime <- readRDS(file.path(in_dir, "leaf_pseudotime_model.rds"))
pseudotime_df <- leaf_pseudotime$pseudotime_df

leaf_metadata_ps <- leaf_metadata %>%
  inner_join(photo_og_proportions, by = "sample_id") %>%
  inner_join(pseudotime_df, by = c("sample_id" = "SampleID")) %>%
  arrange(julian_days) %>%
  mutate(sample_id = factor(sample_id, levels = sample_id)) # Lock factor ordering

# Produce Figure 2D (Photosynthesis Index spanning pseudotime)
photo_ogs_plot <- ggplot(leaf_metadata_ps, aes(x = Pseudotime, y = photo_percent_counts, color = julian_days)) +
  geom_point(size = 1.2, alpha = 0.8, color="black", fill="green", shape=21) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") +
  labs(title = element_blank(), y = "percentage (%) of total transcriptome", x = "pseudotime")

ggsave(file.path(plot_dir, "Figure2D_photo_ogs_v2.tiff"), photo_ogs_plot, width = 2, height = 3, units = "in", dpi = 600, compression = "lzw")

# --- 3. Consolidate Env Data & Model Environmental Correlation --- #
cat("Mapping Daily and Weekly environmental progressions (GDDs)...\n")

calculate_gdd_C <- function(T_max, T_min, T_base = 10, T_upper = 30) {
  pmax(0, ((pmin(T_max, T_upper) + pmax(T_min, T_base))/2 - T_base))
}

process_env_data <- function(file, year) {
  read_tsv(file, locale = locale(encoding = "UTF-8"), show_col_types = FALSE) %>%
    mutate(
      date = mdy(date),
      julian_days = yday(date),
      season = as.character(year),
      avg_soil_moisture = as.numeric(avg_soil_moisture)
    ) %>%
    filter(julian_days >= 200 & julian_days <= 300)
}

env_data <- bind_rows(
  process_env_data(file.path(in_dir, "newa_data_23.txt"), 2023),
  process_env_data(file.path(in_dir, "newa_data_24.txt"), 2024)
) %>%
  mutate(
    max_air_temp_C = (max_air_temp - 32) * 5/9,
    min_air_temp_C = (min_air_temp - 32) * 5/9,
    daily_gdd_C = calculate_gdd_C(max_air_temp_C, min_air_temp_C),
    week = floor((julian_days - 200 + 4) / 7)
  )

# Accumulate Heat to establish total physiological week load
weekly_response_env <- env_data %>%
  group_by(season, week) %>%
  summarize(weekly_gdd_pulse = sum(daily_gdd_C, na.rm = TRUE), .groups = "drop")

plot_response_df <- leaf_metadata_ps %>%
  left_join(weekly_response_env, by = c("season", "week")) %>%
  filter(!is.na(weekly_gdd_pulse))

# Map Dual-Axis plots for each season natively
plot_2023 <- generate_ps_gdd_plot(plot_response_df, "2023", c(1, 14), seq(1, 14, 1))
ggsave(file.path(plot_dir, "photo_gdd_2023_v3.tiff"), plot_2023, width = 3.5, height = 2.5, units = "in", dpi = 600, compression = "lzw")

plot_2024 <- generate_ps_gdd_plot(plot_response_df, "2024", c(5, 14), seq(5, 14, 1))
ggsave(file.path(plot_dir, "photo_gdd_2024_v3.tiff"), plot_2024, width = 2.8, height = 2.5, units = "in", dpi = 600, compression = "lzw")

# Calculate base cross-species daily tracking correlations
cor_input_df_clean <- leaf_metadata_ps %>%
  left_join(env_data %>% select(julian_days, daily_gdd_C, season), by = "julian_days") %>%
  drop_na(photo_percent_counts, daily_gdd_C)

cor_all <- cor(cor_input_df_clean$photo_percent_counts, cor_input_df_clean$daily_gdd_C, method = "pearson")
cat("Pearson correlation between daily GDD and photo OG percentage:", round(cor_all, 3), "\n")

cor_plot <- ggplot(cor_input_df_clean, aes(x = daily_gdd_C, y = photo_percent_counts)) +
  geom_point(alpha = 0.6, size = 1, color = "#5a5a5a") +
  geom_smooth(method = "lm", se = TRUE, color = "purple") +
  theme_classic(base_size = 10) +
  labs(
    x = "Daily GDD (°C)",
    y = "Photo OGs (% of total counts)"
  ) + annotate("text", x = min(cor_input_df_clean$daily_gdd_C), y = max(cor_input_df_clean$photo_percent_counts), label = paste0("r = ", round(cor_all,3)), hjust=0, vjust=1)

ggsave(file.path(plot_dir, "photo_percentage_gdd_correlation.tiff"), cor_plot, width = 2.5, height = 2.5, units="in", dpi=600, compression="lzw")

cat("Photo Index evaluation script completed securely.\n")
