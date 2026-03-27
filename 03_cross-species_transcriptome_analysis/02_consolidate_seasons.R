# Command-line arguments for input/output directories
args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) > 0) args[1] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")
out_dir <- if (length(args) > 1) args[2] else file.path(getwd(), "03_cross-species_transcriptome_analysis")

# Load required libraries
library(tidyverse)  # includes readr and other needed packages
library(lubridate)
library(patchwork)

# Function to process environmental data

process_env_data <- function(file, year) {
  read_tsv(file, locale = locale(encoding = "UTF-8")) %>%
    mutate(
      date = mdy(date),
      julian_days = yday(date),
      season = as.character(year),
      avg_soil_moisture = as.numeric(avg_soil_moisture)
    ) %>%
    filter(julian_days >= 200 & julian_days <= 300)
}

# Function to calculate GDD in Celsius
calculate_gdd_C <- function(T_max, T_min, T_base = 10, T_upper = 30) {
  pmax(0, ((pmin(T_max, T_upper) + pmax(T_min, T_base))/2 - T_base))
}

# Import and process both seasons' data
env_data <- bind_rows(
  process_env_data(file.path(in_dir, "newa_data_23.txt"), 2023),
  process_env_data(file.path(in_dir, "newa_data_24.txt"), 2024)
) %>%
  mutate(
    # Convert temperatures to Celsius
    max_air_temp_C = (max_air_temp - 32) * 5/9,
    min_air_temp_C = (min_air_temp - 32) * 5/9,
    avg_air_temp_C = (avg_air_temp - 32) * 5/9,
    # Calculate GDD
    daily_gdd_C = calculate_gdd_C(max_air_temp_C, min_air_temp_C)
  )

# Optional: verify the processed data
print(head(env_data))
print(summary(env_data))

# Adjust week calculation with Sundays as reference, starting at Julian day 200 (week 0)
# Calculate week numbers first (if not already in env_data)
env_data <- env_data %>%
  mutate(week = floor((julian_days - 200 + 4) / 7))

# Summarize weekly metrics more efficiently
weekly_summary <- env_data %>%
  group_by(week, season) %>%
  summarise(
    across(
      c(daily_gdd_C, solar_radiation, avg_air_temp_C),
      list(
        avg = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE)
      ),
      .names = "weekly_{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# Function to process metadata
process_metadata <- function(file) {
  readRDS(file) %>%
    mutate(
      date_formatted = mdy(date_formatted),
      julian_days = yday(date_formatted),
      week = floor((day_number + 2) / 7) + 1  # Adjust for Sundays
    )
}

# Import and process both years' metadata
metadata <- bind_rows(
  process_metadata(file.path(in_dir, "filtered_metadata_2023.rds")) %>% mutate(season = "2023"),
  process_metadata(file.path(in_dir, "filtered_metadata_2024.rds")) %>% mutate(season = "2024")
)

# Optional: inspect the processed data
print(head(metadata))

#fix julian bin assigment for week 13 (6-day interval sampling week)
# Fix date-formatted week bin assignment
metadata <- metadata %>%
  mutate(
    week = case_when(
      date_formatted == as.Date("2024-10-15") ~ 13,
      TRUE ~ week  # keep original otherwise
    )
  )

# Add week calculation to env_data
env_data <- env_data %>%
  mutate(week = floor((julian_days - 200 + 4) / 7))

# Calculate weekly summaries with start and end days
weekly_summary <- env_data %>%
  group_by(week, season) %>%
  summarize(
    julian_day_start = min(julian_days),
    julian_day_end = max(julian_days),
    avg_weekly_gdd = mean(daily_gdd_C, na.rm = TRUE),
    .groups = "drop"
  )

# Function to create composite plot for a given season
create_season_plot <- function(env_data, metadata, weekly_summary, season_year) {
  # Common parameters
  x_min <- 200
  x_max <- 305
  x_breaks <- seq(200, 305, by = 10)

  # Common theme
  common_theme <- theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )

  # Month tiles plot
  p_span_months <- ggplot() +
    geom_tile(
      data = data.frame(julian_days = 200:300) %>%
        mutate(
          date = as.Date(julian_days - 1, origin = paste0(season_year, "-01-01")),
          month = format(date, "%B")
        ),
      aes(x = julian_days, y = 0.5, fill = month),
      height = 0.7
    ) +
    geom_text(
      data = data.frame(julian_days = 200:300) %>%
        mutate(month = format(as.Date(julian_days - 1, origin = paste0(season_year, "-01-01")), "%B")) %>%
        group_by(month) %>%
        summarise(julian_days = mean(julian_days)),
      aes(x = julian_days, y = 0.5, label = month),
      size = 2.5
    ) +
    scale_x_continuous(limits = c(x_min, x_max), breaks = x_breaks) +
    scale_y_continuous(limits = c(0, 1), breaks = NULL) +
      scale_fill_manual(values = c(
        "July" = "#72c266",      # green (photosynthesis peak)
        "August" = "#a6d854",    # yellow-green
        "September" = "#ffd92f", # yellow
        "October" = "#fdae61"   # light orange
      )) +
    common_theme

  # GDD plot with overlaid sample points
  env_season <- env_data %>% filter(season == season_year)
  metadata_season <- metadata %>% filter(season == season_year)

  # Join metadata with env data to get GDD for sampling days
  sample_points <- metadata_season %>%
    left_join(env_season %>% select(julian_days, daily_gdd_C), by = "julian_days")

  p_gdd <- ggplot(env_season, aes(x = julian_days, y = daily_gdd_C)) +
    geom_line(
      data = weekly_summary %>% filter(season == season_year),
      aes(x = (julian_day_start + julian_day_end)/2, y = avg_weekly_gdd),
      color = "darkslateblue", size = 0.75, alpha = 0.3, inherit.aes = FALSE
    ) +
    geom_point(color = "darkslateblue", alpha = 0.4, size = 2) +
    geom_point(
      data = sample_points,
      aes(x = julian_days, y = daily_gdd_C),
      shape = 25, fill = "goldenrod1", color = "black", size = 2, stroke = 0.5,
      inherit.aes = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, max(env_season$daily_gdd_C, na.rm = TRUE) * 1.1),
      breaks = seq(0, max(env_season$daily_gdd_C, na.rm = TRUE) * 1.1, by = 5)
    ) +
    scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = x_breaks
    ) +
    theme_classic(base_size = 10) +
    labs(
      title = "",
      x = "julian days",
      y = "Growing Degree Days (GDDs)",
      color = "daily GDDs"
    ) +
    theme(
      legend.position = "none",
      legend.title = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  # Combine plots
  combined_plot <- p_span_months / p_gdd +
    plot_layout(heights = c(0.5, 3))

  return(combined_plot)
}

create_season_plot_week <- function(env_data, metadata, weekly_summary, season_year) {
  # Common params
  x_min <- 200
  x_max <- 305
  x_breaks <- seq(200, 305, by = 10)

  # Common theme for the header strips
  common_theme <- theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x  = element_blank(),
      panel.grid   = element_blank(),
      legend.position = "none"
    )

  # --- Build week ranges *from metadata* (this is the important bit) ---
  metadata_season <- metadata %>% filter(season == season_year)
  week_ranges <- metadata_season %>%
    group_by(week) %>%
    summarise(
      julian_day_start = floor(min(julian_days, na.rm = TRUE)),
      julian_day_end   = ceiling(max(julian_days, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(week) %>%
    mutate(week_mid = (julian_day_start + julian_day_end) / 2)

  # --- WEEK strip (top) ---
  p_span_weeks <- ggplot() +
    geom_tile(
      data = week_ranges,
      aes(
        x = week_mid,
        y = 0.5,
        width  = (julian_day_end - julian_day_start + 1),
        height = 0.7,
        fill = factor(week %% 2)  # simple alternating fill
      ),
      color = "grey80", linewidth = 0.2
    ) +
    geom_text(
      data = week_ranges,
      aes(x = week_mid, y = 0.5, label = paste0("W", week)),
      size = 2.6
    ) +
    scale_x_continuous(limits = c(x_min, x_max), breaks = x_breaks) +
    scale_y_continuous(limits = c(0, 1), breaks = NULL) +
    scale_fill_manual(values = c("#ececec", "#d9d9d9")) +
    common_theme

  # --- MONTH strip (unchanged) ---
  p_span_months <- ggplot() +
    geom_tile(
      data = data.frame(julian_days = 200:300) %>%
        mutate(
          date = as.Date(julian_days - 1, origin = paste0(season_year, "-01-01")),
          month = format(date, "%B")
        ),
      aes(x = julian_days, y = 0.5, fill = month),
      height = 0.7
    ) +
    geom_text(
      data = data.frame(julian_days = 200:300) %>%
        mutate(month = format(as.Date(julian_days - 1, origin = paste0(season_year, "-01-01")), "%B")) %>%
        group_by(month) %>%
        summarise(julian_days = mean(julian_days), .groups = "drop"),
      aes(x = julian_days, y = 0.5, label = month),
      size = 2.5
    ) +
    scale_x_continuous(limits = c(x_min, x_max), breaks = x_breaks) +
    scale_y_continuous(limits = c(0, 1), breaks = NULL) +
    scale_fill_manual(values = c(
      "July" = "#72c266",
      "August" = "#a6d854",
      "September" = "#ffd92f",
      "October" = "#fdae61"
    )) +
    common_theme

  # --- GDD plot (unchanged except for env_season definition) ---
  env_season <- env_data %>% filter(season == season_year)

  sample_points <- metadata_season %>%
    left_join(env_season %>% select(julian_days, daily_gdd_C), by = "julian_days")

  p_gdd <- ggplot(env_season, aes(x = julian_days, y = daily_gdd_C)) +
    geom_line(
      data = weekly_summary %>% filter(season == season_year),
      aes(x = (julian_day_start + julian_day_end)/2, y = avg_weekly_gdd),
      color = "darkslateblue", size = 0.75, alpha = 0.3, inherit.aes = FALSE
    ) +
    geom_point(color = "darkslateblue", alpha = 0.4, size = 2) +
    geom_point(
      data = sample_points,
      aes(x = julian_days, y = daily_gdd_C),
      shape = 25, fill = "goldenrod1", color = "black", size = 2, stroke = 0.5,
      inherit.aes = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, max(env_season$daily_gdd_C, na.rm = TRUE) * 1.1),
      breaks = seq(0, max(env_season$daily_gdd_C, na.rm = TRUE) * 1.1, by = 5)
    ) +
    scale_x_continuous(limits = c(x_min, x_max), breaks = x_breaks) +
    theme_classic(base_size = 10) +
    labs(x = "julian days", y = "Growing Degree Days (GDDs)") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  # --- Stack: Weeks (top) -> Months -> GDD ---
  combined_plot <- p_span_weeks / p_span_months / p_gdd +
    plot_layout(heights = c(0.45, 0.45, 3))

  return(combined_plot)
}

# Example:
# plot_2023 <- create_season_plot(env_data, metadata, weekly_summary, "2023")

# Create plots for both seasons
plot_2023 <- create_season_plot_week(env_data, metadata, weekly_summary, "2023")
plot_2024 <- create_season_plot(env_data, metadata, weekly_summary, "2024")

# Display plots
print(plot_2023)
print(plot_2024)


#tiff("gdd_plot_2024.tiff", width = 3, height = 3, units = "in", res = 600)
#print(plot_2024)  # Print the modified plot
#dev.off()  # Close the TIFF device

# Import count data for each season
og_counts_2023 <- readRDS(file.path(in_dir, "filtered_og_counts_2023.rds"))
og_counts_2024 <- readRDS(file.path(in_dir, "filtered_og_counts_2024.rds"))

# Get OGs from 2023
og_2023 <- rownames(og_counts_2023)

# Subset 2024 to keep only OGs present in 2023
og_counts_2024 <- og_counts_2024[og_2023, ]

# Make sure rownames are ordered the same way
og_counts_2023 <- og_counts_2023[sort(rownames(og_counts_2023)), ]
og_counts_2024 <- og_counts_2024[sort(rownames(og_counts_2024)), ]

# Merge the count matrices
og_counts <- cbind(og_counts_2023, og_counts_2024)

# Remove outliers from metadata
#rhizome_outliers
metadata <- metadata %>%
  filter(sample_id != "AC24BRTL000553")

metadata <- metadata %>%
  filter(sample_id != "AC24BRTL004298")

metadata <- metadata %>%
  filter(sample_id != "AC24BRTL004288")

metadata <- metadata %>%
  filter(sample_id != "AC24BRTL003818")

metadata <- metadata %>%
  filter(sample_id != "AC24BRTL002667")

#leaf_outliers
metadata <- metadata %>%
  filter(sample_id != "AC24BRTL004280")

metadata <- metadata %>%
  filter(sample_id != "AC24BRTL004182")

metadata <- metadata %>%
  filter(sample_id != "AC24BRTL003811")

metadata <- metadata %>%
  filter(sample_id != "AC24BRTL002347")#swap

metadata <- metadata %>%
  filter(sample_id != "AC24BRTL004287")#swap

#root_outliers
metadata <- metadata %>%
  filter(sample_id != "AC24BRTL003278")

metadata <- metadata %>%
  filter(sample_id != "AC24BRTL003291")

metadata <- metadata %>%
  filter(sample_id != "AC24BRTL002626")

#stem outlier (remobilizers)
#metadata <- metadata %>%
#  filter(sample_id != "AC24BRTL002869")

# Remove outlier from og_counts
#rhizome outliers
og_counts <- og_counts[, !colnames(og_counts) == "AC24BRTL000553"]
og_counts <- og_counts[, !colnames(og_counts) == "AC24BRTL004298"]
og_counts <- og_counts[, !colnames(og_counts) == "AC24BRTL004288"]
og_counts <- og_counts[, !colnames(og_counts) == "AC24BRTL003818"]
og_counts <- og_counts[, !colnames(og_counts) == "AC24BRTL002667"]

#leaf outliers
og_counts <- og_counts[, !colnames(og_counts) == "AC24BRTL004280"]
og_counts <- og_counts[, !colnames(og_counts) == "AC24BRTL004182"]
og_counts <- og_counts[, !colnames(og_counts) == "AC24BRTL003811"]
og_counts <- og_counts[, !colnames(og_counts) == "AC24BRTL002347"]
og_counts <- og_counts[, !colnames(og_counts) == "AC24BRTL004287"]

#root outliers
og_counts <- og_counts[, !colnames(og_counts) == "AC24BRTL003278"]
og_counts <- og_counts[, !colnames(og_counts) == "AC24BRTL003291"]
og_counts <- og_counts[, !colnames(og_counts) == "AC24BRTL002626"]

#stem outlier (remobilizers); decided to remove it before co-expression
#og_counts <- og_counts[, !colnames(og_counts) == "AC24BRTL002869"]

# Basic check if there are any NAs
any_nas <- anyNA(og_counts)
print(paste("Are there any NAs?", any_nas))

# Verify dimensions
print(dim(og_counts_2023))
print(dim(og_counts_2024))
print(dim(og_counts))

# Modify tissue column and create tissue_type column
metadata <- metadata %>%
  mutate(
    tissue_type = tissue,  # Save original tissue values first
    # Then modify main tissue categories
    tissue = case_when(
      tissue == "leaf" ~ "leaf",
      tissue == "rhizome" ~ "stem",
      tissue %in% c("root", "thick_roots", "fine_roots") ~ "root"
    )
  )

# Update tissue_type for Zea_mays
metadata <- metadata %>%
  mutate(
    tissue_type = case_when(
      species == "Zea_mays" & tissue_type == "rhizome" ~ "stalk",
      TRUE ~ tissue_type
    )
  )

# Verify the changes
print("Check tissue_type by species:")
print(table(metadata$species, metadata$tissue))

# Create a named vector for the replacements
accession_fixes <- c(
  "Pete" = "Td_Pete",
  "Liberty" = "Pv_Liberty",
  "Wisconsin" = "Sn_WI",
  "Iuka" = "Td_Iuka",
  "Niagara" = "Ag_Niagara",
  "PI1441932" = "Zd_momo",
  "Ahring" = "Ss_Ahring",
  "Texas_Sue" = "Td_Texas_Sue",
  "lh287Xlh244" = "Zm_lh287Xlh244",
  "Zea_mays" = "Zm_lh244",
  "Sorghum_bicolor" = "Sb_BTx623"
)

# Update the accession names
metadata <- metadata %>%
  mutate(
    accession = case_when(
      accession %in% names(accession_fixes) ~ accession_fixes[accession],
      TRUE ~ accession
    )
  )

metadata <- metadata %>%
  mutate(
    life_history = case_when(
      species %in% c("Zea_mays", "Sorghum_bicolor") ~ "annual",
      TRUE ~ "perennial"
      )
    )
    
  # Verify the new column
  print("Life history distribution:")
  print(table(metadata$life_history, metadata$species))


#Import counts 
#counts_2023 <- readRDS("filtered_count_data_frames_2023.rds")
#counts_2024 <- readRDS("filtered_count_data_frames_2024.rds")

# Step 2: Identify species groups
#shared_species <- intersect(names(counts_2023), names(counts_2024))
#only_2023_species <- setdiff(names(counts_2023), shared_species)
#only_2024_species <- setdiff(names(counts_2024), shared_species)

# Step 3: Merge overlapping species (full_join on gene rows)
#merged_counts <- map(shared_species, function(sp) {
#  full_join(
#    counts_2023[[sp]] %>% rownames_to_column("gene_id"),
#    counts_2024[[sp]] %>% rownames_to_column("gene_id"),
#    by = "gene_id"
#  ) %>%
#    column_to_rownames("gene_id")
#})

# Assign names to merged list
#names(merged_counts) <- shared_species

# Step 4: Collect all into final counts list
#counts <- c(
#  merged_counts,
#  counts_2023[only_2023_species],
#  counts_2024[only_2024_species]
#)

# Step 5 (optional): Basic check
#cat("✅ Merged species count data summary:\n")
#cat("Total species:", length(counts), "\n")
#cat("Species names:", paste(names(counts), collapse = ", "), "\n")
#cat("Gene and sample dimensions by species:\n")
#map(counts, ~ paste("Genes:", nrow(.), "Samples:", ncol(.))) %>% print()

# Step 1: Extract valid sample IDs from metadata
#valid_samples <- metadata$sample_id

# Step 2: Filter each species-specific count matrix in `counts`
#counts_filtered <- imap(counts, function(df, sp) {
  # Keep only columns (samples) present in metadata
#  filtered_df <- df[, colnames(df) %in% valid_samples, drop = FALSE]
  
  # Optional: Drop species with 0 remaining samples
#  if (ncol(filtered_df) > 0) {
#    return(filtered_df)
#  } else {
#    return(NULL)
#  }
#}) %>% compact()  # Remove NULLs (i.e., species with no valid samples)

# Step 3: Assign back to `counts`
#counts <- counts_filtered

# Step 4: Optional checks
#cat("✅ Filtered counts by metadata sample IDs:\n")
#cat("Remaining species:", length(counts), "\n")
#cat("Species names:", paste(names(counts), collapse = ", "), "\n")
#cat("Sample counts per species:\n")
#map_int(counts, ncol) %>% print()
# Count number of samples (columns) in each species-specific count matrix
#sample_counts <- map_int(counts, ncol)
# Total number of remaining samples
#total_samples <- sum(sample_counts)
# Print result
#cat("✅ Total number of samples remaining after filtering:", total_samples, "\n")
#saveRDS(counts, file = "filtered_species_counts_merged_2023_2024.rds")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(metadata,  file = file.path(out_dir, "metadata.rds"))

write.table(
  metadata,
  file = file.path(out_dir, "metadata_table.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

#saveRDS(og_counts, file = "og_counts.rds")

#Import Orthogroup data frames
metadata <- readRDS(file.path(out_dir, "metadata.rds"))
og_counts <- readRDS(file.path(in_dir, "og_counts.rds"))
OGs <- readRDS(file.path(in_dir, "OGs.rds"))
filteredOG <- readRDS(file.path(in_dir, "filteredOG.rds"))
write.table(filteredOG, 
            file = file.path(out_dir, "Supplementary_Dataset_Orthogroups.txt"), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
gene_og_mapping <- readRDS(file.path(in_dir, "gene_og_mapping.rds"))
partialOGs <- readRDS(file.path(in_dir, "partialOG.rds"))
counts <- readRDS(file.path(in_dir, "filtered_species_counts_merged_2023_2024.rds"))
# Import the dataframes for GO
og_go_df_BP <- readRDS(file.path(in_dir, "og_go_df_BP.rds"))
og_go_df_MF <- readRDS(file.path(in_dir, "og_go_df_MF.rds"))
og_go_df_CC <- readRDS(file.path(in_dir, "og_go_df_CC.rds"))
gmt_data <- readRDS(file.path(in_dir, "gmt_data.rds"))

library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)

# Load total gene counts
total_counts <- read_tsv(file.path(in_dir, "total_gene_counts_named.tsv"))

# Load orthogroup data and compute assigned counts
assigned_counts <- sapply(OGs[, -1], function(x) {
  sum(unlist(lapply(strsplit(x, ", "), function(g) length(g[g != ""]))))
})

assigned_df <- tibble(
  Species = names(assigned_counts),
  Assigned = as.numeric(assigned_counts)
)

# Merge and compute %
merged_df <- total_counts %>%
  inner_join(assigned_df, by = "Species") %>%
  mutate(PercentAssigned = (Assigned / Total) * 100)

ggplot(merged_df, aes(x = reorder(Species, -PercentAssigned), y = PercentAssigned)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%.1f%%", PercentAssigned)), vjust = -0.5, size = 4) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Percentage of Annotated Genes Assigned to Orthogroups",
    x = NULL,
    y = "Assigned Genes (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  ylim(0, 110)

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(forcats)

# Inputs you already have:
# counts      -> list(species_name = data.frame(genes x libraries))
# og_counts   -> data.frame(orthogroups x libraries), assigned counts only

stopifnot(is.list(counts), is.data.frame(og_counts))

all_assignments <- list()

og_libs <- colnames(og_counts)

for (species in names(counts)) {
  raw_counts <- counts[[species]]
  stopifnot(is.data.frame(raw_counts))

  total_per_library <- colSums(raw_counts, na.rm = TRUE)

  # intersect species libraries with the assigned-counts matrix
  common_libraries <- intersect(names(total_per_library), og_libs)

  # friendly diagnostics
  missing_in_og <- setdiff(names(total_per_library), og_libs)
  if (length(missing_in_og) > 0) {
    message(sprintf("[WARN] %s: %d libraries missing from og_counts (e.g., %s)",
                    species, length(missing_in_og), paste(head(missing_in_og, 3), collapse = ", ")))
  }
  if (length(common_libraries) == 0) {
    warning(sprintf("[WARN] %s: no overlapping libraries between raw counts and og_counts; skipping.", species))
    next
  }

  assigned_per_library <- colSums(og_counts[, common_libraries, drop = FALSE], na.rm = TRUE)

  df <- tibble(
    Library   = common_libraries,
    Total     = unname(total_per_library[common_libraries]),
    Assigned  = unname(assigned_per_library),
    Species   = species
  ) %>%
    mutate(
      Unassigned      = pmax(Total - Assigned, 0),
      Assigned_Pct    = if_else(Total > 0, Assigned / Total, NA_real_),
      Unassigned_Pct  = if_else(Total > 0, 1 - Assigned_Pct, NA_real_)
    )

  all_assignments[[species]] <- df
}

all_assignments <- bind_rows(all_assignments)

# Mean % assigned (across libraries)
mean_assigned_pct <- mean(all_assignments$Assigned_Pct, na.rm = TRUE)
message(sprintf("Mean Assigned_Pct across libraries: %.3f", mean_assigned_pct))

# Long format for plotting (Assigned stacked below Unassigned)
plot_data <- all_assignments %>%
  arrange(desc(Assigned_Pct)) %>%
  mutate(Library = factor(Library, levels = unique(Library))) %>%
  select(Library, Species, Assigned_Pct, Unassigned_Pct) %>%
  pivot_longer(c(Assigned_Pct, Unassigned_Pct),
               names_to = "Category", values_to = "Percent") %>%
  mutate(
    Category = recode(Category,
                      Assigned_Pct = "Assigned",
                      Unassigned_Pct = "Unassigned"),
    Category = factor(Category, levels = c("Unassigned", "Assigned"))
  )
print(plot_data)

library(ggplot2)
library(scales)

# expects: plot_data with columns Library, Species, Category (Unassigned/Assigned), Percent
# and mean_assigned_pct computed earlier

n_libs <- length(unique(plot_data$Library))

plot_assigned_per_library <- ggplot(plot_data, aes(x = Library, y = Percent, fill = Category)) +
  geom_col(width = 0.85, color = NA) +
  geom_hline(yintercept = mean_assigned_pct, linetype = "dashed", color = "red") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.05),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = c(Unassigned = "grey80", Assigned = "steelblue")
  ) +
  labs(
    x = paste0("libraries (n = ", n_libs, ")"),
    y = "% of reads assigned to core OGs"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    strip.text = element_text(face = "italic"),
    axis.text.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  ) +
  coord_cartesian(clip = "off")

# Save as TIFF with correct dimensions
tiff(file.path(out_dir, "assigned_percentage.tiff"), width = 2, height = 2, units = "in", res = 600)
print(plot_assigned_per_library)  # Print the modified plot
dev.off()  # Close the TIFF device


# ggsave("assigned_unassigned_by_library.png", p, width = 11, height = 14, dpi = 300)