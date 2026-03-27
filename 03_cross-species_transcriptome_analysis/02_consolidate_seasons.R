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
