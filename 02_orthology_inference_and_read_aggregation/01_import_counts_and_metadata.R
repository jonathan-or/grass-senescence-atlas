# 01_import_counts_and_metadata.R
# This script imports metadata, summary statistics, and gene-level counts.
# It standardizes species names, cleans transcript IDs, and prepares data for orthogroup aggregation.

# Command line arguments for base directory and batch
args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) > 0) args[1] else getwd()
batch <- if (length(args) > 1) args[2] else "2023" # Default to 2023 batch

library(dplyr)
library(tidyr)

# 1. Load Metadata and Summary Stats
metadata_file <- file.path(base_dir, "metadata.txt")
if (!file.exists(metadata_file)) {
    # Try batch-specific subdirectory if available
    metadata_file <- file.path(base_dir, if(batch == "2023") "pot-in-pot" else "remobilization", "metadata.txt")
}
metadata <- read.table(metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Load summary statistics (names vary by batch)
summary_file <- if(batch == "2023") "summary_statistics.txt" else "summary_statistics_all.txt"
summary_path <- file.path(base_dir, if(batch == "2023") "pot-in-pot" else "remobilization", summary_file)

if (file.exists(summary_path)) {
    summary_stats <- read.delim(summary_path, stringsAsFactors = FALSE)
    metadata <- metadata %>% left_join(summary_stats, by = "sample_id")
}

# 2. Import Counts for each species
species_list <- unique(metadata$species)
species_list <- species_list[species_list != "blank" & !is.na(species_list)]
count_data_frames <- list()

for (species in species_list) {
    # Define the file path based on the species name and batch structure
    # Standard format: [base_dir]/[batch_dir]/[species]/[species]_merged_counts.txt
    batch_dir <- if(batch == "2023") "pot-in-pot" else "remobilization"
    file_path <- file.path(base_dir, batch_dir, species, paste0(species, "_merged_counts.txt"))
    
    if (file.exists(file_path)) {
        count_data <- read.delim(file_path, row.names = 1, stringsAsFactors = FALSE)
        count_data_frames[[species]] <- count_data
    } else {
        message("Count file for ", species, " not found at: ", file_path)
    }
}

# 3. Clean Gene IDs
# Zea mays primary transcripts often have suffix "_T0XX" or similar that needs removal for OG matching
if ("Zea_mays" %in% names(count_data_frames)) {
    rownames(count_data_frames$Zea_mays) <- sub("_T0\\d{2}$", "", rownames(count_data_frames$Zea_mays))
}

# 4. Standardize Species and Accession Names
# Apply historical renaming for consistency
metadata <- metadata %>%
  mutate(
    species = case_when(
      species == "Andropogon_virginicus" ~ "Anatherum_virginicum",
      TRUE ~ species
    ),
    accession = case_when(
      accession == "Andropogon_virginicus" ~ "Anatherum_virginicum",
      accession == "Td_FL" ~ "Td_floridanum",
      TRUE ~ accession
    )
  )

# Rename in the count data list as well
names(count_data_frames)[names(count_data_frames) == "Andropogon_virginicus"] <- "Anatherum_virginicum"

# 5. Export for Aggregation
saveRDS(count_data_frames, file = file.path(base_dir, "02_orthology_inference_and_read_aggregation", paste0("count_data_frames_", batch, ".rds")))
saveRDS(metadata, file = file.path(base_dir, "02_orthology_inference_and_read_aggregation", paste0("metadata_", batch, ".rds")))

message("Import completed for batch: ", batch)
