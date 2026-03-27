# 01_import_counts_and_metadata.R
# This script illustrates the processing flow to import metadata, summary statistics, 
# and gene-level counts outputted by the mapping pipeline.
# It expects the following input files in the specified base directory:
#  - metadata.txt
#  - summary_statistics.txt
#  - A directory for each species containing a file named '{species}_merged_counts.txt'

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) > 0) args[1] else getwd()

library(dplyr)
library(tidyr)

# 1. Load Metadata and Summary Stats
metadata_file <- file.path(base_dir, "metadata.txt")
summary_file <- file.path(base_dir, "summary_statistics.txt")

if (!file.exists(metadata_file)) {
  stop("Metadata file not found at: ", metadata_file)
}
metadata <- read.table(metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

if (file.exists(summary_file)) {
    summary_stats <- read.delim(summary_file, stringsAsFactors = FALSE)
    metadata <- metadata %>% left_join(summary_stats, by = "sample_id")
} else {
    warning("Summary statistics file not found at: ", summary_file, ". Proceeding without it.")
}

# Remove blank/NA species rows if any
metadata <- metadata %>% filter(species != "blank" & !is.na(species))

# 2. Import Counts for each species
species_list <- unique(metadata$species)
count_data_frames <- list()

for (species in species_list) {
    # Expected format: base_dir/[species]/[species]_merged_counts.txt
    file_path <- file.path(base_dir, species, paste0(species, "_merged_counts.txt"))
    
    if (file.exists(file_path)) {
        # Read the count file as a data frame, keeping gene names as row names
        count_data <- read.delim(file_path, row.names = 1, stringsAsFactors = FALSE)
        count_data_frames[[species]] <- count_data
    } else {
        warning("Count file for ", species, " not found at: ", file_path)
    }
}

# Verify imported sample IDs against metadata
sample_ids <- metadata$sample_id 
count_columns_all <- unique(unlist(lapply(count_data_frames, colnames)))

missing_in_counts <- setdiff(sample_ids, count_columns_all)
if (length(missing_in_counts) > 0) {
    message("Sample IDs in metadata missing from count data: ", paste(missing_in_counts, collapse = ", "))
} else {
    message("All sample IDs in metadata are present in count data.")
}

# 3. Clean Gene IDs
# Zea mays primary transcripts often have suffix "_T0XX" which needs removal for OrthoFinder matching
if ("Zea_mays" %in% names(count_data_frames)) {
    rownames(count_data_frames$Zea_mays) <- sub("_T0\\d{2}$", "", rownames(count_data_frames$Zea_mays))
}

# 4. Standardize Species and Accession Names (Optional/Project-specific)
# Example of renaming specific taxa correctly historically
# (e.g., A. virginicus changed to A. virginicum)
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

if ("Andropogon_virginicus" %in% names(count_data_frames)) {
    names(count_data_frames)[names(count_data_frames) == "Andropogon_virginicus"] <- "Anatherum_virginicum"
}

# 5. Export for Aggregation
output_dir <- file.path(base_dir, "02_orthology_inference_and_read_aggregation")
dir.create(output_dir, showWarnings = FALSE)

saveRDS(count_data_frames, file = file.path(output_dir, "count_data_frames.rds"))
saveRDS(metadata, file = file.path(output_dir, "metadata.rds"))

message("Import completed.")
