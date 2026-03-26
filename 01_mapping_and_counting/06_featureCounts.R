#!/usr/bin/env Rscript
# Optional: Set the R library path if needed for your environment
# .libPaths("/home/joo29/R/x86_64-pc-linux-gnu-library/4.2")

# Load required package
library(Rsubread)

# Get species and base directory from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
species <- args[1]
# Default base_dir to project root (one level up from this script in mapping_and_counting/)
base_dir <- if (length(args) >= 2) args[2] else ".."

# Set paths based on the species
bam_directory <- file.path(base_dir, species, "alignments")
annotation_file <- file.path(base_dir, species, paste0(species, "_hxr.gff"))

# List all BAM files
bam_files <- list.files(path = bam_directory, pattern = "\\.bam$", full.names = TRUE)

# Initialize a data frame to store count data
count_data <- NULL

# Loop through each BAM file and run featureCounts
for (bam_file in bam_files) {
  counts <- featureCounts(files = bam_file,
                          annot.ext = annotation_file,
                          isGTFAnnotationFile = TRUE,
                          GTF.featureType = "gene",
                          GTF.attrType = "ID",
                          isPairedEnd = FALSE,  # Set to TRUE if paired-end data
                          primaryOnly = TRUE,
                          strandSpecific = 1)

  # Extract the sample ID from the filename (assuming sample_id is right before '_Aligned')
  sample_id <- gsub("_Aligned.*", "", basename(bam_file))

  # Add count data to the dataframe with sample_id as the column name
  if (is.null(count_data)) {
    count_data <- counts$counts
    colnames(count_data) <- sample_id
  } else {
    count_data <- cbind(count_data, counts$counts)
    colnames(count_data)[ncol(count_data)] <- sample_id
  }
}

# Add gene IDs as row names
row.names(count_data) <- counts$annotation$GeneID

# Output path for the count data
output_file <- file.path(base_dir, species, paste0(species, "_counts.txt"))

# Export the count data to a .txt file
write.table(count_data, file = output_file, row.names = TRUE, sep = "\t", quote = FALSE)

cat("Counting completed for species:", species, "\n")
