# 06_merge_all_libraries.R
# Builds the unfiltered dataset of all 3,798 sequenced libraries (2023 + 2024, before
# library QC). Used only for the QC supplementary figures (Supplementary Figures 4-6);
# all other analyses use the 2,685 libraries that passed QC (data/og_counts.rds).
#
# Inputs are the per-season outputs of modules 01-02:
#   <dir_2023>/metadata.txt, summary_statistics.txt, summarized_expression_OG.rds, Orthogroups.tsv
#   <dir_2024>/metadata.txt, summary_statistics_all.txt, <species>/<species>_merged_counts.txt
#
# Usage (from the repository root):
#   Rscript 03_cross-species_transcriptome_analysis/06_merge_all_libraries.R <dir_2023> <dir_2024> [out_dir]
# Outputs (out_dir, default 03_cross-species_transcriptome_analysis/data):
#   metadata_all_libraries.rds/.txt, og_counts_all_libraries.rds

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 06_merge_all_libraries.R <dir_2023> <dir_2024> [out_dir]")
dir_2023 <- args[1]
dir_2024 <- args[2]
out_dir  <- if (length(args) > 2) args[3] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
library(dplyr)
library(readr)
library(tibble)
library(lubridate)
library(purrr)

cat("Starting merge of the whole dataset...\n")

# ==========================================
# 1. PROCESS 2023 (pot-in-pot) METADATA
# ==========================================
cat("Processing 2023 metadata...\n")
meta_2023 <- read.table(file.path(dir_2023, "metadata.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
stats_2023 <- read.delim(file.path(dir_2023, "summary_statistics.txt"), stringsAsFactors = FALSE)

meta_2023 <- meta_2023 %>%
  left_join(stats_2023, by = "sample_id") %>%
  mutate(
    species = ifelse(species == "Andropogon_virginicus", "Anatherum_virginicum", species),
    accesion = case_when(
      accesion == "Andropogon_virginicus" ~ "Anatherum_virginicum",
      accesion == "Td_FL" ~ "Td_floridanum",
      TRUE ~ accesion
    ),
    date = as.Date(date, format = "%b%d%Y"),
    date_formatted = format(date, "%m/%d/%Y")
  ) %>%
  dplyr::rename(accession = accesion)

start_date_2023 <- as.Date("2023-07-27")
meta_2023 <- meta_2023 %>%
  mutate(
    day_number = as.numeric(difftime(date, start_date_2023, units = "days")) + 1,
    week = ceiling(day_number / 7),
    season = "2023"
  )

# ==========================================
# 2. PROCESS 2024 (remobilization) METADATA
# ==========================================
cat("Processing 2024 metadata...\n")
meta_2024 <- read.table(file.path(dir_2024, "metadata.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
stats_2024 <- read.delim(file.path(dir_2024, "summary_statistics_all.txt"), stringsAsFactors = FALSE)

meta_2024 <- meta_2024 %>%
  left_join(stats_2024, by = "sample_id") %>%
  filter(species != "blank" & !is.na(species)) %>%
  mutate(
    date = as.Date(date, format = "%m/%d/%y"),
    date_formatted = format(date, "%m/%d/%Y")
  )

start_date_2024 <- as.Date("2024-07-27")
meta_2024 <- meta_2024 %>%
  mutate(
    day_number = as.numeric(difftime(date, start_date_2024, units = "days")) + 1,
    week = case_when(
      date == as.Date("2024-08-21") ~ 1,
      date == as.Date("2024-08-28") ~ 2,
      date == as.Date("2024-09-04") ~ 3,
      date == as.Date("2024-09-11") ~ 4,
      date == as.Date("2024-09-18") ~ 5,
      date == as.Date("2024-09-25") ~ 6,
      date == as.Date("2024-10-02") ~ 7,
      date == as.Date("2024-10-09") ~ 8,
      date == as.Date("2024-10-15") ~ 9,
      date == as.Date("2024-10-23") ~ 10,
      TRUE ~ ceiling(day_number / 7)
    ),
    season = "2024"
  )

# ==========================================
# 3. MERGE METADATA & APPLY JOINT TRANSFORMS
# ==========================================
cat("Merging metadata sets...\n")
full_metadata <- bind_rows(meta_2023, meta_2024)

# Apply common mutations
full_metadata <- full_metadata %>%
  mutate(
    date_formatted = mdy(date_formatted),
    julian_days = yday(date_formatted),
    week = floor((day_number + 2) / 7) + 1
  ) %>%
  mutate(
    week = case_when(
      date_formatted == as.Date("2024-10-15") ~ 13,
      TRUE ~ week
    )
  )

# Tissue modifications
full_metadata <- full_metadata %>%
  mutate(
    tissue_type = tissue,
    tissue = case_when(
      tissue == "leaf" ~ "leaf",
      tissue == "rhizome" ~ "stem",
      tissue %in% c("root", "thick_roots", "fine_roots") ~ "root"
    )
  ) %>%
  mutate(
    tissue_type = case_when(
      species == "Zea_mays" & tissue_type == "rhizome" ~ "stalk",
      TRUE ~ tissue_type
    )
  )

# Accession fixes
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

full_metadata <- full_metadata %>%
  mutate(
    accession = case_when(
      accession %in% names(accession_fixes) ~ accession_fixes[accession],
      TRUE ~ accession
    ),
    life_history = case_when(
      species %in% c("Zea_mays", "Sorghum_bicolor") ~ "annual",
      TRUE ~ "perennial"
    )
  )

# ==========================================
# 4. COMPILE 2023 (pot-in-pot) OG COUNTS
# ==========================================
cat("Compiling 2023 OG counts...\n")
se1 <- readRDS(file.path(dir_2023, "summarized_expression_OG.rds"))
# Sort rows to align
se1 <- lapply(se1, function(df) df[sort(rownames(df)), , drop = FALSE])
# cbind without names prepending
og_counts_2023 <- do.call(cbind, unname(se1))

# ==========================================
# 5. GENERATE 2024 (remobilization) OG COUNTS
# ==========================================
cat("Generating 2024 OG counts...\n")
# Load Orthogroups list
OGs <- read.table(file = file.path(dir_2023, "Orthogroups.tsv"), sep = '\t', header = TRUE, quote = '', stringsAsFactors = FALSE)
remob_species <- unique(meta_2024$species)
matching_species <- intersect(names(OGs), remob_species)

OGs_filtered <- OGs %>% select(Orthogroup, all_of(matching_species))
speciesNum <- apply(OGs_filtered[, 2:ncol(OGs_filtered)], 1, function(x) length(which(x == '')))
filteredOG_2024 <- OGs_filtered[which(speciesNum == 0), ]
geneNum <- apply(filteredOG_2024[, 2:ncol(filteredOG_2024)], 1, function(x) sum(unlist(lapply(strsplit(x, ', '), length))))
filteredOG_2024 <- filteredOG_2024[which(geneNum <= 75), ]

# Cleaning function
clean_gene_ids <- function(gene_col, species) {
  cleaned_genes <- sapply(gene_col, function(cell) {
    gene_ids <- unlist(strsplit(cell, ", "))
    if (species == "Zea_mays") {
      gene_ids <- gsub("_T0[0-9]+$", "", gene_ids)
    } else {
      gene_ids <- gsub("\\.\\d+$", "", gene_ids)
    }
    paste(gene_ids, collapse = ", ")
  })
  return(cleaned_genes)
}

# Clean IDs for species columns in 2024
for (sp in matching_species) {
  filteredOG_2024[[sp]] <- clean_gene_ids(filteredOG_2024[[sp]], sp)
}

# Mapping table
prepareGeneOGMapping <- function(fOG) {
  fOG %>%
    tidyr::pivot_longer(cols = -Orthogroup, names_to = "Species", values_to = "GeneIDs") %>%
    tidyr::separate_rows(GeneIDs, sep = ",\\s*") %>%
    dplyr::mutate(GeneIDs = trimws(GeneIDs)) %>%
    dplyr::select(Orthogroup, Species, GeneID = GeneIDs)
}
gene_og_mapping <- prepareGeneOGMapping(filteredOG_2024)

# Expression summarization function
sumExp_OG <- function(count_data, gene_og_mapping, all_OGs) {
  count_data_df <- as.data.frame(count_data)
  count_data_df$GeneID <- rownames(count_data)
  merged_data <- merge(count_data_df, gene_og_mapping, by = "GeneID")
  summarized_counts <- merged_data %>%
    dplyr::group_by(Orthogroup) %>%
    dplyr::summarise(across(starts_with("AC"), ~ sum(.x, na.rm = TRUE)))
  summarized_counts <- summarized_counts %>% tibble::column_to_rownames("Orthogroup")
  missing_OGs <- setdiff(all_OGs, rownames(summarized_counts))
  if (length(missing_OGs) > 0) {
    missing_df <- matrix(as.integer(0), nrow = length(missing_OGs), ncol = ncol(summarized_counts))
    rownames(missing_df) <- missing_OGs
    colnames(missing_df) <- colnames(summarized_counts)
    summarized_counts <- rbind(summarized_counts, missing_df)
  }
  return(summarized_counts)
}

# Load count files
count_data_frames_2024 <- list()
for (species in remob_species) {
  file_path <- file.path(dir_2024, species, paste0(species, "_merged_counts.txt"))
  if (file.exists(file_path)) {
    count_data <- read.delim(file_path, row.names = 1, stringsAsFactors = FALSE)
    if (species == "Zea_mays") {
      rownames(count_data) <- sub("_T0\\d{2}$", "", rownames(count_data))
    }
    count_data_frames_2024[[species]] <- count_data
  } else {
    stop("Count file for ", species, " not found at: ", file_path)
  }
}

all_OGs_2024 <- unique(filteredOG_2024$Orthogroup)
se2 <- lapply(names(count_data_frames_2024), function(species_name) {
  count_data <- count_data_frames_2024[[species_name]]
  species_og_mapping <- gene_og_mapping %>% filter(Species == species_name)
  summarized_counts <- sumExp_OG(count_data, species_og_mapping, all_OGs_2024)
  return(summarized_counts)
})
names(se2) <- names(count_data_frames_2024)

# Sort and cbind 2024 without names prepending
se2 <- lapply(se2, function(df) df[sort(rownames(df)), , drop = FALSE])
og_counts_2024 <- do.call(cbind, unname(se2))

# ==========================================
# 6. CONSOLIDATE 2023 AND 2024 OG COUNTS
# ==========================================
cat("Consolidating counts matrices...\n")
# Get intersection of orthogroups
common_ogs <- intersect(rownames(og_counts_2023), rownames(og_counts_2024))
cat("Total overlapping OGs:", length(common_ogs), "\n")

og_counts_2023_sub <- og_counts_2023[common_ogs, ]
og_counts_2024_sub <- og_counts_2024[common_ogs, ]

# Sort rows to be perfectly aligned
og_counts_2023_sub <- og_counts_2023_sub[sort(rownames(og_counts_2023_sub)), ]
og_counts_2024_sub <- og_counts_2024_sub[sort(rownames(og_counts_2024_sub)), ]

og_counts <- cbind(og_counts_2023_sub, og_counts_2024_sub)

# Keep only columns (samples) present in full_metadata
meta_sample_ids <- full_metadata$sample_id
og_counts <- og_counts[, colnames(og_counts) %in% meta_sample_ids]

# ==========================================
# 7. RUN VALIDATION CHECKS
# ==========================================
cat("Running integrity checks...\n")
counts_sample_ids <- colnames(og_counts)

cat("Libraries in full_metadata:", length(meta_sample_ids), "\n")
cat("Libraries in og_counts:", length(counts_sample_ids), "\n")

# Align counts to metadata ordering
if (!all(counts_sample_ids %in% meta_sample_ids)) {
  stop("Some sample IDs in count data are missing from metadata!")
}
if (!all(meta_sample_ids %in% counts_sample_ids)) {
  stop("Some sample IDs in metadata are missing from count data!")
}

# Reorder og_counts columns to match full_metadata$sample_id order
og_counts <- og_counts[, meta_sample_ids]

if (!all(colnames(og_counts) == full_metadata$sample_id)) {
  stop("Alignment check failed!")
}

if (anyNA(og_counts)) {
  stop("NA values found in final count matrix!")
}

cat("Integrity checks passed successfully!\n")

# ==========================================
# 8. SAVE OUTPUTS
# ==========================================
cat("Saving final objects...\n")
saveRDS(full_metadata, file = file.path(out_dir, "metadata_all_libraries.rds"))
saveRDS(og_counts, file = file.path(out_dir, "og_counts_all_libraries.rds"))

write.table(
  full_metadata,
  file = file.path(out_dir, "metadata_all_libraries.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Done:", ncol(og_counts), "libraries saved to", out_dir, "\n")
