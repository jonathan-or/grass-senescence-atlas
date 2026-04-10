# 04a_develop_GO_annotations.R
# This script processes raw PANNZER output files (*_GO.out) to generate species-specific GMT files 
# for Biological Process (BP), Cellular Component (CC), and Molecular Function (MF).
# It also consolidates these annotations dynamically across species into the core orthogroup space.

library(dplyr)
library(readr)
library(stringr)

# Pass directories as arguments or default to assumed project structure
args <- commandArgs(trailingOnly = TRUE)
go_dir <- if (length(args) > 0) args[1] else file.path(getwd(), "raw_annotations")
in_dir <- if (length(args) > 1) args[2] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")
out_dir <- if (length(args) > 2) args[3] else in_dir  # Save .rds objects to data/

if (!dir.exists(go_dir)) {
  message("Warning: Annotation directory not found. Please assure raw annotations are located at: ", go_dir)
  message("Skipping GO GMT development as raw .out files are not populated.")
  quit(save = "no", status = 0)
}

# List all Pannzer output files
go_files <- list.files(go_dir, pattern = "_GO.out$", full.names = TRUE)

if(length(go_files) == 0) {
  message("No _GO.out files found in ", go_dir, ". Expected raw PANNZER outputs. Exiting.")
  quit(save = "no", status = 0)
}

# 1. Process each species' annotation file into .gmt
for (file in go_files) {
  species <- str_extract(basename(file), "^[^_]+")
  pannzer_data <- suppressMessages(read_delim(file, "\t", col_names = TRUE, show_col_types = FALSE))
  
  # Process the qpid column to extract gene IDs
  if (species == "Zm") {
    pannzer_data$gene_id <- str_remove(pannzer_data$qpid, "_P\\d+$")
  } else {
    pannzer_data$gene_id <- str_remove(pannzer_data$qpid, "\\.1$")
  }
  
  pannzer_data$goid <- paste0("GO:", pannzer_data$goid)

  process_ontology <- function(ontology_label) {
    filtered_data <- pannzer_data %>%
      filter(ontology == ontology_label) %>%
      select(gene_id, goid, desc) %>%
      distinct() %>%
      group_by(goid, desc) %>%
      summarise(genes = paste(gene_id, collapse = "\t"), .groups = 'drop')

    if (nrow(filtered_data) > 0) {
      output_file <- file.path(go_dir, paste0("output_", ontology_label, "_", species, ".gmt"))
      writeLines(paste(filtered_data$goid, filtered_data$desc, filtered_data$genes, sep = "\t"), output_file)
      message("Generated file: ", output_file)
    }
  }

  process_ontology("BP")
  process_ontology("CC")
  process_ontology("MF")
}

# * Note: At this stage in the workflow, the Python script `generate_gmt_data.py`
# should be executed to synthesize `gmt_data.txt` describing term identities. *

# 2. Map Species GMT files to Core Orthogroups
library(clusterProfiler)

metadata <- readRDS(file.path(in_dir, "metadata.rds"))
gene_og_mapping <- readRDS(file.path(in_dir, "gene_og_mapping.rds"))

species_mapping <- c(
  "Ag" = "Andropogon_gerardii", "Av" = "Anatherum_virginicum",
  "Bl" = "Bothriochloa_laguroides", "Cz" = "Chrysopogon_zizanioides",
  "Et" = "Elionorus_tripsacoides", "Mg" = "Miscanthus_giganteus",
  "Pv" = "Panicum_virgatum", "Sb" = "Sorghum_bicolor",
  "Sn" = "Sorghastrum_nutans", "Ss" = "Schizachyrium_scoparium",
  "Td" = "Tripsacum_dactyloides", "Tt" = "Themeda_triandra",
  "Zd" = "Zea_diploperennis", "Zm" = "Zea_mays"
)

species_list <- unique(metadata$species)
species_gmt_list_BP <- list()
species_gmt_list_MF <- list()
species_gmt_list_CC <- list()

for (species in species_list) {
  species_code <- names(species_mapping)[which(species_mapping == species)]
  if (length(species_code) == 0) next
  
  species_file_BP <- file.path(go_dir, paste0("output_BP_", species_code, ".gmt"))
  species_file_MF <- file.path(go_dir, paste0("output_MF_", species_code, ".gmt"))
  species_file_CC <- file.path(go_dir, paste0("output_CC_", species_code, ".gmt"))

  if (file.exists(species_file_BP)) species_gmt_list_BP[[paste0(species, "_BP")]] <- read.gmt(species_file_BP)
  if (file.exists(species_file_MF)) species_gmt_list_MF[[paste0(species, "_MF")]] <- read.gmt(species_file_MF)
  if (file.exists(species_file_CC)) species_gmt_list_CC[[paste0(species, "_CC")]] <- read.gmt(species_file_CC)
}

# 3. Consolidation logic
combine_go_annotations <- function(species_gmt_list, suffix) {
  do.call(rbind, lapply(names(species_gmt_list), function(species) {
    gmt_df <- species_gmt_list[[species]]
    gmt_df$Species <- gsub(paste0("_", suffix, "$"), "", species) 
    return(gmt_df)
  }))
}

process_go_annotations <- function(gene_og_mapping, species_gmt_list, suffix) {
  combined_go <- combine_go_annotations(species_gmt_list, suffix) %>%
    mutate(Species = gsub(paste0("_", suffix, "$"), "", Species))
  
  go_og_mapping <- gene_og_mapping %>%
    inner_join(combined_go, by = c("GeneID" = "gene", "Species" = "Species"))
  
  og_go_df <- go_og_mapping %>%
    select(term, gene = Orthogroup) %>%
    distinct(term, gene)
  
  return(og_go_df)
}

og_go_df_BP <- process_go_annotations(gene_og_mapping, species_gmt_list_BP, "BP")
og_go_df_MF <- process_go_annotations(gene_og_mapping, species_gmt_list_MF, "MF")
og_go_df_CC <- process_go_annotations(gene_og_mapping, species_gmt_list_CC, "CC")

og_go_df_BP$term <- as.character(og_go_df_BP$term)
og_go_df_MF$term <- as.character(og_go_df_MF$term)
og_go_df_CC$term <- as.character(og_go_df_CC$term)

# Export consolidated RDS files
saveRDS(og_go_df_BP, file.path(out_dir, "og_go_df_BP.rds"))
saveRDS(og_go_df_MF, file.path(out_dir, "og_go_df_MF.rds"))
saveRDS(og_go_df_CC, file.path(out_dir, "og_go_df_CC.rds"))
message("Generated unified RDS matrices in ", out_dir)
