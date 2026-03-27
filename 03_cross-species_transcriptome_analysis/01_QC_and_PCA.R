args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) > 0) args[1] else file.path(getwd(), "02_orthology_inference_and_read_aggregation")
out_dir <- if (length(args) > 1) args[2] else file.path(getwd(), "03_cross-species_transcriptome_analysis")

# For this script you'll need to load the following objects:
count_data_frames        <- readRDS(file.path(in_dir, "count_data_frames.rds"))
filteredOG               <- readRDS(file.path(in_dir, "filteredOG.rds"))
metadata                 <- readRDS(file.path(in_dir, "metadata.rds"))
summarized_expression_OG <- readRDS(file.path(in_dir, "summarized_expression_OG.rds"))

#load required libraries
library(DESeq2)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(dplyr)

#Run PCA for QC for each species 
run_interactive_pca <- function(species_name, count_data_frames, metadata) {
  # Subset metadata for the selected species
  species_metadata <- metadata %>% filter(species == species_name)
  
  # Subset count data for the selected species
  species_counts <- count_data_frames[[species_name]]
  
  # Filter metadata based on mapping rate and read count criteria
  species_metadata <- species_metadata %>%
    filter(MappedReads > 100000, MappedReads / RawReads > 0.35)
  
  # Filter count data based on filtered metadata sample IDs
  filtered_count_data <- species_counts[, colnames(species_counts) %in% species_metadata$sample_id]
  
  # Filter by library presence threshold
  threshold_libraries <- ceiling(0.3 * ncol(filtered_count_data))
  expressed_in_enough_libraries <- rowSums(filtered_count_data >= 1) >= threshold_libraries
  filtered_count_data <- filtered_count_data[expressed_in_enough_libraries, ]
  
  # Create DESeqDataSet and normalize
  dds <- DESeqDataSetFromMatrix(countData = filtered_count_data, colData = species_metadata, design = ~ day_number)
  dds <- estimateSizeFactors(dds)
  vst_data <- vst(dds, blind = TRUE)
  vst_mat <- assay(vst_data)

  # Perform PCA
  pca_res <- prcomp(t(vst_mat))
  
  # Prepare PCA data frame with metadata
  pca_scores_df <- as.data.frame(pca_res$x)
  pca_scores_df$sample_id <- rownames(pca_scores_df)
  pca_scores_metadata <- left_join(pca_scores_df, species_metadata, by = "sample_id")
  
  # Define color mappings
  continuous_colors <- list(c(0, 'purple'), c(1, 'orange'))
  categorical_colors_tissue <- c("#37b84f", "#405b6c")
  
  # Generate distinct colors for each batch
  unique_batches <- unique(pca_scores_metadata$batch)
  batch_colors <- setNames(RColorBrewer::brewer.pal(min(length(unique_batches), 8), "Dark2"), unique_batches)
  batch_color_vector <- batch_colors[as.character(pca_scores_metadata$batch)]
  
  # Create Plotly interactive plot with dropdown menu for color options
  interactive_plot <- plot_ly(
    data = pca_scores_metadata,
    x = ~PC1, y = ~PC2,
    text = ~paste("Sample ID:", sample_id, "<br>Tissue:", tissue, "<br>Batch:", batch, "<br>Mapping Rate:", MappedReads / RawReads),
    hoverinfo = "text",
    type = 'scatter', mode = 'markers',
    marker = list(size = 12, color = pca_scores_metadata$day_number, colorscale = continuous_colors, showscale = TRUE) # Default to day colors
  ) %>%
    layout(
      title = paste("PCA of VST Transformed Counts for", species_name),
      xaxis = list(title = paste("PC1 -", round(100 * summary(pca_res)$importance[2, 1], 2), "% variance")),
      yaxis = list(title = paste("PC2 -", round(100 * summary(pca_res)$importance[2, 2], 2), "% variance")),
      updatemenus = list(
        list(
          buttons = list(
            list(
              method = "restyle", 
              args = list(
                list(marker = list(color = pca_scores_metadata$day_number, colorscale = continuous_colors, showscale = TRUE, size = 12))
              ),
              label = "Color by Day"
            ),
            list(
              method = "restyle", 
              args = list(
                list(marker = list(color = as.numeric(factor(pca_scores_metadata$tissue)), 
                                   colorscale = list(c(0, categorical_colors_tissue[1]), c(1, categorical_colors_tissue[2])), showscale = FALSE, size = 12))
              ),
              label = "Color by Tissue"
            ),
            list(
              method = "restyle", 
              args = list(
                list(marker = list(color = pca_scores_metadata$MappedReads / pca_scores_metadata$RawReads, colorscale = continuous_colors, showscale = TRUE, size = 12))
              ),
              label = "Color by Mapping Rate"
            ),
            list(
              method = "restyle", 
              args = list(
                list(marker = list(color = batch_color_vector, showscale = FALSE, size = 12))
              ),
              label = "Color by Batch"
            )
          ),
          direction = "down",
          x = 0.1, xanchor = "left", y = 1.1, yanchor = "top"
        )
      )
    )
  
  return(interactive_plot)
}

# Run for a specific species to test
species_plot <- run_interactive_pca("Miscanthus_giganteus", count_data_frames, metadata)
species_plot

# List to store the plots
species_plots <- list()

# Loop through each species in count_data_frames
for (species_name in names(count_data_frames)) {
  species_plots[[species_name]] <- run_interactive_pca(species_name, count_data_frames, metadata)
}

# Export each interactive PCA plot to HTML
for (species_name in names(species_plots)) {
  saveWidget(
    species_plots[[species_name]], 
    file = paste0(species_name, "_PCA_plot.html"),
    selfcontained = TRUE
  )
}

#After visual inspection of the plots for swaps (leaf to roots and vice-versa) and outliers (kept this to a minimum, 
#didn't remove samples in which the outlier call was not evident)

outlier_samples <- c(
  # Miscanthus
  "AC23BRTL001740",  "AC24BRTL000553",
  # Panicum
  "AC24BRTL002335",
  # Tripsacum
  "AC24BRTL000712", "AC23BRTL001808",
  # Andropogon gerardii
  "AC24BRTL001658", "AC23BRTL000794", "AC24BRTL002171", "AC24BRTL000122", 
  "AC24BRTL000539", "AC24BRTL000347", "AC24BRTL000890", "AC24BRTL000827", 
  "AC24BRTL000731", "AC24BRTL000698", "AC23BRTL001754",
  # Andropogon virginicus
  "AC24BRTL000628", "AC24BRTL000501", "AC24BRTL001365",
  # Schizachyrium scoparium
  "AC24BRTL000475", "AC24BRTL002324", "AC24BRTL001523",
  # Sorghastrum nutans
  # Sorghum bicolor
  "AC24BRTL000534", 
  # Panicum virgatum
  "AC24BRTL002335"
  # Zea mays
  # Zea diploperennis
  # Elionorus
  # Chrysopogon
)

# Apply filtering criteria to metadata (below these parameters the leaf/root clustering is lost; it seems these samples contain
#high rRNA/bacterial/fungal or not enough reads to get biological insights from them)
filtered_metadata <- metadata %>%
  filter(MappedReads / RawReads > 0.35, MappedReads > 100000)

# Filter out outliers from metadata
filtered_metadata <- filtered_metadata[!filtered_metadata$sample_id %in% outlier_samples, ]

# Display the remaining libraries after filtering
cat("Remaining libraries after filtering:\n")
print(filtered_metadata$sample_id)

# Convert OG row names to a column in each data frame within summarized_expression
summarized_expression_aligned <- lapply(summarized_expression_OG, function(df) {
  df %>% tibble::rownames_to_column("OG")
})

# Merge all data frames by the "OG" column
merged_og_counts <- Reduce(function(x, y) {
  full_join(x, y, by = "OG")
}, summarized_expression_aligned)

# Set "OG" as rownames and remove it from columns
rownames(merged_og_counts) <- merged_og_counts$OG
merged_og_counts <- merged_og_counts[, -which(names(merged_og_counts) == "OG")]

# Filter merged_og_counts to retain only columns (libraries) in filtered_metadata
filtered_og_counts <- merged_og_counts[, colnames(merged_og_counts) %in% filtered_metadata$sample_id]

# Reorder columns in filtered_og_counts to match the order of sample IDs in filtered_metadata
filtered_og_counts <- filtered_og_counts[, filtered_metadata$sample_id]

# Check if column names in filtered_og_counts and sample IDs in filtered_metadata match
if (!all(colnames(filtered_og_counts) == filtered_metadata$sample_id)) {
  stop("Sample IDs in count data and metadata do not match.")
}

filtered_libraries <- filtered_metadata$sample_id
filtered_count_data_frames <- lapply(count_data_frames, function(df) {
  df[, colnames(df) %in% filtered_libraries, drop = FALSE]
})

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(filtered_metadata,          file = file.path(out_dir, "filtered_metadata.rds"))
saveRDS(filtered_og_counts,         file = file.path(out_dir, "filtered_og_counts.rds"))
saveRDS(filtered_count_data_frames, file = file.path(out_dir, "filtered_count_data_frames.rds"))
