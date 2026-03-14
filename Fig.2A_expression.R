# Load necessary libraries
library(GenomicFeatures)
library(pheatmap)
library(edgeR)
library(svglite)

# Clear environment and set working directory
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load sample name mapping file
sample_mapping <- read.csv("sample_mapping.csv", header = TRUE)

# Define MHC-I gene IDs (Orcu-U1, Orcu-U2, etc.)
mhc1_gene_ids <- c("Orcu-U1", "Orcu-U2", "Orcu-U3", "Orcu-U4", "Orcu-U5", "Orcu-U6", "Orcu-U7", "Orcu-U8", "Orcu-U9")

# Step 1: Create a TxDb object from your GTF file
gtf_file <- "OryCun2.0_corrected_v12.gtf"  # Adjust path if needed
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

# Extract exons by gene and calculate the total exon length per gene
exons_by_gene <- exonsBy(txdb, by = "gene")
gene_lengths <- sum(width(reduce(exons_by_gene)))

# Convert gene lengths to a data frame
gene_lengths_df <- data.frame(gene_id = names(gene_lengths), length = gene_lengths)

# Filter for MHC-I gene lengths
mhc1_lengths <- gene_lengths_df[gene_lengths_df$gene_id %in% mhc1_gene_ids, ]

# Step 2: Load precomputed counts from your sample files (only *_counts.txt files, exclude summary)
sample_files <- list.files(path = 'counts', pattern = "_counts\\.txt$", full.names = TRUE)

# Filter out any summary files
sample_files <- sample_files[!grepl("summary", sample_files, ignore.case = TRUE)]
cat("Filtered sample files found:\n", sample_files, "\n")

# Check if any sample files were found
if (length(sample_files) == 0) {
  stop("No valid *_counts.txt files were found in the current directory.")
}

# Function to ensure MHC-I gene rows exist in all files
add_missing_genes <- function(data, gene_ids) {
  missing_genes <- setdiff(gene_ids, rownames(data))
  if (length(missing_genes) > 0) {
    # Add missing genes with zero counts
    missing_data <- matrix(0, nrow = length(missing_genes), ncol = ncol(data))
    rownames(missing_data) <- missing_genes
    colnames(missing_data) <- colnames(data)  # Match column names
    data <- rbind(data, missing_data)
  }
  return(data[gene_ids, , drop = FALSE])  # Ensure the same order of rows
}

# Initialize a list to store counts from each file
count_list <- lapply(sample_files, function(file) {
  # Read the data
  data <- read.table(file, header = TRUE, sep = "\t", row.names = 1)
  
  # Extract the column with counts (assuming it's the one with ".bam" in its name)
  count_column <- grep(".bam", colnames(data))
  
  # If no matching column found, return NULL
  if (length(count_column) == 0) {
    cat("Warning: No column with .bam found in file:", file, "\n")
    return(NULL)
  }
  
  # Extract the counts column
  numeric_data <- data[, count_column, drop = FALSE]
  
  # Rename the column to just the sample file name without "_counts.txt"
  sample_name <- gsub("_counts.txt", "", basename(file))
  colnames(numeric_data) <- sample_name
  
  # Filter for MHC-I genes and add missing genes if necessary
  numeric_data <- add_missing_genes(numeric_data, mhc1_gene_ids)
  
  return(numeric_data)
})

# Remove NULL entries from count_list (from files without valid data)
count_list <- Filter(Negate(is.null), count_list)

# Check if any valid count data was found
if (length(count_list) == 0) {
  stop("No valid count data was found.")
}

# Combine the count data into a matrix, where rows are genes and columns are samples
mhc1_counts <- do.call(cbind, count_list)

# Ensure gene lengths are aligned with the counts data
mhc1_counts <- mhc1_counts[rownames(mhc1_counts) %in% mhc1_lengths$gene_id, ]
mhc1_lengths <- mhc1_lengths[mhc1_lengths$gene_id %in% rownames(mhc1_counts), ]

# Sort the counts and lengths to match
mhc1_counts <- mhc1_counts[order(rownames(mhc1_counts)), ]
mhc1_lengths <- mhc1_lengths[order(mhc1_lengths$gene_id), ]

# Step 3: Map SRR IDs to sample names and reorder columns based on the provided mapping
ordered_samples <- sample_mapping$Sample  # Desired order from the CSV
mapped_sample_names <- sapply(colnames(mhc1_counts), function(id) {
  mapped_name <- sample_mapping[sample_mapping$SRR == id, "Sample"]
  if (length(mapped_name) == 0) return(id)  # Return original ID if not found
  return(mapped_name)
})

# Reorder mhc1_counts columns based on the sample_mapping order
mhc1_counts <- mhc1_counts[, match(ordered_samples, mapped_sample_names, nomatch = 0)]
colnames(mhc1_counts) <- ordered_samples  # Assign new column names based on the mapping

# Step 4: Calculate FPKM
# FPKM calculation: FPKM = (counts / (gene length * total reads)) * 1e9
total_counts <- colSums(mhc1_counts)  # Total number of reads in each sample
mhc1_fpkm <- t(t(mhc1_counts) * 1e9 / (mhc1_lengths$length * total_counts))

# Replace NA or infinite values with 0
mhc1_fpkm[is.na(mhc1_fpkm)] <- 0
mhc1_fpkm[is.infinite(mhc1_fpkm)] <- 0

# Filter out genes with all zero FPKM values
mhc1_filtered_fpkm <- mhc1_fpkm[rowSums(mhc1_fpkm) > 0, ]

# Step 5: Generate heatmap
if (nrow(mhc1_filtered_fpkm) >= 2 && ncol(mhc1_filtered_fpkm) >= 2) {
  pheatmap(mhc1_filtered_fpkm, 
           cluster_rows = FALSE,  # No clustering for rows
           cluster_cols = FALSE,  # No clustering for columns
           scale = "none", 
           main = "")
}

# Step 6: Save the FPKM values and heatmap
write.csv(mhc1_filtered_fpkm, "mhc1_fpkm_results_ordered.csv")

########################################################################################
# Function to extract tissue name from the column names
extract_tissue <- function(column_name) {
  return(strsplit(column_name, " ")[[1]][1])
}

# Apply this function to get tissue names
tissue_names <- sapply(colnames(mhc1_filtered_fpkm), extract_tissue)

# Create a new data frame with mean values for each tissue
mhc1_mean_fpkm <- aggregate(t(mhc1_filtered_fpkm), by=list(Tissue=tissue_names), FUN=mean)

# Transpose the data frame back to original orientation
mhc1_mean_fpkm <- t(mhc1_mean_fpkm[,-1])
colnames(mhc1_mean_fpkm) <- unique(tissue_names)



library(ggplot2)
library(ggplotify)

# Create the main heatmap
p <- pheatmap(log2(mhc1_mean_fpkm + 1),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "none",
         main = "",
         fontsize = 8,
         fontsize_row = 8,
         fontsize_col = 8,
         angle_col = 0,
         fontsize_legend = 8,
         legend = TRUE)

heatmap_ggplot <- as.ggplot(p)

ggsave("MHC1_FPKM_heatmap_ordered.svg", p, width = 6, height = 2)

