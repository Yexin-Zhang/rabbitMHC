rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(seqinr)
library(ape)
library(pegas)
library(ggplot2)
library(tidyr)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("Nucleotide_Diversity.Rdata")

# Assign all 9 genes and generate the new nomenclature
target_genes <- c("10972b","11016","10972a","13304","06557","10984","08993b","08993","08993c")
new_names <- paste0("Orcu-U", as.character(1:9))
populations <- c("UK", "AU")

# These are 1-indexed relative to the start of the a1 domain
contact_positions <- c(5, 9, 24, 25, 34, 45, 63, 66, 67, 69, 70, 73, 
                       74, 76, 77, 80, 81, 84, 111, 123, 124, 143, 
                       146, 147, 150, 152, 156, 160, 163, 167)

results_split <- data.frame()

# Loop through populations and genes to isolate sequences
for (pop in populations) {
  for (gene in target_genes) {
    print(gene)
    sequences_matrix <- gene_sequence_matrices[[gene]]
    coordinate <- coordinates[[gene]]
    
    # Filter sequences for the specific population
    is_rabbit_pop <- which(grepl(pop, sequences_matrix[,1]))
    
    if (length(is_rabbit_pop) > 0) {
      rabbit_dna_strings <- Biostrings::DNAStringSet(unlist(sequences_matrix[is_rabbit_pop,2]))
      full_dnabin <- as.matrix(as.DNAbin(rabbit_dna_strings))
      
      # Identify the exact nucleotide boundaries for the PBS (a1 and a2 domains)
      a1_start <- coordinate$Start[coordinate$Domain == "a1"]
      a2_end <- coordinate$End[coordinate$Domain == "a2"]
      print(a1_start)
      print(a2_end)
      pbs_indices <- a1_start:a2_end
      
      # Calculate global nucleotide indices for contact positions based on the a1 start offset
      contact_nuc_indices <- unlist(lapply(contact_positions, function(x) {
        offset <- (x - 1) * 3
        return(c(a1_start + offset, a1_start + offset + 1, a1_start + offset + 2))
      }))
      
      # Define the three distinct regions ensuring indices do not exceed sequence length
      valid_contact_pbs <- contact_nuc_indices[contact_nuc_indices <= ncol(full_dnabin) & contact_nuc_indices %in% pbs_indices]
      valid_other_pbs <- setdiff(pbs_indices, valid_contact_pbs)
      valid_other_pbs <- valid_other_pbs[valid_other_pbs <= ncol(full_dnabin)]
      
      valid_others <- setdiff(1:ncol(full_dnabin), pbs_indices)
      valid_others <- valid_others[valid_others <= ncol(full_dnabin)]
      
      # Extract sequences for each category
      contact_pbs_seqs <- full_dnabin[, valid_contact_pbs, drop = FALSE]
      other_pbs_seqs <- full_dnabin[, valid_other_pbs, drop = FALSE]
      others_seqs <- full_dnabin[, valid_others, drop = FALSE]
      
      # Calculate metrics for Contact in PBS
      contact_k <- seqinr::kaks(ape::as.alignment(contact_pbs_seqs))
      
      # Calculate metrics for Others in PBS
      other_pbs_k <- seqinr::kaks(ape::as.alignment(other_pbs_seqs))
      
      # Calculate metrics for Others
      others_k <- seqinr::kaks(ape::as.alignment(others_seqs))
      
      # Store results for all three categories
      results_split <- rbind(results_split, data.frame(
        Population = pop,
        Gene = gene,
        Category = c("Contact residues in PBS", "Non-contact residues in PBS", "Others"),
        Mean_dN = c(mean(contact_k$ka, na.rm = TRUE), mean(other_pbs_k$ka, na.rm = TRUE), mean(others_k$ka, na.rm = TRUE)),
        Mean_dS = c(mean(contact_k$ks, na.rm = TRUE), mean(other_pbs_k$ks, na.rm = TRUE), mean(others_k$ks, na.rm = TRUE))
      ))
    }
  }
}

print(results_split)
write.csv(results_split, "pN_pS_3_Categories_9genes.csv", row.names = FALSE)

# Reshape the results dataframe into a long format for faceted plotting
results_long <- results_split %>%
  pivot_longer(
    cols = c("Mean_dN", "Mean_dS"),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(Value = ifelse(is.na(Value) | Value < 0, 0, Value))

# Clean up the metric names for plot labels
results_long$Metric <- factor(results_long$Metric, 
                              levels = c("Mean_dN", "Mean_dS"),
                              labels = c("Non-synonymous", "Synonymous"))

# Rename all 9 genes efficiently using factors
results_long$Gene <- factor(results_long$Gene, levels = target_genes, labels = new_names)
results_long$Population <- recode(results_long$Population, "UK" = "United Kingdom", "AU" = "Australia")

# Reorder Category levels to control facet row hierarchy
results_long$Category <- factor(results_long$Category, levels = c("Contact residues in PBS", "Non-contact residues in PBS", "Others"))
# Determine the maximum y value to set identical y-axis limits across all facets
max_y <- max(results_long$Value * 100, na.rm = TRUE)

# Generate the plot
diversity_plot <- ggplot(results_long, aes(x = Gene, y = Value * 100, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, alpha = 0.85) +
  theme_bw() +
  scale_fill_manual(values = c("Non-synonymous" = "#BC340F", "Synonymous" = "#1F77B4")) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 14),
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    axis.text = element_text(color = "sienna", size = 12),
    axis.text.x = element_text(face = "italic", angle = 50, vjust = 1, hjust = 1),
    strip.text = element_text(size = 12),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(1, "lines"),
    strip.background = element_rect(colour = "black", fill = "white", linewidth = 1),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  labs(
    x = "Gene", 
    y = "Genetic Diversity (%)"
  ) +
  scale_y_continuous(limits = c(0, max_y * 1.05), expand = c(0, 0)) +
  facet_grid(Category ~ Population)

print(diversity_plot)

# Save the plot
ggsave("Diversity_3Categories_9genes.svg", plot = diversity_plot, width = 6, height = 9, device = "svg")