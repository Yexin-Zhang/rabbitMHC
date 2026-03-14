rm(list = ls())
library(seqinr)
library(ape)
library(pegas)
library(ggplot2)
library(reshape2)
library(dplyr)
library(Biostrings)
library(readr)
library(tidyr)
library(gridExtra)
library(cowplot)
# library(pdflite)


start_time <- Sys.time()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# # setwd("/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim1/06_Diversity_PiNpiS/")
# 
# # set basic information
# populations <- c("UK", "AU")
# 
# domains <- c("SP","a1","a2","a3","TM+CYT")[1:4]
# 
# genes <- c("10972b","11016","10972a","13304","06557","10984","08993b","08993","08993c")
# 
# new_names <- c(paste0("Orcu-U", as.character(1:9)))
# 
# 
# 
# fasta_files <- paste0("raw/", genes, ".fasta")
# allele_sequences <- lapply(fasta_files, function(f) readDNAStringSet(f))
# names(allele_sequences) <- genes
# 
# 
# # read coordinates of domains for 9 genes
# coordinate_files <- paste0("coorfile/" ,genes, ".txt")  # Adjust the file names as needed
# 
# coordinates <- lapply(coordinate_files, function(f) {
#   lines <- read_delim(f, delim = "\t", col_names = F,show_col_types = FALSE)
#   names(lines) <- c('Gene','Domain','Start','End','nlen','aalen')
#   return(lines)
# })
# 
# names(coordinates) <- genes
# 
# # read colorful spreadsheet on allele numbers
# allele_info <- read.csv("/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim1/04_PacBio2/FINAL_PacBio_alleles_Nov_2022.csv", 
#                         stringsAsFactors = FALSE)
# 
# allele_info_final <- allele_info %>%
#   select(c("Name", (names(allele_info)[sapply(names(allele_info),
#                                               function(x) any(grepl(paste(genes, collapse = "|"), x)))])))
# # 
# # allele_info_final <- allele_info %>% 
# #   select(Name, matches(paste0("Gene\\.(", paste(genes, collapse = "|"), ")")))
# 
# ################## get allele sequences for all samples ################## 
# 
# # Create an empty list to store sequence matrices for each gene
# gene_sequence_matrices <- list()
# 
# for (gene in genes) {
#   
#   # Create an empty list to store sequences for the current gene
#   gene_allele_sequences <- matrix(ncol=2, nrow=2*nrow(allele_info_final))
#   colnames(gene_allele_sequences) <- c("Name","Sequence")
#   
#   # Find relevant column indices for the current gene
#   gene_sequences <- allele_sequences[[gene]]
#   relevant_columns <- grep(gene, names(allele_info_final), value = TRUE)
#   
#   # Initialize a list to store allele name and sequence pairs
#   for (sample in 1:nrow(allele_info_final)) {
#     
#     sample_name <- allele_info_final$Name[sample]
#     # a1 <- allele_info_final[sample, relevant_columns[1]]
#     # a2 <- allele_info_final[sample, relevant_columns[2]]
#     # Check if the entire sample should be skipped
#     # if (isTRUE(allele_info_final[sample, relevant_columns[1]] != '0') && isTRUE(allele_info_final[sample, relevant_columns[2]] != '0')) {
#     if (allele_info_final[sample, relevant_columns[1]] != '0' && allele_info_final[sample, relevant_columns[2]] != '0') {
#       
#       
#       # Get allele indices
#       allele_index1 <- which(names(gene_sequences) == paste0(gene,'.seq', 
#                                                              as.character(allele_info_final[sample,relevant_columns[1]])))
#       allele_index2 <- which(names(gene_sequences) == paste0(gene,'.seq', 
#                                                              as.character(allele_info_final[sample,relevant_columns[2]])))
#       
#       allele1 <- toString(gene_sequences[[allele_index1]])
#       allele2 <- toString(gene_sequences[[allele_index2]])
#       
#       # Store alleles in the list for the current gene
#       gene_allele_sequences[sample,2] <- allele1
#       gene_allele_sequences[sample,1] <- paste0(gene, '.', sample_name, '.allele1')
#       
#       
#       gene_allele_sequences[nrow(allele_info_final)+sample,2] <- allele2
#       gene_allele_sequences[nrow(allele_info_final)+sample,1] <- paste0(gene, '.', sample_name, '.allele2')
#       
#       
#     } else {
#       gene_allele_sequences[sample,] <- NA
#       gene_allele_sequences[nrow(allele_info_final)+sample,] <- NA
#     }
#     
#   }
#   
#   # Remove the rows with no allele recognition - basically the '0's in the colorful spreadsheet
#   gene_allele_sequences <- gene_allele_sequences[!is.na(gene_allele_sequences[, 1]), ]
#   
#   # Store the gene_sequence_matrix in the gene_sequence_matrices list
#   gene_sequence_matrices[[gene]] <- gene_allele_sequences
# }
# 
# 
# ################## Initialize data frames to store results ################## 
# 
# 
# all_genes_results_df1 <- data.frame(Gene = character(), Population = character(), Domain = character(),
#                                     domain_pi_diversity = numeric(), domain_pi_diversity_CI1 = numeric(), omain_pi_diversity_CI2 = numeric(),
#                                     domain_non = numeric(), domain_non_CI1 = numeric(), domain_non_CI2 = numeric(),
#                                     domain_syn = numeric(), domain_syn_CI1 = numeric(), domain_syn_CI2 = numeric(),
#                                     stringsAsFactors = TRUE)
# 
# all_genes_results_df2 <- data.frame(Gene = character(), Population = character(),
#                                     pi_diversity = numeric(), pi_diversity_CI1 = numeric(), pi_diversity_CI2 = numeric(),
#                                     theta = numeric(), theta_CI1 = numeric(), theta_CI2 = numeric(),
#                                     hap_div = numeric(), hap_div_CI1 = numeric(), hap_div_CI2 = numeric(),
#                                     tajima = numeric(), tajima_CI1 = numeric(), tajima_CI2 = numeric(),
#                                     stringsAsFactors = TRUE)
# 
# 
# ################## generate metrics for each domain ################## 
# 
# for(i in 1:length(populations)){
#   # Loop through each gene and calculate genetic diversity for rabbits with names starting with "M_UK"
#   for (gene in genes) {
#     
#     # retrieve the allele matrix and coordinates for this gene
#     sequences_matrix <- gene_sequence_matrices[[gene]]
#     coordinate <- coordinates[[gene]]
#     
#     # Separate rabbits based on sample names. UK or AU
#     is_rabbit_pop <- which(grepl(populations[i], sequences_matrix[,1]))
#     
#     # Check if samples for this gene in this populations are found
#     if (any(is_rabbit_pop)) {
#       
#       # Filter sequences for specific population
#       rabbit_dna_strings <- DNAStringSet(unlist(sequences_matrix[is_rabbit_pop,2]))
#       rabbit_sequences <- as.matrix(as.DNAbin(rabbit_dna_strings))
#       
#       # Calculate Nucleotide Diversity, Haplotype Diversity, ThetaS and TajimaD for all sites
#       
#       pi_diversity <- nuc.div(rabbit_sequences)
#       
#       s <- length(seg.sites(rabbit_sequences))
#       n <- nrow(rabbit_sequences)
#       theta <- theta.s(s, n)/(length(rabbit_sequences)/n)
#       
#       hap_div <- hap.div(rabbit_sequences)
#       
#       tajima <- tajima.test(rabbit_sequences)$D
# 
#       ################## bootstrap ##################
#       
#       # Number of bootstrap iterations
#       num_bootstraps <- 1000
#       # Initialize a list to store the bootstrap samples
#       bootstrap_samples_strings <- vector("list", num_bootstraps)
#       bootstrap_samples_dnabin <- vector("list", num_bootstraps)
#       
#       for (bb in 1:num_bootstraps) {
#         # Sampling indices
#         sampled_dna_strings <- sample(rabbit_dna_strings, length(rabbit_dna_strings), replace = TRUE)
#         sampled_dnabin <- as.matrix(as.DNAbin(sampled_dna_strings))
# 
#         # Subset and store the sampled sequences
#         bootstrap_samples_strings[[bb]] <- sampled_dna_strings
#         bootstrap_samples_dnabin[[bb]] <- sampled_dnabin
#       }
# 
# 
#       bootstrap_pi_diversity <- sapply(bootstrap_samples_dnabin, nuc.div)
#       pi_diversity_CI1 <- as.numeric(quantile(bootstrap_pi_diversity,probs=0.025))
#       pi_diversity_CI2 <- as.numeric(quantile(bootstrap_pi_diversity,probs=0.975))
# 
#       bootstrap_theta <- sapply(bootstrap_samples_dnabin, function(f){
#         s <- length(seg.sites(f))
#         n <- nrow(f)
#         thetas <- theta.s(s, n)/(length(f)/n)
#         return(thetas)
#       })
# 
#       theta_CI1 <- as.numeric(quantile(bootstrap_theta,probs=0.025))
#       theta_CI2 <- as.numeric(quantile(bootstrap_theta,probs=0.975))
# 
#       bootstrap_hap_div <- sapply(bootstrap_samples_dnabin, hap.div)
#       hap_div_CI1 <- as.numeric(quantile(bootstrap_hap_div,probs=0.025))
#       hap_div_CI2 <- as.numeric(quantile(bootstrap_hap_div,probs=0.975))
# 
#       bootstrap_tajima <- sapply(bootstrap_samples_dnabin, function(f){
#         tajimad <- tajima.test(f)$D
#         return(tajimad)
#       })
# 
#       tajima_CI1 <- as.numeric(quantile(bootstrap_tajima,probs=0.025,na.rm = TRUE))
#       tajima_CI2 <- as.numeric(quantile(bootstrap_tajima,probs=0.975,na.rm = TRUE))
# 
# 
#       all_genes_results_df2 <- rbind(all_genes_results_df2,
#                                      data.frame(Gene = gene, Population = populations[i],
#                                                 pi_diversity = pi_diversity, pi_diversity_CI1 = pi_diversity_CI1, pi_diversity_CI2 = pi_diversity_CI2,
#                                                 theta = theta, theta_CI1 = theta_CI1, theta_CI2 = theta_CI2,
#                                                 hap_div = hap_div, hap_div_CI1 = hap_div_CI1, hap_div_CI2 = hap_div_CI2,
#                                                 tajima = tajima, tajima_CI1 = tajima_CI1, tajima_CI2 = tajima_CI2,
#                                                 stringsAsFactors = TRUE))
# 
# 
#       for (d in 1:length(domains)){
# 
# 
#         domain_name <- domains[d]
#         
#         domain_seq <- as.matrix(as.DNAbin(subseq(rabbit_dna_strings,
#                                                  start = coordinate$Start[d] , end = coordinate$End[d])))
# 
#         domain_pi_diversity <- nuc.div(domain_seq)
# 
#         k <- seqinr::kaks(ape::as.alignment(domain_seq))
#         domain_non <- mean(k$ka,na.rm = TRUE)
#         domain_syn <- mean(k$ks,na.rm = TRUE)
# 
# 
#         bootstrap_samples <- lapply(bootstrap_samples_strings, function(f) {
# 
#           domain_sequences <- subseq(f, start = coordinate$Start[d], end = coordinate$End[d])
# 
#           domain_dnabin <- as.matrix(as.DNAbin(domain_sequences))
# 
#           return(domain_dnabin)
#         })
# 
# 
#         bootstrap_domain_pi_diversity <- sapply(bootstrap_samples, nuc.div)
#         domain_pi_diversity_CI1 <- as.numeric(quantile(bootstrap_domain_pi_diversity,probs=0.025,na.rm = TRUE))
#         domain_pi_diversity_CI2 <- as.numeric(quantile(bootstrap_domain_pi_diversity,probs=0.975,na.rm = TRUE))
# 
#         bootstrap_domain_non <- sapply(bootstrap_samples,
#                                        function(f) mean(seqinr::kaks((ape::as.alignment(f)))$ka))
#         
#         domain_non_CI1 <- as.numeric(quantile(bootstrap_domain_non,probs=0.025))
#         domain_non_CI2 <- as.numeric(quantile(bootstrap_domain_non,probs=0.975))
# 
#         bootstrap_domain_syn <- sapply(bootstrap_samples,
#                                        function(f) mean(seqinr::kaks((ape::as.alignment(f)))$ks))
#         domain_syn_CI1 <- as.numeric(quantile(bootstrap_domain_syn,probs=0.025))
#         domain_syn_CI2 <- as.numeric(quantile(bootstrap_domain_syn,probs=0.975))
# 
# 
#         all_genes_results_df1 <- rbind(all_genes_results_df1,
#                                        data.frame(Gene = gene, Population = populations[i], Domain = domain_name,
#                                                   domain_pi_diversity = domain_pi_diversity,
#                                                   domain_pi_diversity_CI1 = domain_pi_diversity_CI1, domain_pi_diversity_CI2 = domain_pi_diversity_CI2,
#                                                   domain_non = domain_non, domain_non_CI1 = domain_non_CI1, domain_non_CI2 = domain_non_CI2,
#                                                   domain_syn = domain_syn, domain_syn_CI1 = domain_syn_CI1, domain_syn_CI2 = domain_syn_CI2,
#                                                   stringsAsFactors = TRUE))
#         
#       }
#       
#     } else
#     {
#       cat("No samples found in", gene,populations[i], "\n")
#     }
#   }
# }
# 
# end_time = Sys.time()
# 
# all_genes_results_df2_converted <- reshape2::melt(all_genes_results_df2,id.vars=c("Gene","Population"))
# 
# # Transform the data
# all_genes_results_df2_converted <- all_genes_results_df2_converted %>%
#   # Extract the measurement type and CI type into separate columns
#   mutate(
#     Category = gsub("(_CI1|_CI2)$", "", variable),
#     CI_Type = ifelse(grepl("_CI", variable), gsub(".*_", "", variable), "Value")
#   ) %>%
#   # Spread the CI_Type into separate columns
#   pivot_wider(
#     id_cols = c(Gene, Population, Category),
#     names_from = CI_Type,
#     values_from = value,
#     values_fill = list(value = NA)  # Fill missing values with NA
#   ) 
# 
# 
# 
# all_genes_results_df1$Gene <- factor(all_genes_results_df1$Gene, levels = genes)
# levels(all_genes_results_df1$Gene) <- new_names
# levels(all_genes_results_df1$Population) <- c("United Kingdom", "Australia")
# 
# all_genes_results_df2_converted$Gene <- factor(all_genes_results_df2_converted$Gene, levels = genes)
# levels(all_genes_results_df2_converted$Gene) <- new_names
# levels(all_genes_results_df2_converted$Population) <- c("United Kingdom", "Australia")
# all_genes_results_df2_converted$Category <- as.factor(all_genes_results_df2_converted$Category)
# levels(all_genes_results_df2_converted$Category) <- c("Haplotype Diversity", "Nucleotide Diversity",
#                                                       "Tajiama's D", "Watterson's Theta")
# 
# 
# # save.image(file = "Nucleotide_Diversity.Rdata")
load("Nucleotide_Diversity.Rdata")

avg_line_data_UK = 0.2961867
avg_line_data_AU = 0.2697530


###############################################################################################################
########################################## Plot Nucleotide Diversity  #########################################
###############################################################################################################
nu_diversity_plot <- ggplot(all_genes_results_df1, aes(x = factor(Domain, level=c("SP","a1","a2","a3","TM+CYT")), y = domain_pi_diversity*100)) +
  geom_bar(fill = "black", stat="identity", position = "dodge", alpha=0.6,show.legend=F) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size = 14)) +
  theme(axis.title.x = element_text(color = "black", size = 14),
        axis.title.y = element_text(color = "black", size = 14),
        axis.text = element_text(color = "sienna", size = 12),
        axis.text.x = element_text(face = "italic",angle = 50, vjust = 1, hjust = 1),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(),
        panel.spacing.y = unit(1, "lines")) +
  labs(x="Domain",y="Nucleotide Diversity (%)") +
  geom_errorbar(aes(ymin=domain_pi_diversity_CI1*100, ymax=domain_pi_diversity_CI2*100),color="brown", width=0.3, alpha=0.9) +
  facet_grid(cols = vars(Gene),rows = vars(Population)) +
  theme(strip.background = element_rect(colour="black", fill="white", linewidth = 1), legend.position = "none") +
  # scale_x_discrete(labels = c("SP" = "SP",
  #                             "a1" = bquote(alpha[1]),
  #                             "a2" = bquote(alpha[2]),
  #                             "a3" = bquote(alpha[3]),
  #                             "TM+CYT" = "TM+CYT")) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10), expand = c(0, 0))



nu_diversity_plot

# ggsave("Nucleotide_diversity_domain.pdf", plot = nu_diversity_plot, height = 4, width = 8, device = "pdf")
ggsave("Nucleotide_diversity_domain.pdf", nu_diversity_plot, width = 8, height = 6, device = "pdf")


all_genes_results_df1_kaks <- subset(pivot_longer(all_genes_results_df1, cols = c("domain_non", "domain_syn"),
                                                  names_to = "Category", values_to = "value_kaks"),
                                     select = c("Gene","Population","Domain","Category","value_kaks"))

all_genes_results_df1_CI1 <- subset(pivot_longer(all_genes_results_df1, cols = c("domain_non_CI1", "domain_syn_CI1"),
                                                 names_to = "Category", values_to = "value_CI1"),
                                    select = "value_CI1")

all_genes_results_df1_CI2 <- subset(pivot_longer(all_genes_results_df1, cols = c("domain_non_CI2", "domain_syn_CI2"),
                                                 names_to = "Category", values_to = "value_CI2"),
                                    select = "value_CI2")

all_genes_results_df1_converted <- cbind(all_genes_results_df1_kaks, all_genes_results_df1_CI1, all_genes_results_df1_CI2)

all_genes_results_df1_converted$Category <- as.factor(all_genes_results_df1_converted$Category)

levels(all_genes_results_df1_converted$Category) <- c("Non-synonymous", "Synonymous")

all_genes_results_df1_converted <- all_genes_results_df1_converted %>%
  mutate(value_kaks = ifelse(value_kaks < 0, 0, value_kaks)) %>%
  mutate(value_CI1 = ifelse(value_kaks <0, 0, value_CI1)) %>%
  mutate(value_CI2 = ifelse(value_kaks <0, 0, value_CI2))


################################################################################################################
########################################### Plot Genetic Diversity  ############################################
################################################################################################################
kaks_plot <- ggplot(all_genes_results_df1_converted, aes(x = factor(Domain, level=c("SP","a1","a2","a3","TM+CYT")))) +
  geom_bar(aes(y = value_kaks * 100, fill = Category), stat = "identity", position = position_dodge(width = 0.8), alpha = 0.85) +
  theme_bw() +
  scale_fill_manual(values = c("#BC340F","#1F77B4")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14)) +
  theme(axis.title.x = element_text(color = "black", size = 14),
        axis.title.y = element_text(color = "black", size = 14),
        axis.text = element_text(color = "sienna", size = 12),
        axis.text.x = element_text(face = "italic", angle = 50, vjust = 1, hjust = 1),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(),
        panel.spacing.y = unit(1, "lines"),
        legend.title = element_blank()) +
  labs(x = "Domain", y = "Nucleotide Diversity (%)") +
  geom_errorbar(aes(ymin=value_CI1*100, ymax=value_CI2*100, color = Category), width=0.2, alpha=0.7, position = position_dodge(width = 0.8)) +
  facet_grid(cols = vars(Gene), rows = vars(Population)) +
  theme(strip.background = element_rect(colour = "black", fill = "white", linewidth = 1), legend.position = "bottom") +
  # scale_x_discrete(labels = c("SP" = "SP",
  #                             "a1" = bquote(alpha[1]),
  #                             "a2" = bquote(alpha[2]),
  #                             "a3" = bquote(alpha[3]),
  #                             "TM+CYT" = "TM+CYT")) +
  scale_y_continuous(breaks = seq(0, 12, by = 4), limits = c(0, 12), expand = c(0, 0))

kaks_plot

ggsave("Syn_nonsyn_domain.pdf", kaks_plot, width = 8, height = 6, device = "pdf")

################################################################################################################
################################### Haplotype Nucelotide Diversity #############################################
################################################################################################################
hap_div <- all_genes_results_df2_converted %>%
  dplyr::filter(Category == 'Nucleotide Diversity')

allsites_hap_div_plot <- ggplot(hap_div, aes(x = factor(Gene))) +
  geom_bar(aes(y = Value * 100), fill = "black", stat = "identity",
           position = "dodge", alpha = 0.6, show.legend = FALSE) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14),
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    axis.text = element_text(color = "sienna", size = 12),
    axis.text.x = element_text(face = "italic", angle = 50, vjust = 1, hjust = 1),
    strip.text = element_text(size = 11),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.spacing.y = unit(1, "lines"),  # Adjust spacing between facet rows
    strip.background = element_rect(colour = "black", fill = "white", linewidth = 1),
    legend.position = "none"
  ) +
  labs(x = "Gene", y = "Nucleotide Diversity (%)") +
  geom_errorbar(aes(ymin = CI1 * 100, ymax = CI2 * 100), color = "brown", width = 0.3, alpha = 0.9) +
  facet_wrap(~ Population, nrow = 2) +  # Facet by population into two rows
  theme(strip.background = element_rect(colour = "black", fill = "white", linewidth = 1)) +
  scale_y_continuous(breaks = seq(0, 5, by = 1), limits = c(0, 5), expand = c(0, 0))

allsites_hap_div_plot


ggsave("allsites_hap_div_plot.pdf", allsites_hap_div_plot, width = 4, height = 8, device = "pdf")

# Create a grid layout with allsites_hap_div_plot on the left, taking 2 rows, and the other two plots on the right
left_side <- plot_grid(allsites_hap_div_plot, ncol = 1)
right_side <- plot_grid(nu_diversity_plot, kaks_plot, ncol = 1, labels = c("B", "C"))

# Combine left and right sides with specified relative widths
combined_plot <- plot_grid(left_side, right_side, labels = c("A", ""), ncol = 2, rel_widths = c(1/3, 2/3))
combined_plot

# Save the combined plot
ggsave("Figure3_Genetic_Diversity.pdf", combined_plot, width = 12, height = 10, device = "pdf")



################################################################################################################
################################### Haplotype metrics (Supplementary) ##########################################
################################################################################################################
allsites_plot <- ggplot(data = filter(all_genes_results_df2_converted, Category != 'Nucleotide Diversity'), aes(x = factor(Gene))) +
  # geom_bar(aes(y = Value,fill = Category), stat = "identity", position = position_dodge(width = 0.8), alpha = 0.85) +
  geom_bar(aes(y = Value), fill = "black", stat="identity", position = "dodge", alpha=0.6, show.legend=F) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14),
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    axis.text = element_text(color = "sienna", size = 12),
    axis.text.x = element_text(face = "italic", angle = 50, vjust = 1, hjust = 1),
    strip.text = element_text(size = 11),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.spacing.y = unit(1, "lines"),  # Adjust spacing between facet rows
    strip.background = element_rect(colour = "black", fill = "white", linewidth = 1),
    legend.position = "none"
  ) +
  labs(x = "Gene") +
  # geom_errorbar(aes(ymin=CI1, ymax=CI2, color = Category), width=0.2, alpha=0.7, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin=CI1, ymax=CI2),color="brown", width=0.3, alpha=0.9) +
  facet_grid(cols = vars(Population), rows = vars(Category), scales = "free_y") +
  # scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme(strip.background = element_rect(colour = "black", fill = "white", linewidth = 1), legend.position = "none")
  # + scale_fill_manual(values = fill_colors) +  # Apply custom fill colors
  # scale_color_manual(values = line_colors)

allsites_plot

# ggsave("General_diversity_allsites_Supp.pdf", plot = allsites_plot, height = 8, width = 6, device = "pdf")
ggsave("General_diversity_allsites_Supp.pdf", allsites_plot, width = 6, height = 8, device = "pdf")




################################################################################################################
################################### NEW FIGURE3 B & SUPP FIGURE  #############################################
################################################################################################################
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare nucleotide diversity data
nu_diversity_data_AU <- all_genes_results_df1 %>%
  filter(Population == "Australia") %>%
  dplyr::select(Gene, Domain, domain_pi_diversity, domain_pi_diversity_CI1, domain_pi_diversity_CI2) %>%
  mutate(Category = "Nucleotide Diversity",
         SubCategory = "Nucleotide Diversity")

# Prepare synonymous/non-synonymous diversity data
kaks_data_AU <- all_genes_results_df1_converted %>%
  dplyr::filter(Population == "Australia") %>%
  dplyr::rename(domain_pi_diversity = value_kaks,
         domain_pi_diversity_CI1 = value_CI1,
         domain_pi_diversity_CI2 = value_CI2) %>%
  dplyr::select(Gene, Domain, Category, domain_pi_diversity, domain_pi_diversity_CI1, domain_pi_diversity_CI2) %>%
  mutate(SubCategory = "Synonymous/Non-synonymous")

# Combine the datasets
combined_data_AU <- bind_rows(nu_diversity_data_AU, kaks_data_AU)


hap_div_AU <- all_genes_results_df2_converted %>%
  dplyr::filter(Category == 'Nucleotide Diversity', Population == 'Australia')


#######################################################################
# Prepare nucleotide diversity data
nu_diversity_data_UK <- all_genes_results_df1 %>%
  filter(Population == "United Kingdom") %>%
  dplyr::select(Gene, Domain, domain_pi_diversity, domain_pi_diversity_CI1, domain_pi_diversity_CI2) %>%
  mutate(Category = "Nucleotide Diversity",
         SubCategory = "Nucleotide Diversity")

# Prepare synonymous/non-synonymous diversity data
kaks_data_UK <- all_genes_results_df1_converted %>%
  dplyr::filter(Population == "United Kingdom") %>%
  dplyr::rename(domain_pi_diversity = value_kaks,
         domain_pi_diversity_CI1 = value_CI1,
         domain_pi_diversity_CI2 = value_CI2) %>%
  dplyr::select(Gene, Domain, Category, domain_pi_diversity, domain_pi_diversity_CI1, domain_pi_diversity_CI2) %>%
  mutate(SubCategory = "Synonymous/Non-synonymous")

# Combine the datasets
combined_data_UK <- bind_rows(nu_diversity_data_UK, kaks_data_UK)


hap_div_UK <- all_genes_results_df2_converted %>%
  dplyr::filter(Category == 'Nucleotide Diversity', Population == 'United Kingdom')


#######################################################################
#######################################################################
#######################################################################

# Create the plot
# ADDED: Create a helper data frame for the average line
avg_line_df_AU <- data.frame(
  avg_value = avg_line_data_AU,
  SubCategory = "Nucleotide Diversity" 
)

combined_plot_AU <- ggplot(combined_data_AU, aes(x = factor(Domain, level=c("SP","a1","a2","a3","TM+CYT")))) +
  geom_bar(data = subset(combined_data_AU, Category == "Nucleotide Diversity"),
           aes(y = domain_pi_diversity * 100), 
           fill = "black", stat = "identity", alpha = 0.6) + 
  geom_bar(data = subset(combined_data_AU, Category != "Nucleotide Diversity"),
           aes(y = domain_pi_diversity * 100, fill = Category), 
           stat = "identity", position = position_dodge(width = 0.8), alpha = 0.85) + 
  geom_hline(data = avg_line_df_AU, aes(yintercept = avg_value), 
             linetype = "dashed", color = "black", size = 0.4) +
  # Add error bars with appropriate colors
  geom_errorbar(data = subset(combined_data_AU, Category == "Nucleotide Diversity"),
                aes(ymin = domain_pi_diversity_CI1 * 100, 
                    ymax = domain_pi_diversity_CI2 * 100),
                color = "brown", width = 0.2, alpha = 0.7) +
  geom_errorbar(data = subset(combined_data_AU, Category != "Nucleotide Diversity"),
                aes(ymin = domain_pi_diversity_CI1 * 100, 
                    ymax = domain_pi_diversity_CI2 * 100,
                    group = Category,
                    color = Category),
                width = 0.2, alpha = 0.7, 
                position = position_dodge(width = 0.8)) +
  facet_grid(SubCategory ~ Gene, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = c("Non-synonymous" = "#BC340F", 
                               "Synonymous" = "#1F77B4")) +
  scale_x_discrete(labels = c("SP" = "SP", 
                              "a1" = expression(alpha[1]), 
                              "a2" = expression(alpha[2]), 
                              "a3" = expression(alpha[3]), 
                              "TM+CYT" = "TM+CYT")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14),
        axis.title.x = element_text(color = "black", size = 14),
        axis.title.y = element_text(color = "black", size = 14),
        axis.text = element_text(color = "sienna", size = 12),
        axis.text.x = element_text(face = "italic", angle = 50, vjust = 1, hjust = 1),
        strip.text.x = element_text(size = 12, face = "bold.italic"),
        strip.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.y = unit(0.5, "lines"),
        panel.spacing.x = unit(0.8, "lines"),
        # panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        legend.title = element_blank(),
        legend.position = c(0.99, 0.4),
        legend.justification = c(1, 1),
        legend.box.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt")) +
  labs(x = "Domain", y = "Nucleotide Diversity (%)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(strip.background = element_rect(colour = "black", fill = "white", linewidth = 0))


# Display the plot
print(combined_plot_AU)

# Create the plot
# ADDED: Create a helper data frame for the average line
avg_line_df_UK <- data.frame(
  avg_value = avg_line_data_UK,
  SubCategory = "Nucleotide Diversity"
)

combined_plot_UK <- ggplot(combined_data_UK, aes(x = factor(Domain, level=c("SP","a1","a2","a3","TM+CYT")))) +
  geom_bar(data = subset(combined_data_UK, Category == "Nucleotide Diversity"),
           aes(y = domain_pi_diversity * 100), 
           fill = "black", stat = "identity", alpha = 0.6) + 
  geom_bar(data = subset(combined_data_UK, Category != "Nucleotide Diversity"),
           aes(y = domain_pi_diversity * 100, fill = Category), 
           stat = "identity", position = position_dodge(width = 0.8), alpha = 0.85) + 
  geom_hline(data = avg_line_df_UK, aes(yintercept = avg_value), 
             linetype = "dashed", color = "black", size = 0.4) +
  # Add error bars with appropriate colors
  geom_errorbar(data = subset(combined_data_UK, Category == "Nucleotide Diversity"),
                aes(ymin = domain_pi_diversity_CI1 * 100, 
                    ymax = domain_pi_diversity_CI2 * 100),
                color = "brown", width = 0.2, alpha = 0.7) +
  geom_errorbar(data = subset(combined_data_UK, Category != "Nucleotide Diversity"),
                aes(ymin = domain_pi_diversity_CI1 * 100, 
                    ymax = domain_pi_diversity_CI2 * 100,
                    group = Category,
                    color = Category),
                width = 0.2, alpha = 0.6, 
                position = position_dodge(width = 0.8)) +
  facet_grid(SubCategory ~ Gene, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = c("Non-synonymous" = "#BC340F", 
                               "Synonymous" = "#1F77B4")) +
  scale_x_discrete(labels = c("SP" = "SP", 
                              "a1" = expression(alpha[1]), 
                              "a2" = expression(alpha[2]), 
                              "a3" = expression(alpha[3]), 
                              "TM+CYT" = "TM+CYT")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14),
        axis.title.x = element_text(color = "black", size = 14),
        axis.title.y = element_text(color = "black", size = 14),
        axis.text = element_text(color = "sienna", size = 12),
        axis.text.x = element_text(face = "italic", angle = 50, vjust = 1, hjust = 1),
        strip.text.x = element_text(size = 12, face = "bold.italic"),
        strip.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.y = unit(0.5, "lines"),
        panel.spacing.x = unit(0.8, "lines"),
        # panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        legend.title = element_blank(),
        legend.position = c(0.99, 0.4),
        legend.justification = c(1, 1),
        legend.box.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt")) +
  labs(x = "Domain", y = "Nucleotide Diversity (%)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(strip.background = element_rect(colour = "black", fill = "white", linewidth = 0))


# Display the plot
print(combined_plot_UK)

# Save the plot
ggsave("Combined_diversity_plot_AU.pdf", combined_plot_AU, width = 10, height = 6, device = "pdf")
ggsave("Combined_diversity_plot_UK.pdf", combined_plot_UK, width = 10, height = 6, device = "pdf")



# Create Australia plot
allsites_hap_div_plot_AU <- ggplot(hap_div_AU, aes(x = factor(Gene))) +
  geom_bar(aes(y = Value * 100), fill = "black", stat = "identity",
           position = "dodge", alpha = 0.6, show.legend = FALSE) +
  geom_hline(data = avg_line_df_AU, aes(yintercept = avg_value), 
             linetype = "dashed", color = "black", size = 0.4) +
  theme_bw() +
  theme(
    # plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14),
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    axis.text = element_text(color = "sienna", size = 12),
    axis.text.x = element_text(face = "italic", angle = 50, vjust = 1, hjust = 1),
    strip.text = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(1, "lines"),
    strip.background = element_rect(colour = "black", fill = "white", linewidth = 1),
    legend.position = "none"
  ) +
  labs(x = "Gene", y = "Nucleotide Diversity (%)") +
  geom_errorbar(aes(ymin = CI1 * 100, ymax = CI2 * 100), color = "brown", width = 0.3, alpha = 0.9) +
  scale_y_continuous(breaks = seq(0, 5, by = 1), limits = c(0, 5), expand = c(0, 0))


allsites_hap_div_plot_UK <- ggplot(hap_div_UK, aes(x = factor(Gene))) +
  geom_bar(aes(y = Value * 100), fill = "black", stat = "identity",
           position = "dodge", alpha = 0.6, show.legend = FALSE) +
  geom_hline(data = avg_line_df_UK, aes(yintercept = avg_value), 
             linetype = "dashed", color = "black", size = 0.4) +
  theme_bw() +
  theme(
    # plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14),
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    axis.text = element_text(color = "sienna", size = 12),
    axis.text.x = element_text(face = "italic", angle = 50, vjust = 1, hjust = 1),
    strip.text = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(1, "lines"),
    strip.background = element_rect(colour = "black", fill = "white", linewidth = 1),
    legend.position = "none"
  ) +
  labs(x = "Gene", y = "Nucleotide Diversity (%)") +
  geom_errorbar(aes(ymin = CI1 * 100, ymax = CI2 * 100), color = "brown", width = 0.3, alpha = 0.9) +
  scale_y_continuous(breaks = seq(0, 5, by = 1), limits = c(0, 5), expand = c(0, 0))

allsites_hap_div_plot_AU
allsites_hap_div_plot_UK
################################################################################################################
################################### NEW TSE ####################################################################
################################################################################################################
# Load necessary libraries
library(GenomicFeatures)
library(pheatmap)
library(edgeR)
library(ggplot2)
library(ggplotify)

# Clear environment and set working directory
setwd("/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim1/07_TSE")

# Load sample name mapping file
sample_mapping <- read.csv("RNAseq_sample_mapping.csv", header = TRUE)

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
  mapped_name <- sample_mapping[sample_mapping$Sample.Code == id, "Sample"]
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

# Step 5: Save the FPKM values and heatmap
write.csv(mhc1_filtered_fpkm, "mhc1_fpkm_results_ordered.csv")



#calculate the mean FPKM of Orcu-U1 and U2
#Orcu-U1
mean(mhc1_filtered_fpkm[1,])

#Orcu-U2:
mean(mhc1_filtered_fpkm[2,])

#Ratio
mean(mhc1_filtered_fpkm[1,])/mean(mhc1_filtered_fpkm[2,])






# Step 6: plot
# Get column names (tissue types)
col_names <- colnames(mhc1_filtered_fpkm)

# Create a function to find middle position for duplicated names
get_label_positions <- function(names) {
  positions <- list()
  for (name in unique(names)) {
    pos <- which(names == name)
    if (length(pos) > 1) {
      mid_pos <- mean(pos)  # Use 'round' to ensure exact index
      # Add 5 spaces before the name if there are exactly 2 positions
      if(length(pos) == 2) {
        positions[[paste0("    ", name)]] <- mid_pos  # Add 5 spaces
      } else {
        positions[[name]] <- mid_pos
      }
    } else {
      positions[[name]] <- pos
    }
  }
  return(positions)
}

# Get positions for labels
label_positions <- get_label_positions(col_names)

# Create a new labels vector
new_labels <- rep("", length(col_names))
for (name in names(label_positions)) {
  new_labels[round(label_positions[[name]])] <- name
}

# Convert new_labels to expressions with bold formatting
new_labels_bold <- lapply(new_labels, function(x) {
  if(x != "") {
    parse(text = paste0("bold('", x, "')"))[[1]]
  } else {
    parse(text = "''")[[1]]
  }
})


# Find positions where gaps should be inserted between different tissue types
tissue_changes <- which(col_names[-1] != col_names[-length(col_names)])


# Create row names with bold italic formatting
newnames.row <- lapply(
  rownames(mhc1_filtered_fpkm),
  function(x) bquote(italic(.(x))))



# Create the heatmap
p <- pheatmap(log2(mhc1_filtered_fpkm + 1),
              show_rownames = TRUE,
              show_colnames = TRUE,
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              scale = "none",
              main = "",
              fontsize = 10,
              angle_col = 45,
              fontsize_legend = 8,
              legend = TRUE,
              labels_col = as.expression(new_labels_bold),
              labels_row = as.expression(newnames.row),
              gaps_col = tissue_changes)

heatmap_ggplot <- as.ggplot(p)
tse <- heatmap_ggplot + 
  annotate("text", 
           x = Inf, 
           y = Inf, 
           parse = TRUE,    # Add this to parse the expression
           label = "log[2](FPKM+1)",  # Remove expression() and use string
           hjust = 1,    
           vjust = 1.1,    
           size = 3.4,     
           fontface = 'bold') +
  theme(plot.margin = margin(t = 0, r = 20, b = 0, l = 40))  # Adds larger left margin


tse

# Save the plot as an pdf file
ggsave("MHC1_FPKM_heatmap_ordered.pdf", plot = tse, width = 6, height = 2)




################################################################################################################
############################################ NEW FIGURE MAIN#######################################################
################################################################################################################



setwd(dirname(rstudioapi::getSourceEditorContext()$path))





bottom_side <- plot_grid(combined_plot_UK)
up_side <- plot_grid(tse, allsites_hap_div_plot_UK, ncol = 2, rel_widths = c(3/4, 1/4), labels = c("A", "B"))

figure_in_the_main <- plot_grid(up_side, bottom_side, labels = c("", "C"), ncol = 1, rel_heights = c(2/5, 3/5))
figure_in_the_main



# Save the combined plot
ggsave("Figure3_Genetic_Diversity_TSE.pdf", figure_in_the_main, width = 12, height = 10, device = "pdf")





################################################################################################################
############################################ NEW FIGURE SUPP #######################################################
################################################################################################################

au_diversity_supp_plot <- plot_grid(allsites_hap_div_plot_AU, combined_plot_AU, labels = c("a", "b"), ncol = 2, rel_widths = c(1/4, 3/4))


# Save the combined plot
ggsave("Supp_AU_Genetic_Diversity.pdf", au_diversity_supp_plot, width = 12, height = 6, device = "pdf")






