
library(ggplot2)
library(pould)
library(tidyr)
library(dplyr)
library(genetics)
library(RColorBrewer)
library(scales)
library(base)
library(patchwork)


rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
path <- dirname(rstudioapi::getSourceEditorContext()$path)


final_allele_table <- read.csv('/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim2/04_museum_typing/mm97_b0010_combined_paired/AF_common_alleles/final_allele_table.csv') %>%
  rowwise() %>%
  # filter(!any(c_across(starts_with("X")) %in% c('*', '-'))) %>%
  mutate(across(starts_with("X"), ~ if_else(. %in% c('*', '-'), NA_character_, .))) %>%
  ungroup()
final_allele_table <- as.data.frame(final_allele_table)

# colnames(final_allele_table) <- c('Sample','Orcu-U2', 'Orcu-U2', 'Orcu-U1', 'Orcu-U1')

final_allele_table$Population <- ifelse(grepl("AU", final_allele_table$Sample), "AU",
                                        ifelse(grepl("FR", final_allele_table$Sample), "FR",
                                               ifelse(grepl("UK", final_allele_table$Sample), "UK", NA)))

new_nomenclature_file <- read.csv("/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim1/01_New_nomenclature/nomenclature_names.csv")
name_lookup <- setNames(new_nomenclature_file$ShortName, new_nomenclature_file$OriginalName)



populations <- c('FR', 'UK', 'AU')

# cutoff <- 0.1
cutoff_list <- c(0.1, 0.05, 0.03, 0.01)
################################################################################################

for (cutoff in cutoff_list){
  

plot_ld_heatmap <- function(population) {
  allele_table <- filter(final_allele_table, Population == population)[2:5]
  colnames(allele_table) <- c('Orcu.U2.1', 'Orcu.U2.2', 'Orcu.U1.1', 'Orcu.U1.2')
  
  
  # Calculate frequencies for Orcu-U1
  haplotypes_U1 <- allele_table %>%
    dplyr::select(Orcu.U1.1, Orcu.U1.2) %>%
    unlist() %>%
    table() %>%
    prop.table()
  
  
  
  # Calculate frequencies for Orcu-U2
  haplotypes_U2 <- allele_table %>%
    dplyr::select(Orcu.U2.1, Orcu.U2.2) %>%
    unlist() %>%
    table() %>%
    prop.table()
  
  # Filter haplotypes with frequency > 0.1 for Orcu-U1
  filtered_hap_U1 <- names(haplotypes_U1[haplotypes_U1 > cutoff])
  
  # Filter haplotypes with frequency > 0.1 for Orcu-U2
  filtered_hap_U2 <- names(haplotypes_U2[haplotypes_U2 > cutoff])
  
  # Convert to numeric 
  hap_U1 <- as.numeric(filtered_hap_U1) %>% sort()
  hap_U2 <- as.numeric(filtered_hap_U2) %>% sort()
  
  colnames(allele_table) <- c('Orcu.U2.1', 'Orcu.U2.2', 'Orcu.U1.1', 'Orcu.U1.2')
  
  
  ld_results_rsquare <- matrix(NA, nrow = length(hap_U1), ncol = length(hap_U2))
  rownames(ld_results_rsquare) <- hap_U1
  colnames(ld_results_rsquare) <- hap_U2
  
  ld_results_Dprime <- matrix(NA, nrow = length(hap_U1), ncol = length(hap_U2))
  rownames(ld_results_Dprime) <- hap_U1
  colnames(ld_results_Dprime) <- hap_U2
  
  
  # Continuing the loop through each pair of haplotypes
  for(i in 1:length(hap_U1)){
    for(j in 1:length(hap_U2)){
      
      print(paste0("Orcu-U1 allele: ",hap_U1[i]))
      print(paste0("Orcu-U2 allele: ",hap_U2[j]))
      
      
      # Create a temporary dataframe for this pair
      allele_temp <- allele_table
      
      # Replace other haplotypes with "non" + haplotype number for Orcu-U1
      allele_temp[!(allele_temp$Orcu.U1.1 %in% c(hap_U1[i])), "Orcu.U1.1"] <- paste("non", hap_U1[i], sep="")
      allele_temp[!(allele_temp$Orcu.U1.2 %in% c(hap_U1[i])), "Orcu.U1.2"] <- paste("non", hap_U1[i], sep="")
      
      # Replace other haplotypes with "non" + haplotype number for Orcu-U2
      allele_temp[!(allele_temp$Orcu.U2.1 %in% c(hap_U2[j])), "Orcu.U2.1"] <- paste("non", hap_U2[j], sep="")
      allele_temp[!(allele_temp$Orcu.U2.2 %in% c(hap_U2[j])), "Orcu.U2.2"] <- paste("non", hap_U2[j], sep="")
      
      
      # Function to convert allele pairs to genotype format, keeping 'non'
      convert_to_genotype <- function(allele1, allele2) {
        # Directly create genotype string from allele names
        return(paste(allele1, "/", allele2, sep=""))
      }
      
      # Apply the function to each pair of alleles
      genotypes <- makeGenotypes(data.frame(
        OrcuU2 = mapply(convert_to_genotype, allele_temp$Orcu.U2.1, allele_temp$Orcu.U2.2),
        OrcuU1 = mapply(convert_to_genotype, allele_temp$Orcu.U1.1, allele_temp$Orcu.U1.2)
      ))
      
      genetics_output <- genetics::LD(genotypes)
      rsquare <- genetics_output$`R^2`[1,2]
      
      
      ld_results_rsquare[i, j] <- rsquare
    }
  }
  
  
  # Convert 'ld_results' matrix to a data frame for plotting
  ld_df <- as.data.frame(as.table(ld_results_rsquare))
  colpop <- c(rep(population, nrow(ld_df)))
  ld_df <- cbind(ld_df, colpop)
  names(ld_df) <- c("Orcu_U1", "Orcu_U2", "Rsquare", "Population")
  
  return(ld_df)
  
}


final_ld_df <- rbind(plot_ld_heatmap("FR"),plot_ld_heatmap("UK"),plot_ld_heatmap("AU"))

final_ld_df$Population <- factor(final_ld_df$Population, levels = c("FR", "UK", "AU"), labels = c("France", "United Kingdom", "Australia"))


my_palette <- colorRampPalette(brewer.pal(9, "Blues"))(1000)


# final_ld_df$Orcu_U1 <- sapply(final_ld_df$Orcu_U1, function(x) paste0(strsplit(name_lookup[paste0("10972b.seq", x)], "\\*")[[1]][2],
#                                                                       "(", x, ")"))
# final_ld_df$Orcu_U2 <- sapply(final_ld_df$Orcu_U2, function(x) paste0(strsplit(name_lookup[paste0("11016.seq", x)], "\\*")[[1]][2],
#                                                                       "(", x, ")"))


final_ld_df$Orcu_U1 <- sapply(final_ld_df$Orcu_U1, function(x) {
  full_name <- name_lookup[paste0("10972b.seq", x)]
  clean_name <- strsplit(full_name, "\\*")[[1]][2]
  return(clean_name)
})

final_ld_df$Orcu_U2 <- sapply(final_ld_df$Orcu_U2, function(x) {
  full_name <- name_lookup[paste0("11016.seq", x)]
  clean_name <- strsplit(full_name, "\\*")[[1]][2]
  return(clean_name)
})



final_ld_df$Orcu_U1 <- factor(final_ld_df$Orcu_U1, levels = sort(unique(final_ld_df$Orcu_U1)))
final_ld_df$Orcu_U2 <- factor(final_ld_df$Orcu_U2, levels = sort(unique(final_ld_df$Orcu_U2)))


p <- ggplot(final_ld_df, aes(x = Orcu_U1, y = Orcu_U2, fill = Rsquare)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = my_palette, limits = c(0, 1), na.value = NA) +
  facet_wrap(~ Population, ncol = 3, scales = "free") + 
  # coord_fixed() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, size = 6),
        axis.text.y = element_text(angle = 0, hjust = 0, size = 6),
        axis.title = element_text(size = 10,  face = 'bold.italic'),
        panel.grid.major = element_blank(), # Keep major grid lines blank if not needed
        panel.grid.minor = element_blank(), # Keep minor grid lines blank if not needed
        axis.ticks = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 12, face = "bold"),
        panel.background = element_rect(fill = "white", colour = NA), # Light grey background for each facet
        panel.spacing = unit(2, "lines")) +  # Adjust spacing between facets
  labs(fill = expression(r^2), x = "Orcu-U1", y = "Orcu-U2")


print(p)


ggsave(filename = paste0('LD_heatmap', cutoff, '.svg'), plot = p, width = 10, height = 3, device = "svg")
         



}




