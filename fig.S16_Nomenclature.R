library(Biostrings)
library(ape)
library(phangorn)
library(DECIPHER)
library(dplyr)

rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


TOTAL_NUMBER_OF_11016 <- 17
TOTAL_NUMBER_OF_10972b <- 29 

#################### translation function ####################
# Define a custom translation function that handles 'N's
translateHandleN <- function(dnaSeq) {
  placeholder <- "X" # Use 'X' for ambiguous amino acids
  # Function to translate a single codon
  translateCodon <- function(codon) {
    if("N" %in% strsplit(codon, "")[[1]]) {
      return(placeholder)
    } else {
      # Translate using Biostrings' translate function and catch errors
      tryCatch({
        return(as.character(Biostrings::translate(DNAString(codon))[1]))
      }, error = function(e) {
        return(placeholder) # Return placeholder if translation fails
      })
    }
  }
  
  # Split the sequence into codons and translate each
  codons <- strsplit(as.character(dnaSeq), split = "(?<=.{3})", perl = TRUE)[[1]]
  aaSeq <- sapply(codons, translateCodon)
  
  return(paste(aaSeq, collapse = ""))
}



#################### Distance matrix calculation ####################
calculateCustomDistance <- function(format,seq1, seq2) {
  # Ensure sequences are of the same length
  if (nchar(as.character(seq1)) != nchar(as.character(seq2))) {
    stop("Sequences must be of the same length.")
  }
  
  # Convert sequences to character vectors for iteration
  seq1_chars <- strsplit(as.character(seq1), "")[[1]]
  seq2_chars <- strsplit(as.character(seq2), "")[[1]]
  
  # Initialize distance
  distance <- 0
  ambiguity <- 'X'
  if (format == 'DNA'){
    ambiguity <- 'N'
  }
  # Iterate over aligned positions
  for (i in 1:length(seq1_chars)) {
    if (seq1_chars[i] != seq2_chars[i] && seq1_chars[i] != ambiguity && seq2_chars[i] != ambiguity) {
      distance <- distance + 1
    }
  }
  
  return(distance)
}


#################### Generate distance matrix ####################
GenerateDistanceMatrix <- function(format, seqs){
  n <- length(seqs)
  distance_matrix <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      distance_matrix[i, j] <- calculateCustomDistance(format, seqs[[i]], seqs[[j]])
      distance_matrix[j, i] <- calculateCustomDistance(format, seqs[[i]], seqs[[j]])
    }
  }
  diag(distance_matrix) <- 0  # Set diagonal to 0
  
  colnames(distance_matrix) <- names(seqs)
  rownames(distance_matrix) <- names(seqs)
  
  return(distance_matrix)
}

#############################################################


setwd(dirname(rstudioapi::getSourceEditorContext()$path))

pacbio_allele_file_folder <- '/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim1/05_PacBio_MHC_paralogs/final_pacbio_alleles_from_9genes'
miseq_allele_file <- '/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim1/01_New_nomenclature/All_seqs_2genes.fasta'


GENE_LIST <- list("10972b" = 1, "11016" = 2, "10972a" = 3, "13304" = 4, 
                  "06557" = 5, "10984" = 6, "08993b" = 7, "08993" = 8, "08993c" = 9)


nomenclature_df <- data.frame(OriginalName = character(), NomenclatureName = character(), stringsAsFactors = FALSE)



#################### AUTO nomenclature ####################
fasta_files <- list.files(path = pacbio_allele_file_folder, full.names = TRUE)

# Define gene on reverse strand
gene_reverse <- c('10972a','10972b','13304','10984')

# Read sequences from each file
pacbio_alleles <- lapply(fasta_files, readDNAStringSet)

# Initialize an empty list to hold processed sequences
final_alleles <- list()
final_aa_pbs <- list()

for (i in seq_along(fasta_files)) {
  
  gene_name <- sub("\\.hap\\.txt$", "", basename(fasta_files[i]))
  
  
  if (gene_name %in% gene_reverse) {
    # Reverse complement the sequences
    final_alleles [[gene_name]] <- reverseComplement(pacbio_alleles[[i]])
  } else {
    
    final_alleles [[gene_name]] <- pacbio_alleles[[i]]
  }
}

miseq_alleles <- readDNAStringSet(miseq_allele_file)

# Generate the list of desired sequence names from 11016.seq10/ 10972b.seq13
desired_names <- c(paste0("11016.seq", 10:TOTAL_NUMBER_OF_11016), paste0("10972b.seq", 13:TOTAL_NUMBER_OF_10972b))

filtered_miseq_alleles <- miseq_alleles[names(miseq_alleles) %in% desired_names]



for (gene in c('11016','10972b')) {
  indices <- grepl(gene, names(filtered_miseq_alleles))
  gene_sequences <- filtered_miseq_alleles[indices]
  final_alleles[[gene]] <- c(final_alleles[[gene]],gene_sequences)
}




for (gene in names(final_alleles)){
  # Sort the seqs with the numeric order of names 
  sequence_numbers <- as.numeric(sub(".*\\.seq", "", names(final_alleles[[gene]])))
  
  # Get the order of the sequence numbers
  order_indices <- order(sequence_numbers)
  
  # Subset the DNAStringSet based on the order
  sequences <- final_alleles[[gene]][order_indices]
  
  final_alleles[[gene]] <- sequences
  
  
  
  aligned_dna_sequences <- AlignSeqs(sequences, verbose = FALSE)
  
  start_matches <- c(start(matchPattern("GGCTCGCAC", aligned_dna_sequences[[1]])), 
                     start(matchPattern("GGTTCCCAC", aligned_dna_sequences[[1]])), 
                     start(matchPattern("GGCTTGCAC", aligned_dna_sequences[[1]])),
                     start(matchPattern("GGCTCACAC", aligned_dna_sequences[[1]])),
                     start(matchPattern("GCCTCCCAC", aligned_dna_sequences[[1]])),
                     start(matchPattern("GGATCGCAC", aligned_dna_sequences[[1]])))
  
  
  
  end_matches <- c(end(matchPattern("GAGATGGGGAAG", aligned_dna_sequences[[1]])),
                   end(matchPattern("CAAGACGGGAAG", aligned_dna_sequences[[1]])),
                   end(matchPattern("TACTTGACAAAG", aligned_dna_sequences[[1]])))
  
  # Nucleotide sequences of PBS region
  chunked_dna_sequences <- narrow(aligned_dna_sequences, start=start_matches, end=end_matches) 
  
  # Amino acid sequences of PBS region
  chunked_aa_sequences <- AAStringSet(unlist(lapply(chunked_dna_sequences, translateHandleN)))
  final_aa_pbs[[gene]] <- chunked_aa_sequences #for phylogeny
  
  # Nucleotide sequences OUT of PBS region
  dna_other1 <- as.character(narrow(aligned_dna_sequences,start=1, end=start_matches))
  dna_other2 <- as.character(narrow(aligned_dna_sequences,start=end_matches, end=width(aligned_dna_sequences)))
  other_dna_sequences <- DNAStringSet(paste0(dna_other1,dna_other2))
  names(other_dna_sequences) <- names(chunked_dna_sequences)
  
  final_nomenclature <- matrix(nrow=length(final_alleles[[gene]]),
                               ncol=4,)
  colnames(final_nomenclature) <- c("Group", "Subgroup", "Synonymous", "Outside")
  rownames(final_nomenclature) <- names(final_alleles[[gene]])
  
  
  ###################################### Distance matrix for three different seq sets ##########################################
  distances_aa <- GenerateDistanceMatrix('AA', chunked_aa_sequences)
  hclust_res <- hclust(as.dist(distances_aa), method = "complete")
  # Dynamic threshold for clustering based on distance <= 8
  groups <- cutree(hclust_res, h = 8)
  # Dynamic threshold for clustering based on distance == 0
  unique_groups <- unique(groups)
  
  distances_dna <- GenerateDistanceMatrix('DNA', chunked_dna_sequences)
  distances_otherdna <- GenerateDistanceMatrix('DNA', other_dna_sequences)
  
  ###################################################################################################
  
  
  for (group in unique_groups){
    
    group_mems <- names(groups[groups == group])
    final_nomenclature[group_mems,c(1)] = group
    
    if (length(group_mems) == 1){
      final_nomenclature[group_mems[1],c(2,3,4)] = c(1,1,1)
      next
      
    }else{
      subgroup_distances <- distances_aa[group_mems, group_mems]
      subgroups <- cutree(hclust(as.dist(subgroup_distances), method = "complete"), h = 0)
      unique_subgroups <- unique(subgroups)
    }
    
    for (sg in unique_subgroups){
      subgroup_mems <- names(subgroups[subgroups == sg])
      final_nomenclature[subgroup_mems,c(2)] = sg
      
      if (length(subgroup_mems) == 1){
        final_nomenclature[subgroup_mems[1],c(3,4)] = c(1,1)
        next
        
      }else{
        synonymous_distances <- distances_dna[subgroup_mems, subgroup_mems]
        synonymous_groups <- cutree(hclust(as.dist(synonymous_distances), method = "complete"), h = 0)
        unique_synonymous_groups <- unique(synonymous_groups)
      }
      for (syng in unique_synonymous_groups){
        synonymous_group_mems <- names(synonymous_groups[synonymous_groups == syng])
        final_nomenclature[synonymous_group_mems,c(3)] = syng
        
        if (length(synonymous_group_mems) == 1){
          final_nomenclature[synonymous_group_mems[1],c(4)] = c(1)
          next
        }else{
          
          otherdna_distances <- distances_otherdna[synonymous_group_mems, synonymous_group_mems]
          tail_groups <- cutree(hclust(as.dist(otherdna_distances), method = "complete"), h = 0)
          
          for (tt in names(tail_groups)){
            final_nomenclature[tt,c(4)] = tail_groups[[tt]]
          }
        }
        
      }
      
    }  
    
  }
  
  for (i in 1:nrow(final_nomenclature)){
    
    original_name <- rownames(final_nomenclature)[i]
    digits <- sprintf("%02d:%02d:%02d:%02d", 
                      final_nomenclature[i, 1], 
                      final_nomenclature[i, 2], 
                      final_nomenclature[i, 3], 
                      final_nomenclature[i, 4])
    
    new_gene_name <- paste0('U',GENE_LIST[[gene]])
    
    nomenclature_name <- c(paste0('Orcu-', new_gene_name, '*', digits))
    nomenclature_df <- rbind(nomenclature_df, data.frame(OriginalName = original_name, NomenclatureName = nomenclature_name, stringsAsFactors = FALSE))
    
  }
  
}

nomenclature_df <- nomenclature_df[order(nomenclature_df$NomenclatureName), ]

# If you wish to reset the row names (optional)
rownames(nomenclature_df) <- NULL




all_names <- nomenclature_df$NomenclatureName

# Function to trim the NomenclatureName
trim_name <- function(name) {
  parts <- unlist(strsplit(name, ":"))
  for (i in 2:length(parts)) {
    trimmed_name <- paste(parts[1:i], collapse = ":")
    trimmed_all <- sapply(all_names, function(x) paste(unlist(strsplit(x, ":"))[1:i], collapse = ":"))
    if (sum(trimmed_name == trimmed_all) == 1) {
      return(trimmed_name)
    }
  }
  return(name)
}

# Apply the function to create the ShortName column

nomenclature_df_with_short_names <- nomenclature_df %>%
  rowwise() %>%
  mutate(ShortName = trim_name(NomenclatureName))



print(nomenclature_df_with_short_names)
write.csv(nomenclature_df_with_short_names, "nomenclature_names.csv", row.names = FALSE)



############### Write ALL FASTA with new names ###############
# Create a new list to store sequences with updated names
renamed_sequences <- list()

for (gene in names(final_alleles)) {
  # Get the sequences for this gene
  gene_sequences <- final_alleles[[gene]]
  
  # Create a new DNAStringSet for this gene
  new_sequences <- DNAStringSet(gene_sequences)
  
  # Update the names using the nomenclature mapping
  for (i in seq_along(gene_sequences)) {
    old_name <- names(gene_sequences)[i]
    # Find the corresponding new name in nomenclature_df
    new_name <- nomenclature_df_with_short_names$NomenclatureName[nomenclature_df$OriginalName == old_name]
    
    # Check if we found a matching name
    if (length(new_name) == 0) {
      warning(sprintf("No matching new name found for %s in gene %s", old_name, gene))
      new_name <- old_name  # Keep original name if no match found
    }
    
    names(new_sequences)[i] <- new_name
  }
  
  # Add to our collection
  renamed_sequences[[gene]] <- new_sequences
}

# Combine all sequences into a single DNAStringSet
all_renamed_sequences <- do.call(c, renamed_sequences)

final_seqs <- DNAStringSet()
for (gene_seqs in all_renamed_sequences) {
  final_seqs <- c(final_seqs, gene_seqs)  
}

# Write to FASTA file
writeXStringSet(final_seqs, "9genes_ALL_sequences_with_new_names.fasta")



# ############### Write NCBI submission files ###############
# # Initialize empty DNAStringSet objects for pacbio and miseq
# # Collect and sort sequences
# pacbio_seqs <- DNAStringSet()
# miseq_seqs <- DNAStringSet()
# 
# # First pass - collect sequences with names
# for (gene in names(final_alleles)) {
#   gene_sequences <- final_alleles[[gene]]
# 
#   for (i in seq_along(gene_sequences)) {
#     orig_name <- names(gene_sequences)[i]
#     new_name <- nomenclature_df$NomenclatureName[nomenclature_df$OriginalName == orig_name]
#     seq_to_add <- gene_sequences[i]
# 
#     if (gene == "11016") {
#       seq_num <- as.numeric(sub(".*\\.seq", "", orig_name))
#       if (seq_num >= 10) {
#         miseq_seqs <- c(miseq_seqs, seq_to_add)
#         names(miseq_seqs)[length(miseq_seqs)] <- new_name
#       } else {
#         pacbio_seqs <- c(pacbio_seqs, seq_to_add)
#         names(pacbio_seqs)[length(pacbio_seqs)] <- new_name
#       }
#     } else if (gene == "10972b") {
#       seq_num <- as.numeric(sub(".*\\.seq", "", orig_name))
#       if (seq_num >= 13) {
#         miseq_seqs <- c(miseq_seqs, seq_to_add)
#         names(miseq_seqs)[length(miseq_seqs)] <- new_name
#       } else {
#         pacbio_seqs <- c(pacbio_seqs, seq_to_add)
#         names(pacbio_seqs)[length(pacbio_seqs)] <- new_name
#       }
#     } else {
#       pacbio_seqs <- c(pacbio_seqs, seq_to_add)
#       names(pacbio_seqs)[length(pacbio_seqs)] <- new_name
#     }
#   }
# }
# 
# # Sort sequences by name
# pacbio_seqs <- pacbio_seqs[order(names(pacbio_seqs))]
# miseq_seqs <- miseq_seqs[order(names(miseq_seqs))]
# 
# # Format descriptions with BioProject
# format_pacbio_desc <- function(name) {
#   sprintf("%s [organism=Oryctolagus cuniculus] [BioProject=PRJNA1027777] %s, PacBio long-reads assembled transcript",
#           name, name)
# }
# 
# format_miseq_desc <- function(name) {
#   sprintf("%s [organism=Oryctolagus cuniculus] [BioProject=PRJNA1027777] %s, MiSeq amplicon assembled exon 2-3 and partial ex",
#           name, name)
# }
# 
# names(pacbio_seqs) <- sapply(names(pacbio_seqs), format_pacbio_desc)
# names(miseq_seqs) <- sapply(names(miseq_seqs), format_miseq_desc)
# 
# writeXStringSet(pacbio_seqs, "pacbio_submission.fasta")
# writeXStringSet(miseq_seqs, "miseq_submission.fasta")
# 


############### Write NCBI submission files ###############

#################### Sequence Processing Functions ####################

# Function to process gene sequences with position-wise consensus and N-filling
processGeneSequences <- function(seqs) {
  # First align all sequences
  aligned_seqs <- AlignSeqs(seqs, verbose = FALSE)
  
  # Convert to character matrix for easier processing
  seq_matrix <- do.call(rbind, strsplit(as.character(aligned_seqs), ""))
  
  # Process each sequence
  processed_seqs <- vector("list", length(seqs))
  names(processed_seqs) <- names(seqs)
  
  # For each sequence
  for (i in 1:nrow(seq_matrix)) {
    current_seq <- seq_matrix[i,]
    
    # For each position that has N in current sequence
    for (j in 1:ncol(seq_matrix)) {
      if (current_seq[j] == "N") {
        # Get all bases at this position from other sequences
        other_bases <- seq_matrix[-i, j]
        # Only consider non-N bases
        valid_bases <- other_bases[other_bases != "N"]
        
        if (length(valid_bases) > 0) {
          # Replace N with most common valid base
          base_counts <- table(valid_bases)
          current_seq[j] <- names(which.max(base_counts))
        }
      }
    }
    
    # Convert back to sequence
    processed_seqs[[i]] <- DNAString(paste(current_seq[current_seq != "-"], collapse=""))
  }
  
  return(DNAStringSet(processed_seqs))
}

# Initialize empty DNAStringSet objects for pacbio and miseq
pacbio_seqs <- DNAStringSet()
miseq_seqs <- DNAStringSet()

# Process each gene
for (gene in names(final_alleles)) {
  # Process all sequences in this gene
  processed_gene_seqs <- processGeneSequences(final_alleles[[gene]])
  
  # Add sequences to appropriate output set
  for (i in seq_along(processed_gene_seqs)) {
    orig_name <- names(processed_gene_seqs)[i]
    new_name <- nomenclature_df$NomenclatureName[nomenclature_df$OriginalName == orig_name]
    seq_to_add <- processed_gene_seqs[[i]]
    
    if (gene == "11016") {
      seq_num <- as.numeric(sub(".*\\.seq", "", orig_name))
      if (seq_num >= 10) {
        miseq_seqs <- c(miseq_seqs, DNAStringSet(seq_to_add))
        names(miseq_seqs)[length(miseq_seqs)] <- new_name
      } else {
        pacbio_seqs <- c(pacbio_seqs, DNAStringSet(seq_to_add))
        names(pacbio_seqs)[length(pacbio_seqs)] <- new_name
      }
    } else if (gene == "10972b") {
      seq_num <- as.numeric(sub(".*\\.seq", "", orig_name))
      if (seq_num >= 13) {
        miseq_seqs <- c(miseq_seqs, DNAStringSet(seq_to_add))
        names(miseq_seqs)[length(miseq_seqs)] <- new_name
      } else {
        pacbio_seqs <- c(pacbio_seqs, DNAStringSet(seq_to_add))
        names(pacbio_seqs)[length(pacbio_seqs)] <- new_name
      }
    } else {
      pacbio_seqs <- c(pacbio_seqs, DNAStringSet(seq_to_add))
      names(pacbio_seqs)[length(pacbio_seqs)] <- new_name
    }
  }
}

# Sort sequences by name
pacbio_seqs <- pacbio_seqs[order(names(pacbio_seqs))]
miseq_seqs <- miseq_seqs[order(names(miseq_seqs))]

# Format descriptions with BioProject
format_pacbio_desc <- function(name) {
  sprintf("%s [organism=Oryctolagus cuniculus] [BioProject=PRJNA1027777] %s, PacBio long-reads assembled transcript",
          name, name)
}

format_miseq_desc <- function(name) {
  sprintf("%s [organism=Oryctolagus cuniculus] [BioProject=PRJNA1027777] %s, MiSeq amplicon assembled exon 2-3 and partial exon 1",
          name, name)
}

# Update sequence names with formatted descriptions
names(pacbio_seqs) <- sapply(names(pacbio_seqs), format_pacbio_desc)
names(miseq_seqs) <- sapply(names(miseq_seqs), format_miseq_desc)

# Write sequences to files
writeXStringSet(pacbio_seqs, "pacbio_submission_v2.fasta")
writeXStringSet(miseq_seqs, "miseq_submission_v2.fasta")


# 
# HLA_A_seq <- "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQIMYGCDVGSDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHEAEQLRAYLDGTCVEWLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWELSSQPTIPIVGIIAGLVLLGAVITGAVVAAVMWRRKSSDRKGGSYTQAASSDSAQGSDVSLTACKV"
# HLA_A_aaseq <- AAStringSet(substr(HLA_A_seq, 25, 200))
# names(HLA_A_aaseq) <- 'HLA-A'
# 
# pbs_9genes <- unlist(AAStringSetList(final_aa_pbs))
# 
# all_pbs <- c(pbs_9genes, HLA_A_aaseq)
# 
# 
# 
# aligned_all_pbs <- AlignSeqs(all_pbs, verbose = FALSE)
# aln <- as.AAbin(aligned_all_pbs)
# 
# 
# names(aln) <- sub("(.*)\\..*\\.", "\\1.", names(aln))
# new_nomenclature_file <- read.csv("nomenclature_names.csv")
# 
# # Update names to new nomenclature
# names(aln) <- sapply(names(aln), function(name) {
#   new_name <- new_nomenclature_file$NomenclatureName[match(name, new_nomenclature_file$OriginalName)]
#   if (!is.na(new_name)) return(new_name) else return(name)
# })
# 
# 
# # phy_dat <- phangorn::phyDat(aln, type = "AA")
# # tree <- phangorn::pratchet(phy_dat)
# # rooted_tree <- ape::root(tree, outgroup="HLA-A", resolve.root = TRUE)
# 
# phy_dat <- phangorn::phyDat(aln, type = "AA")
# dist.mat <- dist.aa(as.AAbin(phy_dat))
# treeNJ <- NJ(dist.mat)
# rooted_tree <- ape::root(treeNJ, outgroup="HLA-A", resolve.root = TRUE)
# 
# 
# 
# 
# 
# gene_color_map <- c("Orcu-U1" = "#E41A1C", "Orcu-U2" = "#FF7F00", "Orcu-U3" = "gray",
#                     "Orcu-U4" = "#4DAF4A", "Orcu-U5" = "cyan", "Orcu-U6" = "#377EB8",
#                     "Orcu-U7" = "#984EA3", "Orcu-U8" = "#A65628", "Orcu-U9" = "#F781BF") # Color map for genes
# 
# 
# 
# tip_colors <- sapply(rooted_tree$tip.label, function(label) {
#   gene_prefix <- sub("\\*.*", "", label)
#   if (gene_prefix %in% names(gene_color_map)) {
#     return(gene_color_map[gene_prefix])
#   } else {
#     return("black")
#   }
# })
# 
# 
# ggt <- ggtree(rooted_tree) +
#   # geom_tree(layout = "rectangle") +
#   geom_tiplab(aes(color = label), size=2, align = FALSE, linesize = 2, hjust = 0.01) +
#   scale_color_manual(values = setNames(tip_colors, rooted_tree$tip.label)) +
#   # ggtitle(domain) +
#   theme_tree2() +
#   # theme(plot.title = element_text(size=10, face="bold")) +
#   theme(legend.position = "none")  # Hide the legend if desired
# 
# ggt
# ggsave("9genes_phylogeny_pbs_domains.pdf", ggt, width = 15, height = 8)




