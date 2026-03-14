# Revised Analysis Script
# Panel A: Expected Heterozygosity (existing)
# Panel B: Expected Mean Allele Divergence (NEW - bars like Panel A)
# Panel C: Distance to Selected Alleles (SIMPLIFIED - one line per panel, pooling timepoints)

######################################################################################################
######################################################################################################
# Load required libraries
library(ggplot2)
library(dplyr)
library(grantham)
library(Biostrings)
library(seqinr)
library(tidyr)
library(patchwork)

# Set working directory
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
cwd <- dirname(rstudioapi::getSourceEditorContext()$path)

cat("=== REVISED COMBINED ANALYSIS ===\n")
cat("Panel A: Expected Heterozygosity (existing)\n")
cat("Panel B: Expected Mean Allele Divergence (NEW)\n")
cat("Panel C: Distance to Selected Alleles (pooled timepoints)\n\n")

### PART 1: LOAD AND PROCESS RAW DATA ###

# Your file paths
wd <- '/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim2/04_museum_typing/mm97_b0010_combined_paired/'
norm_file <- '/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim2/04_museum_typing/norm_cov_combined.txt'
miseq_file <- '/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim1/11_Miseq/FINAL_Miseq_alleles.csv'
nomenclature_file <- '/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim1/01_New_nomenclature/nomenclature_names.csv'

setwd(paste0(wd,'/coverage_metrics'))

# Load nomenclature mapping
cat("Loading nomenclature mapping...\n")
nomenclature <- read.csv(nomenclature_file, stringsAsFactors = FALSE)
lookup <- setNames(nomenclature$ShortName, nomenclature$OriginalName)
cat(sprintf("Loaded %d nomenclature mappings\n", length(lookup)))

# Load MiSeq data
miseq_table <- read.csv(file = miseq_file, header = T)

# Load historical samples
historical_samplelist <- read.table(norm_file, header = T)[,3]
norm_coverage <- read.table(norm_file, header = T)[,6]

genes <- c('11016','10972b')
sample_alleles <- matrix(nrow = length(historical_samplelist), ncol = 5)
colnames(sample_alleles) <- c('Sample','11016.A1', '11016.A2', '10972b.A1', '10972b.A2')

# Prepare modern table
modern_table <- miseq_table %>%
  dplyr::filter(Sequence.Group == 'Miseq (in Sci)' & Category != 'TOO MANY') %>%
  dplyr::select(Nomenclature, Orcu.2.allele1, Orcu.2.allele2, Orcu.1.allele1, Orcu.1.allele2)

# Function to revise alleles
revise_alleles <- function(a1, a2, bad_values) {
  if(a1 %in% bad_values & a2 %in% bad_values) {
    return(c("*", "*"))
  } else if(a1 %in% bad_values) {
    return(c(a2, a2))
  } else if(a2 %in% bad_values) {
    return(c(a1, a1))
  } else {
    return(c(a1, a2))
  }
}

# Revise modern table alleles
for(i in 1:nrow(modern_table)) {
  alleles <- revise_alleles(modern_table[i, "Orcu.2.allele1"], modern_table[i, "Orcu.2.allele2"], c(7))
  modern_table[i, c("Orcu.2.allele1", "Orcu.2.allele2")] <- alleles
  
  alleles <- revise_alleles(modern_table[i, "Orcu.1.allele1"], modern_table[i, "Orcu.1.allele2"], c(11))
  modern_table[i, c("Orcu.1.allele1", "Orcu.1.allele2")] <- alleles
}

colnames(modern_table) <- c('Sample','11016.A1', '11016.A2', '10972b.A1', '10972b.A2')

# Load historical allele calls
fs_withpseudo <- list.files('withpseudo/')
fs_nonpseudo <- list.files("nonpseudo/")

for (i in seq_along(historical_samplelist)){
  sample_alleles[i,1] <- historical_samplelist[i]
  for (g in genes){
    if (g == '11016'){
      x <- 2
      psu_allele <- 7
    } else if (g == '10972b') {
      x <- 4
      psu_allele <- 11
    }
    
    search_string <- paste(historical_samplelist[i],'_result.tsv.',g,sep = '')
    selected_allele_f <- fs_withpseudo[grep(paste0("^", search_string), fs_withpseudo)]
    
    if (any(grepl(psu_allele, selected_allele_f))){
      selected_allele_f <- fs_nonpseudo[grep(paste0("^", search_string), fs_nonpseudo)]
    }
    
    if (length(selected_allele_f) > 0){
      if (length(selected_allele_f) == 2) {
        a <- as.integer(sub('seq', '', unlist(strsplit(selected_allele_f[1], '\\.'))[4]))
        b <- as.integer(sub('seq', '', unlist(strsplit(selected_allele_f[2], '\\.'))[4]))
        allele_1 <- min(a,b)
        allele_2 <- max(a,b)
      } else if (length(selected_allele_f) == 1) {
        allele_1 <- as.integer(sub('seq', '', unlist(strsplit(selected_allele_f[1], '\\.'))[4]))
        allele_2 <- allele_1
      }
      sample_alleles[i,x] <- as.integer(allele_1)
      sample_alleles[i,x+1] <- as.integer(allele_2)
    }
  }
}

# Create historical table and filter by coverage
historical_table <- cbind(sample_alleles, norm_coverage)
historical_table <- as.data.frame(subset(historical_table, norm_coverage >= 3)[,-6])

# Add population information
add_population_info <- function(data_table) {
  data_table %>%
    mutate(
      Population = case_when(
        grepl("AU", Sample) ~ "Australia",
        grepl("FR", Sample) ~ "France", 
        grepl("UK", Sample) ~ "Britain",
        TRUE ~ "Unknown"
      )
    ) %>%
    filter(Population != "Unknown")
}

historical_data <- historical_table %>%
  add_population_info() %>%
  mutate(Time_Period = "Historical")

modern_data <- modern_table %>%
  add_population_info() %>%
  mutate(Time_Period = "Modern")

# Combine data
combined_data <- rbind(
  historical_data %>% dplyr::select(Sample, Population, Time_Period, `11016.A1`, `11016.A2`, `10972b.A1`, `10972b.A2`),
  modern_data %>% dplyr::select(Sample, Population, Time_Period, `11016.A1`, `11016.A2`, `10972b.A1`, `10972b.A2`)
)

# Remove samples with '*' in any allele column
cat(sprintf("Before filtering: %d samples\n", nrow(combined_data)))
combined_data_clean <- combined_data %>%
  filter(
    `11016.A1` != "*" & `11016.A2` != "*" & 
      `10972b.A1` != "*" & `10972b.A2` != "*" &
      `11016.A1` != "-" & `11016.A2` != "-" & 
      `10972b.A1` != "-" & `10972b.A2` != "-" &
      !is.na(`11016.A1`) & !is.na(`11016.A2`) & 
      !is.na(`10972b.A1`) & !is.na(`10972b.A2`)
  )
cat(sprintf("After removing samples with '*' or NA: %d samples\n", nrow(combined_data_clean)))

### PART 2: CALCULATE ALLELE FREQUENCIES ###

calculate_allele_frequencies <- function(data_table, gene_col1, gene_col2, gene_name) {
  freq_data <- data.frame()
  
  for (pop in unique(data_table$Population)) {
    for (time in unique(data_table$Time_Period)) {
      subset_data <- data_table %>% filter(Population == pop, Time_Period == time)
      
      if (nrow(subset_data) > 0) {
        subset_alleles <- c(subset_data[[gene_col1]], subset_data[[gene_col2]])
        subset_alleles <- subset_alleles[!is.na(subset_alleles) & subset_alleles != "*"]
        
        if (length(subset_alleles) > 0) {
          allele_counts <- table(subset_alleles)
          total_alleles <- sum(allele_counts)
          
          for (allele in names(allele_counts)) {
            allele_num <- suppressWarnings(as.numeric(allele))
            
            if (!is.na(allele_num)) {
              original_name <- paste0(ifelse(gene_name == "Orcu-U1", "10972b", "11016"), ".seq", allele_num)
              formatted_allele <- lookup[original_name]
              
              if (is.na(formatted_allele)) {
                formatted_allele <- sprintf("%s*%02d:01", gene_name, allele_num)
                warning(paste("No nomenclature found for:", original_name, "- using fallback"))
              }
              
              freq_data <- rbind(freq_data, data.frame(
                Gene = gene_name,
                Allele = formatted_allele,
                Allele_Number = allele_num,
                Population = pop,
                Time_Period = time,
                Count = as.numeric(allele_counts[allele]),
                Total_Alleles = total_alleles,
                Frequency = as.numeric(allele_counts[allele]) / total_alleles,
                stringsAsFactors = FALSE
              ))
            }
          }
        }
      }
    }
  }
  return(freq_data)
}

# Calculate frequencies
u1_frequencies <- calculate_allele_frequencies(combined_data_clean, "10972b.A1", "10972b.A2", "Orcu-U1")
u2_frequencies <- calculate_allele_frequencies(combined_data_clean, "11016.A1", "11016.A2", "Orcu-U2")
all_frequencies <- rbind(u1_frequencies, u2_frequencies)

cat(sprintf("Total allele frequency records: %d\n", nrow(all_frequencies)))

### PART 3: LOAD SEQUENCE DATA ###

# Load FASTA sequences
fasta_sequences <- readDNAStringSet("/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim1/01_New_nomenclature/All_seqs_2genes.fasta")
sequence_names <- names(fasta_sequences)

# Extract coding region
extract_coding_region <- function(dna_seq) {
  tryCatch({
    seq_char <- as.character(dna_seq)
    if (nchar(seq_char) != 561) return("")
    
    coding_seq <- substr(seq_char, 67, 402)
    if (nchar(coding_seq) != 336) return("")
    
    nt_vector <- unlist(strsplit(coding_seq, ""))
    aa_seq <- seqinr::translate(nt_vector, frame = 0, sens = "F", NAstring = "X")
    aa_string <- paste(aa_seq, collapse = "")
    
    return(aa_string)
  }, error = function(e) {
    return("")
  })
}

# Process sequences
all_sequences <- list()
for (i in 1:length(fasta_sequences)) {
  seq_name <- as.character(sequence_names[i])
  
  if (grepl("10972b", seq_name) || grepl("11016", seq_name)) {
    aa_seq <- extract_coding_region(fasta_sequences[[i]])
    if (nchar(aa_seq) > 100) {
      final_name <- lookup[seq_name]
      
      if (is.na(final_name)) {
        if (grepl("10972b", seq_name)) {
          seq_num <- as.numeric(gsub(".*seq", "", seq_name))
          final_name <- sprintf("Orcu-U1*%02d:01", seq_num)
        } else {
          seq_num <- as.numeric(gsub(".*seq", "", seq_name))
          final_name <- sprintf("Orcu-U2*%02d:01", seq_num)
        }
        warning(paste("No nomenclature found for:", seq_name, "- using fallback:", final_name))
      }
      
      all_sequences[[final_name]] <- aa_seq
    }
  }
}

cat(sprintf("Processed %d amino acid sequences\n", length(all_sequences)))

### PART 4: GRANTHAM DISTANCE CALCULATION ###

calculate_grantham_distance <- function(seq1, seq2) {
  tryCatch({
    seq1_chars <- unlist(strsplit(seq1, ""))
    seq2_chars <- unlist(strsplit(seq2, ""))
    min_length <- min(length(seq1_chars), length(seq2_chars))
    seq1_chars <- seq1_chars[1:min_length]
    seq2_chars <- seq2_chars[1:min_length]
    
    valid_aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                  "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
    valid_positions <- seq1_chars %in% valid_aa & seq2_chars %in% valid_aa
    
    if (sum(valid_positions) < 30) return(NA)
    
    seq1_filtered <- seq1_chars[valid_positions]
    seq2_filtered <- seq2_chars[valid_positions]
    
    tryCatch({
      seq1_three <- as_three_letter(seq1_filtered)
      seq2_three <- as_three_letter(seq2_filtered)
      
      distances_result <- grantham_distance(seq1_three, seq2_three)
      
      if (is.data.frame(distances_result)) {
        if ("distance" %in% colnames(distances_result)) {
          dist_values <- distances_result$distance
        } else if ("d" %in% colnames(distances_result)) {
          dist_values <- distances_result$d
        } else if (ncol(distances_result) >= 1) {
          dist_values <- distances_result[, ncol(distances_result)]
        } else {
          return(NA)
        }
      } else if (is.vector(distances_result) && is.numeric(distances_result)) {
        dist_values <- distances_result
      } else if (is.matrix(distances_result)) {
        dist_values <- as.vector(distances_result)
      } else {
        return(NA)
      }
      
      if (is.numeric(dist_values) && length(dist_values) > 0) {
        mean_dist <- mean(dist_values, na.rm = TRUE)
        return(mean_dist)
      } else {
        return(NA)
      }
    }, error = function(e2) {
      return(NA)
    })
  }, error = function(e) {
    return(NA)
  })
}


######################################################################################################
######################################################################################################

# Panel B: Expected Mean Allele Divergence - EXACT Panel A Style
# Following the exact same approach as Panel A
# Use the same data variable names as Panel A

# Use the same combined_clean dataset from Panel A
# If your data is named differently, replace 'combined_clean' with your actual variable name
# For example, if it's 'combined_data_clean', change all instances below

### 1. CALCULATE REAL EXPECTED MEAN DIVERGENCE VALUES ###

# Function to calculate expected mean divergence from real data
calculate_real_expected_mean_divergence <- function(allele_vector, sequences, gene_name) {
  # Remove missing values
  alleles <- allele_vector[!is.na(allele_vector) & allele_vector != "*" & allele_vector != "-"]
  
  if(length(alleles) < 4) return(NA)  # Need minimum sample size
  
  # Calculate allele frequencies
  allele_counts <- table(alleles)
  allele_frequencies <- allele_counts / length(alleles)
  
  # Convert allele numbers to nomenclature names
  allele_names <- c()
  for (allele in names(allele_frequencies)) {
    allele_num <- suppressWarnings(as.numeric(allele))
    if (!is.na(allele_num)) {
      original_name <- paste0(ifelse(gene_name == "Orcu-U1", "10972b", "11016"), ".seq", allele_num)
      formatted_allele <- lookup[original_name]
      
      if (is.na(formatted_allele)) {
        formatted_allele <- sprintf("%s*%02d:01", gene_name, allele_num)
      }
      allele_names <- c(allele_names, formatted_allele)
    }
  }
  names(allele_frequencies) <- allele_names
  
  # Filter to alleles present in sequence data
  available_alleles <- allele_names[allele_names %in% names(sequences)]
  if (length(available_alleles) < 2) return(NA)
  
  allele_frequencies <- allele_frequencies[available_alleles]
  
  # Calculate expected mean divergence using Hardy-Weinberg weights
  total_weighted_divergence <- 0
  total_weight <- 0
  
  for (i in 1:(length(available_alleles)-1)) {
    for (j in (i+1):length(available_alleles)) {
      allele_i <- available_alleles[i]
      allele_j <- available_alleles[j]
      
      # Calculate Grantham distance between alleles
      divergence <- calculate_grantham_distance(sequences[[allele_i]], sequences[[allele_j]])
      
      if (!is.na(divergence)) {
        # Weight by probability of this heterozygote under Hardy-Weinberg
        weight <- 2 * allele_frequencies[allele_i] * allele_frequencies[allele_j]
        
        total_weighted_divergence <- total_weighted_divergence + (weight * divergence)
        total_weight <- total_weight + weight
      }
    }
  }
  
  if (total_weight > 0) {
    return(total_weighted_divergence / total_weight)
  } else {
    return(NA)
  }
}

# Calculate REAL Expected Mean Divergence values for each population-time-gene combination
cat("\n=== STEP 1: CALCULATING REAL EXPECTED MEAN DIVERGENCE VALUES ===\n")

real_divergence_results <- data.frame()

for (pop in unique(combined_data_clean$Population)) {
  for (time in unique(combined_data_clean$Time_Period)) {
    for (gene_info in list(c("10972b.A1", "10972b.A2", "Orcu-U1"),
                           c("11016.A1", "11016.A2", "Orcu-U2"))) {
      
      subset_data <- combined_data_clean %>% 
        filter(Population == pop, Time_Period == time) %>%
        filter(!is.na(!!sym(gene_info[1])) & !is.na(!!sym(gene_info[2])) &
                 !!sym(gene_info[1]) != "*" & !!sym(gene_info[2]) != "*" &
                 !!sym(gene_info[1]) != "-" & !!sym(gene_info[2]) != "-")
      
      if (nrow(subset_data) >= 5) {
        # Get all alleles from the real data
        all_alleles <- c(subset_data[[gene_info[1]]], subset_data[[gene_info[2]]])
        
        # Calculate real expected mean divergence
        real_div <- calculate_real_expected_mean_divergence(all_alleles, all_sequences, gene_info[3])
        
        if (!is.na(real_div)) {
          cat(sprintf("Real Expected Mean Divergence: %s %s %s = %.4f (n=%d individuals)\n",
                      pop, time, gene_info[3], real_div, nrow(subset_data)))
          
          real_divergence_results <- rbind(real_divergence_results, data.frame(
            Population = pop,
            Time_Period = time,
            Gene = gene_info[3],
            Real_Mean_Divergence = real_div,
            Sample_size = nrow(subset_data),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}

### 2. BOOTSTRAP FOR STANDARD ERRORS ONLY ###
cat("\n=== STEP 2: BOOTSTRAP FOR STANDARD ERRORS ONLY ===\n")

# Function to calculate bootstrap CIs for divergence
calculate_divergence_bootstrap_ci <- function(individuals_data, gene_cols, sequences, gene_name, n_bootstrap = 1000) {
  cat(sprintf("Bootstrapping %s for 95%% confidence intervals...\n", gene_cols[3]))
  
  bootstrap_div_values <- numeric(n_bootstrap)
  
  for (b in 1:n_bootstrap) {
    boot_sample <- individuals_data[sample(nrow(individuals_data), replace = TRUE), ]
    boot_alleles <- c(boot_sample[[gene_cols[1]]], boot_sample[[gene_cols[2]]])
    bootstrap_div_values[b] <- calculate_real_expected_mean_divergence(boot_alleles, sequences, gene_name)
  }
  
  # Remove NA values
  bootstrap_div_values <- bootstrap_div_values[!is.na(bootstrap_div_values)]
  
  if (length(bootstrap_div_values) > 0) {
    ci <- quantile(bootstrap_div_values, probs = c(0.025, 0.975), na.rm = TRUE)
    return(list(lower = ci[1], upper = ci[2]))
  } else {
    return(list(lower = NA, upper = NA))
  }
}

# Calculate bootstrap CIs for each group
bootstrap_div_ci_results <- data.frame()

for (pop in unique(combined_data_clean$Population)) {
  for (time in unique(combined_data_clean$Time_Period)) {
    for (gene_info in list(c("10972b.A1", "10972b.A2", "Orcu-U1"),
                           c("11016.A1", "11016.A2", "Orcu-U2"))) {
      
      subset_data <- combined_data_clean %>% 
        filter(Population == pop, Time_Period == time) %>%
        filter(!is.na(!!sym(gene_info[1])) & !is.na(!!sym(gene_info[2])) &
                 !!sym(gene_info[1]) != "*" & !!sym(gene_info[2]) != "*" &
                 !!sym(gene_info[1]) != "-" & !!sym(gene_info[2]) != "-")
      
      if (nrow(subset_data) >= 5) {
        bootstrap_ci <- calculate_divergence_bootstrap_ci(subset_data, gene_info, all_sequences, gene_info[3], n_bootstrap = 1000)
        
        if (!is.na(bootstrap_ci$lower) && !is.na(bootstrap_ci$upper)) {
          cat(sprintf("Bootstrap 95%% CI: %s %s %s = [%.4f, %.4f]\n",
                      pop, time, gene_info[3], bootstrap_ci$lower, bootstrap_ci$upper))
          
          bootstrap_div_ci_results <- rbind(bootstrap_div_ci_results, data.frame(
            Population = pop,
            Time_Period = time,
            Gene = gene_info[3],
            CI_lower = bootstrap_ci$lower,
            CI_upper = bootstrap_ci$upper,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}

### 3. STATISTICAL TESTING ON REAL INDIVIDUAL DATA ###
cat("\n=== STEP 3: STATISTICAL TESTING ON REAL INDIVIDUAL DATA ===\n")

# Function to perform proper statistical tests on individual-level data for divergence
perform_divergence_statistical_tests <- function(individual_data, sequences) {
  results <- data.frame()
  
  for (pop in unique(individual_data$Population)) {
    for (gene_info in list(c("10972b.A1", "10972b.A2", "Orcu-U1"),
                           c("11016.A1", "11016.A2", "Orcu-U2"))) {
      
      # Get historical and modern data for this population-gene
      hist_data <- individual_data %>% 
        filter(Population == pop, Time_Period == "Historical") %>%
        filter(!is.na(!!sym(gene_info[1])) & !is.na(!!sym(gene_info[2])) &
                 !!sym(gene_info[1]) != "*" & !!sym(gene_info[2]) != "*" &
                 !!sym(gene_info[1]) != "-" & !!sym(gene_info[2]) != "-")
      
      mod_data <- individual_data %>% 
        filter(Population == pop, Time_Period == "Modern") %>%
        filter(!is.na(!!sym(gene_info[1])) & !is.na(!!sym(gene_info[2])) &
                 !!sym(gene_info[1]) != "*" & !!sym(gene_info[2]) != "*" &
                 !!sym(gene_info[1]) != "-" & !!sym(gene_info[2]) != "-")
      
      if (nrow(hist_data) >= 5 && nrow(mod_data) >= 5) {
        
        # Calculate divergence for each group
        hist_alleles <- c(hist_data[[gene_info[1]]], hist_data[[gene_info[2]]])
        mod_alleles <- c(mod_data[[gene_info[1]]], mod_data[[gene_info[2]]])
        
        hist_div <- calculate_real_expected_mean_divergence(hist_alleles, sequences, gene_info[3])
        mod_div <- calculate_real_expected_mean_divergence(mod_alleles, sequences, gene_info[3])
        
        if (!is.na(hist_div) && !is.na(mod_div)) {
          
          observed_diff <- mod_div - hist_div
          
          # PERMUTATION TEST on individual data
          # Combine all individuals from both time periods
          all_individuals <- rbind(
            hist_data %>% mutate(Original_Time = "Historical"),
            mod_data %>% mutate(Original_Time = "Modern")
          )
          
          n_hist <- nrow(hist_data)
          n_mod <- nrow(mod_data)
          n_permutations <- 1000
          
          cat(sprintf("Testing %s %s: n_hist=%d, n_mod=%d\n",
                      pop, gene_info[3], n_hist, n_mod))
          
          # Permutation test (1000 iterations to match Panel A preference)
          perm_diffs <- replicate(n_permutations, {
            # Randomly reassign time periods to individuals
            shuffled_indices <- sample(nrow(all_individuals))
            perm_hist_indices <- shuffled_indices[1:n_hist]
            perm_mod_indices <- shuffled_indices[(n_hist + 1):(n_hist + n_mod)]
            
            # Calculate divergence for permuted groups
            perm_hist_alleles <- c(all_individuals[perm_hist_indices, gene_info[1]],
                                   all_individuals[perm_hist_indices, gene_info[2]])
            perm_mod_alleles <- c(all_individuals[perm_mod_indices, gene_info[1]],
                                  all_individuals[perm_mod_indices, gene_info[2]])
            
            perm_hist_div <- calculate_real_expected_mean_divergence(perm_hist_alleles, sequences, gene_info[3])
            perm_mod_div <- calculate_real_expected_mean_divergence(perm_mod_alleles, sequences, gene_info[3])
            
            if (!is.na(perm_hist_div) && !is.na(perm_mod_div)) {
              return(perm_mod_div - perm_hist_div)
            } else {
              return(NA)
            }
          })
          
          # Remove NA values from permutation results
          perm_diffs <- perm_diffs[!is.na(perm_diffs)]
          
          if (length(perm_diffs) > 100) {  # Need sufficient permutations (lowered threshold)
            # Calculate p-value (two-tailed)
            p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
            
            # Ensure p-value is not exactly 0
            if (p_value == 0) {
              p_value <- 1 / length(perm_diffs)
            }
            
            # Effect size (using SD of permutation distribution)
            effect_size <- observed_diff / sd(perm_diffs)
            
            cat(sprintf("  Observed diff: %.4f, p-value: %.6f, effect size: %.3f\n",
                        observed_diff, p_value, effect_size))
            
            results <- rbind(results, data.frame(
              Population = pop,
              Gene = gene_info[3],
              Historical_Divergence = round(hist_div, 4),
              Modern_Divergence = round(mod_div, 4),
              Difference = round(observed_diff, 4),
              p_value = round(p_value, 6),
              Significant = p_value < 0.05,
              Effect_size = round(effect_size, 3),
              n_historical_individuals = n_hist,
              n_modern_individuals = n_mod,
              n_permutations = length(perm_diffs),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  return(results)
}

# Perform proper statistical tests
proper_divergence_test_results <- perform_divergence_statistical_tests(combined_data_clean, all_sequences)

cat("\n=== PROPER DIVERGENCE STATISTICAL TEST RESULTS ===\n")
print(proper_divergence_test_results)

### 4. COMBINE RESULTS FOR PLOTTING ###

# Merge real divergence values with bootstrap CIs and statistical results
final_divergence_results <- real_divergence_results %>%
  left_join(bootstrap_div_ci_results, by = c("Population", "Time_Period", "Gene")) %>%
  left_join(proper_divergence_test_results %>% dplyr::select(Population, Gene, Significant, p_value),
            by = c("Population", "Gene"))

cat("\n=== FINAL COMBINED DIVERGENCE RESULTS ===\n")
print(final_divergence_results)

### 5. CREATE PROPER VISUALIZATION - EXACT SAME STYLE AS PANEL A ###

create_proper_divergence_plot <- function(plot_data, test_results) {
  
  # Prepare data for plotting - consistent with the he_plot
  plotting_data <- plot_data %>%
    mutate(
      # REVISION 1: Use bolditalic() to make gene names bold and italic
      Gene = case_when(
        Gene == "Orcu-U1" ~ "italic(Orcu)*'-'*italic(U1)",
        Gene == "Orcu-U2" ~ "italic(Orcu)*'-'*italic(U2)"
      ),
      # Use "Time" as the column name for consistency
      Time = factor(Time_Period, levels = c("Historical", "Modern")),
      Population = factor(Population, levels = c("France", "Britain", "Australia"))
    )
  
  # Add significance annotations - consistent with the he_plot
  if (nrow(test_results) > 0) {
    significance_data <- test_results %>%
      mutate(
        significance = case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01  ~ "**",
          p_value < 0.05  ~ "*",
          TRUE ~ ""
        ),
        # REVISION 1 (repeated): Use bolditalic()
        Gene = case_when(
          Gene == "Orcu-U1" ~ "italic(Orcu)*'-'*italic(U1)",
          Gene == "Orcu-U2" ~ "italic(Orcu)*'-'*italic(U2)"
        ),
        Population = factor(Population, levels = c("France", "Britain", "Australia"))
      ) %>%
      filter(significance != "")
    
    # Calculate bracket positions (values adjusted for divergence scale)
    if (nrow(significance_data) > 0) {
      max_heights <- plotting_data %>%
        group_by(Population, Gene) %>%
        summarise(max_height = max(CI_upper, na.rm = TRUE), .groups = "drop")
      
      significance_data <- significance_data %>%
        left_join(max_heights, by = c("Population", "Gene")) %>%
        mutate(
          bracket_y = max_height + 5, # Adjusted for divergence data scale
          text_y = bracket_y + 2    # Adjusted for divergence data scale
        )
    }
  } else {
    significance_data <- data.frame()
  }
  
  # Create the plot using styles from the he_plot
  p <- ggplot(plotting_data, aes(x = Time, y = Real_Mean_Divergence, fill = Time)) +
    geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                  width = 0.2, position = position_dodge(0.7)) +
    geom_text(aes(y = 0, label = paste0("n=", Sample_size)),
              vjust = 1.2, 
              size = 4, # REVISION 2: Matched sample size font to he_plot
              position = position_dodge(0.7)) +
    
    # Add significance brackets if any
    {if (nrow(significance_data) > 0) {
      list(
        geom_segment(data = significance_data,
                     aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y),
                     inherit.aes = FALSE, color = "black", linewidth = 0.5),
        geom_segment(data = significance_data,
                     aes(x = 1, xend = 1, y = bracket_y, yend = bracket_y - 0.5), # Adjusted tick size
                     inherit.aes = FALSE, color = "black", linewidth = 0.5),
        geom_segment(data = significance_data,
                     aes(x = 2, xend = 2, y = bracket_y, yend = bracket_y - 0.5), # Adjusted tick size
                     inherit.aes = FALSE, color = "black", linewidth = 0.5),
        geom_text(data = significance_data,
                  aes(x = 1.5, y = text_y, label = significance),
                  inherit.aes = FALSE, 
                  size = 6, # REVISION 3: Matched significance star size to he_plot
                  vjust = 0.5)
      )
    }} +
    
    facet_grid(Gene ~ Population,
               labeller = labeller(Gene = label_parsed, Population = label_value),
               scales = "free") +
    scale_fill_manual(values = c("Historical" = "#d73027", "Modern" = "#4575b4")) +
    labs(
      x = "Time Period",
      y = "Expected Mean Allele Divergence", # Updated y-axis title
      fill = "Time Period"
    ) +
    # REVISION 4: COPIED ENTIRE THEME BLOCK FROM he_plot FOR CONSISTENCY
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x = element_text(face = "bold", size = 12),
      strip.text.y = element_text(size = 13, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(
        size = 14,
        margin = margin(r = 15)
      ),
      legend.position = "none",
      panel.spacing.y = unit(1.05, "lines"),
      panel.spacing.x = unit(0.5, "lines"),
      plot.caption = element_text(size = 12, color = "gray50", hjust = 0.5)
    ) +
    # REVISION 5: Updated y-axis scale to start at 0 and add padding at the top
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1)))
  
  return(p)
}

# Create the proper plot - EXACT same style as Panel A
proper_divergence_plot <- create_proper_divergence_plot(final_divergence_results, proper_divergence_test_results)

print(proper_divergence_plot)

# BIOLOGICAL INTERPRETATION - EXACT same style as Panel A
cat("\n=== BIOLOGICAL INTERPRETATION ===\n")
significant_div_results <- proper_divergence_test_results %>% filter(Significant == TRUE)

if (nrow(significant_div_results) > 0) {
  cat("SIGNIFICANT CHANGES IN EXPECTED MEAN DIVERGENCE:\n")
  
  decreases <- sum(significant_div_results$Difference < 0)
  increases <- sum(significant_div_results$Difference > 0)
  
  for (i in 1:nrow(significant_div_results)) {
    row <- significant_div_results[i,]
    direction <- ifelse(row$Difference < 0, "DECREASE", "INCREASE")
    
    cat(sprintf("%s in %s: %.4f → %.4f (REAL %s of %.4f, p=%.6f, ES=%.3f)\n",
                row$Gene, row$Population,
                row$Historical_Divergence, row$Modern_Divergence,
                direction, abs(row$Difference), row$p_value, row$Effect_size))
  }
  
  cat(sprintf("\nOVERALL PATTERN: %d decreases, %d increases\n", decreases, increases))
} else {
  cat("No significant changes detected in expected mean divergence.\n")
}

# Save results - EXACT same as Panel A
panel_b_output_dir <- paste0(cwd, "/Divergent Alleles")
if (!dir.exists(panel_b_output_dir)) {
  dir.create(panel_b_output_dir, recursive = TRUE)
}

# Save results
write.csv(proper_divergence_test_results, paste0(panel_b_output_dir, "/proper_divergence_statistical_tests.csv"), row.names = FALSE)
write.csv(final_divergence_results, paste0(panel_b_output_dir, "/final_divergence_results.csv"), row.names = FALSE)

# Save plot
ggsave(paste0(panel_b_output_dir, "/proper_divergence_plot.pdf"), proper_divergence_plot, width = 12, height = 8)
ggsave(paste0(panel_b_output_dir, "/proper_divergence_plot.png"), proper_divergence_plot, width = 12, height = 8, dpi = 300)



cat("\n=== PANEL B ANALYSIS COMPLETE ===\n")


######################################################################################################
######################################################################################################

### PANEL C: ARE SELECTED ALLELES ESPECIALLY DIFFERENT FROM COMMON ALLELES? ###

cat("=== PANEL C: TESTING IF SELECTED ALLELES ARE ESPECIALLY DIFFERENT ===\n")

# Define the selected alleles
positively_selected_alleles <- c("Orcu-U1*02:01", "Orcu-U2*01:01")
cat(sprintf("Selected alleles: %s\n", paste(positively_selected_alleles, collapse = ", ")))

# Pool modern and historical data
pooled_frequencies <- all_frequencies %>%
  group_by(Population, Gene, Allele, Allele_Number) %>%
  summarise(
    Total_Count = sum(Count),
    Total_Alleles = sum(Total_Alleles),
    Pooled_Frequency = Total_Count / Total_Alleles,
    .groups = "drop"
  ) %>%
  filter(Allele %in% names(all_sequences)) %>%
  filter(Pooled_Frequency >= 0.01)  # Include alleles ≥1% frequency

cat(sprintf("Pooled data: %d alleles across populations\n", nrow(pooled_frequencies)))

# Function to calculate average distance to common alleles
calculate_distance_to_common_alleles <- function(data, sequences, common_threshold = 0.05) {
  data$Avg_Distance_To_Common <- NA
  
  for (i in 1:nrow(data)) {
    row <- data[i, ]
    
    # Get common alleles for this population-gene combination (≥5% frequency)
    common_alleles <- data %>%
      filter(Population == row$Population, Gene == row$Gene) %>%
      filter(Pooled_Frequency >= common_threshold) %>%
      filter(Allele != row$Allele) %>%  # Exclude the target allele itself
      pull(Allele)
    
    if (length(common_alleles) > 0) {
      distances <- c()
      
      # Calculate distance to each common allele
      for (common_allele in common_alleles) {
        if (row$Allele %in% names(sequences) && common_allele %in% names(sequences)) {
          dist <- calculate_grantham_distance(
            sequences[[row$Allele]], 
            sequences[[common_allele]]
          )
          
          if (!is.na(dist)) {
            distances <- c(distances, dist)
          }
        }
      }
      
      if (length(distances) > 0) {
        data$Avg_Distance_To_Common[i] <- mean(distances)
      }
    }
  }
  
  return(data)
}

# Calculate distance to common alleles
analysis_data <- calculate_distance_to_common_alleles(pooled_frequencies, all_sequences)

# Filter complete data and add log frequency
analysis_complete <- analysis_data %>%
  filter(!is.na(Avg_Distance_To_Common)) %>%
  mutate(
    Log_Frequency = log10(Pooled_Frequency),
    Is_Selected = Allele %in% positively_selected_alleles
  )

cat(sprintf("Analysis complete: %d alleles with distance data\n", nrow(analysis_complete)))

# Calculate statistics for each population-gene combination
calculate_stats <- function(data) {
  data %>%
    group_by(Population, Gene) %>%
    summarise(
      n = n(),
      correlation = ifelse(n >= 3, cor(Avg_Distance_To_Common, Log_Frequency, 
                                       method = "spearman", use = "complete.obs"), NA),
      p_value = ifelse(n >= 3, {
        test_result <- try(cor.test(Avg_Distance_To_Common, Log_Frequency, 
                                    method = "spearman"), silent = TRUE)
        if (class(test_result) == "try-error") NA else test_result$p.value
      }, NA),
      .groups = "drop"
    ) %>%
    mutate(
      significance = case_when(
        is.na(p_value) ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**", 
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      stat_label = case_when(
        is.na(correlation) ~ paste0("n=", n),
        n < 3 ~ paste0("n=", n),
        TRUE ~ paste0("ρ=", round(correlation, 3), significance, "\nn=", n)
      )

    )
}

stats_results <- calculate_stats(analysis_complete)

# Set factor levels for proper ordering
analysis_complete$Population <- factor(analysis_complete$Population, 
                                       levels = c("France", "Britain", "Australia"))
analysis_complete$Gene <- factor(analysis_complete$Gene, 
                                 levels = c("Orcu-U1", "Orcu-U2"))
stats_results$Population <- factor(stats_results$Population, 
                                   levels = c("France", "Britain", "Australia"))
stats_results$Gene <- factor(stats_results$Gene, 
                             levels = c("Orcu-U1", "Orcu-U2"))

# Create the plot - 2x3 facet grid (genes as rows, populations as columns)
panel_c_plot <- ggplot(analysis_complete, aes(x = Avg_Distance_To_Common, y = Log_Frequency)) +
  geom_point(aes(color = Is_Selected), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, linewidth = 1.2, color = "black") +
  # Highlight selected alleles with larger points
  geom_point(data = analysis_complete %>% filter(Is_Selected == TRUE),
             aes(color = Is_Selected), size = 2, stroke = 2) +
  # Add correlation statistics for each facet
  geom_text(data = stats_results,
            aes(x = Inf, y = Inf, label = stat_label),
            hjust = 1.05, vjust = 1.2,
            size = 3.5, fontface = "bold", color = "black",
            inherit.aes = FALSE) +
  facet_grid(Gene ~ Population, scales = "free") +
  scale_color_manual(
    values = c("TRUE" = "#FF6B6B", "FALSE" = "gray60"),
    name = "Allele Type",
    labels = c("TRUE" = "Selected", "FALSE" = "Other")
  ) +
  labs(
    x = "Average Grantham pairwise divergence",
    y = "Log10(Allele Frequency)",
    # title = "Are Selected Alleles Especially Different from Common Alleles?",
    # caption = paste0("Selected alleles: ", paste(positively_selected_alleles, collapse = ", "), 
    #                  " (highlighted in red)\nSpearman correlation coefficients (ρ) shown for each gene-population combination\n",
    #                  "*p<0.05, **p<0.01, ***p<0.001 | Positive correlation = more divergent alleles have higher frequency")
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.text.x = element_text(size = 12, face = "bold"),  # Population labels
    strip.text.y = element_text(size = 12, face = "bold.italic"),  # Gene labels  
    strip.background = element_rect(fill = "grey90", color = "black"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 0.5),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 15, 10, 10)
  )

# Display the plot
print(panel_c_plot)

# Save the plot
ggsave(paste0(panel_b_output_dir,"/Panel_C_Selected_vs_Common_6Facets.pdf"), panel_c_plot,
       device = "pdf", width = 15, height = 10)

# Print detailed results for selected alleles by gene and population
cat("\n=== SELECTED ALLELES DETAILED RESULTS ===\n")
selected_results <- analysis_complete %>%
  filter(Is_Selected == TRUE) %>%
  arrange(Gene, Population)

if (nrow(selected_results) > 0) {
  print(selected_results %>% 
          select(Gene, Population, Allele, Pooled_Frequency, Avg_Distance_To_Common))
  
  cat("\nSummary by Gene:\n")
  summary_stats <- selected_results %>%
    group_by(Gene) %>%
    summarise(
      Mean_Distance = mean(Avg_Distance_To_Common, na.rm = TRUE),
      SD_Distance = sd(Avg_Distance_To_Common, na.rm = TRUE),
      Min_Distance = min(Avg_Distance_To_Common, na.rm = TRUE),
      Max_Distance = max(Avg_Distance_To_Common, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_stats)
}

# Print statistics summary with gene-population breakdown
cat("\n=== CORRELATION STATISTICS BY GENE AND POPULATION ===\n")
print(stats_results %>% select(Gene, Population, n, correlation, p_value, significance))

# Test: Are selected alleles significantly more divergent than others? (by gene-population)
cat("\n=== STATISTICAL TESTS: SELECTED vs OTHER ALLELES (BY GENE-POPULATION) ===\n")

comparison_results <- data.frame()

for (gene in unique(analysis_complete$Gene)) {
  for (pop in unique(analysis_complete$Population)) {
    subset_data <- analysis_complete %>% filter(Gene == gene, Population == pop)
    
    if (nrow(subset_data) > 0) {
      selected_divergence <- subset_data %>% filter(Is_Selected == TRUE) %>% pull(Avg_Distance_To_Common)
      other_divergence <- subset_data %>% filter(Is_Selected == FALSE) %>% pull(Avg_Distance_To_Common)
      
      if (length(selected_divergence) > 0 && length(other_divergence) > 1) {
        # Wilcoxon test (non-parametric)
        test_result <- try(wilcox.test(selected_divergence, other_divergence, alternative = "greater"), silent = TRUE)
        
        if (class(test_result) != "try-error") {
          comparison_results <- rbind(comparison_results, data.frame(
            Gene = gene,
            Population = pop,
            Selected_Mean = mean(selected_divergence, na.rm = TRUE),
            Other_Mean = mean(other_divergence, na.rm = TRUE),
            Difference = mean(selected_divergence, na.rm = TRUE) - mean(other_divergence, na.rm = TRUE),
            P_Value = test_result$p.value,
            N_Selected = length(selected_divergence),
            N_Other = length(other_divergence),
            Test_Statistic = test_result$statistic,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}

if (nrow(comparison_results) > 0) {
  comparison_results <- comparison_results %>%
    mutate(
      Significance = case_when(
        P_Value < 0.001 ~ "***",
        P_Value < 0.01 ~ "**",
        P_Value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  cat("Wilcoxon test results (testing if selected alleles are MORE divergent):\n")
  print(comparison_results %>% arrange(Gene, Population))
}

# Save results
write.csv(analysis_complete, paste0(panel_b_output_dir, "/selected_alleles_divergence_analysis.csv"), row.names = FALSE)
write.csv(stats_results, paste0(panel_b_output_dir, "/correlation_statistics.csv"), row.names = FALSE)



######################################################################################################
######################################################################################################

# Load the heterozygosity plot
load(paste0(cwd, "/Expected Heterozygosity/proper_he_plot.RData"))


# Load required library for panel creation
library(patchwork)

# # Create right side panel first
# upper_panel <- proper_he_plot | proper_divergence_plot
# 
# # Create the full panel
# panel_plot <- upper_panel / panel_c_plot
# panel_plot <- panel_plot + 
#   plot_layout(heights = c(2, 1)) + 
#   plot_annotation(tag_levels = c('a', 'b', 'c')) &
#   theme(plot.tag = element_text(face = "bold", size = 16))

supp_plot <- proper_divergence_plot / panel_c_plot +
  plot_layout(heights = c(1.5, 1)) + 
  plot_annotation(tag_levels = c('a', 'b')) &
  theme(plot.tag = element_text(face = "bold", size = 16))

# Display the panel
print(panel_plot)

# Save the panel plot
ggsave(paste0(cwd, "/Combined_Analysis_Panel_main.pdf"), proper_he_plot,
       device = cairo_pdf, width = 6, height = 8)

ggsave(paste0(cwd, "/Combined_Analysis_Panel_supp.pdf"), supp_plot,
       device = cairo_pdf, width = 12, height = 8)


# 
# ######################################################################################################
# ######################################################################################################
# ########################################################################
# ### PHYLOGENY ANALYSIS: Orcu-U1 AMINO ACID SEQUENCES ###
# ########################################################################
# cat("\n=== PHYLOGENY ANALYSIS: Orcu-U1 AMINO ACID SEQUENCES ===\n")
# 
# # Load required libraries for phylogeny
# if (!require(ape)) install.packages("ape")
# if (!require(phangorn)) install.packages("phangorn")
# if (!require(ggtree)) {
#   if (!require(BiocManager)) install.packages("BiocManager")
#   BiocManager::install("ggtree")
# }
# 
# library(ape)
# library(phangorn)
# library(ggtree)
# 
# # Extract Orcu-U1 sequences only
# orcu_u1_sequences <- all_sequences[grepl("Orcu-U1", names(all_sequences))]
# cat(sprintf("Found %d Orcu-U1 amino acid sequences for phylogeny\n", length(orcu_u1_sequences)))
# 
# # Check if target allele exists
# target_allele <- "Orcu-U1*20:01"
# if (!target_allele %in% names(orcu_u1_sequences)) {
#   cat("Warning: Target allele Orcu-U1*20:01 not found in sequences!\n")
#   cat("Available Orcu-U1 alleles:\n")
#   print(names(orcu_u1_sequences))
# } else {
#   cat(sprintf("Target allele %s found in dataset\n", target_allele))
# }
# 
# # Create alignment if we have enough sequences
# if (length(orcu_u1_sequences) >= 3) {
#   
#   # Load Biostrings first
#   if (!require(Biostrings)) {
#     if (!require(BiocManager)) install.packages("BiocManager")
#     BiocManager::install("Biostrings")
#   }
#   library(Biostrings)
#   
#   # Create AAStringSet properly
#   aa_sequences_vector <- unlist(orcu_u1_sequences)
#   aa_stringset <- AAStringSet(aa_sequences_vector)
#   names(aa_stringset) <- names(orcu_u1_sequences)
#   
#   # Check sequence lengths and clean if necessary
#   cat("Checking sequence integrity...\n")
#   seq_lengths <- nchar(aa_stringset)
#   cat(sprintf("Sequence lengths: min=%d, max=%d\n", min(seq_lengths), max(seq_lengths)))
#   
#   # Remove sequences that are too short or contain too many X's
#   valid_seqs <- seq_lengths > 50 & 
#     sapply(as.character(aa_stringset), function(x) {
#       x_count <- sum(unlist(strsplit(x, "")) == "X")
#       return(x_count < (nchar(x) * 0.1))
#     })
#   
#   if (sum(valid_seqs) < length(aa_stringset)) {
#     cat(sprintf("Filtering out %d sequences with length < 50 or >10%% X residues\n", 
#                 sum(!valid_seqs)))
#     aa_stringset <- aa_stringset[valid_seqs]
#     orcu_u1_sequences <- orcu_u1_sequences[valid_seqs]  # Keep both in sync
#   }
#   
#   # Load msa package for alignment
#   if (!require(msa)) {
#     if (!require(BiocManager)) install.packages("BiocManager")
#     BiocManager::install("msa")
#   }
#   library(msa)
#   
#   cat("Performing multiple sequence alignment...\n")
#   
#   # Try ClustalW alignment with error handling
#   alignment <- tryCatch({
#     msaClustalW(aa_stringset)
#   }, error = function(e) {
#     cat("ClustalW failed, trying ClustalOmega...\n")
#     tryCatch({
#       msaClustalOmega(aa_stringset)
#     }, error = function(e2) {
#       cat("ClustalOmega failed, trying Muscle...\n")
#       tryCatch({
#         msaMuscle(aa_stringset)
#       }, error = function(e3) {
#         cat("All alignment methods failed. Using simple concatenation...\n")
#         return(NULL)
#       })
#     })
#   })
#   
#   # Handle alignment result and create distance matrix
#   if (!is.null(alignment)) {
#     # Convert alignment to matrix
#     alignment_matrix <- as.matrix(alignment)
#     
#     # Convert to phyDat format for phylogenetic analysis
#     phyDat_seqs <- phyDat(alignment_matrix, type = "AA")
#     
#     # Calculate distance matrix
#     cat("Calculating amino acid distance matrix...\n")
#     dist_matrix <- dist.ml(phyDat_seqs, model = "JTT")  # JTT model good for amino acids
#     
#   } else {
#     # Fallback: calculate pairwise distances without alignment
#     cat("Using pairwise distance calculation without formal alignment...\n")
#     
#     # Use the original character sequences for distance calculation
#     seq_names <- names(orcu_u1_sequences)
#     n_seqs <- length(orcu_u1_sequences)
#     dist_matrix_raw <- matrix(0, nrow = n_seqs, ncol = n_seqs)
#     rownames(dist_matrix_raw) <- colnames(dist_matrix_raw) <- seq_names
#     
#     for (i in 1:(n_seqs-1)) {
#       for (j in (i+1):n_seqs) {
#         # Calculate simple amino acid distance using original sequences
#         seq1 <- orcu_u1_sequences[[i]]
#         seq2 <- orcu_u1_sequences[[j]]
#         
#         # Use your existing Grantham distance function
#         grantham_dist <- calculate_grantham_distance(seq1, seq2)
#         
#         if (!is.na(grantham_dist)) {
#           dist_matrix_raw[i, j] <- dist_matrix_raw[j, i] <- grantham_dist
#         } else {
#           # Fallback: simple amino acid difference proportion
#           min_len <- min(nchar(seq1), nchar(seq2))
#           if (min_len > 0) {
#             seq1_chars <- substr(seq1, 1, min_len)
#             seq2_chars <- substr(seq2, 1, min_len)
#             differences <- sum(unlist(strsplit(seq1_chars, "")) != unlist(strsplit(seq2_chars, "")))
#             dist_matrix_raw[i, j] <- dist_matrix_raw[j, i] <- differences / min_len
#           } else {
#             dist_matrix_raw[i, j] <- dist_matrix_raw[j, i] <- 1.0
#           }
#         }
#       }
#     }
#     
#     # Convert to dist object
#     dist_matrix <- as.dist(dist_matrix_raw)
#   }
#   
#   # Build neighbor-joining tree with bootstrap support
# 
#   cat("Building neighbor-joining tree with 1000 bootstrap replicates...\n")
#   nj_tree <- nj(dist_matrix)
#   
#   # Calculate bootstrap support
#   cat("Calculating bootstrap support (this may take a few minutes)...\n")
#   # Convert alignment to matrix format for boot.phylo
#   alignment_matrix_for_boot <- as.matrix(alignment)
#   
#   bootstrap_trees <- boot.phylo(nj_tree, alignment_matrix_for_boot,
#                                 function(x) {
#                                   # Convert matrix back to phyDat and calculate tree
#                                   phyDat_boot <- phyDat(x, type = "AA")
#                                   dist_boot <- dist.ml(phyDat_boot, model = "JTT")
#                                   return(nj(dist_boot))
#                                 },
#                                 B = 1000, quiet = TRUE)
#   
#   # Add bootstrap values to tree (convert raw counts to percentages)
#   nj_tree$node.label <- round((bootstrap_trees / 10))  # Convert to percentages
#   
#   # Add bootstrap values to tree
#   # Debug bootstrap values first
#   cat("Bootstrap values range:", range(bootstrap_trees, na.rm = TRUE), "\n")
#   
#   
#   
#   # Root the tree (using midpoint rooting)
#   nj_tree <- midpoint(nj_tree)
#   # Create phylogeny plot with ggtree - rectangular layout only
#   cat("Creating phylogeny visualization...\n")
#   
#   phylo_plot_rect <- ggtree(nj_tree) +
#     geom_tiplab(aes(color = label == target_allele),
#                 size = 3.5,
#                 hjust = -0.1,
#                 fontface = "italic") +
#     geom_tippoint(aes(color = label == target_allele,
#                       size = label == target_allele),
#                   alpha = 0.8) +
#     geom_nodelab(aes(label = ifelse(as.numeric(label) >= 40, label, "")),
#                  size = 2.5, hjust = 1.4, vjust = -0.3, color = "black") +  # Back to ≥70
#     scale_color_manual(values = c("FALSE" = "black", "TRUE" = "#FF0000"),
#                        name = "Allele Type",
#                        labels = c("FALSE" = "Other", "TRUE" = target_allele)) +
#     scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4),
#                       guide = "none") +
#     theme_tree() +
#     theme(
#       legend.position = "none",
#       legend.title = element_text(size = 11, face = "bold"),
#       plot.title = element_text(size = 14, face = "bold"),
#       plot.margin = margin(20, 20, 20, 20)
#     ) +
#     xlim(0, max(node.depth.edgelength(nj_tree)) * 1.4)
#   
#   # Display plot
#   print(phylo_plot_rect)
#   
#   # Calculate some basic phylogenetic statistics
#   cat("\n=== PHYLOGENETIC STATISTICS ===\n")
#   cat(sprintf("Number of sequences: %d\n", length(nj_tree$tip.label)))
#   cat(sprintf("Tree length: %.4f\n", sum(nj_tree$edge.length)))
#   
#   # Find closest relatives to target allele
#   if (target_allele %in% nj_tree$tip.label) {
#     target_distances <- cophenetic(nj_tree)[target_allele, ]
#     target_distances <- target_distances[target_distances > 0]  # Remove self
#     target_distances <- sort(target_distances)
#     
#     cat(sprintf("\nClosest relatives to %s:\n", target_allele))
#     closest_5 <- head(target_distances, 5)
#     for (i in 1:length(closest_5)) {
#       cat(sprintf("  %d. %s (distance: %.4f)\n", i, names(closest_5)[i], closest_5[i]))
#     }
#     
#     # Calculate average distance to all other alleles
#     avg_distance <- mean(target_distances)
#     cat(sprintf("\nAverage distance from %s to all other alleles: %.4f\n", 
#                 target_allele, avg_distance))
#   }
#   
#   # Save phylogeny results
#   phylo_output_dir <- paste0(cwd, "/Phylogeny_Analysis")
#   if (!dir.exists(phylo_output_dir)) {
#     dir.create(phylo_output_dir, recursive = TRUE)
#   }
#   
#   # Save single rectangular plot
#   ggsave(paste0(phylo_output_dir, "/Orcu_U1_phylogeny.pdf"), 
#          phylo_plot_rect, width = 12, height = 8, device = "pdf")
#   
#   # Save tree in Newick format
#   write.tree(nj_tree, paste0(phylo_output_dir, "/Orcu_U1_tree.newick"))
#   
#   # Save distance matrix
#   write.csv(as.matrix(dist_matrix), paste0(phylo_output_dir, "/Orcu_U1_distance_matrix.csv"))
#   
#   # Save alignment only if it exists
#   if (!is.null(alignment)) {
#     tryCatch({
#       writeXStringSet(unmasked(alignment), paste0(phylo_output_dir, "/Orcu_U1_alignment.fasta"))
#     }, error = function(e) {
#       cat("Could not save alignment, saving original sequences instead\n")
#       writeXStringSet(aa_stringset, paste0(phylo_output_dir, "/Orcu_U1_sequences.fasta"))
#     })
#   } else {
#     writeXStringSet(aa_stringset, paste0(phylo_output_dir, "/Orcu_U1_sequences.fasta"))
#   }
#   
#   cat(sprintf("\nPhylogeny analysis complete. Results saved to: %s\n", phylo_output_dir))
#   
#   # Store phylo_plot for potential inclusion in combined figure
#   phylo_plot_for_panel <- phylo_plot_rect
#   
# } else {
#   cat("Insufficient Orcu-U1 sequences for phylogenetic analysis (need ≥3)\n")
#   phylo_plot_for_panel <- NULL
# }
# 
# cat("\n=== PHYLOGENY ANALYSIS SECTION COMPLETE ===\n")


