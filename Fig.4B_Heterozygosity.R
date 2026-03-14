# PROPER Bootstrap Expected Heterozygosity Analysis
# CORRECT APPROACH:
# - Bar heights = Real He values from observed data
# - Bootstrap = Standard errors only
# - Statistical testing = On REAL individual data (not bootstrap samples)

library(ggplot2)
library(tidyr)
library(dplyr)

# Set working directory (modify as needed)
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
cwd <- dirname(rstudioapi::getSourceEditorContext()$path)

# File paths - UPDATE THESE TO YOUR ACTUAL PATHS
wd <- '/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim2/04_museum_typing/mm97_b0010_combined_paired/'
norm_file <- '/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim2/04_museum_typing/norm_cov_combined.txt'
miseq_file <- '/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim1/11_Miseq/FINAL_Miseq_alleles.csv'

cat("=== PROPER BOOTSTRAP EXPECTED HETEROZYGOSITY ANALYSIS ===\n")
cat("CORRECT METHOD:\n")
cat("- Bar heights: Real He from observed data\n")
cat("- Error bars: Bootstrap standard errors\n")
cat("- Statistical tests: On real individual data\n\n")

### DATA LOADING SECTION (same as before) ###
setwd(paste0(wd, '/coverage_metrics'))

# Load data (using your existing data loading code - abbreviated for clarity)
miseq_table <- read.csv(file = miseq_file, header = T)
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
  # For 11016 (Orcu-U2)
  alleles <- revise_alleles(modern_table[i, "Orcu.2.allele1"], modern_table[i, "Orcu.2.allele2"], c(7))
  modern_table[i, c("Orcu.2.allele1", "Orcu.2.allele2")] <- alleles

  # For 10972b (Orcu-U1)
  alleles <- revise_alleles(modern_table[i, "Orcu.1.allele1"], modern_table[i, "Orcu.1.allele2"], c(11))
  modern_table[i, c("Orcu.1.allele1", "Orcu.1.allele2")] <- alleles
}

# Rename columns to match historical table
colnames(modern_table) <- c('Sample','11016.A1', '11016.A2', '10972b.A1', '10972b.A2')

# Load historical allele calls (abbreviated - use your existing code)
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

# Prepare datasets with population info
historical_data <- historical_table %>%
  add_population_info() %>%
  mutate(Time = "Historical")

modern_data <- modern_table %>%
  add_population_info() %>%
  mutate(Time = "Modern")

# Combine datasets
combined_data <- rbind(
  historical_data %>% select(Sample, Population, Time, `11016.A1`, `11016.A2`, `10972b.A1`, `10972b.A2`),
  modern_data %>% select(Sample, Population, Time, `11016.A1`, `11016.A2`, `10972b.A1`, `10972b.A2`)
)

# Clean data - remove missing values
combined_clean <- combined_data %>%
  filter(
    `11016.A1` != "*" & `11016.A2` != "*" &
      `10972b.A1` != "*" & `10972b.A2` != "*" &
      `11016.A1` != "-" & `11016.A2` != "-" &
      `10972b.A1` != "-" & `10972b.A2` != "-" &
      !is.na(`11016.A1`) & !is.na(`11016.A2`) &
      !is.na(`10972b.A1`) & !is.na(`10972b.A2`)
  )

cat(sprintf("Clean dataset: %d individuals\n", nrow(combined_clean)))

# Print sample sizes per group
sample_sizes <- combined_clean %>%
  group_by(Population, Time) %>%
  summarise(n = n(), .groups = 'drop')
cat("\nSample sizes per group:\n")
print(sample_sizes)

### 1. CALCULATE REAL EXPECTED HETEROZYGOSITY VALUES ###

# Function to calculate expected heterozygosity (He) from real data
calculate_real_expected_heterozygosity <- function(allele_vector) {
  # Remove missing values
  alleles <- allele_vector[!is.na(allele_vector) & allele_vector != "*" & allele_vector != "-"]

  if(length(alleles) < 4) return(NA)  # Need minimum sample size

  # Calculate allele frequencies
  allele_freqs <- table(alleles) / length(alleles)

  # Calculate expected heterozygosity: He = 1 - sum(pi^2)
  he <- 1 - sum(allele_freqs^2)

  return(he)
}

# Calculate REAL He values for each population-time-gene combination
cat("\n=== STEP 1: CALCULATING REAL He VALUES ===\n")

real_he_results <- data.frame()

for (pop in unique(combined_clean$Population)) {
  for (time in unique(combined_clean$Time)) {
    for (gene_info in list(c("10972b.A1", "10972b.A2", "Orcu-U1"),
                           c("11016.A1", "11016.A2", "Orcu-U2"))) {

      subset_data <- combined_clean %>%
        filter(Population == pop, Time == time) %>%
        filter(!is.na(!!sym(gene_info[1])) & !is.na(!!sym(gene_info[2])) &
                 !!sym(gene_info[1]) != "*" & !!sym(gene_info[2]) != "*" &
                 !!sym(gene_info[1]) != "-" & !!sym(gene_info[2]) != "-")

      if (nrow(subset_data) >= 5) {
        # Get all alleles from the real data
        all_alleles <- c(subset_data[[gene_info[1]]], subset_data[[gene_info[2]]])

        # Calculate real He
        real_he <- calculate_real_expected_heterozygosity(all_alleles)

        if (!is.na(real_he)) {
          cat(sprintf("Real He: %s %s %s = %.4f (n=%d individuals)\n",
                      pop, time, gene_info[3], real_he, nrow(subset_data)))

          real_he_results <- rbind(real_he_results, data.frame(
            Population = pop,
            Time = time,
            Gene = gene_info[3],
            Real_He = real_he,
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

# Function to calculate bootstrap CIs
calculate_bootstrap_ci <- function(individuals_data, gene_cols, n_bootstrap = 1000) {
  cat(sprintf("Bootstrapping %s for 95%% confidence intervals...\n", gene_cols[3]))

  bootstrap_he_values <- numeric(n_bootstrap)

  for (b in 1:n_bootstrap) {
    boot_sample <- individuals_data[sample(nrow(individuals_data), replace = TRUE), ]
    boot_alleles <- c(boot_sample[[gene_cols[1]]], boot_sample[[gene_cols[2]]])
    bootstrap_he_values[b] <- calculate_real_expected_heterozygosity(boot_alleles)
  }

  # Remove NA values
  bootstrap_he_values <- bootstrap_he_values[!is.na(bootstrap_he_values)]

  if (length(bootstrap_he_values) > 0) {
    ci <- quantile(bootstrap_he_values, probs = c(0.025, 0.975), na.rm = TRUE)
    return(list(lower = ci[1], upper = ci[2]))
  } else {
    return(list(lower = NA, upper = NA))
  }
}


# Calculate bootstrap CIs for each group
bootstrap_ci_results <- data.frame()
for (pop in unique(combined_clean$Population)) {
  for (time in unique(combined_clean$Time)) {
    for (gene_info in list(c("10972b.A1", "10972b.A2", "Orcu-U1"),
                           c("11016.A1", "11016.A2", "Orcu-U2"))) {

      subset_data <- combined_clean %>%
        filter(Population == pop, Time == time) %>%
        filter(!is.na(!!sym(gene_info[1])) & !is.na(!!sym(gene_info[2])) &
                 !!sym(gene_info[1]) != "*" & !!sym(gene_info[2]) != "*" &
                 !!sym(gene_info[1]) != "-" & !!sym(gene_info[2]) != "-")

      if (nrow(subset_data) >= 5) {
        bootstrap_ci <- calculate_bootstrap_ci(subset_data, gene_info, n_bootstrap = 1000)

        if (!is.na(bootstrap_ci$lower) && !is.na(bootstrap_ci$upper)) {
          cat(sprintf("Bootstrap 95%% CI: %s %s %s = [%.4f, %.4f]\n",
                      pop, time, gene_info[3], bootstrap_ci$lower, bootstrap_ci$upper))

          bootstrap_ci_results <- rbind(bootstrap_ci_results, data.frame(
            Population = pop,
            Time = time,
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

# Function to perform proper statistical tests on individual-level data
perform_proper_statistical_tests <- function(individual_data) {
  results <- data.frame()

  for (pop in unique(individual_data$Population)) {
    for (gene_info in list(c("10972b.A1", "10972b.A2", "Orcu-U1"),
                           c("11016.A1", "11016.A2", "Orcu-U2"))) {

      # Get historical and modern data for this population-gene
      hist_data <- individual_data %>%
        filter(Population == pop, Time == "Historical") %>%
        filter(!is.na(!!sym(gene_info[1])) & !is.na(!!sym(gene_info[2])) &
                 !!sym(gene_info[1]) != "*" & !!sym(gene_info[2]) != "*" &
                 !!sym(gene_info[1]) != "-" & !!sym(gene_info[2]) != "-")

      mod_data <- individual_data %>%
        filter(Population == pop, Time == "Modern") %>%
        filter(!is.na(!!sym(gene_info[1])) & !is.na(!!sym(gene_info[2])) &
                 !!sym(gene_info[1]) != "*" & !!sym(gene_info[2]) != "*" &
                 !!sym(gene_info[1]) != "-" & !!sym(gene_info[2]) != "-")

      if (nrow(hist_data) >= 5 && nrow(mod_data) >= 5) {

        # Calculate He for each group
        hist_alleles <- c(hist_data[[gene_info[1]]], hist_data[[gene_info[2]]])
        mod_alleles <- c(mod_data[[gene_info[1]]], mod_data[[gene_info[2]]])

        hist_he <- calculate_real_expected_heterozygosity(hist_alleles)
        mod_he <- calculate_real_expected_heterozygosity(mod_alleles)

        if (!is.na(hist_he) && !is.na(mod_he)) {

          observed_diff <- mod_he - hist_he

          # PERMUTATION TEST on individual data
          # Combine all individuals from both time periods
          all_individuals <- rbind(
            hist_data %>% mutate(Original_Time = "Historical"),
            mod_data %>% mutate(Original_Time = "Modern")
          )

          n_hist <- nrow(hist_data)
          n_mod <- nrow(mod_data)
          n_permutations <- 10000

          cat(sprintf("Testing %s %s: n_hist=%d, n_mod=%d\n",
                      pop, gene_info[3], n_hist, n_mod))

          # Permutation test
          perm_diffs <- replicate(n_permutations, {
            # Randomly reassign time periods to individuals
            shuffled_indices <- sample(nrow(all_individuals))
            perm_hist_indices <- shuffled_indices[1:n_hist]
            perm_mod_indices <- shuffled_indices[(n_hist + 1):(n_hist + n_mod)]

            # Calculate He for permuted groups
            perm_hist_alleles <- c(all_individuals[perm_hist_indices, gene_info[1]],
                                   all_individuals[perm_hist_indices, gene_info[2]])
            perm_mod_alleles <- c(all_individuals[perm_mod_indices, gene_info[1]],
                                  all_individuals[perm_mod_indices, gene_info[2]])

            perm_hist_he <- calculate_real_expected_heterozygosity(perm_hist_alleles)
            perm_mod_he <- calculate_real_expected_heterozygosity(perm_mod_alleles)

            if (!is.na(perm_hist_he) && !is.na(perm_mod_he)) {
              return(perm_mod_he - perm_hist_he)
            } else {
              return(NA)
            }
          })

          # Remove NA values from permutation results
          perm_diffs <- perm_diffs[!is.na(perm_diffs)]

          if (length(perm_diffs) > 1000) {  # Need sufficient permutations
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
              Historical_He = round(hist_he, 4),
              Modern_He = round(mod_he, 4),
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
proper_test_results <- perform_proper_statistical_tests(combined_clean)

cat("\n=== PROPER STATISTICAL TEST RESULTS ===\n")
print(proper_test_results)

### 4. COMBINE RESULTS FOR PLOTTING ###

# Merge real He values with bootstrap SEs and statistical results
final_results <- real_he_results %>%
  left_join(bootstrap_ci_results, by = c("Population", "Time", "Gene")) %>%
  left_join(proper_test_results %>% select(Population, Gene, Significant, p_value),
            by = c("Population", "Gene"))

cat("\n=== FINAL COMBINED RESULTS ===\n")
print(final_results)



### 4. PLOTTING ###
create_proper_he_plot <- function(plot_data, test_results) {
  
  # Prepare data for plotting
  plotting_data <- plot_data %>%
    mutate(
      Gene = case_when(
        Gene == "Orcu-U1" ~ "italic(Orcu)*'-'*italic(U1)",
        Gene == "Orcu-U2" ~ "italic(Orcu)*'-'*italic(U2)"
      ),
      Time = factor(Time, levels = c("Historical", "Modern")),
      Population = factor(Population, levels = c("France", "Britain", "Australia"))
    )
  
  # Add significance annotations
  if (nrow(test_results) > 0) {
    significance_data <- test_results %>%
      mutate(
        significance = case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01  ~ "**",
          p_value < 0.05  ~ "*",
          TRUE ~ ""
        ),
        Gene = case_when(
          Gene == "Orcu-U1" ~ "italic(Orcu)*'-'*italic(U1)",
          Gene == "Orcu-U2" ~ "italic(Orcu)*'-'*italic(U2)"
        ),
        Population = factor(Population, levels = c("France", "Britain", "Australia"))
      ) %>%
      filter(significance != "")
    
    # Calculate bracket positions
    if (nrow(significance_data) > 0) {
      max_heights <- plotting_data %>%
        group_by(Population, Gene) %>%
        summarise(max_height = max(CI_upper, na.rm = TRUE), .groups = "drop")
      
      significance_data <- significance_data %>%
        left_join(max_heights, by = c("Population", "Gene")) %>%
        mutate(
          bracket_y = pmax(max_height + 0.05, 0.85),
          text_y = bracket_y + 0.03
        )
    }
  } else {
    significance_data <- data.frame()
  }
  
  # Create the plot
  p <- ggplot(plotting_data, aes(x = Time, y = Real_He, fill = Time)) +
    geom_col(position = position_dodge(width = 0.4), alpha = 0.8, width = 1) +
    geom_errorbar(aes(ymin = pmax(CI_lower, 0),
                      ymax = pmin(CI_upper, 1)),
                  width = 0.2, position = position_dodge(0.5)) +
    geom_text(aes(y = 0.03, label = Sample_size),
              vjust = 0.5, size = 3.5, color = "black",
              position = position_dodge(0.5)) +
    
    # Add significance brackets if any
    {if (nrow(significance_data) > 0) {
      list(
        geom_segment(data = significance_data,
                     aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y),
                     inherit.aes = FALSE, color = "black", linewidth = 0.5),
        geom_segment(data = significance_data,
                     aes(x = 1, xend = 1, y = bracket_y, yend = bracket_y - 0.01),
                     inherit.aes = FALSE, color = "black", linewidth = 0.5),
        geom_segment(data = significance_data,
                     aes(x = 2, xend = 2, y = bracket_y, yend = bracket_y - 0.01),
                     inherit.aes = FALSE, color = "black", linewidth = 0.5),
        geom_text(data = significance_data,
                  aes(x = 1.5, y = text_y, label = significance),
                  inherit.aes = FALSE, 
                  size = 6,
                  vjust = 0.5)
      )
    }} +
    
    facet_grid(Gene ~ Population,
               labeller = labeller(Gene = label_parsed, Population = label_value),
               scales = "free") +
    scale_fill_manual(values = c("Historical" = "#d73027", "Modern" = "#4575b4")) +
    labs(
      y = "Expected Frequency of Heterozygotes",
      fill = NULL
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x = element_text(face = "bold", size = 14),
      strip.text.y = element_text(size = 13, face = "bold"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 12), 
      axis.title.x = element_blank(),
      axis.title.y = element_text(
        size = 13,
        margin = margin(r = 15)
      ),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      panel.spacing.y = unit(1.05, "lines"),
      # CHANGED: Set horizontal panel spacing to 0 to remove the gap
      panel.spacing.x = unit(1, "lines"),
      plot.caption = element_text(size = 12, color = "gray50", hjust = 0.5)
    ) +
    scale_y_continuous(limits = c(0, 1.05), labels = scales::percent_format())
  
  return(p)
}




# Create the proper plot
proper_he_plot <- create_proper_he_plot(final_results, proper_test_results)
print(proper_he_plot)

### 6. BIOLOGICAL INTERPRETATION ###
cat("\n=== BIOLOGICAL INTERPRETATION ===\n")
significant_results <- proper_test_results %>% filter(Significant == TRUE)

if (nrow(significant_results) > 0) {
  cat("SIGNIFICANT CHANGES IN EXPECTED HETEROZYGOSITY:\n")

  decreases <- sum(significant_results$Difference < 0)
  increases <- sum(significant_results$Difference > 0)

  for (i in 1:nrow(significant_results)) {
    row <- significant_results[i,]
    direction <- ifelse(row$Difference < 0, "DECREASE", "INCREASE")

    cat(sprintf("%s in %s: %.4f → %.4f (REAL %s of %.4f, p=%.6f, ES=%.3f)\n",
                row$Gene, row$Population,
                row$Historical_He, row$Modern_He,
                direction, abs(row$Difference), row$p_value, row$Effect_size))
  }

  cat(sprintf("\nOVERALL PATTERN: %d decreases, %d increases\n", decreases, increases))

} else {
  cat("No significant changes detected with proper statistical testing.\n")
}

### SAVE RESULTS ###
output_dir <- paste0(cwd, "/Expected Heterozygosity")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save results
write.csv(proper_test_results, paste0(output_dir, "/proper_statistical_tests.csv"), row.names = FALSE)
write.csv(final_results, paste0(output_dir, "/final_results_with_real_he.csv"), row.names = FALSE)

# Save plot
ggsave(paste0(output_dir, "/proper_bootstrap_he_plot.pdf"), proper_he_plot, width = 4, height = 7)
ggsave(paste0(output_dir, "/proper_bootstrap_he_plot.png"), proper_he_plot, width = 4, height = 7, dpi = 300)

cat("\n=== PROPER ANALYSIS COMPLETE ===\n")


# Save the plot object
save(proper_he_plot, file = paste0(cwd, "/Expected Heterozygosity/proper_he_plot.RData"))

