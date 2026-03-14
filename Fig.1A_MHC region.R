library(ggplot2)
library(dplyr)
library(ggrepel)
library(readxl)
library(rtracklayer)
library(cowplot)
library(geneviewer)


rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

###################################################################################################
########################### Gene Map ##############################################################
###################################################################################################

# Read the GTF file
gtf_file <- '/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim1/04_PacBio2/GTF/Oryctolagus_cuniculus.OryCun2.0.103.gtf'
gtf_data <- rtracklayer::import(gtf_file)
gtf_df <- as.data.frame(gtf_data)

mhc1_genes_on_ensembl <- c("ENSOCUG00000010972",
                           "ENSOCUG00000011016",
                           "ENSOCUG00000013304",
                           "ENSOCUG00000006557",
                           "ENSOCUG00000010984",
                           "ENSOCUG00000008993")

# Filter for chromosome 12 and the region of interest
filtered_gtf <- gtf_df %>%
  filter(seqnames == "12" &
           type == "gene" &
           start >= 20000000 &
           end <= 24000000 &
           !gene_id %in% mhc1_genes_on_ensembl &
           # is.na(gene_name) == F &
           !grepl("HLA", gene_name)) %>%
  mutate(., Category = 'Unknown') %>%
  select(gene_name, strand, start, end, Category, gene_id, source)

colnames(filtered_gtf) <- c("Gene", "Strand", "Start", "End", "Category", "ID", "Source")



mhc1_genes <- c("10972b","11016","10972a","13304","6557","10984","8993b", "8993","8993c")
mhc1_genes_on_ensembl <- c("ENSOCUG00000010972",
                           "ENSOCUG00000011016",
                           "ENSOCUG00000013304",
                           "ENSOCUG00000006557",
                           "ENSOCUG00000010984",
                           "ENSOCUG00000008993")
mhc1_new_names <- c(paste0("Orcu-U", as.character(1:9)))
mhc1_name_correlation <- data.frame(old = mhc1_genes, new = mhc1_new_names)



mhc1_data <- as.data.frame(read.csv("MHC1_genemodel.csv")) %>%
  dplyr::select(Start,End, Strand, Gene) %>%
  left_join(mhc1_name_correlation, by = c("Gene" = "old")) %>%
  # If you want to remove the old names and keep only the new names
  select(-Gene) %>%
  dplyr::rename(Gene = new) %>%
  group_by(Gene, Strand) %>%
  summarize(Start = min(Start), End = max(End), .groups = 'drop') %>%
  mutate(., Category = "MHC-I") %>%
  mutate(., ID = NA) %>%
  mutate(., Source = "Curated") %>%
  ungroup()


mhc2_data <- as.data.frame(read.csv("MHC2_genemodel.csv")) %>%
  dplyr::select(Start,End, Strand, Gene) %>%
  filter(., Gene != '' & (!grepl("Pseudo",  Gene))) %>%
  group_by(Gene, Strand) %>%
  summarize(Start = min(Start), End = max(End), .groups = 'drop') %>%
  mutate(., Category = "MHC-II") %>%
  mutate(., ID = NA) %>%
  mutate(., Source = "Curated") %>%
  ungroup()

mhc3_gene <- c("TNF", "C2", "C4-A")
peptide_loading_gene <- c("TAP1", "TAP2", "Tapasin", "MOG", "TRIM10", "TRIM15",
                          "TRIM26", "TRIM39", "TRIM40", "GABBR1", "TCF19","RNF39", "ZFP57",
                          "PPP1R11", "IER3")

### Annotate the MHC-III from the GTF, and add the predicted genes
filtered_gtf <- filtered_gtf %>%
  mutate(Category = ifelse(Gene %in% mhc3_gene, "MHC-III", Category))

new_mhc3 <- setdiff(mhc3_gene, filtered_gtf$Gene)

new_mhc3_data <- as.data.frame(read.csv(file = "MHC3_genemodel.csv")) %>%
  filter(Gene %in% new_mhc3) %>%
  mutate(Category = "MHC-III", Source = "NCBI prediction") %>%
  select(Gene, Strand, Start, End, Category, ID, Source)


### Annotate peptide loading genes from the GTF, and add the predicted genes
filtered_gtf <- filtered_gtf %>%
  mutate(Category = ifelse(Gene %in% peptide_loading_gene, "Others", Category))

new_peptide_loading <- setdiff(peptide_loading_gene, filtered_gtf$Gene)

new_peptide_loading_data <- as.data.frame(read.csv(file = "Others_genemodel.csv")) %>%
  filter(Gene %in% new_peptide_loading) %>%
  mutate(Category = "Others", Source = "NCBI prediction") %>%
  select(Gene, Strand, Start, End, Category, ID, Source)



combined_data <- rbind(filtered_gtf, mhc1_data, mhc2_data,
                       new_mhc3_data, new_peptide_loading_data) %>%
  filter(Start >2e+7)

write.csv(combined_data, "All_genemodel.csv")



combined_data$Category <- as.factor(combined_data$Category)

classical_genes <- c("Orcu-U1", "Orcu-U2",
                     "DPB", "DPA", "DQB", "DQA", "DRB1", "DRB2", "DRB3", "DRA")


# Create the plot
chrmap <- ggplot(combined_data, aes(x = Start/1e+6, xend = (Start+5000)/1e+6, y = 0, color = Category)) +
  geom_segment(data = combined_data %>% filter(Category == "Unknown"),
               linewidth = 12, show.legend = FALSE, aes(color = 'lightgray')) +
  geom_segment(data = combined_data %>% filter(Category != "Unknown"),
               linewidth = 15) +
  scale_color_manual(values = c("MHC-I" = "#FF0000", "MHC-II" = "#0070C0", "MHC-III" = "#00802b", "Others" = "black")) +

  geom_text_repel(data = combined_data %>% filter(Category != "Unknown"),
                  aes(label = sub("Orcu-", "", Gene),
                      x = (Start + End) / (2*1e+6),
                      # y = 0,  # Conditionally adjust y based on Category
                      vjust = ifelse(Category == "Others", 8, -8),
                      fontface = ifelse(Gene %in% classical_genes, "bold", "plain")),
                  size = 3.5,
                  angle = 0,
                  hjust = 0.2,
                  direction = "x",
                  max.overlaps = Inf,
                  segment.color = "lightgray",
                  nudge_x = 0.2,
                  segment.size = 0.5) +
  geom_hline(yintercept = 0, linewidth = 1, color = "black") +
  theme_minimal() +
  xlab("Chromosome 12 (Mb)") +
  guides(color = guide_legend(
    title = NULL,  # Remove legend title
    override.aes = list(linewidth = 4),  # Make legend segments smaller
    label.theme = element_text(size = 12)  # Increase legend font size
  )) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(face = 'bold', size = 12, hjust = 0.5, vjust = 0),
        axis.text.x = element_text(face = 'italic', size = 10),
        legend.title = element_text(face = 'bold', size = 12),
        legend.position = "right")


chr_map <- plot_grid(chrmap, ncol = 1, labels = "A")
chr_map
ggsave("chr_map_mhc1mhc2.svg", plot = chr_map, width = 13, height = 3, device = "svg")


# ###################################################################################################
# ########################### Mummerplot ############################################################
# ###################################################################################################
# 
# 
# # Read the delta file
# # Mummer_path <- '/Users/monica/Library/CloudStorage/OneDrive-UniversityofCambridge/Aim1/10_Assembly_matches_plot/'
# mums_file <- "ory2.MHC1.mums"
# 
# mb_formatter <- function(x) {
#   sprintf("%.2f", x / 1e6)
# }
# 
# # Define a function to read and parse the delta data
# readMums <- function(mumsfile) {
#   lines <- readLines(mumsfile)
#   parsed_data <- data.frame()
#   current_header <- NULL
# 
#   # Process lines to extract data
#   for (line in lines) {
#     if (startsWith(line, ">")) {
#       # Handle the header line
#       current_header <- line
#     } else {
#       # Process data lines
#       elements <- as.numeric(unlist(strsplit(trimws(line), "\\s+")))
#       if (length(elements) == 3) {
#         parsed_data <- rbind(parsed_data, c(current_header, elements))
#       }
#     }
#   }
# 
#   # Convert to data frame and clean up
#   colnames(parsed_data) <- c("Header", "Ref_Start", "Query_Start", "Length")
#   parsed_data$Ref_Start <- as.numeric(parsed_data$Ref_Start)
#   parsed_data$Query_Start <- as.numeric(parsed_data$Query_Start)
#   parsed_data$Length <- as.numeric(parsed_data$Length)
# 
#   return(parsed_data)
# }
# 
# # Load your delta file
# mums <- readMums(mums_file)
# 
# # Adjust the positions based on the provided information
# mums_data <- mums %>%
#   mutate(Ref_End = Ref_Start + Length - 1,
#          Query_End = Query_Start + Length - 1,
#          Ref_Start = Ref_Start + 20180000,
#          Ref_End = Ref_End + 20180000,
#          Query_Start = Query_Start + 20180000,
#          Query_End = Query_End + 20180000)
# 
# mums_data$Header <- factor(mums_data$Header, levels = c('> 12:20180000-20319999', '> 12:20180000-20319999 Reverse'))
# levels(mums_data$Header) <- c("Forward", "Reverse")
# 
# 
# mummerplot <- ggplot(mums_data, aes(x = Ref_Start, xend = Ref_End, y = Query_Start, yend = Query_End, color = Header)) +
#   geom_segment(show.legend = FALSE) +  # Draws segments for each alignment
#   geom_point(alpha = 0.8, show.legend = FALSE) +  # Adds points with some transparency
#   coord_fixed(ratio = 1) +
#   theme_bw() +  # Uses a white background theme
#   theme(strip.text.y = element_blank(),  # Rotates and sizes strip text
#         # legend.position.inside = c(0.99, 0.01),  # Places legend at bottom right
#         # legend.justification = c(1, 0),  # Adjusts legend position
#         # strip.background = element_blank(),  # Removes background of strip
#         axis.text= element_text(size = 12, face = "bold"),  # Matches y-axis text size with x-axis
#         axis.ticks.y = element_blank(),
#         axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # Moves y-axis title further away
#         axis.title.x = element_blank()) +  # Removes y-axis ticks
#   ylab('OryCun 2.0 Chr12 (Mb)') +  # Labels y-axis
#   scale_x_continuous(limit = c(20180000, 20320000), breaks = seq(20180000, 20320000, by = 20000),
#                      expand = c(0,0), labels = mb_formatter) +  # Set x-axis breaks and labels
#   scale_y_continuous(breaks = seq(20180000, max(mums_data$Query_End, na.rm = TRUE) + 10000, by = 20000),
#                      expand = c(0,0), labels = mb_formatter) +  # Set y-axis breaks and labels
#   scale_color_manual(values = c("Forward" = "#9400D4", "Reverse" = "#57b5e8"))
# 
# 
# mummer_plot <- plot_grid(mummerplot, ncol = 1, labels = "B")
# ggsave("Mummer_dotplot.png", plot = mummer_plot, width = 10, height = 8, device = 'png')



# 
# 
# 
# ###################################################################################################
# ########################### MHC-I cluster 1 zoom in ###############################################
# ###################################################################################################
# 
# mhc1_data_cluster1 <- as.data.frame(mhc1_data) %>%
#   filter(Start / 1e+6 >= 20.18 & Start / 1e+6 <= 20.32 & End / 1e+6 >= 20.18 & End / 1e+6 <= 20.32)
# 
# # Function to identify continuous segments
# find_continuous_segments <- function(data, start_col, end_col) {
#   data <- data %>%
#     arrange(!!sym(start_col)) %>%
#     mutate(diff_ref = c(NA, diff(!!sym(start_col))),
#            diff_query = c(NA, diff(!!sym(end_col))),
#            group = cumsum((diff_ref != diff_query) | is.na(diff_ref)))
#   
#   continuous_segments <- data %>%
#     group_by(group) %>%
#     summarize(start = min(!!sym(start_col)),
#               end = max(!!sym(end_col)),
#               .groups = 'drop') %>%
#     filter(end - start > 0)  # Filter out segments that aren't continuous
#   
#   return(continuous_segments)
# }
# 
# continuous_segments <- find_continuous_segments(mums_data, 'Ref_Start', 'Query_Start')
# 
# 
# mhc1_map <- ggplot(mhc1_data_cluster1, aes(x = Start/1e+6, xend = End/1e+6, y = 0)) +
#   geom_segment(linewidth = 15, color = '#FF0000') +
#   geom_text(aes(label = sub("Orcu-", "", Gene),
#                       x = (Start + End) / (2 * 1e+6),
#                       y = 0),
#                   fontface = "bold",
#                   size = 5,
#             vjust = 4) +
#   geom_hline(yintercept = 0, linewidth = 1, color = "black") +
#   theme_minimal() +
#   xlab("Chromosome 12 (Mb)") +
#   theme(axis.title = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.x = element_text(face = 'italic', size = 14),
#         # legend.title = element_text(face = 'bold', size = 10),
#         panel.grid = element_blank(),
#         legend.position = "none") +
#   scale_x_continuous(limit = c(20.17, 20.33), breaks = seq(20.18, 20.32, by = 0.04), expand = c(0,0))
# 
# 
# 
# mhc1_map
# ggsave("mhc1_region_map.svg", mhc1_map, width = 12, height = 2, device = "svg")
# 
# 

###################################################################################################
########################### MHC-I cluster 1 gene duplication ###############################################
###################################################################################################
############################################################################################################################
## MHC-I cluster 1 gene duplication with geneviewer
############################################################################################################################


# Simple plot with MHC-I genes and duplication arrows in separate clusters
library(geneviewer)

# Prepare MHC-I cluster data
mhc1_genes_clean <- as.data.frame(mhc1_data) %>%
  filter(Start / 1e+6 >= 20.18 & Start / 1e+6 <= 20.32 & End / 1e+6 >= 20.18 & End / 1e+6 <= 20.32) %>%
  mutate(
    cluster = "MHC-I",
    start = Start / 1e+6,  # Convert to Mb
    end = End / 1e+6,      # Convert to Mb
    strand = ifelse(Strand == "+", 1, -1),
    gene = sub("Orcu-", "", Gene),
    group = "MHC-I"
  ) %>%
  select(cluster, start, end, strand, gene, group)

# Define duplication pairs using your coordinates
duplicated_pairs <- c(
  20.207, 20.213, "inverted1", 
  20.241, 20.246, "inverted1", 
  20.209, 20.213, "direct1",
  20.282, 20.286, "direct1",
  20.218, 20.225, "direct2",
  20.267, 20.275, "direct2",
  20.226, 20.253, "inverted2",
  20.273, 20.299, "inverted2"
)

# Convert to data frame
duplication_matrix <- matrix(duplicated_pairs, ncol = 3, byrow = TRUE)
duplication_df <- data.frame(
  start_pos = as.numeric(duplication_matrix[, 1]),
  end_pos = as.numeric(duplication_matrix[, 2]),
  dup_type = duplication_matrix[, 3],
  stringsAsFactors = FALSE
)

# Create duplication arrows - all in one cluster
duplication_arrows <- duplication_df %>%
  mutate(
    start = start_pos,
    end = end_pos,
    strand = ifelse(grepl("direct", dup_type), 1, 
                    ifelse(duplicated(dup_type), -1, 1)),
    gene = paste0(substr(dup_type, 1, 3), "_", 
                  ifelse(duplicated(dup_type), "2", "1")),
    group = dup_type  # Each pair gets its own group for different colors
  ) %>%
  # Add row numbers to determine which cluster
  mutate(row_num = row_number()) %>%
  mutate(
    # First 4 duplications go to cluster 1, next 4 go to cluster 2
    cluster = ifelse(row_num == 1 | row_num == 2 | row_num == 7 | row_num == 8, 
                     "Duplications 1", 
                     "Duplications 2")) %>%
  select(cluster, start, end, strand, gene, group)



# Combine all data
combined_data <- rbind(
  mhc1_genes_clean,
  duplication_arrows
)

# Create the visualization with two clusters
p_genes_with_duplications <- GC_chart(combined_data, 
                                      cluster = "cluster", 
                                      group = "group",
                                      height = 350,
                                      ) %>%
  # Style MHC-I cluster as red boxes
  GC_genes(
    cluster = "MHC-I",
    marker = "rbox",
    marker_size = "large",
    markerHeight = 20
  ) %>%
  # Style Duplications cluster as arrows with overlap prevention
  GC_genes(
    cluster = "Duplications 1", 
    marker = "arrow",
    marker_size = "small",
    markerHeight = 10,
    itemStyle = list(list(index = 0, stroke = "black",  strokeWidth = 3),
                     list(index = 2, stroke = "black",  strokeWidth = 3)
  )) %>%
  
  GC_genes(
    cluster = "Duplications 2", 
    marker = "arrow",
    marker_size = "small",
    markerHeight = 10,
    itemStyle = list(list(index = 0, stroke = "black",  strokeWidth = 3),
                     list(index = 3, stroke = "black",  strokeWidth = 3)
    )) %>%
  
  # Prevent overlap in the duplications cluster
  # GC_cluster(
  #   prevent_gene_overlap = TRUE,
  #   overlap_spacing = 10
  # ) %>%
  GC_labels("gene", 
            cluster = "MHC-I",  # Only show labels for MHC-I
            fontSize = "14px", 
            color = "white",
            y=14,
            fontWeight = "bold") %>%
  GC_scale(
    cluster = "MHC-I",
    # hidden = TRUE,
    start = 20.20,
    end = 20.32,
    tickValues = c(20.20, 20.24, 20.28, 20.32),
    ticksFormat = ".2f",
    textStyle = list("font-size" = "12px")
  ) %>%
  GC_scale(
    cluster = "Duplications 1",
    start = 20.20,
    end = 20.32,
    tickValues = c(20.20, 20.24, 20.28, 20.32),
    ticksFormat = ".2f",
    textStyle = list("font-size" = "12px"),
    axis_position = "bottom",
    y=10
  ) %>%
  GC_scale(
    cluster = "Duplications 2",
    start = 20.20,
    end = 20.32,
    tickValues = c(20.20, 20.24, 20.28, 20.32),
    ticksFormat = ".2f",
    textStyle = list("font-size" = "12px"),
    axis_position = "bottom",
    y=10
  ) %>%
  
  # Set colors for all groups
  GC_color(
    customColors = list(
      "MHC-I" = "#FF0000",      # Red for MHC-I
      "direct1" = "#9400D4",    # Green for direct1
      "direct2" = "#9400D4",    # Orange for direct2  
      "inverted1" = "#57b5e8",  # Blue for inverted1
      "inverted2" = "#57b5e8"   # Purple for inverted2
    )
  ) %>%
  GC_annotation(
    cluster = "MHC-I",
    type = "promoter", # terminator
    x = 20.21176,
    y = 55, 
    direction = "revserse", # reverse
    style = list(
      fill = "none",
      stroke = "black",
      strokeWidth = 1.5
      # Any other CSS style
    ),
    rotation = 0,
    scale = 1
  ) %>% 
  GC_annotation(
    cluster = "MHC-I",
    type = "promoter", # terminator
    x = 20.24234,
    y = 55, 
    direction = "forward", # reverse
    style = list(
      fill = "none",
      stroke = "black",
      strokeWidth = 1.5
      # Any other CSS style
    ),
    rotation = 0,
    scale = 1
  ) %>% 
  GC_annotation(
    cluster = "MHC-I",
    type = "promoter", # terminator
    x = 20.28429,
    y = 55, 
    direction = "revserse", # reverse
    style = list(
      fill = "none",
      stroke = "black",
      strokeWidth = 1.5
      # Any other CSS style
    ),
    rotation = 0,
    scale = 1
  ) %>% 
  GC_annotation(
    cluster = "MHC-I",
    type = "promoter", # terminator
    x = 20.30245,
    y = 55, 
    direction = "revserse", # reverse
    style = list(
      fill = "none",
      stroke = "black",
      strokeWidth = 1.5
      # Any other CSS style
    ),
    rotation = 0,
    scale = 1
  ) %>% 
  
  
  
  GC_legend(FALSE)

# Print summary
cat("Data structure:\n")
print(table(combined_data$cluster, combined_data$group))

p_genes_with_duplications

