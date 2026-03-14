rm(list = ls()) 
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(ggplot2)
library(dplyr)
library(readr)
library(geneviewer)
library(rtracklayer)

################################################################################
### RABBIT DATA PREP
################################################################################

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
           start >= 10000000 &
           end <= 24000000 &
           !gene_id %in% mhc1_genes_on_ensembl &
           # is.na(gene_name) == F &
           !grepl("HLA", gene_name)) %>%
  mutate(., Category = 'Unannotated') %>%
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



rabbit_data <- rbind(filtered_gtf, mhc1_data, mhc2_data,
                       new_mhc3_data, new_peptide_loading_data) %>%
  filter(Start >1e+7)






# Read data files
framework_mapping <- read.csv("Framework_genes.csv", stringsAsFactors = FALSE)
# rabbit_data <- read.csv("All_genemodel.csv", stringsAsFactors = FALSE)

# Read human GTF
gtf <- import("Homo_sapiens.GRCh38.114.chr.gtf")
chr6_genes <- gtf[seqnames(gtf) == "6" & gtf$type == "gene" & 
                    start(gtf) >= 26000000 & end(gtf) <= 34000000]
human_data <- data.frame(
  gene = chr6_genes$gene_name,
  start = start(chr6_genes),
  end = end(chr6_genes),
  stringsAsFactors = FALSE
)
human_data <- human_data[!is.na(human_data$gene), ]

# Function to find gene coordinates
find_gene_coord <- function(gene_name, gene_data, species, coord_type = "start") {
  if(species == "rabbit") {
    match_by_name <- gene_data[!is.na(gene_data$Gene) & gene_data$Gene == gene_name, ]
    if(nrow(match_by_name) == 0) {
      match_by_id <- gene_data[!is.na(gene_data$ID) & gene_data$ID == gene_name, ]
      if(nrow(match_by_id) > 0) {
        if(coord_type == "start") return(match_by_id$Start[1])
        if(coord_type == "end") return(match_by_id$End[1])
      }
    } else {
      if(coord_type == "start") return(match_by_name$Start[1])
      if(coord_type == "end") return(match_by_name$End[1])
    }
  } else {
    match <- gene_data[!is.na(gene_data$gene) & gene_data$gene == gene_name, ]
    if(nrow(match) > 0) {
      if(coord_type == "start") return(match$start[1])
      if(coord_type == "end") return(match$end[1])
    }
  }
  return(NA)
}

# Create framework data - FIXED scaling
create_framework_data <- function() {
  human_frameworks <- data.frame()
  rabbit_frameworks <- data.frame()
  
  for(i in 1:nrow(framework_mapping)) {
    fw <- framework_mapping[i, ]
    
    # Process Human framework
    if(fw$Start_Human != "-" && fw$End_Human != "-") {
      human_start <- find_gene_coord(fw$Start_Human, human_data, "human", "start")
      human_end <- find_gene_coord(fw$End_Human, human_data, "human", "end")
      
      if(!is.na(human_start) && !is.na(human_end)) {
        # Convert to Mb - no additional scaling needed
        human_start_mb <- human_start / 1000000
        human_end_mb <- human_end / 1000000
        
        human_frameworks <- rbind(human_frameworks, data.frame(
          cluster = "Human Chr6",
          Gene = fw$Framework,
          start = min(human_start_mb, human_end_mb),
          end = max(human_start_mb, human_end_mb),
          strand = 1,
          region = fw$Region,
          framework_name = fw$Framework,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Process Rabbit framework
    if(fw$Start_Rabbit != "-" && fw$End_Rabbit != "-") {
      rabbit_start <- find_gene_coord(fw$Start_Rabbit, rabbit_data, "rabbit", "start")
      rabbit_end <- find_gene_coord(fw$End_Rabbit, rabbit_data, "rabbit", "end")
      
      if(!is.na(rabbit_start) && !is.na(rabbit_end)) {
        # Convert to Mb - no additional scaling needed
        rabbit_start_mb <- rabbit_start / 1000000
        rabbit_end_mb <- rabbit_end / 1000000
        
        rabbit_frameworks <- rbind(rabbit_frameworks, data.frame(
          cluster = "Rabbit Chr12",
          Gene = fw$Framework,
          start = min(rabbit_start_mb, rabbit_end_mb),
          end = max(rabbit_start_mb, rabbit_end_mb),
          strand = 1,
          region = fw$Region,
          framework_name = fw$Framework,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Sort by genomic position
  human_frameworks <- human_frameworks[order(human_frameworks$start), ]
  rabbit_frameworks <- rabbit_frameworks[order(rabbit_frameworks$start), ]
  
  return(rbind(human_frameworks, rabbit_frameworks))
}

# Create the framework data
framework_data <- create_framework_data()


# Process FW2 modifications - keep human blocks separate with distinct names for linking
# For Human: rename FW2-1 and FW2-2 but keep separate for proper linking
human_fw2_1_rows <- framework_data$cluster == "Human Chr6" & framework_data$framework_name == "FW2-1"
human_fw2_2_rows <- framework_data$cluster == "Human Chr6" & framework_data$framework_name == "FW2-2"

# Change display labels to "FW2" but keep distinct framework names for linking
framework_data[human_fw2_1_rows, "Gene"] <- "FW2"
framework_data[human_fw2_2_rows, "Gene"] <- "FW2"
# Keep framework_name as FW2-1 and FW2-2 for linking

# For Rabbit: change labels but use distinct framework names that match human
rabbit_fw2_1_rows <- framework_data$cluster == "Rabbit Chr12" & framework_data$framework_name == "FW2-1"
rabbit_fw2_2_rows <- framework_data$cluster == "Rabbit Chr12" & framework_data$framework_name == "FW2-2"

framework_data[rabbit_fw2_1_rows, "Gene"] <- "FW2"
framework_data[rabbit_fw2_2_rows, "Gene"] <- "FW2"
# Swap the framework names to show cross-correspondence
framework_data[rabbit_fw2_1_rows, "framework_name"] <- "FW2-1"  # rabbit first → human second
framework_data[rabbit_fw2_2_rows, "framework_name"] <- "FW2-2"  # rabbit second → human first

framework_data <- framework_data %>%
  mutate(Gene = if_else(Gene %in% c("Extend MHC-I-1", "Extend MHC-I-2"), "Extend MHC-I", Gene))






framework_data$region <- factor(framework_data$region,
                                levels = c("MHC-I", "MHC-II", "MHC-III"))

# Create manual links data
create_framework_links <- function() {
  links_data <- data.frame(
    framework = c("Extend MHC-I-1", "Extend MHC-I-2", "Alpha", "FW1", "Kappa", "FW2-1", "FW2-2",
                  "Beta", "MHC-III", "MHC-II", "Extend MHC-II"),
    relationship = c("conserved", "conserved", "lost", "conserved", "conserved", "translocated", "conserved",
                     "translocated", "translocated", "conserved", "conserved"),
    stringsAsFactors = FALSE
  )
  return(links_data)
}

# Add link relationships
links_info <- create_framework_links()
framework_data <- framework_data %>%
  left_join(links_info, by = c("framework_name" = "framework")) %>%
  mutate(
    relationship = ifelse(is.na(relationship), "conserved", relationship),
    link_group = paste(framework_name, relationship, sep = "_")
  )

# Add identity column
framework_data$identity <- case_when(
  framework_data$relationship == "conserved" ~ 35,
  framework_data$relationship == "translocated" ~ 75,
  framework_data$relationship == "lost" ~ 100,
  TRUE ~ 0
)

# Custom colors
custom_colors <- c("#FF0000", "#00802B", "#0070C0")  # red, green, blue


# Create the synteny plot with FIXED positioning
mhc_synteny_plot <- GC_chart(
  framework_data,
  cluster = "cluster",
  strand = "strand", 
  group = "region",
  height = 430
) %>%
  # FIXED: Move labels further from segments (increased y distance)
  GC_labels("Gene", cluster = "Human Chr6", y = 55, fontSize = "15px",
            itemStyle = 
              list(
                list(index = 5, fill = "", fontWeight = "bold"),
                list(index = 6, fill = "", fontWeight = "bold"),
                list(index = 7, fill = "", fontWeight = "bold")
  )) %>%
  GC_labels("Gene", cluster = "Rabbit Chr12", y = 55, fontSize = "15px",
            itemStyle = 
              list(
                list(index = 1, fill = "", fontWeight = "bold"),
                list(index = 2, fill = "", fontWeight = "bold"),
                list(index = 3, fill = "", fontWeight = "bold")
  )) %>%
  GC_genes(
    marker_size = "medium",
    marker = "rbox",
    markerHeight = 13,
    customColors = custom_colors, # A vector of color name
  ) %>%
  GC_links(
    group = "framework_name",
    measure = "identity",
    label = FALSE,
    color_bar = FALSE,
    curve = 0.4
  ) %>%
  GC_clusterLabel(title = c("Human chr6 (Mb)", "Rabbit chr12 (Mb)"),
    fontSize = "16px") %>%
  # GC_color(customColors = custom_colors) %>%

  GC_scale(
    cluster = "Human Chr6",
    start = 26,
    end = 34,
    textStyle = list("font-size" = "14px"),
    ticksCount = 3,
    y=45
  ) %>%
  # GC_scale(
  #   cluster = "Rabbit Chr12",
  #   start = 10,
  #   end = 24,
  #   textStyle = list("font-size" = "14px"),
  #   ticksCount = 5,
  #   y=45
  # ) %>%
  # The updated GC_scale() call for the rabbit track
  GC_scale(
    cluster = "Rabbit Chr12",
    start = 10,
    end = 24,
    textStyle = list("font-size" = "14px"),
    ticksCount =10,
    y = 45,
    breaks = list(
      list(start = 12, end = 17) # Adds a break from 11 Mb to 19.5 Mb
    )
  ) %>%
  GC_legend(FALSE)

# No legend (already removed)
print(mhc_synteny_plot)
# Load the necessary library
library(htmlwidgets)

# Save the plot widget to an HTML file
saveWidget(mhc_synteny_plot, 
           file = "mhc_synteny_plot.html", 
           selfcontained = TRUE)
