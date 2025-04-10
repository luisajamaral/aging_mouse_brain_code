color_region = read.csv("../../Figures/color_scheme/MajorRegion-id.csv")
head(color_region)

color_celltype = read.csv("../../Figures/color_scheme/updated_celltype_palette.csv")
head(color_celltype)



# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Specify the path to the folder containing your CSV files
path <- "."  # Update with the actual path if needed
files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
files <- files[grep("--", files)]

# Thresholds for upregulation/downregulation
logFC_threshold <- 0.5  # Adjust as necessary
p_adj_threshold <- 0.05

# Function to extract clade and region from the filename
extract_metadata <- function(file_name) {
  parts <- strsplit(basename(file_name), "--")[[1]]
  clade <- parts[1]
  region <- sub("\\.csv$", "", parts[2])
  list(clade = clade, region = region)
}

# Initialize a data frame to store gene statistics
all_data <- data.frame()

# Process each file
for (file in files) {
  metadata <- extract_metadata(file)
  clade <- metadata$clade
  region <- metadata$region
  
  data <- read.csv(file)
  colnames(data)[1] <- "Gene"
  
  # Classify genes as upregulated or downregulated
  data <- data %>%
    mutate(Status = case_when(
      avg_log2FC > logFC_threshold & p_val_adj < p_adj_threshold ~ "Upregulated",
      avg_log2FC < -logFC_threshold & p_val_adj < p_adj_threshold ~ "Downregulated",
      TRUE ~ "Neutral"
    )) %>%
    filter(Status != "Neutral") %>%
    mutate(Celltype = clade, Region = region)
  
  all_data <- bind_rows(all_data, data)
}

# Summarize the data by Gene, Celltype, Status, and Region
summary_data <- all_data %>%
  group_by(Gene, Celltype, Status, Region) %>%
  summarize(Count = n(), .groups = "drop")

# Define the clades
clades <- c("NN", "Gaba", "Glut")

# Assign Clade based on Celltype
summary_data <- summary_data %>%
  mutate(Clade = "")  # Initialize Clade column

for (clade in clades) {
  summary_data$Clade[grep(clade, summary_data$Celltype)] <- clade
}

# Filter out rows with unmatched Clades
summary_data <- summary_data %>%
  filter(Clade != "")

# Select top 15 genes for each clade and status without aggregation
top_genes_data <- summary_data %>%
  group_by(Gene, Clade, Status) %>%
  mutate(Unique_DEG_Count = n()) %>%  # Count unique DEGs
  arrange(desc(Unique_DEG_Count)) %>%
  slice_head(n = 15) %>%
  ungroup()

# Load the color keys
color_region <- read.csv("../../Figures/color_scheme/MajorRegion-id.csv")
color_celltype <- read.csv("../../Figures/color_scheme/updated_celltype_palette.csv")
color_celltype$CellType <- gsub(" ", "_", color_celltype$CellType)
color_celltype$CellType <- gsub("/", "-", color_celltype$CellType)

# Add colors to the plotting data
top_genes_data <- top_genes_data %>%
  left_join(color_region, by = c("Region" = "Region")) %>%
  rename(Region_Color = Color) %>%
  left_join(color_celltype, by = c("Celltype" = "CellType")) %>%
  rename(Celltype_Color = Color)

# Plot for each clade
for (clade in clades) {
  clade_data <- top_genes_data %>%
    filter(Clade == clade)
  
  # Skip plotting if there's no data for the clade
  if (nrow(clade_data) == 0) {
    message(paste("No data available for", clade, "- skipping plot"))
    next
  }
  
  # Plot
  plot <- ggplot(clade_data, aes(
    x = reorder(Gene, -Unique_DEG_Count), 
    y = Unique_DEG_Count, 
    fill = Celltype, 
    color = Region
  )) +
    geom_bar(stat = "identity", position = "stack", size = 0.8) +  # `size` for border thickness
    scale_fill_manual(values = setNames(color_celltype$Color, color_celltype$CellType)) +  # Fill colors
    scale_color_manual(values = setNames(color_region$Color, color_region$Region)) +  # Outline colors
    coord_flip() +
    facet_wrap(~ Status, scales = "free", nrow = 2) +
    labs(
      title = paste("Top 15 Upregulated and Downregulated Genes for", clade),
      x = "Gene",
      y = "Unique DEG Count",
      fill = "Cell Type",
      color = "Brain Region"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  print(plot)
}




options(repr.plot.width = 5, repr.plot.height = 7)


# Load required libraries
library(ggplot2)
library(dplyr)

# Thresholds for upregulation/downregulation
logFC_threshold <- 0.25
p_adj_threshold <- 0.05

# Function to extract metadata (cell type and region) from filenames
extract_metadata <- function(file_name) {
  parts <- strsplit(basename(file_name), "--")[[1]]
  clade <- parts[1]
  region <- sub("\\.csv$", "", parts[2])
  list(clade = clade, region = region)
}

# Initialize a data frame to store all DEGs
all_data <- data.frame()

# Process files and classify DEGs
files <- list.files(".", pattern = "\\.csv$", full.names = TRUE)
files <- files[grep("--", files)]  # Ensure valid filenames

for (file in files) {
  metadata <- extract_metadata(file)
  clade <- metadata$clade
  region <- metadata$region
  
  data <- read.csv(file)
  colnames(data)[1] <- "Gene"
  
  data$avg_log2FC = -data$avg_log2FC
  # Classify genes as upregulated or downregulated
  data <- data %>%
    mutate(Status = case_when(
      avg_log2FC > logFC_threshold & p_val_adj < p_adj_threshold ~ "Upregulated",
      avg_log2FC < -logFC_threshold & p_val_adj < p_adj_threshold ~ "Downregulated",
      TRUE ~ "Neutral"
    )) %>%
    filter(Status != "Neutral") %>%
    mutate(Celltype = clade, Region = region)
  
  all_data <- bind_rows(all_data, data)
}

# Define the clades
clades <- c("NN", "Gaba", "Glut")

# Load color keys
color_region <- read.csv("../../Figures/color_scheme/MajorRegion-id.csv")
color_celltype <- read.csv("../../Figures/color_scheme/updated_celltype_palette.csv")
color_celltype$CellType <- gsub(" ", "_", color_celltype$CellType)
color_celltype$CellType <- gsub("/", "-", color_celltype$CellType)

# Generate plots separately for each clade
for (clade in clades) {
  # Filter data for the specific clade
  clade_data <- all_data %>%
    filter(grepl(clade, Celltype))
    
  # Summarize the data to find the top 15 genes separately for each Status
  top_upregulated_genes <- clade_data %>%
    filter(Status == "Upregulated") %>%
    count(Gene) %>%
    arrange(desc(n)) %>%
    slice_head(n = 20) %>%
    pull(Gene)
  
  top_downregulated_genes <- clade_data %>%
    filter(Status == "Downregulated") %>%
    count(Gene) %>%
    arrange(desc(n)) %>%
    slice_head(n = 20) %>%
    pull(Gene)
    
     
  # Filter the original data for the top genes
  filtered_data <- clade_data %>%
    filter((Status == "Upregulated" & Gene %in% top_upregulated_genes) |
             (Status == "Downregulated" & Gene %in% top_downregulated_genes))
  
    
  # Add colors to the data
  filtered_data <- filtered_data %>%
    left_join(color_region, by = c("Region" = "Region")) %>%
    rename(Region_Color = Color) %>%
    left_join(color_celltype, by = c("Celltype" = "CellType")) %>%
    rename(Celltype_Color = Color)
  
  # Plot
  plot <- ggplot(filtered_data, aes(
    x = reorder(Gene, table(Gene)[Gene]), 
    #color = Celltype, 
    fill = Region
  )) +
    geom_bar(stat = "count",size = 0.8) +  # `size` for border thickness
    coord_flip() +
    facet_wrap( ~ Status, scales = "free", nrow = 2) +
    #scale_fill_manual(values = setNames(color_celltype$Color, color_celltype$CellType)) +  # Fill colors
    scale_fill_manual(values = setNames(color_region$Color, color_region$Region)) +  # Outline colors
    labs(
      title = paste("Top 15 Upregulated and Downregulated Genes for", clade),
      x = "Gene",
      y = "Frequency",
      fill = "Cell Type",
      color = "Brain Region"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  print(plot)
}


head(top_genes)

# Load required libraries
library(ggplot2)
library(dplyr)

# Thresholds for upregulation/downregulation
logFC_threshold <- 0.25
p_adj_threshold <- 0.05

# Function to extract metadata (cell type and region) from filenames
extract_metadata <- function(file_name) {
  parts <- strsplit(basename(file_name), "--")[[1]]
  clade <- parts[1]
  region <- sub("\\.csv$", "", parts[2])
  list(clade = clade, region = region)
}

# Initialize a data frame to store all DEGs
all_data <- data.frame()

# Process files and classify DEGs
files <- list.files(".", pattern = "\\.csv$", full.names = TRUE)
files <- files[grep("--", files)]  # Ensure valid filenames

for (file in files) {
  metadata <- extract_metadata(file)
  clade <- metadata$clade
  region <- metadata$region
  
  data <- read.csv(file)
  colnames(data)[1] <- "Gene"
  
  data$avg_log2FC = -data$avg_log2FC
  # Classify genes as upregulated or downregulated
  data <- data %>%
    mutate(Status = case_when(
      avg_log2FC > logFC_threshold & p_val_adj < p_adj_threshold ~ "Upregulated",
      avg_log2FC < -logFC_threshold & p_val_adj < p_adj_threshold ~ "Downregulated",
      TRUE ~ "Neutral"
    )) %>%
    filter(Status != "Neutral") %>%
    mutate(Celltype = clade, Region = region)
  
  all_data <- bind_rows(all_data, data)
}

# Define the clades
clades <- c("NN", "Gaba", "Glut")

# Load color keys
color_region <- read.csv("../../Figures/color_scheme/MajorRegion-id.csv")
color_celltype <- read.csv("../../Figures/color_scheme/updated_celltype_palette.csv")
color_celltype$CellType <- gsub(" ", "_", color_celltype$CellType)
color_celltype$CellType <- gsub("/", "-", color_celltype$CellType)

# Generate plots separately for each clade
for (clade in clades) {
  # Filter data for the specific clade
  clade_data <- all_data %>%
    filter(grepl(clade, Celltype))
    
  # Summarize the data to find the top 15 genes separately for each Status
  top_upregulated_genes <- clade_data %>%
    filter(Status == "Upregulated") %>%
    count(Gene) %>%
    arrange(desc(n)) %>%
    slice_head(n = 20) %>%
    pull(Gene)
  
  top_downregulated_genes <- clade_data %>%
    filter(Status == "Downregulated") %>%
    count(Gene) %>%
    arrange(desc(n)) %>%
    slice_head(n = 20) %>%
    pull(Gene)
    
     
  # Filter the original data for the top genes
  filtered_data <- clade_data %>%
    filter((Status == "Upregulated" & Gene %in% top_upregulated_genes) |
             (Status == "Downregulated" & Gene %in% top_downregulated_genes))
  
    
  # Add colors to the data
  filtered_data <- filtered_data %>%
    left_join(color_region, by = c("Region" = "Region")) %>%
    rename(Region_Color = Color) %>%
    left_join(color_celltype, by = c("Celltype" = "CellType")) %>%
    rename(Celltype_Color = Color)
  
  # Plot
  plot <- ggplot(filtered_data, aes(
    x = reorder(Gene, table(Gene)[Gene]), 
    #color = Celltype, 
    fill = Region
  )) +
    geom_bar(stat = "count",size = 0.8) +  # `size` for border thickness
    coord_flip() +
    facet_wrap( ~ Status, scales = "free", nrow = 2) +
    #scale_fill_manual(values = setNames(color_celltype$Color, color_celltype$CellType)) +  # Fill colors
    scale_fill_manual(values = setNames(color_region$Color, color_region$Region)) +  # Outline colors
    labs(
      title = paste("Top 15 Age-Differential Genes for", clade),
      x = "Gene",
      y = "Frequency",
      #fill = "Cell Type",
      fill = "Brain Region"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  svg(paste(clade,"DEG_freq_plot.svg", sep = ""), width =5, height = 7)
  print(plot)
  dev.off()
}



