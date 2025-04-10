library(ggplot2)
library(pheatmap)
library(edgeR)
library(dplyr)

library(data.table)
meta = fread("~/projects/combined_all/final_meta.csv")

head(meta)

setwd("../histone_volcano/")

setwd("female_AMY/")

library(ggplot2)
library(dplyr)
library(patchwork)

sex <- "male"

# Function to parse row names and extract domain size
parse_domain_size <- function(rowname) {
  parts <- unlist(strsplit(rowname, "[:-]"))
  if (length(parts) >= 3) {
    start_pos <- as.numeric(parts[2])
    end_pos <- as.numeric(parts[3])
    if (!is.na(start_pos) && !is.na(end_pos)) {
      return(abs(end_pos - start_pos))
    }
  }
  return(NA)
}

# Function to add clade information based on cell type
assign_clade <- function(cell_type) {
  if (grepl("Glut", cell_type)) {
    "Glut"
  } else if (grepl("Gaba", cell_type)) {
    "Gaba"
  } else {
    "NN"
  }
}

# Comparisons and regions
comparisons <- c("2mo_vs_18mo", "9mo_vs_18mo", "2mo_vs_9mo")
regions <- c("FC", "ENT", "HCA", "HCP", "AMY", "NAC", "RLP", "CP")

# Store all data for faceting
all_data <- list()

# Loop through comparisons
for (comparison in comparisons) {
  combined_data <- data.frame() # Initialize combined data for this comparison
  
  for (reg in regions) {
    # Dynamically construct file paths for the current region
    region_path <- paste0("../", sex, "_", reg)
    setwd(region_path)
    
    # Load files for the current comparison
    comparison_files <- list.files(pattern = paste0("_", comparison, "_differential_peaks.txt$"))
    
    for (file in comparison_files) {
      cell_type <- sub(paste0("_", comparison, "_differential_peaks.txt"), "", file)
      data <- read.table(file, header = TRUE)
      data$CellType <- cell_type
      data$Region <- reg
      data$Comparison <- comparison
      data$Clade <- assign_clade(cell_type)
      data$DomainSize <- sapply(rownames(data), parse_domain_size)
      combined_data <- rbind(combined_data, data)
    }
  }
  
  # Add significance and filter invalid rows
  combined_data <- combined_data %>%
    filter(!is.na(logFC) & !is.na(PValue) & !is.na(DomainSize)) %>%
    mutate(
      Significance = ifelse(logFC > 1 & FDR < 0.05, "Upregulated",
                            ifelse(logFC < -1 & FDR < 0.05, "Downregulated", "Neutral")),
      Clade = ifelse(FDR > 0.05, "NS", Clade)
    )
  
  # Store the combined data
  all_data[[comparison]] <- combined_data
}

# Combine all comparisons into a single data frame
all_data_combined <- do.call(rbind, all_data)

# Create faceted volcano plots
volcano_plot <- ggplot(all_data_combined, aes(x = logFC, y = -log10(FDR), color = Clade, size = DomainSize)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Glut" = "firebrick1", "Gaba" = "green2", "NN" = "blue3", "NS" = "gray")) +
  theme_minimal() +
  labs(
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  ylim(0, 10.5) + # Fixed y-axis limits
  facet_wrap(~ Region + Comparison, ncol = 2, scales = "free") + # Facet by region and comparison
  theme(
    legend.position = "right", # Single legend
    strip.text = element_text(size = 10, face = "bold")
  )

# Save the plot
#ggsave(filename = "faceted_volcano_plots.pdf", plot = volcano_plot, height = 15, width = 12)

# View the plot
print(volcano_plot)


all_data_combined$Region = factor(all_data_combined$Region,
                                  levels = c("FC", "ENT", "HCA", "HCP", "AMY", "NAC", "RLP", "CP"))


all_data_combined = all_data_combined[-which(all_data_combined$CellType=="Astroependymal_NN"),]

all = all_data_combined

tab = table(meta$region, meta$celltype_final)
colnames(tab) = gsub(" " , "_" , colnames(tab))
colnames(tab) = gsub("-" , "." , colnames(tab))
colnames(tab) = gsub("/" , "_" , colnames(tab))

all$CellType[which(all$CellType=="L2_3_IT_PIR")] = "L2_3_IT_PIR.ENTl_Glut"
all$CellType[which(all$CellType=="L5_6_IT_TPE")] = "L5_6_IT_TPE.ENT_Glut"

unique(all$CellType)[which(!unique(all$CellType)%in%colnames(tab))]


tab <- as.data.frame(tab)
colnames(tab) = c("region", "celltype", "count")
head(tab)



# Create unique key in tab_long
tab$RegionCellType <- paste(tab$region, tab$celltype, sep = "_")

# Create unique key in combined_data
all$RegionCellType <- paste(all$Region, all$CellType, sep = "_")

# Filter combined_data based on the count in tab_long
filtered_combined_data <- all %>%
  filter(RegionCellType %in% tab$RegionCellType[tab$count >= 500])

# Drop the temporary RegionCellType column (optional)
filtered_combined_data <- filtered_combined_data %>%
  select(-RegionCellType)

# View the filtered combined_data
head(filtered_combined_data)


head(filtered_combined_data[which(filtered_combined_data$Comparison=="2mo_vs_18mo"),])

setwd("../histone_volcano/")

volcano_plot <- ggplot(filtered_combined_data[which(filtered_combined_data$Comparison=="2mo_vs_18mo"),], aes(x = logFC, y = -log10(FDR), color = Clade, size = DomainSize)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Glut" = "firebrick1", "Gaba" = "green2", "NN" = "blue3", "NS" = "gray")) +
  theme_minimal() +
  labs(
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  ylim(0, 25) + # Fixed y-axis limits
  xlim(-9, 9) +
  facet_wrap(~ Region, ncol = 1, scales = "free") + # Facet by region and comparison
  theme(
    legend.position = "right", # Single legend
    strip.text = element_text(size = 10, face = "bold")
  )
options(repr.plot.width=3, repr.plot.height=10)

volcano_plot
#ggsave(volcano_plot, height = 6.5, width = 5, file = "Male_2vs18mo_plots.pdf")

volcano_plot <- ggplot(filtered_combined_data[which(filtered_combined_data$Comparison=="2mo_vs_9mo"),], aes(x = logFC, y = -log10(FDR), color = Clade, size = DomainSize)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Glut" = "firebrick1", "Gaba" = "green2", "NN" = "blue3", "NS" = "gray")) +
  theme_minimal() +
  labs(
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  ylim(0, 25) + # Fixed y-axis limits
  xlim(-9, 9) +
  facet_wrap(~ Region, ncol = 2, scales = "free") + # Facet by region and comparison
  theme(
    legend.position = "right", # Single legend
    strip.text = element_text(size = 10, face = "bold")
  )

volcano_plot

volcano_plot <- ggplot(filtered_combined_data[which(filtered_combined_data$Comparison=="9mo_vs_18mo"),], aes(x = logFC, y = -log10(FDR), color = Clade, size = DomainSize)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Glut" = "firebrick1", "Gaba" = "green2", "NN" = "blue3", "NS" = "gray")) +
  theme_minimal() +
  labs(
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  ylim(0, 25) + # Fixed y-axis limits
  xlim(-9, 9) +
  facet_wrap(~ Region, ncol = 2, scales = "free") + # Facet by region and comparison
  theme(
    legend.position = "right", # Single legend
    strip.text = element_text(size = 10, face = "bold")
  )

volcano_plot


volcano_plot <- ggplot(all_data_combined[which(all_data_combined$Comparison=="2mo_vs_18mo"),], aes(x = logFC, y = -log10(FDR), color = Clade, size = DomainSize)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Glut" = "firebrick1", "Gaba" = "green2", "NN" = "blue3", "NS" = "gray")) +
  theme_minimal() +
  labs(
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  ylim(0, 25) + # Fixed y-axis limits
  xlim(-9, 9) +
  facet_wrap(~ Region, ncol = 2, scales = "free") + # Facet by region and comparison
  theme(
    legend.position = "right", # Single legend
    strip.text = element_text(size = 10, face = "bold")
  )

ggsave(volcano_plot, height = 6.5, width = 5, file = "Male_2vs18mo_plots.pdf")

volcano_plot <- ggplot(all_data_combined[which(all_data_combined$Comparison=="9mo_vs_18mo"),], aes(x = logFC, y = -log10(FDR), color = Clade, size = DomainSize)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Glut" = "firebrick1", "Gaba" = "green2", "NN" = "blue3", "NS" = "gray")) +
  theme_minimal() +
  labs(
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  ylim(0, 25) + # Fixed y-axis limits
  xlim(-9, 9) +
  facet_wrap(~ Region, ncol = 2, scales = "free") + # Facet by region and comparison
  theme(
    legend.position = "right", # Single legend
    strip.text = element_text(size = 10, face = "bold")
  )
ggsave(volcano_plot, height = 6.5, width = 5, file = "Male_9vs18mo_plots.pdf")


volcano_plot <- ggplot(all_data_combined[which(all_data_combined$Comparison=="2mo_vs_9mo"),], aes(x = logFC, y = -log10(FDR), color = Clade, size = DomainSize)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Glut" = "firebrick1", "Gaba" = "green2", "NN" = "blue3", "NS" = "gray")) +
  theme_minimal() +
  labs(
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  ylim(0, 25) + # Fixed y-axis limits
  xlim(-9, 9) +
  facet_wrap(~ Region, ncol = 2, scales = "free") + # Facet by region and comparison
  theme(
    legend.position = "right", # Single legend
    strip.text = element_text(size = 10, face = "bold")
  )

ggsave(volcano_plot, height = 6.5, width = 5, file = "Male_2vs9mo_plots.pdf")

cur = all_data_combined[which(all_data_combined$Comparison=="2mo_vs_18mo"),]
cur = cur[order(cur$FDR),]
head(cur,10)

unique(cur[which(cur$Region=="CP"),"CellType"]) 


