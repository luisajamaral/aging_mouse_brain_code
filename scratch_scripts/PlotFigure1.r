library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)

setwd("../../Figures/Figure1-overview/")

age_color = read.csv("../color_scheme/Age.csv")
region_color = read.csv("../color_scheme/MajorRegion.csv")
celltype_color = read.csv("../color_scheme/updated_celltype_palette.csv")
modality_color = read.csv("../color_scheme/Modality.csv") 
sex_color = read.csv("../color_scheme/Sex.csv") 

rownames(age_color) = age_color$Age
rownames(region_color) = region_color$Region
rownames(celltype_color) = celltype_color$CellType
rownames(modality_color) = modality_color$Modality
rownames(sex_color) = sex_color$Sex

setwd("../../h5ads_final/")
meta = read.csv("../final_meta.csv")
meta$celltype_final = gsub("Lymphoid","Lymphoid NN", meta$celltype_final)
meta = meta[-which(meta$celltype_final=="doublet"),]
meta$celltype_final=paste(meta$celltype_final)
meta$celltype_final[which(meta$best_celltype_fixed=="IOL")] = "IOL NN"
unique(meta$celltype_final)[which(! unique(meta$celltype_final) %in% celltype_color$CellType)]

ord <- names(table(meta$celltype_final)[order(table(meta$celltype_final))])
meta$celltype_final <- factor(meta$celltype_final, levels = ord)

cat(unique(meta$celltype_final)[which(! unique(meta$celltype_final) %in% celltype_color$CellType)], sep = "\n")

color_vector <- setNames(celltype_color$Color, celltype_color$CellType)

meta$major = sapply(strsplit(as.character(meta$celltype_final), " "), function(x) tail(x, 1))
meta$major[which(meta$major%in%"Neur")] = "Gaba"
meta$major[which(meta$major%in%"Gaba-Chol")] = "Gaba"
meta$major[which(meta$major%in%"Glut-Sero")] = "Glut"
meta$major[which(meta$major%in%"Dopa")] = "Glut"
meta$major[which(meta$major%in%"IMN")] = "NN"



setwd("../../combined_all/Figures/Figure1-overview/")

options(repr.plot.width=9, repr.plot.height=6.5)
library(scales)
g = ggplot(meta_sub, aes(x = celltype_final, fill = celltype_final)) +
  geom_bar(stat = "count", color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() + 
  geom_text(stat='count', aes(label=label_number(scale_cut = cut_si(""), accuracy = 1)(after_stat(count))), 
            position=position_stack(vjust=1), hjust=-0.1)+ 
  scale_y_continuous(breaks = c(0, 50000, 100000, 150000, 200000), 
                     labels = label_number(scale_cut = cut_si(""), accuracy = 1),
                     limits = c(0, 225000),  # Adjusts the upper limit to 110% of the max value
                     expand = expansion(mult = c(0, 0.05)))   # Adds some space at the top
pdf("celltype_barplot_over5k.pdf",width=9, height=6.5)
print(g)
dev.off()


g = ggplot(meta, aes(x = celltype_final, fill = celltype_final)) +
  geom_bar(stat = "count", color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() + 
  geom_text(stat='count', aes(label=label_number(scale_cut = cut_si(""), accuracy = 1)(after_stat(count))), 
            position=position_stack(vjust=1), hjust=-0.1)+ 
  scale_y_continuous(breaks = c(0, 50000, 100000, 150000, 200000), 
                     labels = label_number(scale_cut = cut_si(""), accuracy = 1),
                     limits = c(0, 225000),  # Adjusts the upper limit to 110% of the max value
                     expand = expansion(mult = c(0, 0.05)))   # Adds some space at the top
g

options(repr.plot.width=12, repr.plot.height=6)
library(scales)
g  = ggplot(meta_sub, aes(x = celltype_final, fill = celltype_final)) +
  geom_bar(stat = "count", color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() + 
  geom_text(stat='count', aes(label=label_number(scale_cut = cut_si(""), accuracy = 1)(after_stat(count))), 
            position=position_stack(vjust=1), hjust=-0.1)+ 
  scale_y_continuous(breaks = c(0, 10000, 50000, 100000, 150000, 200000), 
                     labels = label_number(scale_cut = cut_si(""), accuracy = 1),
                     limits = c(0, 225000),  # Adjusts the upper limit to 110% of the max value
                     expand = expansion(mult = c(0, 0.05))) +  # Adds some space at the top
  facet_wrap(~major, scales = "free")
g

options(repr.plot.width=6, repr.plot.height=6)
sample = meta[sample(1:nrow(meta),100000),]

g=ggplot(sample, aes(x = umap_x, y=umap_y, color = celltype_final)) +
  geom_point(alpha=.8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() +
  theme(legend.position = "none") 
g
pdf("celltype_umap.pdf", width=6, height =6)
print(g)
dev.off()

library(scales)
options(repr.plot.width=12, repr.plot.height=6)


g <- ggplot(meta, aes(x = celltype_final, fill = celltype_final)) +
  geom_bar(stat = "count", color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() + 
  geom_text(stat='count', aes(label=label_number(scale_cut = cut_si(""), accuracy = 1)(after_stat(count))), 
            position=position_stack(vjust=1), hjust=-0.1) + 
  scale_y_continuous(breaks = c(0, 50000, 100000, 150000, 200000), 
                     labels = label_number(scale_cut = cut_si(""), accuracy = 1),
                     limits = c(0, 225000),  # Adjusts the upper limit to 110% of the max value
                     expand = expansion(mult = c(0, 0.05))) +  # Adds some space at the top
  facet_wrap(~major, scales = "free") +
  theme(legend.position = "none")  # Remove the legend entirely

# Display the plot
print(g)


pdf("celltype_barplot.pdf", width=12, height =6)
print(g)
dev.off()

options(repr.plot.width=5.5, repr.plot.height=12)

meta_counts <- meta %>%
  count(celltype_final)
ggplot(meta_counts, aes(x = reorder(celltype_final, n), y = n, color = celltype_final)) +
  geom_point(size = 3) + ylim(c(0,200000))+
  scale_color_manual(values = color_vector) +theme(legend.position = "none") +
  theme_minimal() +
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() +
  geom_text(aes(label = n), hjust = -0.15)+theme(legend.position = "none") 


options(repr.plot.width=10, repr.plot.height=14)

# Use a constant y-value for all points
meta_counts$y_position <- 1

legend_plot <- ggplot(meta_counts, aes(x = reorder(celltype_final, n), y = y_position, color = celltype_final)) +
  geom_point(size = 5) + # Create points for the legend
  scale_color_manual(values = color_vector) +
  theme_void() + # Remove axes and background
  theme(legend.position = "none") + # Hide default legend
  geom_text(aes(label = paste(celltype_final, " (", n, " cells)", sep = "")),
            size = 5, hjust = -.05, nudge_x = 0,nudge_y = 0) + # Add text labels for each cell type and count
  coord_flip() + # Flip coordinates for better readability
  labs(title = "Cell Type Legend") # Add a title for the legend plot

# Display the legend plot
print(legend_plot)

rownames(corrs)
clustering <- hclust(corrs)
ordered_celltypes <- clustering$labels[clustering$order]
ordered_celltypes

hclust(as.data.frame(corrs))

clustering

options(repr.plot.width=10, repr.plot.height=14)

# Use a constant y-value for all points
meta_counts$y_position <- 1

legend_plot <- ggplot(meta_counts, aes(x = reorder(celltype_final, n), y = y_position, color = celltype_final)) +
  geom_point(size = 5) + # Create points for the legend
  scale_color_manual(values = color_vector) +
  theme_void() + # Remove axes and background
  theme(legend.position = "none") + # Hide default legend
  geom_text(aes(label = paste(celltype_final, " (", n, " cells)", sep = "")),
            size = 5, hjust = -.05, nudge_x = 0,nudge_y = 0) + # Add text labels for each cell type and count
  coord_flip() + # Flip coordinates for better readability
  labs(title = "Cell Type Legend") # Add a title for the legend plot

# Display the legend plot
print(legend_plot)

options(repr.plot.width=6, repr.plot.height=5)
meta$Modality = "ATAC"
meta$Modality[which(meta$batch=="Female")] = "Multiome"
color_vector <- setNames(modality_color$Color, modality_color$Modality)

sample = meta[sample(1:nrow(meta),100000),]

ggplot(sample, aes(x = umap_x, y=umap_y, color = Modality)) +
  geom_point(alpha=.8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() 

options(repr.plot.width=6, repr.plot.height=4.8)
meta$Sex = "Male"
meta$Sex[which(meta$batch=="Female")] = "Female"
color_vector <- setNames(sex_color$Color, sex_color$Sex)

sample = meta[sample(1:nrow(meta),100000),]

g = ggplot(sample, aes(x = umap_x, y=umap_y, color = Sex)) +
  geom_point(alpha=.8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() 
g
pdf("sex_umap.pdf", width=6, height =5.5)
print(g)
dev.off()

dev.off()

color_vector <- setNames(sex_color$Color, sex_color$Sex)


options(repr.plot.width=4.5, repr.plot.height=12)

ggplot(meta, aes(x = celltype_final, fill = Sex)) +
  geom_bar(stat = "count", position = "fill",color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal()  + 
  labs(title = "Cells per Sex",
       x = "Cell Type",
       y = "Cells") +
  coord_flip() 


color_vector <- setNames(region_color$Color, region_color$Region)


options(repr.plot.width=5.5, repr.plot.height=12)

ggplot(meta, aes(x = celltype_final, fill = region_name)) +
  geom_bar(stat = "count", position = "fill",color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal()  + 
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() 


options(repr.plot.width=6, repr.plot.height=5)
color_vector <- setNames(age_color$Color, age_color$Age)
meta$age = factor(meta$age, levels = c("2mo", "9mo", "18mo"))
sample = meta[sample(1:nrow(meta),80000),]

ggplot(sample, aes(x = umap_x, y=umap_y, color = age)) +
  geom_point(alpha=.8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() 

# Split the data by age
meta_split <- split(meta, meta$age)

# Sample 50,000 cells from each age group
sampled_meta <- do.call(rbind, lapply(meta_split, function(x) x[sample(1:nrow(x), min(20000, nrow(x))),]))

# Set the factor levels for age
sampled_meta$age <- factor(sampled_meta$age, levels = c("2mo", "9mo", "18mo"))
sampled_meta = sampled_meta[sample(1:nrow(sampled_meta), nrow(sampled_meta)),]
# Plot using ggplot2
g = ggplot(sampled_meta, aes(x = umap_x, y = umap_y, color = age)) +
  geom_point(alpha = 0.8, size = 0.2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic()+
  theme(legend.position = "none") 
                                      
pdf("age_UMAP.pdf",width=6, height=6)
print(g)
dev.off()


meta$is_prog = meta$major
meta$is_prog[grep("IMN", meta$celltype_final)] = paste(meta$celltype_final[grep("IMN", meta$celltype_final)])
meta$is_prog[grep("IOL", meta$celltype_final)] = "IOL NN"
table(meta$is_prog)

meta_sub = meta[which(meta$celltype_final %in% c("DG-PIR Ex IMN","OB-STR-CTX Inh IMN",
                                                 "IOL NN", "OPC NN", "Oligo NN", 
                                                 "Astro-TE NN", "DG Glut")),]

options(repr.plot.width=6, repr.plot.height=4)

ggplot(meta_sub, aes(x = celltype_final, fill = age)) +
  geom_bar(stat = "count", position = "fill",color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() 

head(meta_sub)

over5k = names(table(meta$celltype_final)[which(table(meta$celltype_final)>5000)])

meta_sub = meta[which(meta$celltype_final %in% over5k),]

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)  # Load the purrr package for map functions
library(broom)  # broom is useful for tidying up statistical test results

# Calculate the proportions within each sample
meta_percentage <- meta_sub %>%
  group_by(sample, major, celltype_final, age, Sex) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(sample) %>%
  mutate(percentage = n / sum(n))

# Calculate weights for each age group
sample_weights <- meta_percentage %>%
  group_by(age) %>%
  summarise(weight = 1 / n_distinct(sample))

# Merge weights into the percentage data
meta_percentage <- meta_percentage %>%
  left_join(sample_weights, by = "age") %>%
  mutate(weighted_percentage = percentage * weight)


# Perform a Wilcoxon rank-sum test for each cell type
significance_tests <- meta_percentage %>%
  filter(age %in% c("2mo", "18mo")) %>%
  group_by(celltype_final) %>%
  summarise(p_value = wilcox.test(weighted_percentage[age == "2mo"],
                                  weighted_percentage[age == "18mo"])$p.value,
            direction = ifelse(mean(weighted_percentage[age == "18mo"]) >
                               mean(weighted_percentage[age == "2mo"]),
                               "Increase", "Decrease"),
            .groups = 'drop') %>%
  arrange(p_value, direction)
significance_tests$sig = -log(significance_tests$p_value )
significance_tests$sig[which(significance_tests$direction=="Decrease")] = -significance_tests$sig[which(significance_tests$direction=="Decrease")]
significance_tests = significance_tests[order(-significance_tests$sig),]
significance_tests
# Apply the significance-based ordering
meta_percentage <- meta_percentage %>%
  mutate(celltype_final = factor(celltype_final, levels = significance_tests$celltype_final))

# Plot the weighted percentage of each cell type within each sample
g <- ggplot(meta_percentage, aes(x = celltype_final, y = weighted_percentage, fill = age)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Celltype Age Proportion",
       x = "Cell Type",
       y = "Weighted Proportion of Cells") +
  coord_flip() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

# Display the plot
print(g)

significance_tests$celltype_final = factor(significance_tests$celltype_final,
                                           levels = significance_tests$celltype_final[order(significance_tests$sig, decreasing = T)])
g2 = ggplot(significance_tests, aes(x = celltype_final, y = sig, fill = direction)) +
  geom_bar(stat = "identity") + coord_flip() + 
  labs(title = "Celltype proportion changes in age",
       x = "Cell Type",
       y = "-log(P-value)*direction")

g2

svg("Celltype_age_plot.svg", height = 6, width = 5)
print(g)
dev.off()

options(repr.plot.width=12, repr.plot.height=6)

library(patchwork)
combined_plot <- g + g2

# Display the combined plot
print(combined_plot)

pdf("celltype_age_proportion.pdf", width=12, height =6)
print(combined_plot)
dev.off()



head(significance_tests)
head(meta_percentage)

library(dplyr)

# Calculate the total number of cells per sample
sample_totals <- meta %>%
  group_by(sample) %>%
  summarise(total_cells = n())

# Calculate the proportion of each cell type within each sample
proportions <- meta %>%
  group_by(sample, age, celltype_final) %>%
  summarise(count = n(), .groups = 'drop') %>%
  left_join(sample_totals, by = "sample") %>%
  mutate(proportion = count / total_cells)


head(proportions)

# List to store the results
results <- list()

# Loop through each cell type
celltypes <- unique(proportions$celltype_final)
for (ct in celltypes) {
  # Filter the data for the current cell type
  ct_data <- proportions %>%
    filter(celltype_final == ct, age %in% c("2mo", "18mo"))
  
  # Perform the T-Test
  test_result <- t.test(proportion ~ age, data = ct_data)
  
  # Store the result
  results[[ct]] <- data.frame(
    celltype = ct,
    p_value = test_result$p.value,
    statistic = test_result$statistic,
    direction = ifelse(mean(ct_data$proportion[ct_data$age == "18mo"]) > 
                       mean(ct_data$proportion[ct_data$age == "2mo"]),
                       "Increase", "Decrease")
  )
}

# Combine results into a single data frame
significance_results <- do.call(rbind, results)


# List to store the results
results <- list()

# Loop through each cell type
celltypes <- unique(proportions$celltype_final)
for (ct in celltypes) {
    cat(ct)
  # Filter the data for the current cell type
  ct_data <- proportions %>%
    filter(celltype_final == ct, age %in% c("2mo", "18mo"))
  
  # Perform the Wilcoxon Rank-Sum Test
  test_result <- wilcox.test(proportion ~ age, data = ct_data)
  
  # Store the result
  results[[ct]] <- data.frame(
    celltype = ct,
    p_value = test_result$p.value,
    statistic = test_result$statistic,
    direction = ifelse(mean(ct_data$proportion[ct_data$age == "18mo"]) > 
                       mean(ct_data$proportion[ct_data$age == "2mo"]),
                       "Increase", "Decrease")
  )
}

# Combine results into a single data frame
significance_results <- do.call(rbind, results)


head(significance_results[order(significance_results$p_value),])

significance_results = (significance_results[order(significance_results$p_value),])
head(significance_results)

meta_percentage$celltype_final=  factor(meta_percentage$celltype_final, levels = significance_results$celltype)

head(meta_percentage)

g

significance_results = significance_results[order(significance_results$statistic),]
meta_percentage$celltype_final=  factor(meta_percentage$celltype_final, levels = significance_results$celltype)
g <- ggplot(meta_percentage, aes(x = celltype_final, y = weighted_percentage, fill = age)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Celltype Age Proportion",
       x = "Cell Type",
       y = "Weighted Proportion of Cells") +
  coord_flip() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

# Display the plot
print(g)

significance_results$celltype = factor(significance_results$celltype,
                                           levels = significance_results$celltype)
g2 = ggplot(significance_results, aes(x = celltype, y = -statistic, fill = direction)) +
  geom_bar(stat = "identity") + coord_flip() + 
  labs(title = "Celltype proportion changes in age",
       x = "Cell Type",
       y = "-log(P-value)*direction")
g2


options(repr.plot.width=4.5, repr.plot.height=12)

ggplot(meta[which(meta$batch=="Male"),], aes(x = celltype_final, fill = age)) +
  geom_bar(stat = "count", position = "fill",color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal()  + 
  labs(title = "Male",
       x = "Cell Type",
       y = "Cells") +
  coord_flip() 


ggplot(meta[which(meta$batch=="Female"),], aes(x = celltype_final, fill = age)) +
  geom_bar(stat = "count", position = "fill",color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal()  + 
  labs(title = "Female",
       x = "Cell Type",
       y = "Cells") +
  coord_flip() 




options(repr.plot.width=18, repr.plot.height=10)
color_vector <- setNames(celltype_color$Color, celltype_color$CellType)
sample = meta[sample(1:nrow(meta),100000),]

ggplot(sample, aes(x = umap_x, y=umap_y, color = celltype_final)) +
  geom_point(alpha=.8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() 

meta$region_name = ""
meta$region_name[which(meta$region=="HCA")] = "Anterior_Hippocampus"
meta$region_name[which(meta$region=="HCP")] = "Posterior_Hippocampus"
meta$region_name[which(meta$region=="ENT")] = "Entorhinal_Cortex"
meta$region_name[which(meta$region=="AMY")] = "Amygdala"
meta$region_name[which(meta$region=="FC")] = "Frontal_Cortex"
meta$region_name[which(meta$region=="NAC")] = "Nucleus_accumbens"
meta$region_name[which(meta$region=="CP")] = "Caudate_Putamen"
meta$region_name[which(meta$region=="RLP")] = "PAG-PCG"


options(repr.plot.width=7, repr.plot.height=5)
color_vector <- setNames(region_color$Color, region_color$Region)
sample = meta[sample(1:nrow(meta),100000),]

g = ggplot(sample, aes(x = umap_x, y=umap_y, color = region_name)) +
  geom_point(alpha=.8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() 

options(repr.plot.width=7, repr.plot.height=5)
color_vector <- setNames(region_color$Color, region_color$Region)

g1= ggplot(sample, aes(x = umap_x, y=umap_y, color = region_name)) +
  geom_point(alpha=.8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() 
g1

ggplot(meta, aes(x = region_name, fill = region_name)) +
  geom_bar(stat = "count", color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + ylim(c(0,175000)) + 
  labs(title = "Cells per Region",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() +
  geom_text(stat='count', aes(label=after_stat(count)), position=position_stack(vjust=1), hjust=-0.1)

library(scales)
options(repr.plot.width=7
        , repr.plot.height=8)
meta$region_name=factor(meta$region_name, levels = sort(unique(meta$region_name), decreasing =T))
g=ggplot(meta, aes(x = region_name, fill = region_name)) +
  geom_bar(stat = "count", color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + ylim(c(0,195000)) + 
  labs(title = "Cells per Region",
       x = "Brain Region",
       y = "Number of Cells", fill = "Brain Region") +
  coord_flip() +
    geom_text(stat='count', aes(label=label_number(scale_cut = cut_si(""), accuracy = 1)(after_stat(count))), 
            position=position_stack(vjust=1), hjust=-0.1) +
  theme(legend.pos = "bottom",
    axis.title = element_text(size = 14),    # Increase axis title size
    axis.text = element_text(size = 12)      # Increase axis text size
  )+ 
  scale_y_continuous(breaks = c(0, 50000, 100000, 150000), 
                     labels = label_number(scale_cut = cut_si(""), accuracy = 1),
                     limits = c(0, 195000),  # Adjusts the upper limit to 110% of the max value
                     expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~Sex, nrow = 2)# Adds some space at the top
g



pdf("num_cell.pdf", width=6, height =5)
print(g)
dev.off()



pdf("region_umap.pdf", width=6.5, height =5)
print(g1)
dev.off()

g

dev.off()


ameta = meta
ameta$age_rep = paste(ameta$age, ameta$rep,ameta$region, ameta$batch)

head(ameta)

tab = table(ameta$celltype_final,ameta$age_rep)
#head(tab)
tab = sweep(tab,2,colSums(tab),'/')
#head(tab)
melted = melt(tab)
#head(melted)

melted$rep = sapply(strsplit(as.character(melted$Var2), " "), "[[", 2)
melted$age = sapply(strsplit(as.character(melted$Var2), " "), "[[", 1)
melted$region = sapply(strsplit(as.character(melted$Var2), " "), "[[", 3)
melted$age = factor(melted$age, levels = c("2mo", "9mo", "18mo"))
melted$sex = sapply(strsplit(as.character(melted$Var2), " "), "[[", 4)


unique(melted$Var1)

color_vector

color_vector <- setNames(age_color$Color, age_color$Age)

g2 = ggplot(melted[which(melted$Var1 == "OB-STR-CTX Inh IMN" ),]) + 
  geom_col(aes(factor(region),value,fill=age, color = rep ),position="dodge") + 
  xlab("Brain region") + ylab("Fraction of all cells") + 
  guides(fill=guide_legend(title="Age Group"), color = guide_legend(title="Rep"))  +
  theme_bw() +ggtitle("IOL NN") + facet_wrap(~sex, nrow = 2)+ scale_color_manual(values = c("black", "black", "black"))+ scale_fill_manual(values = color_vector) 
#svg("OB_IMN_celltype_frac.svg", height =5, width =5)
#g2
#dev.off()
g2

color_vector <- setNames(age_color$Color, age_color$Age)

g2 = ggplot(melted[which(melted$Var1 == "OPC NN" ),]) + 
  geom_col(aes(factor(region),value,fill=age, color = rep ),position="dodge") + 
  xlab("Brain region") + ylab("Fraction of all cells") + 
  guides(fill=guide_legend(title="Age Group"), color = guide_legend(title="Rep"))  +
  theme_bw() +ggtitle("IOL NN") + facet_wrap(~sex, nrow = 2)+ scale_color_manual(values = c("black", "black", "black"))+ scale_fill_manual(values = color_vector) 
#svg("OB_IMN_celltype_frac.svg", height =5, width =5)
#g2
#dev.off()
g2

head(melted)



cl = "IOL NN"

pdf("cell_fraction.pdf")
for (cl in unique(melted$Var1)) {

g2 = ggplot(melted[which(melted$Var1 == cl ),]) + 
  geom_col(aes(factor(region),value,fill=age, color =rep ),position="dodge") + 
  xlab("Cell Types") + ylab("Fraction") + 
  guides(fill=guide_legend(title="Age Group"), color = guide_legend(title="Rep"))  +
  theme_bw() +ggtitle(cl) + facet_wrap(~sex, nrow = 2)+ scale_fill_manual(values = color_vector) +scale_color_manual(values = c("black", "black", "black"))

print(g2)
    }

dev.off()


