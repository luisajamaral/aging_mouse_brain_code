library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)

setwd("~/projects/combined_all/Figures/Figure1-overview//")

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

mer_meta = read.csv("../Figure2-IMN-IOL/MERFISH_meta.csv")

head(mer_meta)

mer_meta$Modality = "MERFISH"

nrow(mer_meta)

table(mer_meta$sample)

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

meta$major = sapply(strsplit(as.character(meta$celltype_final), " "), function(x) tail(x, 1))
meta$major[which(meta$major%in%"Neur")] = "Gaba"
meta$major[which(meta$major%in%"Gaba-Chol")] = "Gaba"
meta$major[which(meta$major%in%"Glut-Sero")] = "Glut"
meta$major[which(meta$major%in%"Dopa")] = "Glut"
meta$major[which(meta$major%in%"IMN")] = "NN"



head(meta)


unique(mer_meta$celltype_final)[which(unique(mer_meta$celltype_final) %in% unique(meta$celltype_final) )]

unique(mer_meta$celltype_final)[which(!unique(mer_meta$celltype_final) %in% unique(meta$celltype_final) )]

head(mer_meta$batch
    )

meta$Modality = "scATAC"
meta$Modality[which(meta$batch=="Female")] = "Multiome"
table(meta$Modality)

mer_meta = mer_meta[-which(mer_meta$batch == "2moRep1"),]

mer_meta$rep[which(mer_meta$batch == "2moRep3")] = "Rep1"

mer_meta = mer_meta[-which(mer_meta$celltype_final %in% names(table(mer_meta$celltype_final)[which(table(mer_meta$celltype_final)<200)])),]

names(table(mer_meta$celltype_final)[which(table(mer_meta$celltype_final)<200)])

table(mer_meta$celltype_final)

mer_meta$sample = mer_meta$batch

head(mer_meta)

mer_meta$region = "all"

mer = mer_meta[,c('celltype_final','age', 'rep', 'Modality',"sample","region")]



mt = meta[,c('celltype_final','age', 'rep', 'Modality', "sample", "region")]


all = rbind(mt, mer)

head(all)

library(stringr)
all$rep = str_to_title(all$rep)

library(matrixStats) # rowMins()
good = names(rowMins(table(all$celltype_final, all$Modality))[which(rowMins(table(all$celltype_final, all$Modality))>225)])
good = good[-grep("CA1-ProS Glut", good)]

#write.table(all, file = "all_meta.csv", sep = ",", row.names = F) 

all = all[which(all$celltype_final %in% unique(mer$celltype_final)),]
all = all[which(all$celltype_final %in% good),]

head(all)

all = as.data.frame(all)

# Check column names and structure of `all`
colnames(all)
str(all)

# Ensure the columns are in the correct format
all <- all %>%
  mutate(
    celltype_final = as.character(celltype_final),
    age = as.character(age),
    Modality = as.character(Modality),
    rep = as.character(rep)
  )

# Count cells for each combination of celltype, age, modality, and rep
cell_counts <- all %>%
  group_by(celltype_final, age, Modality, rep) %>%
  summarise(n = n(), .groups = 'drop')

# Display the resulting table
head(cell_counts)


library(dplyr)

# Step 1: Calculate total cell count per replicate
celltype_fractions <- cell_counts %>%
  group_by(rep, age, Modality) %>%
  mutate(total_cells = sum(n)) %>%  # Total cell count within each rep, age, and modality
  ungroup()

# Step 2: Calculate fraction of each cell type
celltype_fractions <- celltype_fractions %>%
  mutate(fraction = n / total_cells)

# Step 3: Calculate weights for each age group and modality
sample_weights <- celltype_fractions %>%
  group_by(age, Modality) %>%
  summarise(weight = 1 / n_distinct(rep), .groups = 'drop')

# Step 4: Merge weights into the cell type fractions and calculate weighted fractions
weighted_fractions <- celltype_fractions %>%
  left_join(sample_weights, by = c("age", "Modality")) %>%
  mutate(weighted_fraction = fraction * weight)


# Perform Wilcoxon rank-sum test for each cell type across modalities using the `fraction` column
significance_tests <- weighted_fractions %>%
  filter(age %in% c("2mo", "18mo")) %>%
  group_by(celltype_final) %>%
  summarise(
    # Wilcoxon test using `fraction`
    p_value = wilcox.test(fraction[age == "2mo"],
                          fraction[age == "18mo"])$p.value,
    direction = ifelse(mean(fraction[age == "18mo"]) >
                       mean(fraction[age == "2mo"]),
                       "Increase", "Decrease"),
    .groups = 'drop'
  ) %>%
  mutate(
    sig = -log(p_value),
    sig = ifelse(direction == "Decrease", -sig, sig)
  ) %>%
  arrange( -sig)

# Display the results
print(significance_tests)



# Join the significance tests with the weighted data
plot_data <- weighted_data %>%
  left_join(significance_tests, by = c("celltype_final"))

plot_data$age = factor(plot_data$age , levels = rev(c("2mo", "9mo", "18mo")))

# Plotting
library(ggplot2)

plot_data$celltype_final = factor(plot_data$celltype_final,
                                           levels = significance_tests$celltype_final[order(significance_tests$sig, decreasing = T)])

plot <- ggplot(plot_data, aes(x = celltype_final, y = weighted_percentage, fill = age)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Weighted Cell Type Proportions by Age",
       x = "Cell Type",
       y = "Weighted Proportion of Cells") +
  coord_flip() + 
  facet_grid(~Modality) +
  theme(legend.position = "right") 

# Display the plot
print(plot)


significance_tests$celltype_final = factor(significance_tests$celltype_final,
                                           levels = significance_tests$celltype_final[order(significance_tests$sig, decreasing = T)])
g2 = ggplot(significance_tests, aes(x = celltype_final, y = sig, fill = direction)) +
  geom_bar(stat = "identity") + coord_flip() + 
  scale_fill_manual(values = c("royalblue","tomato")) +
  labs(title = "Significance",
       x = "Cell Type",
       y = "-log(P-value)*direction")+
     geom_text(data = significance_tests, aes(x = celltype_final, y = (sig * 1.05)+0.2,
                                           label = ifelse(p_value < 0.05, "*", "")),
            color = "black", size = 5, hjust = 1, vjust = .75)+theme_minimal() +theme(axis.title.y = element_blank(), axis.text.y = element_blank())


g2


options(repr.plot.width = 15, repr.plot.height = 7)
library(scales)


# Adjust the width of the individual plots within the combined plot
combined_plot <- (plot | g2 | count) + 
  plot_layout(widths = c(15, 6, 5), guides = "collect") & 
  theme(legend.position = "right", text = element_text(color = "black"))

# Display the combined plot
print(combined_plot)


setwd("~/projects/combined_all/Figures/Figure1-overview/")

pdf("updated_cell_composition_plot.pdf", height = 7, width = 15)
print(combined_plot)
dev.off()

dev.off()

modality_vector = c(modality_vector)

modality_vector["Multiome"] = '#7D4195'
modality_vector["MERFISH"] = "#13ba1e"
modality_vector["scATAC"] = '#EF7D1A'




modality_vector

head(all)

options(repr.plot.width = 10, repr.plot.height = 8)

all <- all %>%
  mutate(celltype_final = factor(celltype_final, levels = significance_tests$celltype_final))


color_vector <- setNames(celltype_color$Color, celltype_color$CellType)
modality_vector <- setNames(modality_color$Color, modality_color$Modality)

modality_vector = c(modality_vector)

modality_vector["Multiome"] = '#7D4195'
modality_vector["MERFISH"] = "#13ba1e"
modality_vector["scATAC"] = '#EF7D1A'



count = ggplot(all, aes(x = celltype_final,  fill = celltype_final, color = Modality)) +
  geom_bar(stat = "count", position = "stack") +
  scale_fill_manual(values = color_vector) +
  scale_color_manual(values = modality_vector) +

  theme_minimal() + 
  labs(title = "Cell Counts",
       x = "Cell Type",
       y = "Count") +
  coord_flip() + 
 # facet_grid(~Modality) +
  theme(legend.position = "none") +theme(axis.title.y = element_blank(), axis.text.y = element_blank())+
  guides(fill = "none")+ scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()))


count

# Calculate log2 fold change between 18mo and 2mo for each cell type and modality
fold_change_df <- weighted_data %>%
  filter(age %in% c("2mo", "18mo")) %>%
  group_by(celltype_final, Modality) %>%
  summarise(
    mean_2mo = mean(weighted_percentage[age == "2mo"], na.rm = TRUE),
    mean_18mo = mean(weighted_percentage[age == "18mo"], na.rm = TRUE),
    log2_fold_change = log2(mean_18mo / mean_2mo),
    .groups = 'drop'
  ) %>%
  arrange(Modality, desc(log2_fold_change))

# Define a function to create a fold change plot for a given modality
create_fc_plot <- function(data, modality) {
  ggplot(data, aes(x = celltype_final, y = log2_fold_change, fill = log2_fold_change > 0)) +
    geom_bar(stat = "identity", aes(color = log2_fold_change > 0), size = 0.8) +
    scale_fill_manual(values = c("red", "blue"), labels = c("Decrease", "Increase")) +
    scale_color_manual(values = c("red", "blue"), guide = "none") +  # Outline for bars
    theme_minimal() +
    labs(title = paste("Log2 Fold Change in", modality),
         y = "Log2 Fold Change",
         fill = "Direction") +
    coord_flip() +
    theme(legend.position = "right", axis.title.y = element_blank(), axis.text.y = element_blank())
}

# Create fold change plots for each modality
fc_plot_scATAC <- create_fc_plot(fold_change_df %>% filter(Modality == "scATAC"), "scATAC")
fc_plot_Multiome <- create_fc_plot(fold_change_df %>% filter(Modality == "Multiome"), "Multiome")
fc_plot_MERFISH <- create_fc_plot(fold_change_df %>% filter(Modality == "MERFISH"), "MERFISH")

# Combine the plots using patchwork
library(patchwork)

# Place the fold change plots next to the weighted proportion plot
combined_plot <- plot | (fc_plot_scATAC | fc_plot_Multiome | fc_plot_MERFISH) +  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

# Display the combined plot
print(combined_plot)


# Calculate Cohen's d for effect size between 18mo and 2mo for each cell type and modality
effect_size_df <- weighted_data %>%
  filter(age %in% c("2mo", "18mo")) %>%
  group_by(celltype_final, Modality) %>%
  summarise(
    mean_2mo = mean(weighted_percentage[age == "2mo"], na.rm = TRUE),
    mean_18mo = mean(weighted_percentage[age == "18mo"], na.rm = TRUE),
    sd_2mo = sd(weighted_percentage[age == "2mo"], na.rm = TRUE),
    sd_18mo = sd(weighted_percentage[age == "18mo"], na.rm = TRUE),
    n_2mo = sum(age == "2mo"),
    n_18mo = sum(age == "18mo"),
    # Calculate pooled standard deviation
    pooled_sd = sqrt(((n_2mo - 1) * sd_2mo^2 + (n_18mo - 1) * sd_18mo^2) / (n_2mo + n_18mo - 2)),
    # Cohen's d
    effect_size = (mean_18mo - mean_2mo) / pooled_sd,
    .groups = 'drop'
  ) %>%
  arrange(Modality, desc(effect_size))

# Define a function to create effect size plot for a given modality
create_effect_size_plot <- function(data, modality) {
  ggplot(data, aes(x = celltype_final, y = effect_size, fill = effect_size > 0)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("red", "blue"), labels = c("Decrease", "Increase")) +
    theme_minimal() +
    labs(title = paste("Effect Size (Cohen's d) in", modality),
         y = "Effect Size (Cohen's d)",
         fill = "Direction") +
    coord_flip() +
    theme(legend.position = "right", axis.title.y = element_blank(), axis.text.y = element_blank())
}

# Create effect size plots for each modality
effect_size_plot_scATAC <- create_effect_size_plot(effect_size_df %>% filter(Modality == "scATAC"), "scATAC")
effect_size_plot_Multiome <- create_effect_size_plot(effect_size_df %>% filter(Modality == "Multiome"), "Multiome")
effect_size_plot_MERFISH <- create_effect_size_plot(effect_size_df %>% filter(Modality == "MERFISH"), "MERFISH")

# Combine the effect size plots with the main plot
combined_plot <- plot | (effect_size_plot_scATAC | effect_size_plot_Multiome |effect_size_plot_MERFISH) +  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

# Display the combined plot
print(combined_plot)


all$age = factor(all$age, levels= c("2mo", "9mo", "18mo"))

color_vector <- setNames(age_color$Color, age_color$Age)

ggplot(all, aes(x = celltype_final, fill = age)) +
  geom_bar(stat = "count", position = "fill",color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() + facet_grid(~Modality)

# Calculate the proportions within each sample and modality
meta_percentage <- all %>%
  group_by(region, rep, celltype_final, age, Modality) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(region, rep, age, Modality) %>%
  mutate(percentage = n / sum(n))

# Calculate weights for each age group and modality
sample_weights <- meta_percentage %>%
  group_by(age, Modality) %>%
  summarise(weight = 1 / n_distinct(region, rep), .groups = 'drop')

# Merge weights into the percentage data
meta_percentage <- meta_percentage %>%
  left_join(sample_weights, by = c("age", "Modality")) %>%
  mutate(weighted_percentage = percentage * weight)

# Perform a Wilcoxon rank-sum test for each cell type across modalities
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

# Apply significance transformation and sorting
significance_tests$sig <- -log(significance_tests$p_value)
significance_tests$sig[significance_tests$direction == "Decrease"] <- 
  -significance_tests$sig[significance_tests$direction == "Decrease"]
significance_tests <- significance_tests[order(-significance_tests$sig), ]
significance_tests_all = significance_tests



meta_percentage$Modality <- factor(meta_percentage$Modality, levels = c("scATAC", "Multiome", "MERFISH"))
options(repr.plot.width = 10, repr.plot.height = 8)

g <- ggplot(meta_percentage, aes(x = celltype_final, y = weighted_percentage, fill = age)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Celltype Age Proportion",
       x = "Cell Type",
       y = "Weighted Proportion of Cells") +
  coord_flip() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(~Modality)

# Display the plot
print(g)

library(purrr)  # Load the purrr package for map functions
library(broom)  # broom is useful for tidying up statistical test results
library(dplyr)
library(ggplot2)

# Calculate the proportions within each sample and modality
meta_percentage <- all %>%
  group_by(sample, celltype_final, age, Modality) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(sample, Modality) %>%
  mutate(percentage = n / sum(n))

# Calculate weights for each age group and modality
sample_weights <- meta_percentage %>%
  group_by(age, Modality) %>%
  summarise(weight = 1 / n_distinct(sample), .groups = 'drop')

# Merge weights into the percentage data
meta_percentage <- meta_percentage %>%
  left_join(sample_weights, by = c("age", "Modality")) %>%
  mutate(weighted_percentage = percentage * weight)

# Perform a Wilcoxon rank-sum test for each cell type across modalities
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

# Apply significance transformation and sorting
significance_tests$sig <- -log(significance_tests$p_value)
significance_tests$sig[significance_tests$direction == "Decrease"] <- 
  -significance_tests$sig[significance_tests$direction == "Decrease"]
significance_tests <- significance_tests[order(-significance_tests$sig), ]

# Apply the significance-based ordering to cell types
meta_percentage <- meta_percentage %>%
  mutate(celltype_final = factor(celltype_final, levels = significance_tests_all$celltype_final))


# Plot the results
meta_percentage$Modality <- factor(meta_percentage$Modality, levels = c("scATAC", "Multiome", "MERFISH"))
options(repr.plot.width = 10, repr.plot.height = 8)

g <- ggplot(meta_percentage, aes(x = celltype_final, y = weighted_percentage, fill = age)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Celltype Age Proportion",
       x = "Cell Type",
       y = "Weighted Proportion of Cells") +
  coord_flip() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(~Modality)

# Display the plot
print(g)


color_vector <- setNames(age_color$Color, age_color$Age)

# Apply the significance-based ordering to cell types
meta_percentage <- meta_percentage %>%
  mutate(celltype_final = factor(celltype_final, levels = ord))


# Plot the results
meta_percentage$Modality <- factor(meta_percentage$Modality, levels = c("scATAC", "Multiome", "MERFISH"))
options(repr.plot.width = 10, repr.plot.height = 8)

g <- ggplot(meta_percentage, aes(x = celltype_final, y = weighted_percentage, fill = age)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Celltype Age Proportion",
       x = "Cell Type",
       y = "Weighted Proportion of Cells") +
  coord_flip() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(~Modality)

# Display the plot
print(g)


options(repr.plot.width = 10, repr.plot.height = 8)

all$Modality <- factor(all$Modality, levels = c("scATAC", "Multiome", "MERFISH"))
all = all[which(all$celltype_final%in% meta_percentage$celltype_final),]
color_vector <- setNames(celltype_color$Color, celltype_color$CellType)
all <- all %>%
  mutate(celltype_final = factor(celltype_final, levels = ord))


g2 = ggplot(all, aes(x = celltype_final,  fill = celltype_final)) +
  geom_bar(stat = "count", position = "stack") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Cell Counts",
       x = "Cell Type",
       y = "Count") +
  coord_flip() + 
  facet_grid(~Modality) +
  theme(legend.position = "none") 

svg("Cell_counts_by_modality.svg",width = 10, height = 8)
g2
dev.off()

g2

meta_percentage <- meta_percentage %>%
  mutate(celltype_final = factor(celltype_final, levels = ord))
color_vector <- setNames(age_color$Color, age_color$Age)

# Plot the results
meta_percentage$Modality <- factor(meta_percentage$Modality, levels = c("scATAC", "Multiome", "MERFISH"))
options(repr.plot.width = 10, repr.plot.height = 8)

g <- ggplot(meta_percentage, aes(x = celltype_final, y = weighted_percentage, fill = age)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Celltype Age Proportion",
       x = "Cell Type",
       y = "Weighted Proportion of Cells") +
  coord_flip() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(~Modality)

# Display the plot
print(g)
svg("Cell_proportions_by_modality.svg",width = 10, height = 8)
g
dev.off()


options(repr.plot.width = 10, repr.plot.height = 8)

all <- all %>%
  mutate(celltype_final = factor(celltype_final, levels = significance_tests_all$celltype_final))


color_vector <- setNames(celltype_color$Color, celltype_color$CellType)


g2 = ggplot(all, aes(x = celltype_final,  fill = celltype_final)) +
  geom_bar(stat = "count", position = "stack") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Cell Counts",
       x = "Cell Type",
       y = "Count") +
  coord_flip() + 
  facet_grid(~Modality) +
  theme(legend.position = "none") 

g2

color_vector

# Apply the significance-based ordering
color_vector <- setNames(age_color$Color, age_color$Age)

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

options(repr.plot.width=10, repr.plot.height=8)

meta_percentage$Modality = factor(meta_percentage$Modality, levels=c("scATAC", "Multiome", "MERFISH"))
g <- ggplot(meta_percentage, aes(x = celltype_final, y = weighted_percentage, fill = age)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Celltype Age Proportion",
       x = "Cell Type",
       y = "Weighted Proportion of Cells") +
  coord_flip() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
facet_grid(~Modality)

# Display the plot
print(g)

head(meta_percentage)

meta_percentage

head(significance_tests)

significance_tests$celltype_final = factor(significance_tests$celltype_final,
                                           levels = significance_tests$celltype_final[order(significance_tests$sig, decreasing = T)])
g2 = ggplot(significance_tests, aes(x = celltype_final, y = sig, fill = direction)) +
  geom_bar(stat = "identity") + coord_flip() + 
  labs(title = "Celltype proportion changes in age",
       x = "Cell Type",
       y = "-log(P-value)*direction")

g2

head(meta_percentage)

# Calculate effect size (Cohen's d) between 2mo and 18mo for each cell type
effect_size_df <- meta_percentage %>%
  filter(age %in% c("2mo", "18mo")) %>%
  group_by(celltype_final) %>%
  summarise(
    mean_2mo = mean(weighted_percentage[age == "2mo"], na.rm = TRUE),
    mean_18mo = mean(weighted_percentage[age == "18mo"], na.rm = TRUE),
    sd_2mo = sd(weighted_percentage[age == "2mo"], na.rm = TRUE),
    sd_18mo = sd(weighted_percentage[age == "18mo"], na.rm = TRUE),
    n_2mo = sum(age == "2mo"),
    n_18mo = sum(age == "18mo"),
    # Calculate pooled standard deviation
    pooled_sd = sqrt(((n_2mo - 1) * sd_2mo^2 + (n_18mo - 1) * sd_18mo^2) / (n_2mo + n_18mo - 2)),
    # Cohen's d
    effect_size = (mean_18mo - mean_2mo) / pooled_sd
  ) %>%
  arrange(desc(effect_size))

# Display the effect size for each cell type
print(effect_size_df)

library(ggplot2)

# Plot the effect sizes for each cell type
p=ggplot(effect_size_df, aes(x = celltype_final, y = effect_size, fill = effect_size > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red", "blue"), labels = c("Decrease", "Increase")) +
  theme_minimal() +
  labs(title = "Effect Sizes of Cell Type Proportions by Age",
       x = "Cell Type",
       y = "Effect Size (Cohen's d)",
       fill = "Direction") +
  coord_flip() +
  theme(legend.position = "top")
p






options(repr.plot.width=5, repr.plot.height=8)

# Calculate log2 fold change between 18mo and 2mo for each cell type
fold_change_df <- meta_percentage %>%
  filter(age %in% c("2mo", "18mo")) %>%
  group_by(celltype_final) %>%
  summarise(
    mean_2mo = mean(weighted_percentage[age == "2mo"], na.rm = TRUE),
    mean_18mo = mean(weighted_percentage[age == "18mo"], na.rm = TRUE),
    # Calculate log2 fold change
    log2_fold_change = log2(mean_18mo / mean_2mo)
  ) %>%
  arrange(desc(log2_fold_change))

# Plot the log2 fold changes for each cell type
ggplot(fold_change_df, aes(x = celltype_final, y = log2_fold_change, fill = log2_fold_change > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red", "blue"), labels = c("Decrease", "Increase")) +
  theme_minimal() +
  labs(title = "Log2 Fold Change of Cell Type Proportions by Age",
       x = "Cell Type",
       y = "Log2 Fold Change",
       fill = "Direction") +
  coord_flip() +
  theme(legend.position = "top")


options(repr.plot.width=10, repr.plot.height=8)

meta_percentage$Modality = factor(meta_percentage$Modality, levels=c("scATAC", "Multiome", "MERFISH"))
meta_percentage$celltype_final= factor(meta_percentage$celltype_final , levels = fold_change_df$celltype_final)
g <- ggplot(meta_percentage, aes(x = celltype_final, y = weighted_percentage, fill = age)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Celltype Age Proportion",
       x = "Cell Type",
       y = "Weighted Proportion of Cells") +
  coord_flip() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
facet_grid(~Modality)

# Display the plot
print(g)

# Load necessary libraries
library(ggplot2)
library(patchwork)

# Calculate log2 fold change between 18mo and 2mo for each cell type
fold_change_df <- meta_percentage %>%
  filter(age %in% c("2mo", "18mo")) %>%
  group_by(celltype_final) %>%
  summarise(
    mean_2mo = mean(weighted_percentage[age == "2mo"], na.rm = TRUE),
    mean_18mo = mean(weighted_percentage[age == "18mo"], na.rm = TRUE),
    # Calculate log2 fold change
    log2_fold_change = log2(mean_18mo / mean_2mo)
  ) %>%
  arrange(desc(log2_fold_change))

# Create the log2 fold change plot
fc_plot <- ggplot(fold_change_df, aes(x = celltype_final, y = log2_fold_change, fill = log2_fold_change > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red", "blue"), labels = c("Decrease", "Increase")) +
  theme_minimal() +
  labs(title = "Log2 Fold Change of Cell Type Proportions by Age",
       y = "Log2 Fold Change",
       fill = "Direction") +
  coord_flip() +
  theme(legend.position = "top", axis.title.y = element_blank(), axis.text.y = element_blank())

# Create the cell type age proportion plot
meta_percentage$Modality <- factor(meta_percentage$Modality, levels = c("scATAC", "Multiome", "MERFISH"))
meta_percentage$celltype_final <- factor(meta_percentage$celltype_final, levels = fold_change_df$celltype_final)

proportion_plot <- ggplot(meta_percentage, aes(x = celltype_final, y = weighted_percentage, fill = age)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Cell Type Age Proportion",
       x = "Cell Type",
       y = "Weighted Proportion of Cells") +
  coord_flip() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(~Modality)

# Combine the plots
combined_plot <- proportion_plot + fc_plot + plot_layout(ncol = 2, widths = c(2, 1))
print(combined_plot)


# Load necessary libraries
library(ggplot2)
library(patchwork)

# Set a significance threshold
significance_threshold <- 0.05  # Adjust as appropriate for your analysis

# Calculate log2 fold change and add significance indicator
fold_change_df <- meta_percentage %>%
  filter(age %in% c("2mo", "18mo")) %>%
  group_by(celltype_final) %>%
  summarise(
    mean_2mo = mean(weighted_percentage[age == "2mo"], na.rm = TRUE),
    mean_18mo = mean(weighted_percentage[age == "18mo"], na.rm = TRUE),
    # Calculate log2 fold change
    log2_fold_change = log2(mean_18mo / mean_2mo),
    # Placeholder for p-value, replace with actual test if available
    p_value = wilcox.test(weighted_percentage[age == "2mo"],
                          weighted_percentage[age == "18mo"])$p.value
  ) %>%
  mutate(significant = p_value < significance_threshold) %>%
  arrange(desc(log2_fold_change))

# Create the log2 fold change plot with highlighted significant cell types
fc_plot <- ggplot(fold_change_df, aes(x = celltype_final, y = log2_fold_change, fill = log2_fold_change > 0)) +
  geom_bar(stat = "identity", aes(color = significant), size = 0.8, show.legend = TRUE) +
  scale_fill_manual(values = c("red", "blue"), labels = c("Decrease", "Increase")) +
  scale_color_manual(values = c("black", NA), guide = "none") +  # Outline for significant bars
  theme_minimal() +
  labs(title = "Log2 Fold Change of Cell Type Proportions by Age",
       y = "Log2 Fold Change",
       fill = "Direction") +
  coord_flip() +
  theme(legend.position = "right", axis.title.y = element_blank(), axis.text.y = element_blank())

# Create the cell type age proportion plot
meta_percentage$Modality <- factor(meta_percentage$Modality, levels = c("scATAC", "Multiome", "MERFISH"))
meta_percentage$celltype_final <- factor(meta_percentage$celltype_final, levels = fold_change_df$celltype_final)

proportion_plot <- ggplot(meta_percentage, aes(x = celltype_final, y = weighted_percentage, fill = age)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Cell Type Age Proportion",
       x = "Cell Type",
       y = "Weighted Proportion of Cells") +
  coord_flip() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(~Modality) +
  theme(legend.position = "right")

# Combine the plots with legends on the right
combined_plot <- proportion_plot + fc_plot + plot_layout(ncol = 2, widths = c(2, 1))
print(combined_plot)


fold_change_df


# Load necessary libraries
library(ggplot2)
library(patchwork)

# Set a significance threshold
significance_threshold <- 0.05  # Adjust as appropriate for your analysis

# Calculate log2 fold change and add significance indicator, grouped by Modality
fold_change_df <- meta_percentage %>%
  filter(age %in% c("2mo", "18mo")) %>%
  group_by(celltype_final, Modality) %>%
  summarise(
    mean_2mo = mean(weighted_percentage[age == "2mo"], na.rm = TRUE),
    mean_18mo = mean(weighted_percentage[age == "18mo"], na.rm = TRUE),
    # Calculate log2 fold change
    log2_fold_change = log2(mean_18mo / mean_2mo),
    # Placeholder for p-value, replace with actual test if available
    p_value = wilcox.test(weighted_percentage[age == "2mo"],
                          weighted_percentage[age == "18mo"])$p.value,
    .groups = 'drop'
  ) %>%
  mutate(significant = p_value < significance_threshold) %>%
  arrange(desc(log2_fold_change))

# Ensure `celltype_final` is a factor with levels ordered by log2 fold change across all modalities
celltype_order <- unique(fold_change_df$celltype_final[order(-fold_change_df$log2_fold_change)])
fold_change_df$celltype_final <- factor(fold_change_df$celltype_final, levels = celltype_order)
meta_percentage$celltype_final <- factor(meta_percentage$celltype_final, levels = celltype_order)

# Define a function to create the log2 fold change plot for a given modality
create_fc_plot <- function(data, modality) {
  ggplot(data, aes(x = celltype_final, y = log2_fold_change, fill = log2_fold_change > 0)) +
    geom_bar(stat = "identity", aes(color = significant), size = 0.8, show.legend = modality == "scATAC") +
    scale_fill_manual(values = c("red", "blue"), labels = c("Decrease", "Increase")) +
    scale_color_manual(values = c("black", NA), guide = "none") +  # Outline for significant bars
    theme_minimal() +
    labs(title = paste("Log2 Fold Change in", modality),
         y = "Log2 Fold Change",
         fill = "Direction") +
    coord_flip() +
    theme(legend.position = "right", axis.title.y = element_blank(), axis.text.y = element_blank())
}

# Create fold change plots for each modality
fc_plot_scATAC <- create_fc_plot(fold_change_df %>% filter(Modality == "scATAC"), "scATAC")
fc_plot_Multiome <- create_fc_plot(fold_change_df %>% filter(Modality == "Multiome"), "Multiome")
fc_plot_MERFISH <- create_fc_plot(fold_change_df %>% filter(Modality == "MERFISH"), "MERFISH")

# Create the cell type age proportion plot
meta_percentage$Modality <- factor(meta_percentage$Modality, levels = c("scATAC", "Multiome", "MERFISH"))

proportion_plot <- ggplot(meta_percentage, aes(x = celltype_final, y = weighted_percentage, fill = age)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Cell Type Age Proportion",
       x = "Cell Type",
       y = "Weighted Proportion of Cells") +
  coord_flip() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(~Modality) +
  theme(legend.position = "right")

# Combine the proportion plot with the three fold change plots
combined_plot <- proportion_plot | (fc_plot_scATAC | fc_plot_Multiome | fc_plot_MERFISH)
print(combined_plot)


(meta_percentage)

combined_plot <- (proportion_plot | (fc_plot_scATAC | fc_plot_Multiome | fc_plot_MERFISH)) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

print(combined_plot)


