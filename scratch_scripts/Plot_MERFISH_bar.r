library(Seurat)
library(ggplot2)
library(dplyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)  # Load the purrr package for map functions
library(broom)  # broom is useful for tidying up statistical test results
library(ggpubr) # For adding statistical comparisons


obj = readRDS("~/projects/merfish/filtered.RDS")

head(obj@meta.data)
 
    metaf = obj@meta.data
    #metaf = metaf[which(!is.na(sub$predictions_cca)),]
    metaf = metaf[which(metaf$predictions_cca_max_score>0.25),]

    predictions <- table(metaf$sub_leiden,metaf$predictions_cca)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    dplyr::select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    #new_df = new_df[which(new_df$Freq>0.8),]
    #rownames(new_df) = new_df$Var1

    mat = match(obj@meta.data$sub_leiden, new_df$Var1, nomatch = NA)
    obj@meta.data$prediction = new_df[mat, "Var2"]
    obj@meta.data$prediction = as.character(obj@meta.data$prediction)
    obj@meta.data$prediction = new_df[mat, "Var2"]

obj@meta.data$age = obj@meta.data$batch
obj@meta.data$age = gsub("Rep1", "", obj@meta.data$age)
obj@meta.data$age = gsub("Rep2", "", obj@meta.data$age)
obj@meta.data$age = gsub("Rep3", "", obj@meta.data$age)
obj$age = factor(obj$age, levels = c("2mo", "9mo", "18mo"))
obj@meta.data$rep = obj@meta.data$batch
obj@meta.data$rep = gsub("2mo", "", obj@meta.data$rep)
obj@meta.data$rep = gsub("9mo", "", obj@meta.data$rep)
obj@meta.data$rep = gsub("18mo", "", obj@meta.data$rep)
obj$rep = factor(obj$rep, levels = c("Rep1", "Rep2", "Rep3"))

obj@meta.data$umap_x = obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$umap_y = obj@reductions$umap@cell.embeddings[,2]
obj@meta.data$celltype_final = paste(obj@meta.data$prediction)
obj@meta.data$celltype_final[which(obj@meta.data$celltype_final=="IOL")] ="IOL NN"
obj@meta.data$celltype_final[which(obj@meta.data$celltype_final=="RGL")] ="OB-STR-CTX Inh IMN"

meta= obj@meta.data

#key[which( key$CellType %in%meta$celltype_final ),]
#setwd("~/projects/combined_all/Figures/Figure2-IMN-IOL/")
key_file = "../color_scheme/updated_celltype_palette.csv"
key = read.csv(key_file)
head(key)
color_vector <- setNames(key$Color, key$CellType)
unique(meta$celltype_final)[which(! unique(meta$celltype_final) %in% key$CellType  )]

options(repr.plot.width=12, repr.plot.height=5)

FeaturePlot(obj, split.by = "age", "Meg3")


options(repr.plot.width=12, repr.plot.height=5)

FeaturePlot(obj, split.by = "age", "Pax6")


rownames(obj)[grep("Neat", rownames(obj))] 

VlnPlot(subset(obj, prediction == "OB-STR-CTX Inh IMN") , 
        split.by = "age", 
        group.by = "prediction","Pax6", pt.size = 0)

VlnPlot(subset(obj, celltype_final == "Oligo NN") , 
        split.by = "age", 
        group.by = "prediction","Neat1", pt.size = 0)


VlnPlot(subset(obj, celltype_final == "DG Glut") , split.by = "age", group.by = "prediction","Nipbl", pt.size = 0)

VlnPlot(subset(obj, celltype_final == "DG Glut") , split.by = "age", group.by = "prediction","Meg3", pt.size = 0)

obj$Mal = obj@assays$SCT$scale.data["Mal",]

obj$Il33 = obj@assays$SCT$scale.data["Il33",]

obj$Nrg1 = obj@assays$SCT$scale.data["Nrg1",]

obj$Apoe = obj@assays$SCT$scale.data["Apoe",]

obj$Cdkn2a = obj@assays$SCT$scale.data["Cdkn2a",]

obj$Sox10 = obj@assays$SCT$scale.data["Sox10",]

obj$Sox6 = obj@assays$SCT$scale.data["Sox6",]


obj$Neurod1 = obj@assays$SCT$scale.data["Neurod1",]


obj$Robo1 = obj@assays$SCT$scale.data["Robo1",]


obj$Meg3 = obj@assays$SCT$scale.data["Meg3",]


obj$Nrg1 = obj@assays$SCT$scale.data["Nrg1",]


"Nrg1" %in% rownames(obj)

obj$Pax6 = obj@assays$SCT$scale.data["Pax6",]

# Subset the data for the specified cell type
subset_dg <- obj@meta.data %>%
  filter(celltype_final == "OB-STR-CTX Inh IMN")

gene = "Pax6" 


options(repr.plot.width=4.5, repr.plot.height=4)
subset_dg <- obj@meta.data %>%
  filter(celltype_final == "DG Glut")

gene = "Nrg1" 

# Create a violin plot with boxplot-like statistics and statistical comparisons
g = ggplot(subset_dg, aes(x = age, y = Nrg1, fill = age, 
                           )) +
  # Add violin plot
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  # Add boxplot overlay for summary statistics (median, IQR)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), alpha = 0.6, outlier.shape = NA) +
  # Add mean as a point for extra clarity
  stat_summary(fun = "mean", geom = "point", shape = 21, size = 3, fill = "white", color = "black") +
  # Add statistical comparisons (3 tests between age groups)
  stat_compare_means(comparisons = list(c("2mo", "9mo"), c("2mo", "18mo"), c("9mo", "18mo")), 
                     method = "t.test", label = "p.format") +
  # Customize the theme
  theme_classic() +#ylim(c(0,3))+
  labs(
    x = "Age Group",
    y = paste("Expression of ",gene),
    fill = "Age Group",
    title = paste("Oligo NN ",gene," Expression (MERFISH)")
  ) +
 # scale_color_manual(values = c("Rep1" = "black", "Rep2" = "black","Rep3" = "black")) +

  scale_fill_manual(values = c("2mo" = '#480080', "9mo" = '#e23c5d', "18mo" = '#ffb42c')) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16)
  )#+
  #scale_fill_manual(values = color_vector) 
g

options(repr.plot.width=4.5, repr.plot.height=4)
subset_dg <- obj@meta.data %>%
  filter(celltype_final == "Oligo NN")

gene = "Mal" 

# Create a violin plot with boxplot-like statistics and statistical comparisons
g = ggplot(subset_dg, aes(x = age, y = Mal, fill = age, 
                           )) +
  # Add violin plot
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  # Add boxplot overlay for summary statistics (median, IQR)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), alpha = 0.6, outlier.shape = NA) +
  # Add mean as a point for extra clarity
  stat_summary(fun = "mean", geom = "point", shape = 21, size = 3, fill = "white", color = "black") +
  # Add statistical comparisons (3 tests between age groups)
  stat_compare_means(comparisons = list(c("2mo", "9mo"), c("2mo", "18mo"), c("9mo", "18mo")), 
                     method = "t.test", label = "p.format") +
  # Customize the theme
  theme_classic() +#ylim(c(0,3))+
  labs(
    x = "Age Group",
    y = paste("Expression of ",gene),
    fill = "Age Group",
    title = paste("Oligo NN ",gene," Expression (MERFISH)")
  ) +
 # scale_color_manual(values = c("Rep1" = "black", "Rep2" = "black","Rep3" = "black")) +

  scale_fill_manual(values = c("2mo" = '#480080', "9mo" = '#e23c5d', "18mo" = '#ffb42c')) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16)
  )#+
  #scale_fill_manual(values = color_vector) 
g

pdf(g, file = "../Figure3-oligos//MAL_MERFISH.pdf", height = 4, width = 5)
g
dev.off()

options(repr.plot.width=4.5, repr.plot.height=4)
subset_dg <- obj@meta.data %>%
  filter(celltype_final == "OB-STR-CTX Inh IMN")

gene = "Pax6" 

# Create a violin plot with boxplot-like statistics and statistical comparisons
g = ggplot(subset_dg, aes(x = age, y = Pax6, fill = age, 
                           )) +
  # Add violin plot
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  # Add boxplot overlay for summary statistics (median, IQR)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), alpha = 0.6, outlier.shape = NA) +
  # Add mean as a point for extra clarity
  stat_summary(fun = "mean", geom = "point", shape = 21, size = 3, fill = "white", color = "black") +
  # Add statistical comparisons (3 tests between age groups)
  stat_compare_means(comparisons = list(c("2mo", "9mo"), c("2mo", "18mo"), c("9mo", "18mo")), 
                     method = "t.test", label = "p.format") +
  # Customize the theme
  theme_classic() +#ylim(c(0,3))+
  labs(
    x = "Age Group",
    y = paste("Expression of ",gene),
    fill = "Age Group",
    title = paste("OB-STR-CTX Inh IMN ",gene," Expression (MERFISH)")
  ) +
 # scale_color_manual(values = c("Rep1" = "black", "Rep2" = "black","Rep3" = "black")) +

  scale_fill_manual(values = c("2mo" = '#480080', "9mo" = '#e23c5d', "18mo" = '#ffb42c')) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16)
  )#+
  #scale_fill_manual(values = color_vector) 
g

pdf(g, file = "../Figure2-IMN-IOL/Pax6_MERFISH.pdf", height = 4, width = 5)
g
dev.off()

# Subset the data for the specified cell type
subset_dg <- obj@meta.data %>%
  filter(celltype_final == "DG Glut")

gene = "Meg3" 
options(repr.plot.width=4.5, repr.plot.height=4)
# Create a violin plot with boxplot-like statistics and statistical comparisons
g = ggplot(subset_dg, aes(x = age, y = Meg3, fill = age, color = rep
                           )) +
  # Add violin plot
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  # Add boxplot overlay for summary statistics (median, IQR)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), alpha = 0.6, outlier.shape = NA) +
  # Add mean as a point for extra clarity
  stat_summary(fun = "median", geom = "point", shape = 21, size = 3, fill = "white", color = "black") +
  # Add statistical comparisons (3 tests between age groups)
  stat_compare_means(comparisons = list(c("2mo", "9mo"), c("2mo", "18mo"), c("9mo", "18mo")), 
                     method = "t.test", label = "p.format") +
  # Customize the theme
  theme_classic() +ylim(c(0,3))+
  labs(
    x = "Age Group",
    y = paste("Expression of ",gene),
    fill = "Age Group",
    title = paste("DG ",gene," Expression (MERFISH)")
  ) +
  scale_color_manual(values = c("Rep1" = "black", "Rep2" = "black","Rep3" = "black")) +

  scale_fill_manual(values = c("2mo" = "#E69F00", "9mo" = "#56B4E9", "18mo" = "#009E73")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16)
  )+
  scale_fill_manual(values = color_vector) 
g

# Subset the data for the specified cell type
subset_dg <- obj@meta.data %>%
  filter(celltype_final == "DG Glut")

gene = "Il33" 
options(repr.plot.width=4.5, repr.plot.height=4)
# Create a violin plot with boxplot-like statistics and statistical comparisons
g = ggplot(subset_dg, aes(x = age, y = Robo1, fill = age
                           )) +
  # Add violin plot
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  # Add boxplot overlay for summary statistics (median, IQR)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), alpha = 0.6, outlier.shape = NA) +
  # Add mean as a point for extra clarity
  stat_summary(fun = "median", geom = "point", shape = 21, size = 3, fill = "white", color = "black") +
  # Add statistical comparisons (3 tests between age groups)
  stat_compare_means(comparisons = list(c("2mo", "9mo"), c("2mo", "18mo"), c("9mo", "18mo")), 
                     method = "t.test", label = "p.format") +
  # Customize the theme
  theme_classic() +#ylim(c(0,5))+
  labs(
    x = "Age Group",
    y = paste("Expression of ",gene),
    fill = "Age Group",
    title = paste("DG ",gene," Expression (MERFISH)")
  ) +
  scale_color_manual(values = c("Rep1" = "black", "Rep2" = "black","Rep3" = "black")) +

  scale_fill_manual(values = c("2mo" = "#E69F00", "9mo" = "#56B4E9", "18mo" = "#009E73")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16)
  )+
  scale_fill_manual(values = color_vector) 
g

options(repr.plot.width=4.5, repr.plot.height=4)
# Create a violin plot with boxplot-like statistics and statistical comparisons
g = ggplot(subset_dg, aes(x = age, y = Meg3, fill = age)) +
  # Add violin plot
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  # Add boxplot overlay for summary statistics (median, IQR)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), alpha = 0.6, outlier.shape = NA) +
  # Add mean as a point for extra clarity
  stat_summary(fun = "mean", geom = "point", shape = 21, size = 3, fill = "white", color = "black") +
  # Add statistical comparisons (3 tests between age groups)
  stat_compare_means(comparisons = list(c("2mo", "9mo"), c("2mo", "18mo"), c("9mo", "18mo")), 
                     method = "wilcox.test", label = "p.signif") +
  # Customize the theme
  theme_classic() +
  labs(
    x = "Age Group",
    y = "Expression of Mal",
    fill = "Age Group",
    title = "Oligo Mal Expression (MERFISH)"
  ) +
  scale_fill_manual(values = c("2mo" = "#E69F00", "9mo" = "#56B4E9", "18mo" = "#009E73")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16)
  )+ylim(c(0,10))+
  scale_fill_manual(values = color_vector) 

g

# Subset the data for the specified cell type
subset_data <- obj@meta.data %>%
  filter(celltype_final == "Oligo NN")


options(repr.plot.width=4.5, repr.plot.height=4)
# Create a violin plot with boxplot-like statistics and statistical comparisons
g = ggplot(subset_data, aes(x = age, y = Mal, fill = age, color = rep
                           )) +
  # Add violin plot
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  # Add boxplot overlay for summary statistics (median, IQR)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), alpha = 0.6, outlier.shape = NA) +
  # Add mean as a point for extra clarity
  stat_summary(fun = "median", geom = "point", shape = 21, size = 3, fill = "white", color = "black") +
  # Add statistical comparisons (3 tests between age groups)
  stat_compare_means(comparisons = list(c("2mo", "9mo"), c("2mo", "18mo"), c("9mo", "18mo")), 
                     method = "t.test", label = "p.format") +
  # Customize the theme
  theme_classic() +ylim(c(0,15))+
  labs(
    x = "Age Group",
    y = "Expression of Il33",
    fill = "Age Group",
    title = "Oligo Il33 Expression (MERFISH)"
  ) +
  scale_color_manual(values = c("Rep1" = "black", "Rep2" = "black","Rep3" = "black")) +

  scale_fill_manual(values = c("2mo" = "#E69F00", "9mo" = "#56B4E9", "18mo" = "#009E73")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16)
  )+
  scale_fill_manual(values = color_vector) 
g

pseudobulk_data <- subset_data %>%
  group_by(age, rep) %>%
  summarize(Il33_mean = mean(Il33), .groups = "drop")

t.test(Il33_mean ~ age, data = pseudobulk_data[-which(pseudobulk_data$age == "9mo"),])


# Subset the data for the specified cell type
subset_data <- obj@meta.data %>%
  filter(celltype_final == "Oligo NN")


options(repr.plot.width=4.5, repr.plot.height=4)
# Create a violin plot with boxplot-like statistics and statistical comparisons
g = ggplot(subset_data, aes(x = age, y = Mal, fill = age)) +
  # Add violin plot
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  # Add boxplot overlay for summary statistics (median, IQR)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), alpha = 0.6, outlier.shape = NA) +
  # Add mean as a point for extra clarity
  stat_summary(fun = "mean", geom = "point", shape = 21, size = 3, fill = "white", color = "black") +
  # Add statistical comparisons (3 tests between age groups)
  stat_compare_means(comparisons = list(c("2mo", "9mo"), c("2mo", "18mo"), c("9mo", "18mo")), 
                     method = "wilcox.test", label = "p.signif") +
  # Customize the theme
  theme_classic() +
  labs(
    x = "Age Group",
    y = "Expression of Mal",
    fill = "Age Group",
    title = "Oligo Mal Expression (MERFISH)"
  ) +
  scale_fill_manual(values = c("2mo" = "#E69F00", "9mo" = "#56B4E9", "18mo" = "#009E73")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16)
  )+
  scale_fill_manual(values = color_vector) 

g

svg("Oligo_MAL_MERFISH_Exp.svg", height = 4, width =4)
g
dev.off()

unique(meta$celltype_final)

library(reshape2)
options(repr.plot.width=7, repr.plot.height=5)
ameta = meta
ameta$age_rep = paste(ameta$age, ameta$rep)
tab = table(ameta$celltype_final,ameta$age_rep)
#head(tab)
tab = sweep(tab,2,colSums(tab),'/')
#head(tab)
melted = melt(tab)
#head(melted)
melted$rep = sapply(strsplit(as.character(melted$Var2), " "), "[[", 2)
melted$age = sapply(strsplit(as.character(melted$Var2), " "), "[[", 1)
melted$age = factor(melted$age, levels = c("2mo", "9mo", "18mo"))

cl = "DG Glut"
key_file = "../color_scheme/Age.csv"
key = read.csv(key_file)
head(key)
color_vector <- setNames(key$Color, key$Age)
g2 = ggplot(melted[which(melted$Var1 == cl ),]) + 
  geom_col(aes(factor(Var1),value,fill=age, color =rep ),position="dodge") + 
  xlab("Cell Types") + ylab("Fraction") + 
  scale_fill_manual(values = color_vector) +

  guides(fill=guide_legend(title="Age Group"), color = guide_legend(title="Rep"))  +
  theme_bw() +ggtitle(cl) 

g2

options(repr.plot.width=11, repr.plot.height=9)

obj@meta.data$is= "no"
obj@meta.data$is[which(obj$prediction==cl)]= cl

obj@meta.data$is = factor(obj@meta.data$is, levels = c(cl, "no"))

obj@meta.data %>%
  arrange(desc(is)) %>%
  ggplot(aes(x = global_x, y = global_y, color = is)) +
  geom_point(size=.5) + facet_wrap(~age+rep,nrow = 3, scales="free", drop=F) + theme_bw()+
  scale_color_manual(values = c("tomato", "grey"))

subset_melted = melted[which(melted$Var1 == cl ),]
g2 <- ggplot(subset_melted, aes(x = age, y = value, fill = age)) +
  geom_col(aes(color = rep), position = "dodge", alpha = 0.8) +
  xlab("Age Group") +
  ylab("Fraction") +
  scale_fill_manual(values = color_vector) +
  guides(fill = guide_legend(title = "Age Group"), color = guide_legend(title = "Rep")) +
  theme_classic() +
  ggtitle(cl) +#ylim(c(0,.05))+
  # Add statistical comparisons between ages
  stat_compare_means(comparisons = list(c("2mo", "9mo"), c("2mo", "18mo"), c("9mo", "18mo")),
                     method = "t.test", 
                     label = "p.signif",
                     vjust = 0) # Adjust vertical position of significance labels

g2

svg("../Figure5-DG//DG_cell_prop_MERFISH.svg", height =3.5, width =3.5)
g2
dev.off()

svg("../Figure2-IMN-IOL/DG-PIR_Ex_IMN_cell_prop_MERFISH.svg", height =3.5, width =3.5)
g2
dev.off()

svg("../Figure3-oligos/Oligo_cell_prop_MERFISH.svg", height =3.5, width =3.5)
g2
dev.off()

svg("../Figure3-oligos/Oligo_cell_prop_MERFISH.svg", height =3, width =3.5)
g2
dev.off()

library(dplyr)
library(ggplot2)

# Filter for points where "is" equals the specified class
filtered_data <- obj@meta.data %>%
  filter(is == "Oligo NN")

# Plot density of positive "is" values
filtered_data %>%
  ggplot(aes(x = global_x, y = global_y)) +
  # Add density layer
  stat_density_2d(aes(fill = after_stat(nlevel), alpha = after_stat(nlevel)), 
                  geom = "polygon", contour = TRUE) +
  # Facet by age and rep
  facet_wrap(~age + rep, nrow = 3, scales = "free", drop = FALSE) +
  theme_bw() +
  scale_fill_viridis_c(option = "plasma") +  # Use viridis for density visualization
  scale_alpha(range = c(0.4, 0.9), guide = "none") +  # Adjust transparency
  xlab("Global X") + ylab("Global Y") +
  labs(fill = "Density") +
  ggtitle(paste("Density of", cl, "Points")) +
  theme(legend.position = "right")

library(dplyr)
library(ggplot2)
library(ggrepel)

# Define the rotation function
rotate_coordinates <- function(data, angle) {
  radians <- angle * pi / 180
  data %>%
    mutate(
      global_x_rotated = cos(radians) * global_x - sin(radians) * global_y,
      global_y_rotated = sin(radians) * global_x + cos(radians) * global_y
    )
}

# Apply specific rotations and flip
apply_rotations <- function(data) {
  data %>%
    mutate(
      global_x = case_when(
        age == "9mo" & rep == "Rep2" ~ -rotate_coordinates(., 120)$global_x_rotated,
        age == "18mo" & rep == "Rep1" ~ -rotate_coordinates(., 0)$global_x_rotated,
        TRUE ~ global_x
      ),
      global_y = case_when(
        age == "9mo" & rep == "Rep2" ~ rotate_coordinates(., 120)$global_y_rotated, # Flip along y-axis
        age == "18mo" & rep == "Rep1" ~ rotate_coordinates(., 0)$global_y_rotated,
        TRUE ~ global_y
      )
    )
}

# Step 1: Filter the data for specific replicates
filtered_meta <- obj@meta.data %>%
  filter((age == "2mo" & rep == "Rep2") |
         (age == "9mo" & rep == "Rep2") |
         (age == "18mo" & rep == "Rep1"))

# Step 2: Apply rotations and flips
filtered_meta <- apply_rotations(filtered_meta)

# Step 3: Subsample the data (up to 3000 cells per cell type)
sampled_meta <- filtered_meta %>%
  group_by(celltype_final) %>%
  sample_n(size = min(5000, n())) %>%
  ungroup()

# Step 4: Update `is` column based on specific cell types
sampled_meta$is <- "other"
sampled_meta$is[which(sampled_meta$celltype_final %in% c("IOL NN", "DG-PIR Ex IMN", "OB-STR-CTX Inh IMN"))] <- 
  sampled_meta$celltype_final[which(sampled_meta$celltype_final %in% c("IOL NN", "DG-PIR Ex IMN", "OB-STR-CTX Inh IMN"))]

sampled_meta$is <- factor(sampled_meta$is, levels = c("DG-PIR Ex IMN","OB-STR-CTX Inh IMN", "IOL NN",   "other"))
color_vector["other"] = "lightgray"


# Step 5: Calculate centroids for the three cell types to be labeled
centroids <- sampled_meta %>%
  filter(is != "NA") %>%
  group_by(is) %>%
  summarise(global_x = median(global_x), global_y = median(global_y))

# Step 6: Plot the final data with labels
options(repr.plot.width = 12, repr.plot.height = 5)

g1 = sampled_meta %>%
  arrange(desc(is)) %>%
  ggplot(aes(x = global_x, y = global_y, color = is)) +
  geom_point(size = .5) +
  facet_wrap(~age , nrow = 1, scales = "free", drop = FALSE) +
  theme_classic() +
  scale_color_manual(values = color_vector) +
  theme(legend.position = "top", text = element_text(size = 14), # Move the legend to the bottom
        legend.title = element_text(size = 0),  # Adjust legend title size
        legend.text = element_text(size = 14)) +  # Adjust legend text size
  guides(color = guide_legend(override.aes = list(size = 3)))  # Increase dot size in the legend


# g2 = sampled_meta %>%
#   arrange(Pax6) %>%
#   ggplot(aes(x = global_x, y = global_y, color = Pax6)) +
#   geom_point(size = .5) +
#   facet_wrap(~age , nrow = 1, scales = "free", drop = FALSE) +
#   theme_classic() + ggtitle("Pax6")+
#   scale_color_gradient(low = "gray", high = "red")+
#   theme(legend.position = "top",  # Move the legend to the bottom
#         legend.title = element_text(size = 0),  # Adjust legend title size
#         legend.text = element_text(size = 8)) +  # Adjust legend text size
#   guides(color = guide_legend(override.aes = list(size = 3)))  # Increase dot size in the legend



# g3 = sampled_meta %>%
#   arrange(Sirt2) %>%
#   ggplot(aes(x = global_x, y = global_y, color = Sirt2)) +
#   geom_point(size = .5) +
#   facet_wrap(~age , nrow = 1, scales = "free", drop = FALSE) +
#   theme_classic() + ggtitle("Sirt2")+
#   scale_color_gradient(low = "gray", high = "red")+
#   theme(legend.position = "top",  # Move the legend to the bottom
#         legend.title = element_text(size = 0),  # Adjust legend title size
#         legend.text = element_text(size = 8)) +  # Adjust legend text size
#   guides(color = guide_legend(override.aes = list(size = 3)))  # Increase dot size in the legend



# sampled_meta %>%
#   arrange(Eif5a) %>%
#   ggplot(aes(x = global_x, y = global_y, color = Eif5a)) +
#   geom_point(size = .5) +
#   facet_wrap(~age , nrow = 1, scales = "free", drop = FALSE) +
#   theme_classic() +
#   scale_color_gradient(low = "gray", high = "red")+
#   theme(legend.position = "top",  # Move the legend to the bottom
#         legend.title = element_text(size = 0),  # Adjust legend title size
#         legend.text = element_text(size = 8)) +  # Adjust legend text size
#   guides(color = guide_legend(override.aes = list(size = 3)))  # Increase dot size in the legend




pdf("../Figure2-IMN-IOL/prog_fov.pdf", height = 4, width = 10)
g1
dev.off()

g1


# Plot the rotated and ordered data
options(repr.plot.width=20, repr.plot.height=10)
celltype_order <- meta %>%
  count(celltype_final) %>%
  arrange(n) %>%
  pull(celltype_final)

# Apply the rotation and cell type ordering
meta <- meta %>%
  mutate(celltype_final = factor(celltype_final, levels = celltype_order))


meta[which(meta$batch == "2moRep2"),] %>%
  ggplot(aes(x = global_x, y = global_y, color = celltype_final)) +
  scale_color_manual(values = color_vector) +
  geom_point(size = .5) 

# Plot the rotated and ordered data
options(repr.plot.width=19, repr.plot.height=26)
celltype_order <- meta %>%
  count(celltype_final) %>%
  arrange(n) %>%
  pull(celltype_final)

# Apply the rotation and cell type ordering
meta <- meta %>%
  mutate(celltype_final = factor(celltype_final, levels = celltype_order))

 p = meta[which(meta$batch == "2moRep2"),] %>%
  ggplot(aes(x = global_x, y = global_y, color = celltype_final)) +
  scale_color_manual(values = color_vector) +
  geom_point(size = .5) +
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  theme_classic()+
  theme(legend.position = "right", 
        legend.box = "vertical", 
        legend.title = element_text(size = 0),  # Adjust legend title size if needed
        legend.text = element_text(size = 10))   # Adjust legend text size if needed


table(obj$celltype_final)
2434/ncol(obj)

table(obj$prediction)


