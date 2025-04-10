library(Seurat)
library(ggplot2)
library(dplyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)  # Load the purrr package for map functions
library(broom)  # broom is useful for tidying up statistical test results


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

options(repr.plot.width=12, repr.plot.height=12)

DimPlot(obj, group.by = "prediction", label =T,pt.size = 0)+NoLegend()


options(repr.plot.width=12, repr.plot.height=5)

FeaturePlot(obj, split.by = "age", "Eif5a")


options(repr.plot.width=12, repr.plot.height=5)

FeaturePlot(obj, split.by = "age", "Il33")


FeaturePlot(obj, split.by = "age", "Neat1")


table(obj@meta.data$batch)

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

write.table(meta, "MERFISH_meta.csv", sep =",") 

nrow(meta[which(meta$batch!="2moRep1"),])
length(unique(meta[which(meta$batch!="2moRep1"), "celltype_final"]))

#key[which( key$CellType %in%meta$celltype_final ),]
setwd("~/projects/combined_all/Figures/Figure2-IMN-IOL/")
key_file = "../color_scheme/updated_celltype_palette.csv"
key = read.csv(key_file)
head(key)
color_vector <- setNames(key$Color, key$CellType)
unique(meta$celltype_final)[which(! unique(meta$celltype_final) %in% key$CellType  )]

unique(meta$celltype_final)[which(! unique(meta$celltype_final) %in% key$CellType  )]

options(repr.plot.width=6, repr.plot.height=6)
sample = meta[sample(1:nrow(meta),100000),]

g=ggplot(sample, aes(x = umap_x, y=umap_y, color = celltype_final)) +
  geom_point(alpha=.8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() +
  theme(legend.position = "none") 
g
#pdf("celltype_umap.pdf", width=6, height =6)
#print(g)
#dev.off()

library(ggrepel)
# Step 1: Calculate the median UMAP coordinates and filter for cell types with more than 3000 cells
centroids <- meta %>%
  group_by(celltype_final) %>%
  filter(n() > 100) %>%  # Keep only cell types with more than 3000 cells
  summarise(umap_x = median(umap_x), umap_y = median(umap_y))

# Step 2: Sample up to 3000 cells per cell type
sample <- meta %>%
  group_by(celltype_final) %>%
  sample_n(size = min(5000, n())) %>%
  ungroup()

# Step 3: Create the plot and add labels using geom_text_repel
g <- ggplot(sample, aes(x = umap_x, y = umap_y, color = celltype_final)) +
  geom_point(alpha = .8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_text_repel(data = centroids, aes(label = celltype_final, color = "black"), 
                  size = 4, 
                  box.padding = 0.5, 
                  point.padding = 0.3, 
                  segment.color = 'grey50', 
                  max.overlaps = Inf)  # Use ggrepel for better label positioning

g

options(repr.plot.width=7, repr.plot.height=7)

# Step 1: Calculate the average UMAP coordinates and filter for cell types with more than 3000 cells
centroids <- meta %>%
  group_by(celltype_final) %>%
  filter(n() > 500) %>%  # Keep only cell types with more than 3000 cells
  summarise(umap_x = median(umap_x), umap_y = median(umap_y))

# Step 2: Sample up to 3000 cells per cell type
sample <- meta %>%
  group_by(celltype_final) %>%
  sample_n(size = min(10000, n())) %>%
  ungroup()

# Step 3: Create the plot and add labels
g <- ggplot(sample, aes(x = umap_x, y = umap_y, color = celltype_final)) +
  geom_point(alpha = .8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_text(data = centroids, aes(label = celltype_final), color = "black", size = 3, check_overlap = TRUE)  # Add labels at average positions

g


library(scales)             # Load the package
options(repr.plot.width=15, repr.plot.height=8)

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
  scale_y_continuous(breaks = c(0, 50000, 100000), 
                     labels = label_number(scale_cut = cut_si(""), accuracy = 1),
                     limits = c(0, 100000),  # Adjusts the upper limit to 110% of the max value
                     expand = expansion(mult = c(0, 0.05)))   # Adds some space at the top
g

key_file = "../color_scheme/Age.csv"
key = read.csv(key_file)
head(key)
color_vector <- setNames(key$Color, key$Age)

options(repr.plot.width=6, repr.plot.height=5)
#color_vector <- setNames(age_color$Color, age_color$Age)
meta$age = factor(meta$age, levels = c("2mo", "9mo", "18mo"))
sample = meta[sample(1:nrow(meta),80000),]

g = ggplot(sample, aes(x = umap_x, y=umap_y, color = age)) +
  geom_point(alpha=.8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() 
g
pdf("age_UMAP.pdf",width=6, height=6)
print(g)
dev.off()

meta_sub = meta[which(meta$celltype_final %in% c("DG-PIR Ex IMN","OB-STR-CTX Inh IMN",
                                                 "IOL NN", "IOL","OPC NN", "Oligo NN", 
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

meta$sample = meta$batch

table(meta$celltype_final)

over1k = names(table(meta$celltype_final)[which(table(meta$celltype_final)>3000)])
over1k = c(over1k, 'DG-PIR Ex IMN')
over1k = c(over1k, 'OB-STR-CTX Inh IMN')


options(repr.plot.width=7, repr.plot.height=10)
meta_sub = meta[which(meta$celltype_final%in%over1k),]
meta_sub = meta_sub[-which(meta_sub$batch == "2moRep1"),]

ggplot(meta_sub, aes(x = celltype_final, fill = age, color = rep)) +
  geom_bar(stat = "count", position = 'fill') +
  theme_minimal() +
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip()+
  scale_fill_manual(values = color_vector) 

over1k

color_vector

options(repr.plot.width=5, repr.plot.height=7)

library(dplyr)
library(ggplot2)

# Filter the data
meta_sub <- meta %>%
  filter(celltype_final %in% over1k & batch != "2moRep1")
meta_sub$rep[which(meta_sub$rep=="Rep3")] = "Rep1"
# Step 1: Calculate the proportion of cells that are 2mo within each cell type
proportion_2mo <- meta_sub %>%
  group_by(celltype_final) %>%
  summarise(proportion_2mo = sum(age == "2mo") / n()) %>%
  arrange(desc(proportion_2mo))

# Step 2: Reorder the celltype_final factor based on the proportion of 2mo cells
meta_sub$celltype_final <- factor(meta_sub$celltype_final, levels = rev(proportion_2mo$celltype_final))

# Step 3: Plot the data
options(repr.plot.width = 5.5, repr.plot.height = 6)
g = ggplot(meta_sub, aes(x = celltype_final, fill = age, color = rep)) +
  geom_bar(stat = "count", position = 'fill') +
  theme_minimal() +
  labs(title = NULL,
       x = "Cell Type",
       y = "Proportion of Cells") +
  coord_flip() +
  scale_fill_manual(values = color_vector)+ theme_bw()
svg("celltype_proportions.svg")
print(g)
dev.off()

options(repr.plot.width=12, repr.plot.height=5)

FeaturePlot(obj, "Gpc5", split.by="age",reduction = "umap")
FeaturePlot(obj, "Crlf2", split.by="age",reduction = "umap")


options(repr.plot.width=12, repr.plot.height=5)

VlnPlot(obj, "Crlf2", group.by="age")


options(repr.plot.width=12, repr.plot.height=12)

DimPlot(obj, group.by = "predicted.id", label =T,pt.size = 0)+NoLegend()
DimPlot(obj, group.by = "predicted_id_ext", label =T,pt.size = 0)+NoLegend()
DimPlot(obj, group.by = "celltype_final", label =T,pt.size = 0)+NoLegend()
DimPlot(obj, group.by = "age", label =T,pt.size = 0)+NoLegend()


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

FeaturePlot(obj, "Crlf2", pt.size = 0)


proj = obj[which(obj$celltype_final %in% c("IOL NN", "DG-PIR Ex IMN","OB-STR-CTX Inh IMN" )),]

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

#cl = "DG Glut"

g2 = ggplot(melted[which(melted$Var1 %in%  c("IOL NN", "DG-PIR Ex IMN","OB-STR-CTX Inh IMN" ) ),]) + 
  geom_col(aes(factor(Var1),value,fill=age, color =rep ),position="dodge") + 
  xlab("Cell Types") + ylab("Fraction") + 
  scale_fill_manual(values = color_vector) +
  guides(fill=guide_legend(title="Age Group"), color = guide_legend(title="Rep"))  +
  theme_bw() 

g2

options(repr.plot.width=11, repr.plot.height=9)

obj@meta.data$is= "NA"
obj@meta.data$is[which(obj$celltype_final %in% c("IOL NN", "DG-PIR Ex IMN","OB-STR-CTX Inh IMN" ))]= obj$celltype_final[which(obj$celltype_final %in% c("IOL NN", "DG-PIR Ex IMN","OB-STR-CTX Inh IMN" ))]

obj@meta.data$is = factor(obj@meta.data$is, levels = c("IOL NN", "OB-STR-CTX Inh IMN", "DG-PIR Ex IMN", "NA"))

obj@meta.data %>%
  arrange(desc(is)) %>%
  ggplot(aes(x = global_x, y = global_y, color = is)) +
  geom_point(size=.5) + facet_wrap(~age+rep,nrow = 3, scales="free", drop=F) + theme_bw()+
  scale_color_manual(values = color_vector) 


options(repr.plot.width=8, repr.plot.height=6)
m = melted[which(melted$Var1 %in%  c("IOL NN", "DG-PIR Ex IMN", "OB-STR-CTX Inh IMN") ),]
m$Var1 = factor(m$Var1, levels = c("IOL NN", "DG-PIR Ex IMN", "OB-STR-CTX Inh IMN"))
g = ggplot(m) + 
  geom_col(aes(factor(Var1), value, fill = age), color = "black", position = "dodge") + 
  xlab("Cell Types") + 
  ylab("Fraction") + 
  scale_fill_manual(values = color_vector) +
  guides(fill = guide_legend(title = "Age Group"), color = "white") + 
  theme_bw() +
  theme(
    text = element_text(size = 14),        # Base font size
    axis.title = element_text(size = 16), # Axis titles
    axis.text = element_text(size = 14),  # Axis tick labels
    legend.title = element_text(size = 14), # Legend title
    legend.text = element_text(size = 12)  # Legend labels
  )

svg("MERFISH_prog_barplot.svg", height =6, width =8)
print(g)
dev.off()
dev.off()

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
        age == "18mo" & rep == "Rep2" ~ rotate_coordinates(., 180)$global_x_rotated,
        TRUE ~ global_x
      ),
      global_y = case_when(
        age == "9mo" & rep == "Rep2" ~ rotate_coordinates(., 120)$global_y_rotated, # Flip along y-axis
        age == "18mo" & rep == "Rep2" ~ rotate_coordinates(., 180)$global_y_rotated,
        TRUE ~ global_y
      )
    )
}

# Step 1: Filter the data for specific replicates
filtered_meta <- meta %>%
  filter((age == "2mo" & rep == "Rep2") |
         (age == "9mo" & rep == "Rep2") |
         (age == "18mo" & rep == "Rep2"))

# Step 2: Apply rotations and flips
filtered_meta <- apply_rotations(filtered_meta)

# Step 3: Subsample the data (up to 3000 cells per cell type)
sampled_meta <- filtered_meta %>%
  group_by(celltype_final) %>%
  sample_n(size = min(100000, n())) %>%
  ungroup()

# Step 4: Update `is` column based on specific cell types
sampled_meta$is <- "other"
sampled_meta$is[which(sampled_meta$celltype_final %in% c("IOL NN", "DG-PIR Ex IMN", "OB-STR-CTX Inh IMN"))] <- 
  sampled_meta$celltype_final[which(sampled_meta$celltype_final %in% c("IOL NN", "DG-PIR Ex IMN", "OB-STR-CTX Inh IMN"))]

sampled_meta$is <- factor(sampled_meta$is, levels = c("IOL NN", "OB-STR-CTX Inh IMN", "DG-PIR Ex IMN", "other"))
color_vector["other"] = "lightgray"
# Step 5: Calculate centroids for the three cell types to be labeled
centroids <- sampled_meta %>%
  filter(is != "NA") %>%
  group_by(is) %>%
  summarise(global_x = median(global_x), global_y = median(global_y))

# Step 6: Plot the final data with labels
options(repr.plot.width = 12, repr.plot.height = 5)

g <- sampled_meta %>%
  arrange(desc(is)) %>%
  ggplot(aes(x = global_x, y = global_y, color = is)) +
  geom_point(size = .5) +
  facet_wrap(~age , nrow = 1, scales = "free", drop = FALSE) +
  theme_classic() +
  scale_color_manual(values = color_vector) +
  theme(legend.position = "top",  # Move the legend to the bottom
        legend.title = element_text(size = 0),  # Adjust legend title size
        legend.text = element_text(size = 8)) +  # Adjust legend text size
  guides(color = guide_legend(override.aes = list(size = 3)))  # Increase dot size in the legend

#pdf("fov_progenitor_plot.pdf", width = 10, height =4)
print(g)
#dev.off()
#print(g)

#svg("fov_progenitor_plot.svg", width = 10, height =4)
print(g)
#dev.off()
#print(g)


proge = sampled_meta[which(sampled_meta$is != "other"),] 
samp =  sampled_meta[which(sampled_meta$is == "other"),] 
cat(nrow(samp))
samp = samp[sample(nrow(samp),75000),]
samp2 = rbind(proge,samp)
g <- samp2 %>%
  arrange(desc(is)) %>%
  ggplot(aes(x = global_x, y = global_y, color = is)) +
  geom_point(size = .5) +
  facet_wrap(~age , nrow = 1, scales = "free", drop = FALSE) +
  theme_classic() +
  scale_color_manual(values = color_vector) +
  theme(legend.position = "top",  # Move the legend to the bottom
        legend.title = element_text(size = 0),  # Adjust legend title size
        legend.text = element_text(size = 8)) +  # Adjust legend text size
  guides(color = guide_legend(override.aes = list(size = 3)))  # Increase dot size in the legend
svg("fov_progenitor_plot.svg", width = 10, height =4)
print(g)
dev.off()
print(g)

library(dplyr)
library(ggplot2)

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
        age == "18mo" & rep == "Rep2" ~ rotate_coordinates(., 180)$global_x_rotated,
        TRUE ~ global_x
      ),
      global_y = case_when(
        age == "9mo" & rep == "Rep2" ~ rotate_coordinates(., 120)$global_y_rotated, # Flip along y-axis
        age == "18mo" & rep == "Rep2" ~ rotate_coordinates(., 180)$global_y_rotated,
        TRUE ~ global_y
      )
    )
}

# Step 1: Filter the data for specific replicates
filtered_meta <- meta %>%
  filter((age == "2mo" & rep == "Rep2") |
         (age == "9mo" & rep == "Rep2") |
         (age == "18mo" & rep == "Rep2"))

# Step 2: Apply rotations and flips
filtered_meta <- apply_rotations(filtered_meta)

# Step 3: Subsample the data (up to 3000 cells per cell type)
sampled_meta <- filtered_meta %>%
  group_by(celltype_final) %>%
  sample_n(size = min(3000, n())) %>%
  ungroup()

# Step 4: Plot the final data
options(repr.plot.width = 11, repr.plot.height = 9)

sampled_meta %>%
  ggplot(aes(x = global_x, y = global_y, color = celltype_final)) +
  geom_point(size = .5) +
  facet_wrap(~age , nrow = 3, scales = "free", drop = FALSE) +
  theme_bw() +
  scale_color_manual(values = color_vector)




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

g2 = ggplot(melted[which(melted$Var1 == cl ),]) + 
  geom_col(aes(factor(Var1),value,fill=age, color =rep ),position="dodge") + 
  xlab("Cell Types") + ylab("Fraction") + 
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

head(meta)

obj$Eif5a = obj@assays$SCT$scale.data["Eif5a",]
obj$Neat1 = obj@assays$SCT$scale.data["Neat1",]

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
        age == "18mo" & rep == "Rep2" ~ rotate_coordinates(., 180)$global_x_rotated,
        TRUE ~ global_x
      ),
      global_y = case_when(
        age == "9mo" & rep == "Rep2" ~ rotate_coordinates(., 120)$global_y_rotated, # Flip along y-axis
        age == "18mo" & rep == "Rep2" ~ rotate_coordinates(., 180)$global_y_rotated,
        TRUE ~ global_y
      )
    )
}

# Step 1: Filter the data for specific replicates
filtered_meta <- obj@meta.data %>%
  filter((age == "2mo" & rep == "Rep2") |
         (age == "9mo" & rep == "Rep2") |
         (age == "18mo" & rep == "Rep2"))

# Step 2: Apply rotations and flips
filtered_meta <- apply_rotations(filtered_meta)

# Step 3: Subsample the data (up to 3000 cells per cell type)
sampled_meta <- filtered_meta %>%
  group_by(celltype_final) %>%
  sample_n(size = min(100000, n())) %>%
  ungroup()

# Step 5: Calculate centroids for the three cell types to be labeled
#centroids <- sampled_meta %>%
#  filter(is != "NA") %>%
#  group_by(is) %>%
#  summarise(global_x = median(global_x), global_y = median(global_y))

# Step 6: Plot the final data with labels
options(repr.plot.width = 12, repr.plot.height = 5)

sampled_meta %>%
  arrange(Neat1) %>%
  ggplot(aes(x = global_x, y = global_y, color = Neat1)) +
  geom_point(size = .5) +
  facet_wrap(~age , nrow = 1, scales = "free", drop = FALSE) +
  theme_classic() +
  scale_color_gradient(low = "gray", high = "red")+
  theme(legend.position = "top",  # Move the legend to the bottom
        legend.title = element_text(size = 0),  # Adjust legend title size
        legend.text = element_text(size = 8)) +  # Adjust legend text size
  guides(color = guide_legend(override.aes = list(size = 3)))  # Increase dot size in the legend



head(sampled_meta)


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
  sample_n(size = min(100000, n())) %>%
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
  theme(legend.position = "top",  # Move the legend to the bottom
        legend.title = element_text(size = 0),  # Adjust legend title size
        legend.text = element_text(size = 8)) +  # Adjust legend text size
  guides(color = guide_legend(override.aes = list(size = 3)))  # Increase dot size in the legend


g2 = sampled_meta %>%
  arrange(Il33) %>%
  ggplot(aes(x = global_x, y = global_y, color = Il33)) +
  geom_point(size = .5) +
  facet_wrap(~age , nrow = 1, scales = "free", drop = FALSE) +
  theme_classic() + ggtitle("Il33")+
  scale_color_gradient(low = "gray", high = "red")+
  theme(legend.position = "top",  # Move the legend to the bottom
        legend.title = element_text(size = 0),  # Adjust legend title size
        legend.text = element_text(size = 8)) +  # Adjust legend text size
  guides(color = guide_legend(override.aes = list(size = 3)))  # Increase dot size in the legend



g3 = sampled_meta %>%
  arrange(Sirt2) %>%
  ggplot(aes(x = global_x, y = global_y, color = Sirt2)) +
  geom_point(size = .5) +
  facet_wrap(~age , nrow = 1, scales = "free", drop = FALSE) +
  theme_classic() + ggtitle("Sirt2")+
  scale_color_gradient(low = "gray", high = "red")+
  theme(legend.position = "top",  # Move the legend to the bottom
        legend.title = element_text(size = 0),  # Adjust legend title size
        legend.text = element_text(size = 8)) +  # Adjust legend text size
  guides(color = guide_legend(override.aes = list(size = 3)))  # Increase dot size in the legend



sampled_meta %>%
  arrange(Eif5a) %>%
  ggplot(aes(x = global_x, y = global_y, color = Eif5a)) +
  geom_point(size = .5) +
  facet_wrap(~age , nrow = 1, scales = "free", drop = FALSE) +
  theme_classic() +
  scale_color_gradient(low = "gray", high = "red")+
  theme(legend.position = "top",  # Move the legend to the bottom
        legend.title = element_text(size = 0),  # Adjust legend title size
        legend.text = element_text(size = 8)) +  # Adjust legend text size
  guides(color = guide_legend(override.aes = list(size = 3)))  # Increase dot size in the legend





sampled_meta %>%
  arrange(Eif5a) %>%
  ggplot(aes(x = global_x, y = global_y, color = Eif5a)) +
  geom_point(size = .5) +
  facet_wrap(~age , nrow = 1, scales = "free", drop = FALSE) +
  theme_classic() +
  scale_color_gradient(low = "gray", high = "red")+
  theme(legend.position = "top",  # Move the legend to the bottom
        legend.title = element_text(size = 0),  # Adjust legend title size
        legend.text = element_text(size = 8)) +  # Adjust legend text size
  guides(color = guide_legend(override.aes = list(size = 3)))  # Increase dot size in the legend



dev.off()


options(repr.plot.width=12, repr.plot.height=4)
pdf("~/projects/combined_all/Figures/Figure2/fov_prog_sirt.pdf", height = 4, width = 12)
print(g1+ggtitle("progenitor celltypes label")+
  theme(legend.position = "none",  # Move the legend to the bottom
        legend.title = element_text(size = 0),  # Adjust legend title size
        legend.text = element_text(size = 8)))  # Adjust legend text size)
print(g3+
  theme(legend.position = "none",  # Move the legend to the bottom
        legend.title = element_text(size = 0),  # Adjust legend title size
        legend.text = element_text(size = 8)))
print(g2+
  theme(legend.position = "none",  # Move the legend to the bottom
        legend.title = element_text(size = 0),  # Adjust legend title size
        legend.text = element_text(size = 8)))


dev.off()


options(repr.plot.width=15, repr.plot.height=15)
obj$Mal = obj@assays$SCT$scale.data["C4b",]
obj@meta.data %>%
  arrange(Mal) %>%  # Arrange the data so that lower values are plotted first
  ggplot(aes(x = global_x, y = global_y, color = Sirt2)) +
  geom_point(size=.5) + 
  facet_wrap(~age+rep, nrow = 3, scales="free", drop=F) + 
  theme_bw() +
  scale_color_gradient(low = "gray", high = "red")


options(repr.plot.width=15, repr.plot.height=15)
obj$Mal = obj@assays$SCT$scale.data["Mal",]
obj@meta.data %>%
  arrange(Mal) %>%  # Arrange the data so that lower values are plotted first
  ggplot(aes(x = global_x, y = global_y, color = Sirt2)) +
  geom_point(size=.5) + 
  facet_wrap(~age+rep, nrow = 3, scales="free", drop=F) + 
  theme_bw() +
  scale_color_gradient(low = "gray", high = "red")




options(reprobj@assays$SCT$scale.data$height=10)


sampled_meta <- meta %>%
  group_by(celltype_final) %>%
  sample_n(size = min(20000, n())) %>%
  ungroup()
sampled_meta = sampled_meta %>%
  mutate(celltype_final = factor(celltype_final, levels = celltype_order))

sampled_meta %>%
  ggplot(aes(x = global_x, y = global_y, color = celltype_final)) +
  scale_color_manual(values = color_vector) +
  geom_point(size=.5) + facet_wrap(~age+rep,nrow = 3, scales="free", drop=F) + theme_bw()


library(dplyr)
library(ggplot2)

# Filter for the specific sample `2mo_rep1`
sample_data <- sampled_meta %>%
  filter(age == "2mo" & rep == "Rep1")

# Apply the 90-degree counterclockwise rotation using the rotation matrix
rotated_sample_data <- sample_data %>%
  mutate(
    global_x_rotated = -global_y,
    global_y_rotated = global_x
  )

# Calculate the frequency of each cell type and order them
celltype_order <- rotated_sample_data %>%
  count(celltype_final) %>%
  arrange(n) %>% # Order by count (smallest to largest)
  pull(celltype_final) # Extract the ordered cell type list


# Replace the original coordinates with the rotated ones in the original dataset
sampled_meta_rotated <- sampled_meta %>%
  mutate(
    global_x = ifelse(age == "2mo" & rep == "Rep1", rotated_sample_data$global_x_rotated, global_x),
    global_y = ifelse(age == "2mo" & rep == "Rep1", rotated_sample_data$global_y_rotated, global_y),
    celltype_final = ifelse(age == "2mo" & rep == "Rep1", as.character(rotated_sample_data$celltype_final), celltype_final)
  )

# Convert the `celltype_final` column into a factor and order by celltype count
sampled_meta_rotated <- sampled_meta_rotated %>%
  mutate(celltype_final = factor(celltype_final, levels = celltype_order))


# Plot the rotated data with rarer cell types on top
sampled_meta_rotated %>%
  ggplot(aes(x = global_x, y = global_y, color = celltype_final)) +
  scale_color_manual(values = color_vector) +
  geom_point(size=.5) + 
  facet_wrap(~age+rep, nrow = 3, scales = "free", drop = FALSE) + 
  theme_bw()


library(dplyr)
library(ggplot2)

# Define the rotation function
rotate_coordinates <- function(data, angle) {
  radians <- angle * pi / 180
  data %>%
    mutate(
      global_x_rotated = cos(radians) * global_x - sin(radians) * global_y,
      global_y_rotated = sin(radians) * global_x + cos(radians) * global_y
    )
}

# Define a function to apply specific rotations and flip
apply_rotations <- function(data) {
  data %>%
    mutate(
      global_x = case_when(
        age == "2mo" & rep == "Rep1" ~ rotate_coordinates(., 180)$global_x_rotated,
        age == "2mo" & rep == "Rep3" ~ rotate_coordinates(., 45)$global_x_rotated,
        age == "9mo" & rep == "Rep1" ~ rotate_coordinates(., 195)$global_x_rotated,
        age == "9mo" & rep == "Rep2" ~ rotate_coordinates(., 275)$global_x_rotated,
        age == "18mo" & rep == "Rep1" ~ -global_x, # Flip across y-axis
        age == "18mo" & rep == "Rep2" ~ rotate_coordinates(., 180)$global_x_rotated,
        TRUE ~ global_x
      ),
      global_y = case_when(
        age == "2mo" & rep == "Rep1" ~ rotate_coordinates(., 180)$global_y_rotated,
        age == "2mo" & rep == "Rep3" ~ rotate_coordinates(., 45)$global_y_rotated,
        age == "9mo" & rep == "Rep1" ~ rotate_coordinates(., 195)$global_y_rotated,
        age == "9mo" & rep == "Rep2" ~ rotate_coordinates(., 275)$global_y_rotated,
        age == "18mo" & rep == "Rep1" ~ global_y, # Flip across y-axis
        age == "18mo" & rep == "Rep2" ~ rotate_coordinates(., 180)$global_y_rotated,
        TRUE ~ global_y
      )
    )
}

sampled_meta <- meta %>%
  group_by(celltype_final) %>%
  sample_n(size = min(25000, n())) %>%
  ungroup()

# Calculate the consistent cell type order across all samples
celltype_order <- sampled_meta %>%
  count(celltype_final) %>%
  arrange(n) %>%
  pull(celltype_final)

# Apply the rotation and cell type ordering
sampled_meta_rotated <- sampled_meta %>%
  mutate(celltype_final = factor(celltype_final, levels = celltype_order)) %>%
  apply_rotations()

# Plot the rotated and ordered data
sampled_meta_rotated %>%
  ggplot(aes(x = global_x, y = global_y, color = celltype_final)) +
  scale_color_manual(values = color_vector) +
  geom_point(size = .5) + 
  facet_wrap(~age + rep, nrow = 3, scales = "free", drop = FALSE) + 
  theme_bw()


head(sampled_meta_rotated$batch)

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


pdf("celltype_fov.pdf", width=9, height =11)
print(p)
dev.off()

av_rna = load("~/projects//combined_all/female_RNA/celltype_age_average_expression.RData")

av_rna = av$SCT

obj$celltype_age = paste(obj$celltype_final, obj$age, sep = "_")
Idents(obj) = "celltype_age"
av = AverageExpression(obj)

av_og = av$originalexp


dim(av_og)

library(pheatmap)
pheatmap(av_og, scale = "row")

unique(obj$celltype_age)

write.csv(fm , file = paste("~/projects/merfish/",ct, "DEGs.csv",sep = ""))

which(unique(obj$celltype_final) == ct)


Idents(obj) = 'celltype_age'
for (ct in unique(obj$celltype_final)[60:100] ) {
       fm = FindMarkers(obj, `ident.1` = paste(ct,"2mo", sep = "_"), `ident.2` = paste(ct,"18mo", sep = "_"))
       head(fm)
        ct = gsub(" ", "-", ct)
        ct = gsub("/", "--", ct)
        write.csv(fm , file = paste("~/projects/merfish/",ct, "_DEGs.csv",sep = ""))
}

obj$yng_old = "yng"
obj$yng_old[which(obj$age != "2mo")] = "old"

obj$celltype_yng_old = paste(obj$celltype_final,obj$yng_old, sep ="_")
unique(obj$celltype_yng_old)

prog_cts = c("IOL NN", 'DG-PIR Ex IMN','OB-STR-CTX Inh IMN' )

Idents(obj) = 'celltype_yng_old'
for (ct in prog_cts) {
       fm = FindMarkers(obj, `ident.1` = paste(ct,"yng", sep = "_"), `ident.2` = paste(ct,"old", sep = "_"))
       head(fm)
        ct = gsub(" ", "-", ct)
        ct = gsub("/", "--", ct)
        write.csv(fm , file = paste("~/projects/merfish/",ct, "_DEGs_yng_old.csv",sep = ""))
}

ct = "IOL NN"
fm = FindMarkers(obj, `ident.1` = paste(ct,"2mo", sep = "_"), `ident.2` = paste(ct,"18mo", sep = "_"), test.use = "MAST")
head(fm)

ct = "IOL NN"
fmd = FindMarkers(obj, `ident.1` = paste(ct,"2mo", sep = "_"), `ident.2` = paste(ct,"18mo", sep = "_"))
head(fm)

down = (rownames(fm[which(fm$p_val_adj<0.05 & fm$avg_log2FC>0),]))
up= (rownames(fm[which(fm$p_val_adj<0.05 & fm$avg_log2FC<0),]))

downd = (rownames(fmd[which(fmd$p_val_adj<0.05 & fmd$avg_log2FC>0),]))
upd= (rownames(fmd[which(fmd$p_val_adj<0.05 & fmd$avg_log2FC<0),]))

pheatmap(av_og[which(rownames(av_og) %in% upd), grep("IOL", colnames(av_og))], scale = "row")

pheatmap(av_og[which(rownames(av_og) %in% downd), grep("IOL", colnames(av_og))], scale = "row")

head(av_og[which(rownames(av_og) %in% up), grep("IOL", colnames(av_og))])

jaccard_similarity()

# Load required libraries
library(dplyr)
library(pheatmap)

# Set the directory containing your CSV files
file_dir <- "~/projects/merfish/"

# List all CSV files in the directory
file_list <- list.files(path = file_dir, pattern = "*.csv", full.names = TRUE)

# Function to read and extract significant genes (e.g., p_val_adj < 0.05)
extract_genes <- function(file) {
  df <- read.csv(file, header = TRUE)
  significant_genes <- df %>% filter(p_val_adj < 0.05 & avg_log2FC>0) %>% pull(X)  # Assuming gene names are in the first column
  return(significant_genes)
}

# Create a list to store gene sets from each file
gene_sets <- lapply(file_list, extract_genes)

# Function to calculate Jaccard similarity
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

# Create an empty matrix to store Jaccard similarities
n <- length(gene_sets)
jaccard_matrix <- matrix(0, nrow = n, ncol = n)
colnames(jaccard_matrix) <- basename(file_list)
rownames(jaccard_matrix) <- basename(file_list)

# Calculate Jaccard similarities
for (i in 1:n) {
  for (j in i:n) {
    jaccard_matrix[i, j] <- jaccard_similarity(gene_sets[[i]], gene_sets[[j]])
    jaccard_matrix[j, i] <- jaccard_matrix[i, j]
  }
}

# Plot heatmap using pheatmap
pheatmap(jaccard_matrix, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         display_numbers = TRUE, 
         main = "Jaccard Similarity Between Differential Genes")


options(repr.plot.width=19, repr.plot.height=26)
pheatmap(jaccard_matrix, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         display_numbers = TRUE, 
         main = "Jaccard Similarity Between Differential Genes")


head(classified_genes)

# Load required libraries
library(dplyr)
library(pheatmap)
library(ggplot2)

# Set the directory containing your CSV files
file_dir <- "~/projects/merfish/"

# List all CSV files in the directory
file_list <- list.files(path = file_dir, pattern = "*.csv", full.names = TRUE)
file_list = c(file_list[grep("Glut", file_list)])

#file_list = c(file_list[grep("IMN", file_list)], file_list[grep("IOL", file_list)])
# Initialize lists to store upregulated and downregulated genes
upregulated_genes <- list()
downregulated_genes <- list()

# Function to classify genes as upregulated or downregulated
classify_genes <- function(file) {
  df <- read.csv(file, header = TRUE)
  
  # Upregulated genes (e.g., avg_log2FC > 0 and p_val_adj < 0.05)
  upregulated <- df %>% filter(avg_log2FC > 0, p_val_adj < 0.05) %>% pull(X)
  
  # Downregulated genes (e.g., avg_log2FC < 0 and p_val_adj < 0.05)
  downregulated <- df %>% filter(avg_log2FC < 0, p_val_adj < 0.05) %>% pull(X)
  
  return(list(upregulated = upregulated, downregulated = downregulated))
}

# Process each file and classify the genes
for (file in file_list) {
  classified_genes <- classify_genes(file)
  upregulated_genes <- c(upregulated_genes, classified_genes$upregulated)
  downregulated_genes <- c(downregulated_genes, classified_genes$downregulated)
}

# Count occurrences of each gene in upregulated and downregulated lists
upregulated_count <- table(unlist(upregulated_genes))
downregulated_count <- table(unlist(downregulated_genes))

# Find the most frequently upregulated and downregulated genes
most_upregulated <- sort(upregulated_count, decreasing = TRUE)
most_downregulated <- sort(downregulated_count, decreasing = TRUE)

# Display the top 10 most frequently upregulated and downregulated genes
cat("Top 10 Most Frequently Downregulated Genes:\n")
print(head(most_upregulated, 15))

cat("\nTop 10 Most Frequently Upregulated Genes:\n")
print(head(most_downregulated, 15))


options(repr.plot.width=5.5, repr.plot.height=4)

up_df <- as.data.frame(head(most_upregulated, 20))
colnames(up_df) <- c("Gene", "Frequency")
down_df <- as.data.frame(head(most_downregulated, 20))
colnames(down_df) <- c("Gene", "Frequency")

up_plot <- ggplot(up_df, aes(x = reorder(Gene, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Most Frequently Downregulated Genes", x = "Gene", y = "Frequency")
up_plot

down_plot <- ggplot(down_df, aes(x = reorder(Gene, Frequency), y = as.integer(Frequency))) +
  geom_bar(stat = "identity", fill = "tomato") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Most Frequently Upregulated Genes", x = "Gene", y = "Frequency")
down_plot

up_plot_gaba= down_plot
down_plot_gaba = up_plot

up_plot_glut= down_plot
down_plot_glut = up_plot

up_plot_nn= down_plot
down_plot_nn = up_plot

library(gridExtra)

options(repr.plot.width=8, repr.plot.height=7)

p=grid.arrange(up_plot_glut+ggtitle("Top 10 Up Glut"), up_plot_gaba+ggtitle("Top 10 Up Gaba"), up_plot_nn+ggtitle("Top 10 Up NN"), 
             down_plot_glut+ggtitle("Top 10 Down Glut"), down_plot_gaba+ggtitle("Top 10 Down Gaba"), down_plot_nn+ggtitle("Top 10 Down NN"),nrow = 2)

svg("../Figure2-IMN-IOL/MERFISH_top_Diff_genes.svg", height = 7, width =9 )
grid.arrange(up_plot_glut+ggtitle("Top 10 Up Glut"), up_plot_gaba+ggtitle("Top 10 Up Gaba"), up_plot_nn+ggtitle("Top 10 Up NN"), 
             down_plot_glut+ggtitle("Top 10 Down Glut"), down_plot_gaba+ggtitle("Top 10 Down Gaba"), down_plot_nn+ggtitle("Top 10 Down NN"),nrow = 2)
dev.off()

print(head(most_upregulated, 100))

up_plot <- ggplot(most_upregulated, aes(x = reorder(Gene, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Most Frequently Upregulated Genes", x = "Gene", y = "Frequency")




file_list

alldf = list()
pval = 0.5

# Set the directory containing your CSV files
file_dir <- "~/projects/merfish/"

# List all CSV files in the directory
file_list <- list.files(path = file_dir, pattern = "*.csv", full.names = TRUE)

for(f in file_list){
    a = read.csv(f)
       if(nrow(a)<1){
        next
    }
    a$dir = "Up"
    a$dir[which(a$avg_log2FC>0)] = "Down"
    ct = gsub("/home/lamaral/projects/merfish//", "", f)
    ct = gsub("_DEGs.csv", "", ct)
    a$ct = ct
    alldf[[ct]]=a
}

alldf = list()
pval = 0.5

# Set the directory containing your CSV files
file_dir <- "~/projects/merfish/"

# List all CSV files in the directory
file_list <- list.files(path = file_dir, pattern = "*.csv", full.names = TRUE)
file_list = file_list[grep("old", file_list)]

for(f in file_list){
    a = read.csv(f)
       if(nrow(a)<1){
        next
    }
    a$dir = "Up"
    a$dir[which(a$avg_log2FC>0)] = "Down"
    ct = gsub("/home/lamaral/projects/merfish//", "", f)
    ct = gsub("_DEGs_yng_old.csv", "", ct)
    a$ct = ct
    alldf[[ct]]=a
}

all = do.call(rbind, alldf)

options(repr.plot.width=7, repr.plot.height=5)
library(RColorBrewer)

gene = "C4b"
cdh8 = all[which(all$X==gene ),]
cdh8 = cdh8[grep("NN", cdh8$ct),]
cdh8 = cdh8[-grep("csv", cdh8$ct),]
cdh8$p_val_adj = cdh8$p_val_adj+1e-100
colors <- (colorRampPalette(brewer.pal(9, "Oranges")[4:9])(100))
g = ggplot(cdh8, aes(x = ct, y = -avg_log2FC, fill = -log10(p_val_adj))) +
  geom_bar(stat = "identity") +  
#  scale_fill_gradientn(colors = colors) + theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = paste("MERFISH", gene,"Age LogFC (18mo/2mo)"), x = "Region", y = "Avg LogFC (18mo/2mo)") +
  scale_fill_gradientn(colors = colors) + theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability
g
#pdf("~/projects/combined_all/Figures/Figure3/eif5a_merfish.pdf", height = 5, width = 5)
#print(g)
#dev.off()

setwd("~/projects/combined_all/female_RNA/DEG_results_.01_.01/")

files <- list.files(path=".", pattern=".csv", full.names=TRUE)

alldf = list()
pval = 0.5
for(f in files){
    a = read.csv(f)
 #   a <- a[which(a$p_val_adj < pval & abs(a$avg_log2FC) > 0.05), ]
    #if(nrow(a)>15000){
    #    a = a[1:15000,]
    #}
    if(nrow(a)<1){
        next
    }
   # a$pos <- floor(as.integer(sub("(.*):(.*)-(.*)", "\\2", a$feature.name)) / win_size)
   # a$sig = log(a$adjusted.p.value+1e-200)*sign(a$log2.fold_change.)
    a$dir = "Up"
    a$dir[which(a$avg_log2FC>0)] = "Down"
    a$ct = gsub("./", "", f)
    a$ct = gsub(".csv", "", a$ct)
    alldf[[a$ct[1]]]=a
}

setwd("~/projects/combined_all/female_RNA/progenitor_DEGs/")

files <- list.files(path=".", pattern="default.csv", full.names=TRUE)

alldf = list()
pval = 0.5
for(f in files){
    a = read.csv(f)
 #   a <- a[which(a$p_val_adj < pval & abs(a$avg_log2FC) > 0.05), ]
    #if(nrow(a)>15000){
    #    a = a[1:15000,]
    #}
    if(nrow(a)<1){
        next
    }
   # a$pos <- floor(as.integer(sub("(.*):(.*)-(.*)", "\\2", a$feature.name)) / win_size)
   # a$sig = log(a$adjusted.p.value+1e-200)*sign(a$log2.fold_change.)
    a$dir = "Up"
    a$dir[which(a$avg_log2FC>0)] = "Down"
    a$ct = gsub("./", "", f)
    a$ct = gsub("2vsOld_default.csv", "", a$ct)
    alldf[[a$ct[1]]]=a
}

prog= do.call(rbind, alldf)

head(prog)

options(repr.plot.width=15, repr.plot.height=5)


gene = "Neat1"
cdh8 = prog[which(prog$X==gene ),]
ggplot(cdh8, aes(x = ct, y = -avg_log2FC, fill = -log10(p_val_adj+1e-100))) +
  geom_bar(stat = "identity") +  
#  scale_fill_gradientn(colors = colors) + theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = paste(gene,"Age LogFC (18mo/2mo)"), x = "Region", y = "Avg LogFC (18mo/2mo)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability


iol_merfish = all[which(all$ct == "Oligo-NN"), ] 

iol_rna = rna[which(rna$ct == "Oligo-NN"), ] 
iol_rna = prog[which(prog$ct == "Oligo-NN"), ] 

unique(rna$ct)

rna= prog

rna$ct = gsub("_", "-", rna$ct)
all$ct = gsub("--", "-", all$ct)

head(all)


# Rename columns of the first table to include 'merfish'
colnames(all) <- paste0(colnames(all), "_merfish")
colnames(all)[colnames(all) == "X_merfish"] <- "gene"
colnames(all)[colnames(all) == "ct_merfish"] <- "ct"

# Rename columns of the second table to include 'rna'
colnames(rna) <- paste0(colnames(rna), "_rna")
colnames(rna)[colnames(rna) == "X_rna"] <- "gene"
colnames(rna)[colnames(rna) == "ct_rna"] <- "ct"

# Merge the two tables by 'gene' and 'ct'
merged_data <- merge(all, rna, by = c("gene", "ct"), all = TRUE)

# View the merged table
head(merged_data)


head(merged_data)

library(ggplot2)

# Filter the data for the gene "Eif5a"
plot_data <- merged_data[which(merged_data$gene == "C4b"),]

# Create the scatter plot
ggplot(plot_data, aes(x = avg_log2FC_rna, y = avg_log2FC_merfish, label = ct)) +
  geom_point() +
  geom_text(vjust = -1) +
  labs(
    title = "Scatter Plot of logFC for RNA and Merfish (Eif5a)",
    x = "logFC (RNA)",
    y = "logFC (Merfish)"
  ) +
  theme_minimal()


options(repr.plot.width=7, repr.plot.height=15)

library(ggplot2)
library(reshape2)

# Filter data for the specific genes
genes_of_interest <- c("Eif5a", "C4b", "Pisd", "Sirt2")
heatmap_data <- merged_data[merged_data$gene %in% genes_of_interest, ]

# Reshape the data into a long format for ggplot
heatmap_long <- melt(heatmap_data, id.vars = c("gene", "ct"), measure.vars = c("avg_log2FC_rna", "avg_log2FC_merfish"))

# Rename the 'variable' column to distinguish between RNA and Merfish
heatmap_long$variable <- gsub("avg_log2FC_", "", heatmap_long$variable)

# Combine gene and method (RNA/Merfish) into one column for the x-axis
heatmap_long$gene_method <- paste(heatmap_long$gene, heatmap_long$variable, sep = "_")

# Create the heatmap
ggplot(heatmap_long, aes(x = gene_method, y = ct, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  labs(title = "Heatmap of logFC for RNA and Merfish across Genes and Cell Types",
       x = "Gene (Method)",
       y = "Cell Type",
       fill = "logFC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


head(heatmap_data)
head(heatmap_matrix)

head(all)

head(rna)
head(all)

plot_rna_volcano("Oligo-NN")
plot_mer_volcano("Oligo-NN")

library(ggplot2)
library(ggrepel)
options(repr.plot.width=10, repr.plot.height=7)


ct = "IOL-NN"
plot_rna_volcano = function(ct){
    # Filter the data for the specific cell type
    #data_filtered <- all[which(all$ct == "STR-D12-Gaba"),]
    data_filtered <- rna[which(rna$ct == ct),]

    # Identify the top 20 most significant genes based on adjusted p-value
    top_genes <- data_filtered[order(data_filtered$p_val_adj_rna), ][1:20, ]

    # Create the volcano plot with top genes labeled
    volcano_rna <- ggplot(data_filtered, aes(x = -avg_log2FC_rna, y = -log10(p_val_adj_rna+1e-200), color = dir_rna)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("Up" = "#E69F00", "Down" = "#56B4E9")) +
      labs(title = paste("Volcano Plot - RNA -",ct),
           x = "Log2 Fold Change (RNA)",
           y = "-Log10 Adjusted P-Value",
           color = "Direction") +
      theme_minimal(base_size = 15) +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)
      ) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")  +# Significance threshold
      geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = Inf) # +# Label the top 20 genes
      #ylim(c(0,250))
    #  coord_cartesian(ylim = c(0, max(-log10(data_filtered$p_val_adj_merfish)) + 2))  # Adjust y-axis to prevent cutoff

    # Print the plot
    return(volcano_rna)
}
plot_rna_volcano(ct)

library(ggplot2)
library(ggrepel)
options(repr.plot.width=10, repr.plot.height=7)


ct = "IOL-NN"
plot_mer_volcano = function(ct){
    # Filter the data for the specific cell type
    #data_filtered <- all[which(all$ct == "STR-D12-Gaba"),]
    data_filtered <- all[which(all$ct == ct),]

    # Identify the top 20 most significant genes based on adjusted p-value
    top_genes <- data_filtered[order(data_filtered$p_val_adj_merfish), ][1:20, ]

    # Create the volcano plot with top genes labeled
    volcano_rna <- ggplot(data_filtered, aes(x = -avg_log2FC_merfish, y = -log10(p_val_adj_merfish+1e-200), color = dir_merfish)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("Up" = "#E69F00", "Down" = "#56B4E9")) +
      labs(title = paste("Volcano Plot - MERFISH -",ct),
           x = "Log2 Fold Change (RNA)",
           y = "-Log10 Adjusted P-Value",
           color = "Direction") +
      theme_minimal(base_size = 15) +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)
      ) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")  +# Significance threshold
      geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = Inf) # +# Label the top 20 genes
      #ylim(c(0,250))
    #  coord_cartesian(ylim = c(0, max(-log10(data_filtered$p_val_adj_merfish)) + 2))  # Adjust y-axis to prevent cutoff

    # Print the plot
    return(volcano_rna)
}
plot_mer_volcano(ct)

plot_mer_volcano("IOL-NN")
plot_rna_volcano("IOL-NN")

plot_mer_volcano("DG-PIR-Ex-IMN")
plot_rna_volcano("DG-PIR-Ex-IMN")

plot_mer_volcano("OB-STR-CTX-Inh-IMN")
plot_rna_volcano("OB-STR-CTX-Inh-IMN")



plot_mer_volcano("IOL-NN")
plot_rna_volcano("IOL-NN")

plot_mer_volcano("DG-PIR-Ex-IMN")
plot_rna_volcano("DG-PIR-Ex-IMN")

plot_mer_volcano("OB-STR-CTX-Inh-IMN")
plot_rna_volcano("OB-STR-CTX-Inh-IMN")



volcano_rna <- ggplot(all[which(all$ct == "Microglia-NN"),], aes(x = avg_log2FC_merfish, y = -log10(p_val_adj_merfish), color = dir_merfish)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Up" = "#E69F00", "Down" = "#56B4E9")) +
  labs(title = "Volcano Plot - RNA",
       x = "Log2 Fold Change (RNA)",
       y = "-Log10 Adjusted P-Value",
       color = "Direction") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")  # Significance threshold
volcano_rna



rna$label = ""
rna$label[which(rna$p_val_adj_rna<0.001 & abs(rna$avg_log2FC_rna) > .1)] = rna$gene[which(rna$p_val_adj_rna<0.001 &abs(rna$avg_log2FC_rna) > 0.1 )] 
selected_points <- rna[which(rna$label !=""), ]
selected_points <- selected_points %>%
  group_by(ct) %>%
  slice(1:10) %>%
  ungroup()
# Create the ggplot with selected points
ggplot(rna[which(abs(rna$avg_log2FC_rna) > 0 & rna$pct.1_rna>0.0),], aes(x = ct, y = -avg_log2FC_rna,
                 color = ifelse(p_val_adj_rna < 0.1, ifelse(avg_log2FC_rna > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential Expression",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_classic() + ylim(c(-10,10))+
  coord_flip() +
  # Add labels for selected points
  ggrepel::geom_text_repel(data = selected_points, aes(label = label), max.overlaps=10000,
            hjust = -0.2, vjust = -0.2, size = 3, color = "black", fontface = "bold")

selected_points



plot(iol$avg_log2FC_merfish, iol$avg_log2FC_rna)

head(iol)

library(ggplot2)
library(dplyr)
library(ggrepel)

iol = merged_data[which(merged_data$ct == "Oligo-NN"),]
iol = iol[which(!is.na(iol$avg_log2FC_merfish)),]
iol = iol[which(!is.na(iol$avg_log2FC_rna)),]
#iol = iol[which(iol$p_val_adj_merfish < 0.5 | iol$p_val_adj_rna < 0.5),]

cor(iol$avg_log2FC_merfish, iol$avg_log2FC_rna)

# Identify the top 5 highest and lowest points in both directions for merfish and RNA
top_bottom_points <- iol %>%
  arrange(desc(avg_log2FC_merfish)) %>% 
  slice(1:5) %>%
  bind_rows(iol %>% arrange(avg_log2FC_merfish) %>% slice(1:5)) %>%
  bind_rows(iol %>% arrange(desc(avg_log2FC_rna)) %>% slice(1:5)) %>%
  bind_rows(iol %>% arrange(avg_log2FC_rna) %>% slice(1:5)) %>%
  distinct()  # Ensure no duplicates in the combined dataset

# Create the scatter plot with ggplot2
scatter_plot <- ggplot(iol, aes(x = -avg_log2FC_merfish, y = -avg_log2FC_rna)) +
  geom_point(alpha = 0.8, color = "blue") +
  geom_text_repel(data = top_bottom_points, aes(label = gene), size = 3, max.overlaps = Inf) +
  labs(title = "Scatter Plot of logFC: Merfish vs RNA",
       x = "Log2 Fold Change (Merfish)",
       y = "Log2 Fold Change (RNA)") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )

# Print the plot
scatter_plot


library(ggplot2)
library(dplyr)
library(ggrepel)

iol = 
[which(merged_data$ct == "OB-STR-CTX-Inh-IMN"),]
iol = iol[which(!is.na(iol$avg_log2FC_merfish)),]
iol = iol[which(!is.na(iol$avg_log2FC_rna)),]
#iol = iol[which(iol$p_val_adj_merfish < 0.5 | iol$p_val_adj_rna < 0.5),]

cor(iol$avg_log2FC_merfish, iol$avg_log2FC_rna)

# Identify the top 5 highest and lowest points in both directions for merfish and RNA
top_bottom_points <- iol %>%
  arrange(desc(avg_log2FC_merfish)) %>% 
  slice(1:5) %>%
  bind_rows(iol %>% arrange(avg_log2FC_merfish) %>% slice(1:5)) %>%
  bind_rows(iol %>% arrange(desc(avg_log2FC_rna)) %>% slice(1:5)) %>%
  bind_rows(iol %>% arrange(avg_log2FC_rna) %>% slice(1:5)) %>%
  distinct()  # Ensure no duplicates in the combined dataset

# Create the scatter plot with ggplot2
scatter_plot <- ggplot(iol, aes(x = -avg_log2FC_merfish, y = -avg_log2FC_rna)) +
  geom_point(alpha = 0.8, color = "blue") +
  geom_text_repel(data = top_bottom_points, aes(label = gene), size = 3, max.overlaps = Inf) +
  labs(title = "Scatter Plot of logFC: Merfish vs RNA",
       x = "Log2 Fold Change (Merfish)",
       y = "Log2 Fold Change (RNA)") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )

# Print the plot
scatter_plot


#"OB-STR-CTX-Inh-IMN"
#"DG-PIR-Ex-IMN"

options(repr.plot.width=7, repr.plot.height=7)
ct = "IOL-NN"
iol = merged_data[which(merged_data$ct == ct),]
iol <- iol[which(!is.na(iol$avg_log2FC_merfish)),]
iol <- iol[which(!is.na(iol$avg_log2FC_rna)),]

# Calculate Direction * Significance for both Merfish and RNA
iol <- iol %>%
  mutate(
    dir_significance_merfish = -sign(avg_log2FC_merfish) * -log10(p_val_merfish),
    dir_significance_rna = -sign(avg_log2FC_rna) * -log10(p_val_rna)
  )

# Identify the top 5 highest and lowest points in both directions for Merfish and RNA
top_bottom_points <- iol %>%
  arrange(desc(dir_significance_merfish)) %>% 
  slice(1:10) %>%
  bind_rows(iol %>% arrange(dir_significance_merfish) %>% slice(1:10)) %>%
  bind_rows(iol %>% arrange(desc(dir_significance_rna)) %>% slice(1:10)) %>%
  bind_rows(iol %>% arrange(dir_significance_rna) %>% slice(1:10)) %>%
  distinct()  # Ensure no duplicates in the combined dataset

# Create the scatter plot with ggplot2
scatter_plot <- ggplot(iol, aes(x = dir_significance_merfish, y = dir_significance_rna)) +
  geom_point(alpha = 0.8, color = "blue") +
  geom_text_repel(data = top_bottom_points, aes(label = gene), size = 3, max.overlaps = Inf) +
  labs(title = paste(ct,"\nMerfish vs RNA"),
       x = "Direction * Significance (Merfish)",
       y = "Direction * Significance (RNA)") +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  ) +ylim(c(-11,11))+xlim(c(-125,125))

# Print the plot
scatter_plot


#"OB-STR-CTX-Inh-IMN"
#"DG-PIR-Ex-IMN"
iol = merged_data[which(merged_data$ct == "IOL-NN"),]
iol <- iol[which(!is.na(iol$avg_log2FC_merfish)),]
iol <- iol[which(!is.na(iol$avg_log2FC_rna)),]

# Calculate Direction * Significance for both Merfish and RNA
iol <- iol %>%
  mutate(
    dir_significance_merfish = -sign(avg_log2FC_merfish) * -log10(p_val_merfish),
    dir_significance_rna = -sign(avg_log2FC_rna) * -log10(p_val_rna)
  )

iol$significance = "ns"
#iol$significance[which(iol$p_val_adj_merfish<0.01)] = "signifiant merfish"
#iol$significance[which(iol$p_val_adj_rna<0.05)] = "signifiant rna"
iol$significance[which(iol$p_val_adj_rna<0.05 & iol$p_val_adj_merfish<0.05)] = "signifiant both"

# Identify the top 5 highest and lowest points in both directions for Merfish and RNA
top_bottom_points <- iol %>%
  arrange(desc(dir_significance_merfish)) %>% 
  slice(1:20) %>%
  bind_rows(iol %>% arrange(dir_significance_merfish) %>% slice(1:20)) %>%
  bind_rows(iol %>% arrange(desc(dir_significance_rna)) %>% slice(1:20)) %>%
  bind_rows(iol %>% arrange(dir_significance_rna) %>% slice(1:20)) %>%
  distinct()  # Ensure no duplicates in the combined dataset

# Create the scatter plot with ggplot2
scatter_plot <- ggplot(iol, aes(x = avg_log2FC_merfish, y = avg_log2FC_rna,color = significance)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(data = top_bottom_points, aes(label = gene), size = 3, max.overlaps = Inf) +
  labs(title = "Scatter Plot of Direction * Significance: Merfish vs RNA",
       x = "Direction * Significance (Merfish)",
       y = "Direction * Significance (RNA)") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )

# Print the plot
scatter_plot


rownames(obj)

pheatmap(av_og[which(rownames(av_og) %in% c(upd, "Cdkn2a")), c(grep("IMN", colnames(av_og)), grep("IOL", colnames(av_og)))], scale = "row",cluster_cols = F)

degs = unique(merged_data$gene[which(merged_data$p_val_adj_merfish<0.01)])

options(repr.plot.width=6, repr.plot.height=30)

prog_av = av_og[degs, c(grep("IMN", colnames(av_og)), grep("IOL", colnames(av_og)))]
prog_av = prog_av[,c(1,4,6,2,3,5,7,8,9)]
ph = pheatmap(prog_av, scale = "row",cluster_cols = F, cutree_rows = 10)

prog_av_rna = av_rna[degs, c(grep("IMN", colnames(av_rna)), grep("IOL", colnames(av_rna)))]
#prog_av_rna = prog_av_rna[which(!(rowSums(prog_av_rna)==0)),]
prog_av_rna = prog_av_rna[,c(2,4,1,3,5,6,7,8,9)]

head(prog_av_rna)
ph = pheatmap(prog_av_rna
              , scale = "row",cluster_cols = F, cutree_rows = 10)



options(repr.plot.width=6, repr.plot.height=30)

prog_av = av_og[degs, c(grep("IMN", colnames(av_og)), grep("IOL", colnames(av_og)))]
prog_av = prog_av[,c(1,4,6,2,3,5,7,8,9)]
ph = pheatmap(prog_av, scale = "row",cluster_cols = F, cutree_rows = 10,cluster_rows = T)



options(repr.plot.width=6, repr.plot.height=10)

library(dplyr)
library(tidyr)
library(tibble)

# Extract the gene of interest (e.g., Sirt2)
gene_of_interest <- "Sirt2"
df <- av_rna[gene_of_interest, , drop = FALSE]

# Adjust the column names to include only cell type and age
# This assumes that age information is always at the end of the string and is in the format "Xmo"
colnames(df) <- gsub("(.*?)([0-9]+mo).*", "\\1 \\2", colnames(df))
head(df)

# Convert the data frame to a long format
df_long <- as.data.frame(t(df)) %>%
  rownames_to_column(var = "Cell_Age") %>%
  separate(Cell_Age, into = c("CellType", "Age"), sep = "\\s+(?=[0-9])") %>%
  mutate(Age = factor(Age, levels = c("2mo", "9mo", "18mo")))  # Ensure age is ordered correctly

# Reshape the data to have ages as columns and cell types as rows
df_wide <- df_long %>%
  spread(key = Age, value = Sirt2)
sirt2_rna = df_wide

rownames(df_wide) = df_wide$CellType
pheatmap(df_wide[,2:4], scale= "none", cluster_cols = F)

m = melt(b)
# Assuming your data is stored in a data frame called 'df'
# Create the scatter plot
ggplot(m, aes(x = variable, y = value, color = CellType)) +
  geom_point(size = 3) +  # Adjust the size of points as needed
  theme_minimal() +
  xlab("Variable") +
  ylab("Value") +
  ggtitle("Scatter Plot of Values by Cell Type and Age") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        legend.position = "right") +
  scale_shape_manual(values = c(16, 17, 18)) + NoLegend() # Customize shapes for different ages


rownames(sirt2_rna) = sirt2_rna$CellType
rownames(il33_rna) = il33_rna$CellType


c(rep(c(rep("Sirt2", 3), rep("Il33", 3)),2))

paste(colnames(b), c(rep(c(rep("Sirt2", 3), rep("Il33", 3)),2)))

options(repr.plot.width=4.5, repr.plot.height=7)

b = cbind(sirt2[,2:4], il33[,2:4], sirt2_rna[sirt2$CellType,2:4], il33_rna[sirt2$CellType,2:4])
colnames(b) = paste(colnames(b), c(rep(c(rep("Sirt2", 3), rep("Il33", 3)),2)))
rownames(b) = gsub("-$", "", sirt2$CellType)
b = b[over1k,]
b = b[order(b[,1],decreasing = T),]

pheatmap(b[,1:3], cluster_rows = F, main = "Sirt2 MERFISH")
pheatmap(b[,4:6], cluster_rows = F, main = "Il33 MERFISH")

pheatmap(b[,7:9], cluster_rows = F, main = "Sirt2 RNA")
pheatmap(b[,10:12], cluster_rows = F, main = "Il33 RNA")


both_av = cbind(scale(prog_av), scale(prog_av_rna))
ph = pheatmap(both_av[,], scale = "row",cluster_cols = F, cluster_rows = T)

head(merged_data)

dg_rna= merged_data[which(merged_data$ct == 'IOL-NN' & !is.na(merged_data$avg_log2FC_rna)), ] 

dg_rna = dg_rna[order(dg_rna$avg_log2FC_rna),]
head(dg_rna)

dg_rna = dg_rna[which(dg_rna$p_val_adj_rna<1),]

library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse gene annotations


# Sort the genes by logFC for GSEA
geneList <- dg_rna$avg_log2FC_rna
names(geneList) <- dg_rna$gene
geneList <- sort(geneList, decreasing = TRUE)


head(geneList)

gsea_results <- gseGO(
  geneList = geneList,
  OrgDb = org.Mm.eg.db,  # Mouse-specific OrgDb
  keyType = "SYMBOL",    # Assuming gene symbols are used, adjust if using different key types
  ont = "ALL",           # Options: "BP", "CC", "MF", or "ALL"
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

options(repr.plot.width=6, repr.plot.height=6)

dotplot(gsea_results)


head(gsea_results)


