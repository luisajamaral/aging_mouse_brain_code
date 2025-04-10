library(ggplot2)
library(data.table)
library(dplyr)

   subcluster="2"
curr = fread(paste(subcluster,"_combined_meta.csv", sep = ""), sep = ",")
#curr$leiden = curr$`leiden-.11` 
#curr = curr[,-c("leiden-.11")]
#curr = curr[,-c("leiden-.1")]
    metaf = curr
    metaf = metaf[which(metaf$subclass_label_v3!= "none"),]
    predictions <- table(metaf$leiden,metaf$subclass_label_v3)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
      group_by(Var1) %>%
      filter(Freq == max(Freq)) %>%
      select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    rownames(new_df) = new_df$Var1

    curr$best_subclass_label_V3 = new_df[paste(curr$leiden), "Var2"]
    curr$best_subclass_label_V3 = as.character(curr$best_subclass_label_V3)

    metaf = curr
    metaf = metaf[which(metaf$L2Annot.rough!= "none"),]
    predictions <- table(metaf$leiden,metaf$L2Annot.rough)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
      group_by(Var1) %>%
      filter(Freq == max(Freq)) %>%
      select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    rownames(new_df) = new_df$Var1

    curr$best_L2Annot.rough = new_df[paste(curr$leiden), "Var2"]
    curr$best_L2Annot.rough = as.character(curr$best_L2Annot.rough)
    curr$leiden_subcluster = paste(subcluster,'-',curr$leiden, sep = "")
    #list_currs[[subcluster]] = curr

meta=read.table("../meta_best_L2_allen2.txt")

head(curr)

head(meta)

rownames(meta) = paste(meta$cell_id)

head(meta)

table(meta[paste(curr$cell_id), "leiden_subcluster"])

meta[paste(curr$cell_id), "leiden_subcluster"] = curr$leiden_subcluster

table(meta$leiden_subcluster)

meta$cluster_for_peakcall  = paste(meta$leiden_subcluster, meta$age)

hist(table(meta$cluster_for_peakcall ), breaks = 50)

table(meta$leiden_subcluster )[which(table(meta$leiden_subcluster) <500)]

write.table(meta , "../meta_subcluster_use.txt", row.names = F)

head(meta[which(sapply(strsplit(as.character(meta$leiden_subcluster), "-"), "[[", 1) == "2"),])

unique(meta2$leiden_subcluster)

smet = meta2[which(sapply(strsplit(as.character(meta2$leiden_subcluster), "-"), "[[", 1) == "2_res0.25"),]

smet = as.data.frame(smet)

mat = match(smet$cell_id, curr$cell_id)

all(smet$cell_id ==curr$cell_id)

smet$cell_id

smet$lei_sub = curr[,"leiden_subcluster"]

mat = match(smet$cell_id,meta2$cell_id)

mat

meta2$subcluster_leiden= meta2$leiden_subcluster
meta2$subcluster_leiden[mat] = smet$lei_sub

table(meta2$subcluster_leiden)

write.table(meta2, file = "meta_fixed.txt", row.names = F)

scatter_plot <- ggplot(data = curr, aes(x = umap_x, y = umap_y, color = as.factor(leiden_subcluster))) +
      geom_point() +
      scale_color_discrete(name = "Cell Type") +  # Customize legend title
      labs(x = "UMAP X", y = "UMAP Y") +  # Customize axis labels
      theme_minimal() +ggtitle("leiden_subcluster")  # Customize the plot theme if needed
# Calculate the center of each group
center_points <- curr %>%
  group_by(leiden_subcluster) %>%
  summarize(center_x = mean(umap_x), center_y = mean(umap_y))

# Add labels for each group
scatter_plot <- scatter_plot +
  geom_text(data = center_points, aes(x = center_x, y = center_y, label = leiden_subcluster), vjust = 0.5, color = "black")

# Display the plot
scatter_plot

table(curr$leiden_subcluster)

head(curr)

scatter_plot <- ggplot(data = curr, aes(x = umap_x, y = umap_y, color = as.factor(leiden))) +
      geom_point() +
      scale_color_discrete(name = "Cell Type") +  # Customize legend title
      labs(x = "UMAP X", y = "UMAP Y") +  # Customize axis labels
      theme_minimal() +ggtitle("leiden")  # Customize the plot theme if needed
# Calculate the center of each group
center_points <- curr %>%
  group_by(leiden) %>%
  summarize(center_x = mean(umap_x), center_y = mean(umap_y))

# Add labels for each group
scatter_plot <- scatter_plot +
  geom_text(data = center_points, aes(x = center_x, y = center_y, label = leiden), vjust = 0.5, color = "black")

# Display the plot
scatter_plot

 ggplot(data = curr, aes(x = umap_x, y = umap_y, color = best_subclass_label_V3)) +
      geom_point() +
      scale_color_discrete(name = "Cell Type") +  # Customize legend title
      labs(x = "UMAP X", y = "UMAP Y") +  # Customize axis labels
      theme_minimal() +ggtitle("age")  + facet_wrap(~region) # Customize the plot theme if needed
#

 ggplot(data = curr, aes(x = umap_x, y = umap_y, color = age)) +
      geom_point() +
      scale_color_discrete(name = "Cell Type") +  # Customize legend title
      labs(x = "UMAP X", y = "UMAP Y") +  # Customize axis labels
      theme_minimal() +ggtitle("age")  # Customize the plot theme if needed
#

scatter_plot <- ggplot(data = curr, aes(x = umap_x, y = umap_y, color = as.factor(region))) +
      geom_point() +
      scale_color_discrete(name = "Cell Type") +  # Customize legend title
      labs(x = "UMAP X", y = "UMAP Y") +  # Customize axis labels
      theme_minimal() +ggtitle("leiden_1")  # Customize the plot theme if needed
scatter_plot

subcluster="2_res0.25"
curr = fread(paste(subcluster,"_combined_meta.csv", sep = ""), sep = ",")
curr$leiden = curr$`leiden-.11` 
curr = curr[,-c("leiden-.11")]
curr = curr[,-c("leiden-.1")]

 head(curr)

list_currs = list()

for(subcluster in 0:18) {
    if(subcluster ==2) {
        subcluster="2_res0.25"
        curr = fread(paste(subcluster,"_combined_meta.csv", sep = ""), sep = ",")
        curr$leiden = curr$`leiden-.11` 
        curr = curr[,-c("leiden-.11")]
        curr = curr[,-c("leiden-.1")]
    } 
    else if (subcluster == 0) {
        subcluster="OPC_OLG"
        curr = fread(paste(subcluster,"_combined_meta.csv", sep = ""), sep = ",")
    } else if (subcluster == 11) {
        next
    }
    else{
        subcluster=as.character(subcluster)
        curr = fread(paste(subcluster,"_combined_meta.csv", sep = ""), sep = ",")
    }
    metaf = curr
    metaf = metaf[which(metaf$subclass_label_v3!= "none"),]
    predictions <- table(metaf$leiden,metaf$subclass_label_v3)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
      group_by(Var1) %>%
      filter(Freq == max(Freq)) %>%
      select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    rownames(new_df) = new_df$Var1

    curr$best_subclass_label_V3 = new_df[paste(curr$leiden), "Var2"]
    curr$best_subclass_label_V3 = as.character(curr$best_subclass_label_V3)

    metaf = curr
    metaf = metaf[which(metaf$L2Annot.rough!= "none"),]
    predictions <- table(metaf$leiden,metaf$L2Annot.rough)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
      group_by(Var1) %>%
      filter(Freq == max(Freq)) %>%
      select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    rownames(new_df) = new_df$Var1

    curr$best_L2Annot.rough = new_df[paste(curr$leiden), "Var2"]
    curr$best_L2Annot.rough = as.character(curr$best_L2Annot.rough)
    curr$leiden_subcluster = paste(subcluster,'-',curr$leiden, sep = "")
    list_currs[[subcluster]] = curr
    scatter_plot <- ggplot(data = curr, aes(x = umap_x, y = umap_y, color = best_subclass_label_V3)) +
      geom_point() +
      scale_color_discrete(name = "Cell Type") +  # Customize legend title
      labs(x = "UMAP X", y = "UMAP Y") +  # Customize axis labels
      theme_minimal() +ggtitle("best_subclass_label_V3")  # Customize the plot theme if needed

    scatter_plotb <- ggplot(data = curr, aes(x = umap_x, y = umap_y, color = best_L2Annot.rough)) +
      geom_point() +
      scale_color_discrete(name = "Cell Type") +  # Customize legend title
      labs(x = "UMAP X", y = "UMAP Y") +  # Customize axis labels
      theme_minimal() +ggtitle("best_L2Annot.rough") # Customize the plot theme if needed

    scatter_plotl <- ggplot(data = curr, aes(x = umap_x, y = umap_y, color = as.factor(leiden))) +
      geom_point() +
      scale_color_discrete(name = "Cell Type") +  # Customize legend title
      labs(x = "UMAP X", y = "UMAP Y") +  # Customize axis labels
      theme_minimal()  # Customize the plot theme if needed
   # png(paste(subcluster, "_UMAP.png", sep = ""), height = 600, width = 700, units = "px", res = 300)
   # print(scatter_plot)
   # print(scatter_plotb)
   # print(scatter_plotl)
   # dev.off()

}

print("out")

comb =  do.call(rbind, list_currs)


nrow(comb)
length(unique(comb$best_L2Annot.rough))

length(unique(comb$best_subclass_label_V3))

head(curr)
options(repr.plot.width=16, repr.plot.height=8)

ggplot(comb, aes(fill = batch , x = factor(best_subclass_label_V3))) +
  geom_bar(position = "fill", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
comb$age = factor(comb$age, levels = c("2mo", "9mo", "18mo"))
ggplot(comb, aes(fill = age , x = factor(best_subclass_label_V3))) +
  geom_bar(position = "fill", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

comb$age = factor(comb$age, levels = c("2mo", "9mo", "18mo"))
ggplot(comb, aes(fill = age , x = factor(best_L2Annot.rough))) +
  geom_bar(position = "fill", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

head(comb)

meta = fread("../meta_best_L2_allen.txt")

mat = match(meta$cell_id, comb$cell_id, nomatch = NA)


meta$best_subclass_label_V3 = comb[mat, "best_subclass_label_V3"]
meta$best_L2Annot.rough = comb[mat, "best_L2Annot.rough"]
meta$leiden_subcluster = comb[mat, "leiden_subcluster"]

meta$best_celltype_fixed = meta$best_subclass_label_V3
meta$best_celltype_fixed[which(meta$best_celltype_fixed%in%c("STR D1 Gaba", "STR D2 Gaba"))] = "STR D12 Gaba"
meta$best_celltype_fixed[which(meta$leiden_subcluster%in%c("2_res0.25-0"))] = "PAG-PCG Neur"
meta$best_celltype_fixed[which(meta$leiden_subcluster%in%c("2_res0.25-1"))] = "AMY Neur"
meta$best_celltype_fixed[which(meta$best_L2Annot.rough%in%c("RGL"))] = "RGL"
meta$best_celltype_fixed[which(meta$best_L2Annot.rough%in%c("IOL"))] = "IOL"


head(meta)
write.table(meta,"../meta_best_L2_allen2.txt", sep ="\t")

unique(meta$best_subclass_label_V3)



unique(meta$n_fragment)

result <- meta %>%
  group_by(leiden_subcluster) %>%
  summarize(max_percentage = max(table(sample) / n()))


result <- meta %>%
  group_by(leiden_subcluster) %>%
  summarize(
    max_percentage = max(table(sample) / n()),
    avg_nfrag = mean(n_fragment),
    avg_tsse = mean(tsse)
  )


result[order(result$max_percentage, decreasing = T),]



meta$best_celltype_fixed = meta$best_subclass_label_V3
meta$best_celltype_fixed[which(meta$best_celltype_fixed%in%c("STR D1 Gaba", "STR D2 Gaba"))] = "STR D12 Gaba"
meta$best_celltype_fixed[which(meta$leiden_1%in%c(2))] = "PAG-PCG Neur"
meta$best_celltype_fixed[which(meta$leiden_1%in%c(13))] = "AMY Neur"
meta$best_celltype_fixed[which(meta$best_L2Annot.rough%in%c("RGL"))] = "RGL"
meta$best_celltype_fixed[which(meta$best_L2Annot.rough%in%c("IOL"))] = "IOL"


table(meta$best_L2Annot.rough)

table(meta$best_celltype_fixed)

length(unique(meta$best_celltype_fixed))

unique(meta$rachel_celltypes)[which(unique(meta$rachel_celltypes) %in% meta$best_subclass_label_V3)]

unique(meta$rachel_celltypes)[which(! unique(meta$rachel_celltypes) %in% meta$best_subclass_label_V3)]

meta = read.table("../meta_best_L2_allen.txt")

meta2 = read.table("../meta_best_L2_allen2.txt")
table(meta2$best_celltype_fixed)

table(meta$best_celltype_fixed)


