library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
setwd("~/ps-renlab2/projects/combined_all/")

rmeta = read.csv("female_RNA/meta.data_final.csv")
head(rmeta$predicted.id)

meta = read.table("meta_celltype_ext.txt")
head(meta)

head(rmeta$barcode)
mat = match(meta$cell_id, paste("Female:",rmeta$barcode,'-1', sep = ""))
length(mat)
meta$predicted.id = rmeta$predicted.id[mat]

length(which(meta$cell_id%in%paste("Female:",rmeta$barcode,'-1', sep = "")))


metaf = meta
metaf = metaf[which(!is.na(metaf$CellType_1127)),]
metaf$CellType_1127[which(metaf$CellType_1127%in%c("STR D1 Gaba", "STR D2 Gaba"))] ="STR D12 Gaba"
predictions <- table(metaf$leiden_subcluster,metaf$CellType_1127)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

new_df <- predictions %>%
group_by(Var1) %>%
filter(Freq == max(Freq)) %>%
select(Var1, Var2, Freq)

new_df = as.data.frame(new_df)
new_df = new_df[which(new_df$Freq>0.5),]
rownames(new_df) = new_df$Var1

meta$CellType_1127_ext = new_df[paste(meta$`leiden_subcluster`), "Var2"]
meta$CellType_1127_ext = as.character(meta$CellType_1127_ext)
meta$CellType_1127_ext[which(is.na(meta$CellType_1127_ext))] = meta$leiden_subcluster[which(is.na(meta$CellType_1127_ext))]


for(i in 0:18) {

    pdf(paste("subclustering3/",i, "_sub_annotation.pdf", sep = ""), height = 6, width = 12)
    # Read data
    m0 <- read.csv(paste("subclustering3/",i,"_combined_meta.csv", sep = ""))

    # Assuming 'meta' is your metadata object
    mat <- match(m0$cell_id, meta$cell_id)
    m0$celltype <- meta$celltype[mat]
    m0$leiden_subcluster <- meta$leiden_subcluster[mat]
    m0$CellType_1127 <- meta$CellType_1127[mat]
    m0$predicted.id = meta$predicted.id[mat]
    m0$predicted_id_sep <- meta$predicted_id_sep[mat]
    m0$present = "No"
    m0$present[which(!is.na(m0$celltype))]="Yes" 

    
    
    #ext predicted_id

    metaf = m0
    metaf = metaf[which(!is.na(m0$predicted.id)),]
    metaf$predicted.id[which(metaf$predicted.id%in%c("STR D1 Gaba", "STR D2 Gaba"))] ="STR D12 Gaba"

    predictions <- table(metaf$leiden,metaf$predicted.id)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    new_df = new_df[which(new_df$Freq>0.45),]
    rownames(new_df) = new_df$Var1

    m0$predicted.id_ext = new_df[paste(m0$`leiden`), "Var2"]
    m0$predicted.id_ext = as.character(m0$predicted.id_ext)
    m0$predicted.id_ext[which(is.na(m0$predicted.id_ext))] = m0$leiden[which(is.na(m0$predicted.id_ext))]
        pr = predictions[which(predictions$Freq>0.1),]
        pr = pr[order(pr$Var1, pr$Freq),]
        write.table(pr, paste("subclustering3/",i, "predicted.id_ext.txt", sep = ""))

    
    #ext predicted_id_sep

    metaf = m0
    metaf = metaf[which(!is.na(m0$predicted_id_sep)),]
    metaf$predicted_id_sep[which(metaf$predicted_id_sep%in%c("STR D1 Gaba", "STR D2 Gaba"))] ="STR D12 Gaba"

    predictions <- table(metaf$leiden,metaf$predicted_id_sep)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    new_df = new_df[which(new_df$Freq>0.45),]
    rownames(new_df) = new_df$Var1

    m0$predicted_id_sep_ext = new_df[paste(m0$`leiden`), "Var2"]
    m0$predicted_id_sep_ext = as.character(m0$predicted_id_sep_ext)
    m0$predicted_id_sep_ext[which(is.na(m0$predicted_id_sep_ext))] = m0$leiden[which(is.na(m0$predicted_id_sep_ext))]
        pr = predictions[which(predictions$Freq>0.1),]
        pr = pr[order(pr$Var1, pr$Freq),]
        write.table(pr, paste("subclustering3/",i, "predicted_id_ext.txt", sep = ""))

    #ext celltype

    metaf = m0
    metaf = metaf[which(!is.na(m0$celltype)),]
    predictions <- table(metaf$leiden,metaf$celltype)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    new_df = new_df[which(new_df$Freq>0.45),]
    rownames(new_df) = new_df$Var1

    m0$celltypes_ext = new_df[paste(m0$`leiden`), "Var2"]
    m0$celltypes_ext = as.character(m0$celltypes_ext)
    m0$celltypes_ext[which(is.na(m0$celltypes_ext))] = m0$leiden[which(is.na(m0$celltypes_ext))]
      pr = predictions[which(predictions$Freq>0.1),]
        pr = pr[order(pr$Var1, pr$Freq),]
        write.table(pr, paste("subclustering3/",i, "celltypes_ext.txt", sep = ""))


    #ext rachel

    metaf = m0
    metaf = metaf[which(!is.na(m0$CellType_1127)),]
    metaf$CellType_1127[which(metaf$CellType_1127%in%c("STR D1 Gaba", "STR D2 Gaba"))] ="STR D12 Gaba"
    predictions <- table(metaf$leiden,metaf$CellType_1127)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    new_df = new_df[which(new_df$Freq>0.45),]
    rownames(new_df) = new_df$Var1

    m0$CellType_1127_ext = new_df[paste(m0$leiden), "Var2"]
    m0$CellType_1127_ext = as.character(m0$CellType_1127_ext)
    m0$CellType_1127_ext[which(is.na(m0$CellType_1127_ext))] = m0$leiden[which(is.na(m0$CellType_1127_ext))]

        pr = predictions[which(predictions$Freq>0.1),]
        pr = pr[order(pr$Var1, pr$Freq),]
        write.table(pr, paste("subclustering3/",i, "CellType_1127_ext.txt", sep = ""))




    # Create UMAP scatter plot with ggplot

    gc = names(table(m0$predicted_id_sep)[which(table(m0$predicted_id_sep)>50)])


    cluster_centers <- m0[which(m0$predicted_id_sep%in%gc),]%>%
      group_by(predicted_id_sep) %>%
      summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

    umap_plot <- ggplot(m0[which(m0$predicted_id_sep%in%gc),], aes(x = umap_x, y = umap_y, color = predicted_id_sep)) +
      geom_point() + theme_bw()+
      # Add labels to the center of each cluster
      geom_text(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = predicted_id_sep),color = "black" ,
                position = position_nudge(y = 0.05),size = 2 )  # Adjust the 'y' v
    print(umap_plot)


    gc = names(table(m0$predicted.id)[which(table(m0$predicted.id)>50)])


    cluster_centers <- m0[which(m0$predicted.id%in%gc),]%>%
      group_by(predicted.id) %>%
      summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

    umap_plot <- ggplot(m0[which(m0$predicted.id%in%gc),], aes(x = umap_x, y = umap_y, color = predicted.id)) +
      geom_point() + theme_bw()+
      # Add labels to the center of each cluster
      geom_text(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = predicted.id),color = "black" ,
                position = position_nudge(y = 0.05),size = 2 )  # Adjust the 'y' v
    print(umap_plot)


    umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = celltype)) +
      geom_point() +
      theme_bw()
    print(umap_plot)

    gc = names(table(m0$celltype)[which(table(m0$celltype)>50)])

    cluster_centers <- m0[which(m0$celltype%in%gc),] %>%
      group_by(celltype) %>%
      summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

    umap_plot <- ggplot(m0[which(m0$celltype%in%gc),], aes(x = umap_x, y = umap_y, color = celltype)) +
      geom_point() +
      theme_bw() +
      # Add labels to the center of each cluster
      geom_text_repel(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = celltype),
                      color = "black", position = position_nudge(y = 0.0), size = 2)
    print(umap_plot)

    cluster_centers <- m0[which(!is.na(m0$predicted.id_ext)),] %>%
      group_by(predicted.id_ext) %>%
      summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

    umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = predicted.id_ext)) +
      geom_point() +
      theme_bw() +
      # Add labels to the center of each cluster
      geom_text_repel(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = predicted.id_ext),
                      color = "black", position = position_nudge(y = 0.0), size = 2)
    print(umap_plot)
    
    
    cluster_centers <- m0[which(!is.na(m0$predicted_id_sep_ext)),] %>%
      group_by(predicted_id_sep_ext) %>%
      summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

    umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = predicted_id_sep_ext)) +
      geom_point() +
      theme_bw() +
      # Add labels to the center of each cluster
      geom_text_repel(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = predicted_id_sep_ext),
                      color = "black", position = position_nudge(y = 0.0), size = 2)
    print(umap_plot)


    cluster_centers <- m0[which(!is.na(m0$celltypes_ext)),] %>%
      group_by(celltypes_ext) %>%
      summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

    umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = celltypes_ext)) +
      geom_point() +
      theme_bw() +
      # Add labels to the center of each cluster
      geom_text_repel(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = celltypes_ext),
                      color = "black", position = position_nudge(y = 0.0), size = 2)
    print(umap_plot)


    cluster_centers <- m0[which(!is.na(m0$CellType_1127_ext)),] %>%
      group_by(CellType_1127_ext) %>%
      summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

    umap_plot <- ggplot(m0[which(!is.na(m0$CellType_1127_ext)),], aes(x = umap_x, y = umap_y, color = CellType_1127_ext)) +
      geom_point() +
      theme_bw() +
      # Add labels to the center of each cluster
      geom_text_repel(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = CellType_1127_ext),
                      color = "black", position = position_nudge(y = 0.0), size = 2)
    print(umap_plot)


    cluster_centers <- m0%>%
      group_by(leiden) %>%
      summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))


    umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = factor(leiden))) +
      geom_point() +
      theme_bw() +
      # Add labels to the center of each cluster
      geom_text_repel(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = leiden),
                      color = "black", position = position_nudge(y = 0.0), size = 2)
    print(umap_plot)


    umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = age)) +
      geom_point() +
      theme_bw()
    print(umap_plot)


    umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = region)) +
      geom_point() +
      theme_bw()
    print(umap_plot)


    umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, alpha = .5,color = doublet_probability)) +
      geom_point() +
      scale_color_gradient(low = "white", high = "red") +  # Customize the color scale
      theme_bw()

    print(umap_plot)

    
    m0 = m0[which(m0$batch == "Female"), ]
    pr_tab = table(m0$leiden, m0$present)/rowSums(table(m0$leiden, m0$present))
    hist(pr_tab[,2])
    pr_tab[which(pr_tab[,2] < .5),, drop = F]
    m0$bad = "no"
    m0$bad[which(m0$leiden %in% rownames(pr_tab[which(pr_tab[,2] < .12),,drop = F]))]= "yes"
    bad_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = bad)) +
      geom_point() +
      theme_bw()
    print(bad_plot)

    rownames(pr_tab[which(pr_tab[,2] < .1),,drop = F])

    umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = as.factor(present), alpha = 0.5)) +
      geom_point() +
      theme_bw()
    print(umap_plot)

    dev.off()
}

pdf("rachel_clusters.pdf", height = 6,width = 7.5)
for(i in 0:18) {
    # Read data
    m0 <- read.csv(paste("subclustering3/",i,"_combined_meta.csv", sep = ""))

    # Assuming 'meta' is your metadata object
    mat <- match(m0$cell_id, meta$cell_id)
    m0$celltype <- meta$celltype[mat]
    m0$leiden_subcluster <- meta$leiden_subcluster[mat]
    m0$CellType_1127 <- meta$CellType_1127[mat]
    m0$predicted.id = meta$predicted.id[mat]
    m0$predicted_id_sep <- meta$predicted_id_sep[mat]
    m0$present = "No"
    m0$present[which(!is.na(m0$celltype))]="Yes" 

     #ext rachel

    metaf = m0
    metaf = metaf[which(!is.na(m0$CellType_1127)),]
    #metaf$CellType_1127[which(metaf$CellType_1127%in%c("STR D1 Gaba", "STR D2 Gaba"))] ="STR D12 Gaba"
    predictions <- table(metaf$leiden,metaf$CellType_1127)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    new_df = new_df[which(new_df$Freq>0.45),]
    rownames(new_df) = new_df$Var1

    m0$CellType_1127_ext = new_df[paste(m0$leiden), "Var2"]
    m0$CellType_1127_ext = as.character(m0$CellType_1127_ext)
    m0$CellType_1127_ext[which(is.na(m0$CellType_1127_ext))] = m0$leiden[which(is.na(m0$CellType_1127_ext))]
    cluster_centers <- m0[which(!is.na(m0$CellType_1127_ext)),] %>%
       group_by(CellType_1127_ext) %>%
       summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

       umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = CellType_1127_ext)) +
       geom_point() +
       theme_bw() +
      # Add labels to the center of each cluster
      geom_text_repel(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = CellType_1127_ext),
                      color = "black", position = position_nudge(y = 0.0), size = 2) + ggtitle(i)
      print(umap_plot)
    
    }
dev.off()


dev.off()


meta$doub = ""
meta$celltype_final = ""

pdf("final_clusters.pdf", height = 6,width = 7.5)
for(i in 0:18) {
    # Read data
    m0 <- read.csv(paste("subclustering3/",i,"_combined_meta.csv", sep = ""))

    # Assuming 'meta' is your metadata object
    mat <- match(m0$cell_id, meta$cell_id)
    m0$celltype <- meta$celltype[mat]
    m0$leiden_subcluster <- meta$leiden_subcluster[mat]
    m0$CellType_1127 <- meta$CellType_1127[mat]
    m0$predicted.id = meta$predicted.id[mat]
    m0$predicted_id_sep <- meta$predicted_id_sep[mat]
    m0$present = "No"
    m0$present[which(!is.na(m0$celltype))]="Yes" 

    #ext predicted_id_sep

    metaf = m0
    metaf = metaf[which(!is.na(m0$predicted_id_sep)),]
    metaf$predicted_id_sep[which(metaf$predicted_id_sep%in%c("STR D1 Gaba", "STR D2 Gaba"))] ="STR D12 Gaba"

    predictions <- table(metaf$leiden,metaf$predicted_id_sep)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    new_df = new_df[which(new_df$Freq>0.45),]
    rownames(new_df) = new_df$Var1

    m0$predicted_id_sep_ext = new_df[paste(m0$`leiden`), "Var2"]
    m0$predicted_id_sep_ext = as.character(m0$predicted_id_sep_ext)
    m0$predicted_id_sep_ext[which(is.na(m0$predicted_id_sep_ext))] = m0$leiden[which(is.na(m0$predicted_id_sep_ext))]
    
    #ext celltype

    metaf = m0
    metaf = metaf[which(!is.na(m0$celltype)),]
    predictions <- table(metaf$leiden,metaf$celltype)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    new_df = new_df[which(new_df$Freq>0.45),]
    rownames(new_df) = new_df$Var1

    m0$celltypes_ext = new_df[paste(m0$`leiden`), "Var2"]
    m0$celltypes_ext = as.character(m0$celltypes_ext)
    m0$celltypes_ext[which(is.na(m0$celltypes_ext))] = m0$leiden[which(is.na(m0$celltypes_ext))]

    m0$doub = "F"
    
    if(i == 0) {
       m0$doub[ which(m0$leiden == 29 ) ] = "T"
       show(table(m0$doub))
    }
    if(i == 1) {
       m0$predicted_id_sep_ext[which(m0$leiden %in% c(22))] = "STR D12 Gaba"

    }
    if(i == 2) {
       m0$doub[ which(m0$leiden == 8 ) ] = "T"
       show(table(m0$doub))
        m0$celltypes_ext[which(m0$leiden %in% c(6))] = "RLP Neur"
        m0$celltypes_ext[which(m0$leiden %in% c(22))] = "CEA-BST Gaba"
        m0$celltypes_ext[which(m0$leiden %in% c(38))] = "CB Granule Glut"
        m0$celltypes_ext[which(m0$leiden %in% c(32))] = "MEA-BST Gaba"




    }
    if(i == 3) {
       m0$doub[ which(m0$leiden == 14 ) ] = "T"
       show(table(m0$doub))
       m0$predicted_id_sep_ext[which(m0$leiden %in% c(12,15))] = "OB-STR-CTX Inh IMN"
       m0$predicted_id_sep_ext[which(m0$leiden %in% c(16))] = "Astroependymal NN"

    }
    if(i == 4) {
       m0$doub[ which(m0$leiden == 17 ) ] = "T"
       show(table(m0$doub))
    }
    if(i == 5) {
       m0$doub[ which(m0$leiden %in% c(12, 22) ) ] = "T"
       show(table(m0$doub))
    }
    if(i == 6) {
       m0$doub[ which(m0$leiden == 7 & m0$doublet_probability>0.25) ] = "T"
       #show(table(m0$doub))
       #m0$celltypes_ext[which(m0$leiden %in% c(7))] = "L2/3 IT ENT Glut"

       m0$celltypes_ext[which(m0$leiden %in% c(10,11))] = "L6 IT CTX Glut"
       m0$celltypes_ext[which(m0$leiden %in% c(20))] = "IT EP-CLA Glut"
    }
    if(i == 7) {
       m0$doub[ which(m0$leiden %in% c(15) ) ] = "T"
       show(table(m0$doub))
       m0$predicted_id_sep_ext[which(m0$leiden %in% c(18))] = "L2/3 IT ENT Glut"
       m0$predicted_id_sep_ext[which(m0$leiden %in% c(3,12))] = "L2/3 IT PPP Glut"
    }
    if(i == 9) {
       m0$doub[ which(m0$leiden %in% c(19) ) ] = "T"
       show(table(m0$doub))
    }
    if(i == 10) {
       m0$doub[ which(m0$leiden %in% c(10) ) ] = "T"
       show(table(m0$doub))
       m0$celltypes_ext[which(m0$leiden %in% c(17))] = "Lymphoid"

    }
    if(i == 11) {
       m0$doub[ which(m0$leiden %in% c(12,13) ) ] = "T"
       show(table(m0$doub))


    }
    if(i == 14) {
       m0$predicted_id_sep_ext[which(m0$leiden %in% c(12))] = "ABC NN"

    }
    if(i == 15) {
       m0$predicted_id_sep_ext[which(m0$leiden %in% c(13))] = "L5 ET CTX Glut"

    }
    if(i == 17) {
       m0$doub[ which(m0$leiden %in% c(6) ) ] = "T"
       show(table(m0$doub))
    }
    if(i %in% c(2,5,6,9,10,12,13,23,16)) {
        m0$celltypes_ext[which(m0$doub == "T")] = "doublet"
        meta$celltype_final[mat] = m0$celltypes_ext
        
       cluster_centers <- m0[which(!is.na(m0$celltypes_ext)),] %>%
       group_by(celltypes_ext) %>%
       summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

       umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = celltypes_ext)) +
       geom_point() +
       theme_bw() +
      # Add labels to the center of each cluster
      geom_text_repel(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = celltypes_ext),
                      color = "black", position = position_nudge(y = 0.0), size = 2) + ggtitle(i)
      print(umap_plot)
    }
    else {
        m0$predicted_id_sep_ext[which(m0$doub == "T")] = "doublet"

        meta$celltype_final[mat] = m0$predicted_id_sep_ext
        
        cluster_centers <- m0[which(!is.na(m0$predicted_id_sep_ext)),] %>%
       group_by(predicted_id_sep_ext) %>%
       summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

       umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = predicted_id_sep_ext)) +
       geom_point() +
       theme_bw() +
      # Add labels to the center of each cluster
      geom_text_repel(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = predicted_id_sep_ext),
                      color = "black", position = position_nudge(y = 0.0), size = 2) + ggtitle(i)
      print(umap_plot)

    }
    meta$doub[mat] = m0$doub
}
dev.off()

dev.off()

  dev.off()

table(meta$celltype_final)



meta = meta[,c(4:ncol(meta))]


head(meta)

length(which(meta$leiden_3 == 74))

length(which(meta$celltype_final == "doublet"))

meta$celltype_final[which(meta$leiden_3 == 74)]="doublet"
meta$doub[which(meta$leiden_3 == 74)]="T"


write.table(meta, "final_meta.csv", sep = ",") 

length(unique(meta$celltype_final))


