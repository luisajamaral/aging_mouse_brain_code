library(data.table)
library(ggplot2)
library(Seurat)
library(pheatmap)
library(enrichR)
websiteLive <- getOption("enrichR.live")
dbs <- c( "GO_Biological_Process_2023")


library(zellkonverter)
h = readH5AD("../h5ads_final/celltype_batch_age_PMAT_RPM1.h5ad")

count_table <-h@assays@data$X

obj = readRDS("../female_RNA/RNA_final_SCT.RDS")
obj@meta.data$age = factor(obj@meta.data$age, levels = c("2mo", "9mo", "18mo"))
obj@meta.data$region = factor(obj@meta.data$region, levels = c("AMY","RLP","HCA","HCP","FC","ENT","CP","NAC"))
obj$celltype_age = paste(obj$celltype_final, obj$age )
obj$celltype_age_region = paste(obj$celltype_final, obj$age ,obj$region)
Idents(obj) = "celltype_age_region"
av = AverageExpression(obj)

plot_heatmap <- function(features, mat, enrich, num_groups = 3, main) {
    options(repr.plot.width=8, repr.plot.height=8)
    features = unique(features)
    tab = mat[features,]
    tab = tab[which(rowSums(tab)>0),]
    age= rep("02mo", ncol(tab))
    age[grep("9mo", colnames(tab))] = "09mo"
    age[grep("18mo", colnames(tab))] = "18mo"

    col_annot = data.frame(age = age)
    rownames(col_annot) = colnames(tab)
    
  options(repr.plot.width=8, repr.plot.height=8)
  tab = tab[,order(age)]
  ph_with_annotation <- pheatmap(tab, cutree_rows = num_groups, scale = "row",show_rownames = (nrow(tab)<=50),
                                 color = colorRampPalette(c("navy", "white", "red"))(40),
                                 cluster_rows = TRUE, main = paste(main , length(features)), cluster_cols = F,
                                  annotation_col = col_annot)
    if(nrow(tab)>50){
        pheatmap(tab[1:50,], cutree_rows = num_groups, scale = "row",show_rownames = T,
                                 color = colorRampPalette(c("navy", "white", "red"))(40),
                                 cluster_rows = TRUE, main = paste(main, "top 50"), cluster_cols = F,
                                 annotation_col = col_annot)
    }
    #options(repr.plot.width=8, repr.plot.height=4)
    
    enriched <- enrichr(rownames(tab), dbs)
    print(plotEnrich(enriched$GO_Biological_Process_2023) + ggtitle(paste("All")))
    
} 
    
    


generate_and_save_plots <- function(genes, mat, enrich, main) {
    library(ggplot2)
    library(pheatmap)
    library(grid)
    library(gridExtra)
    
    # Ensure all devices are closed
    graphics.off()
    
    # Prepare the data
    features <- unique(genes)
    if(length(features)<2){
        return("no genes")
    }
    if(ncol(mat)<2){
        return("no cols")
    }
    tab <- mat[features,]
    tab <- tab[which(rowSums(tab) > 0),]
    if(nrow(tab)<1){
        cat("No rows\n")
        return("")
    }
    age <- rep("02mo", ncol(tab))
    age[grep("9mo", colnames(tab))] <- "09mo"
    age[grep("18mo", colnames(tab))] <- "18mo"
    
    col_annot <- data.frame(age = age)
    rownames(col_annot) <- colnames(tab)
    tab <- tab[, order(age)]
    
    # Generate the heatmaps
    ph_with_annotation <- pheatmap(tab, cutree_rows = 1, scale = "row", show_rownames = (nrow(tab) <= 50),
                                   color = colorRampPalette(c("navy", "white", "red"))(40),
                                   cluster_rows = TRUE, main = paste(main, length(features)), cluster_cols = FALSE,
                                   annotation_col = col_annot)
    
    if (nrow(tab) > 45) {
        ph_top_50 <- pheatmap(tab[1:45,], cutree_rows = 1, scale = "row", show_rownames = TRUE,
                              color = colorRampPalette(c("navy", "white", "red"))(40),
                              cluster_rows = TRUE, main = paste(main, "top 45"), cluster_cols = FALSE,
                              annotation_col = col_annot)
    }
    
    # Save the heatmaps to PDF
    save_pheatmap_pdf(ph_with_annotation, paste(main, "_DE_DAR.pdf", sep = ""))
    if (nrow(tab) > 50) {
        save_pheatmap_pdf(ph_top_50, paste(main, "_DE_DAR_top50.pdf", sep = ""))
    }
    
    # Perform enrichment analysis and save the ggplot
    enriched <- enrichr(rownames(tab), dbs)
    plot_enrich <- plotEnrich(enriched$GO_Biological_Process_2023) + ggtitle(main)
    ggsave(filename = paste(main, "_DE_DAR_enrich.pdf", sep = ""), plot = plot_enrich)
}

# Example call to the function


files = list.files(".","filtered")
head(files)

cts = gsub("filtered.", "", files)
cts = gsub(".8wk.abc_score.csv", "", cts)
cts = gsub(".18mo.abc_score.csv", "", cts)
cts = gsub("_", " " , cts)
cts = unique(cts)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
    stopifnot(inherits(x, 'pheatmap'))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}

cts

meta = read.csv("../ABC/gene_meta.csv")
rownames(meta) = meta$`geneID.1`


cts

head(dar)

ct = "Oligo NN"

 oct = ct
    o2 = read.csv(paste("filtered.",gsub(" ", "_" , ct),".8wk.abc_score.csv",sep =""))
    o18 = read.csv(paste("filtered.",gsub(" ", "_" , ct),".18mo.abc_score.csv",sep =""))
  
    
    oct = gsub("L2", "L2-", oct)
    oct = gsub("L6bCT", "L6b-CT", oct)

    ct = gsub("L2", "L2-", ct)
    ct = gsub("L6bCT", "L6b-CT", ct)

    ct = gsub("D1 G", "D12 G", ct)

    rpm_ct = count_table[,grep(ct, colnames(count_table))]
    colnames(av$RNA) = gsub("/", "-", colnames(av$RNA))
    RNA = av$RNA[,grep(ct, colnames(av$RNA))]

    good_groups = names(table(Idents(obj))[which(table(Idents(obj))>200)])
    good_groups = gsub("/", "-", good_groups)
    RNA= RNA[,which(colnames(RNA)%in%good_groups)]
    colnames(RNA)

    #dar = read.csv(paste("../h5ads_final/combined_diff/diff_csvs/diff_peaks_",gsub(" ", "_", ct),"_2vs18_iter.csv", sep = ""))

    dar = read.csv(paste("../h5ads_final/combined_diff/diff_csvs/diff_peaks_",gsub(" ", "_", ct),"_2vs18_iter.csv", sep = ""))


    dar$log2.fold_change. = - dar$log2.fold_change.
    rownames(dar) = dar$feature.name

    nrow(dar[which(dar$adjusted.p.value<0.05),])
    
    deg = read.csv(paste("../female_RNA/DEG_results_latent_rep_mito_together/",gsub(" ", "_", ct),".csv", sep = ""))
    ct = gsub(" ", ".", ct)
    ct = gsub("-", ".", ct)

    rpm = fread(paste("../h5ads_final/celltype_age_RPM_files/",ct,"_RPM.txt", sep =""))
    rpm = as.data.frame(rpm)
    rownames(rpm) = rpm$V1
    rpm = rpm[,-c(1)]




    rownames(deg) = deg$X
    deg$avg_log2FC = -deg$avg_log2FC

    rownames(o2) = o2$X
    rownames(o18) = o18$X

    merged_df <-merge(o2, o18, by = "X", suffixes = c("_2mo", "_18mo"))
    merged_df$gene = sapply(strsplit(as.character(merged_df$X), "-"), `[`, 3)
    merged_df$loc = paste(sapply(strsplit(as.character(merged_df$X), "-"), `[`, 1), "-", sapply(strsplit(as.character(merged_df$X), "-"), `[`, 2), sep ="")
    merged_df$gene_name = meta[paste(merged_df$gene), "gene_name"]
    rownames(merged_df) = merged_df$X
    merged_df = merged_df[order(merged_df$abc_score_2mo, decreasing = T),]
    merged_df = merged_df[-which(duplicated(merged_df$loc)),]
    #merged_df = merged_df[-which(duplicated(merged_df$gene_name)),]

    merged_df$deg_logfc = deg[paste(merged_df$gene_name), "avg_log2FC"]

    down = deg[which( deg$avg_log2FC< -.15 ), ]
    down = down[1:min(nrow(down),2000),"X"]
    cat(length(down), "Down deg\n")
    up = deg[which( deg$avg_log2FC>.15), ]
    up = up[1:min(nrow(up),2000),"X"]
    cat(length(up), "Up deg\n")


    merged_df$DE = "no"
    merged_df[which(merged_df$gene_name%in% down), "DE"] = "down"
    merged_df[which(merged_df$gene_name%in% up), "DE"] = "up"

    merged_df$dar_logfc = dar[paste(merged_df$loc), "log2.fold_change."]

    down = dar[which( dar$log2.fold_change.< -.15 ), ]
    down = down[1:min(nrow(down),3000),"feature.name"]
    cat(length(down), "Down dar\n")

    up = dar[which( dar$log2.fold_change.>.15), ]
    up = up[1:min(nrow(up),3000),"feature.name"]
    cat(length(up), "Up dar\n")

    merged_df$DAR = "no"
    merged_df[which(merged_df$loc%in% down), "DAR"] = "down"
    merged_df[which(merged_df$loc%in% up), "DAR"] = "up"

    if(oct == "STR D1 Gaba"){
        oct = "STR D12 Gaba"
    }
    merged_df$abc_diff = merged_df$abc_score_18mo-merged_df$abc_score_2mo
    merged_df$rpm_2mo = rpm[paste(merged_df$loc), paste(oct,":2mo", sep = "")]
    merged_df$rpm_9mo = rpm[paste(merged_df$loc), paste(oct,":9mo", sep = "")]
    merged_df$rpm_18mo = rpm[paste(merged_df$loc), paste(oct,":18mo", sep = "")]
    merged_df$rpm_diff = merged_df$rpm_18mo-merged_df$rpm_2mo 
    merged_df$log_rpm_diff = sign(merged_df$rpm_diff) * log2(abs(merged_df$rpm_diff)+1e-10)

    head(merged_df)
    table(merged_df$DAR)
    table(merged_df$DE)


    tab = merged_df[which(merged_df$DAR%in%c("up", "down", "no") & merged_df$DE%in%c("up", "down", "no")),]

    tab = tab[which((tab$abc_score_2mo | tab$abc_score_18mo )>.05),]


    pdf(paste(ct, "_DE_DAR_barplot.pdf", sep = ""))

    options(repr.plot.width=5, repr.plot.height=4)
    g1 = ggplot(tab, aes( x=DAR, fill=DE)) +
    geom_bar(position = "fill", alpha=0.5)+
    #geom_density(alpha=0.6)+ xlim(c(-50,50))+ 
    scale_color_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
    scale_fill_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
    labs(title="abc diff",x="DAR", y = "Density")+

    theme_classic()

    options(repr.plot.width=5, repr.plot.height=4)
    g2 = ggplot(tab, aes( x=DE, fill=DAR)) +
    geom_bar(position = "fill", alpha=0.5)+
    #geom_density(alpha=0.6)+ xlim(c(-50,50))+ 
    scale_color_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
    scale_fill_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
    labs(title="abc diff",x="DEG", y = "Density")+

    theme_classic()
    print(g1)
    print(g2)
    while (!is.null(dev.list())) dev.off()



    tab = tab[order(abs(tab$dar_logfc) , abs(tab$deg_logfc), decreasing = T),]

    genes = tab$gene_name[which(tab$DAR%in%c("up") & tab$DE%in%"up")]
    mat = RNA
    #pdf(paste(ct, "_DE_DAR.pdf", sep = ""))
    #print(plot_heatmap(genes,mat, TRUE,1, paste(ct, "Up")))
    generate_and_save_plots(genes,mat,1, paste(ct, "Up", sep = "_"))

    while (!is.null(dev.list())) dev.off()

    write.table(genes, file = paste(ct, "_up_DE_DAR.txt", sep = ""), sep = "\n", quote = F, row.names = F,col.names = F)

    genes = tab$gene_name[which(tab$DAR%in%c("down")& tab$DE%in%"down")]
    mat = RNA

    #plot_heatmap(genes,mat, TRUE,1, paste(ct, "Down"))
    #while (!is.null(dev.list())) dev.off()
    generate_and_save_plots(genes,mat, 1, paste(ct, "Down", sep = "_"))
    
    
    downlocs = tab$loc[which(tab$DAR%in%c("down")& tab$DE%in%"down")]
    bed_df <- do.call(rbind, strsplit(downlocs, "[:-]"))
    bed_df <- as.data.frame(bed_df, stringsAsFactors = FALSE)
    head(bed_df)

    uplocs = tab$loc[which(tab$DAR%in%c("up")& tab$DE%in%"up")]
    bed_df <- do.call(rbind, strsplit(uplocs, "[:-]"))
    bed_df <- as.data.frame(bed_df, stringsAsFactors = FALSE)
    head(bed_df)



merged_df[grep("Hist1h1", merged_df$gene_name),]

#for(ct in cts[9:length(cts)]) {
for(ct in cts[6:14]) {
    
    oct = ct
    o2 = read.csv(paste("filtered.",gsub(" ", "_" , ct),".8wk.abc_score.csv",sep =""))
    o18 = read.csv(paste("filtered.",gsub(" ", "_" , ct),".18mo.abc_score.csv",sep =""))
  
    
    oct = gsub("L2", "L2-", oct)
    oct = gsub("L6bCT", "L6b-CT", oct)

    ct = gsub("L2", "L2-", ct)
    ct = gsub("L6bCT", "L6b-CT", ct)

    ct = gsub("D1 G", "D12 G", ct)

    rpm_ct = count_table[,grep(ct, colnames(count_table))]
    colnames(av$RNA) = gsub("/", "-", colnames(av$RNA))
    RNA = av$RNA[,grep(ct, colnames(av$RNA))]

    good_groups = names(table(Idents(obj))[which(table(Idents(obj))>200)])
    good_groups = gsub("/", "-", good_groups)
    RNA= RNA[,which(colnames(RNA)%in%good_groups)]
    colnames(RNA)

    #dar = read.csv(paste("../h5ads_final/combined_diff/diff_csvs/diff_peaks_",gsub(" ", "_", ct),"_2vs18_iter.csv", sep = ""))

    dar = read.csv(paste("../h5ads_final/combined_diff/diff_csvs/diff_peaks_",gsub(" ", "_", ct),"_2vs18_iter.csv", sep = ""))


    dar$log2.fold_change. = - dar$log2.fold_change.
    rownames(dar) = dar$feature.name

    nrow(dar[which(dar$adjusted.p.value<0.05),])
    
    deg = read.csv(paste("../female_RNA/DEG_results_latent_rep_mito_together/",gsub(" ", "_", ct),".csv", sep = ""))
    ct = gsub(" ", ".", ct)
    ct = gsub("-", ".", ct)

    rpm = fread(paste("../h5ads_final/celltype_age_RPM_files/",ct,"_RPM.txt", sep =""))
    rpm = as.data.frame(rpm)
    rownames(rpm) = rpm$V1
    rpm = rpm[,-c(1)]




    rownames(deg) = deg$X
    deg$avg_log2FC = -deg$avg_log2FC

    rownames(o2) = o2$X
    rownames(o18) = o18$X

    merged_df <-merge(o2, o18, by = "X", suffixes = c("_2mo", "_18mo"))
    merged_df$gene = sapply(strsplit(as.character(merged_df$X), "-"), `[`, 3)
    merged_df$loc = paste(sapply(strsplit(as.character(merged_df$X), "-"), `[`, 1), "-", sapply(strsplit(as.character(merged_df$X), "-"), `[`, 2), sep ="")
    merged_df$gene_name = meta[paste(merged_df$gene), "gene_name"]
    rownames(merged_df) = merged_df$X
    merged_df = merged_df[order(merged_df$abc_score_2mo, decreasing = T),]
    merged_df = merged_df[-which(duplicated(merged_df$loc)),]
    #merged_df = merged_df[-which(duplicated(merged_df$gene_name)),]

    merged_df$deg_logfc = deg[paste(merged_df$gene_name), "avg_log2FC"]

    down = deg[which( deg$avg_log2FC< -.15 & deg$p_val_adj<0.1), ]
    down = down[1:min(nrow(down),2000),"X"]
    cat(length(down), "Down deg\n")
    up = deg[which( deg$avg_log2FC>.15& deg$p_val_adj<0.1), ]
    up = up[1:min(nrow(up),2000),"X"]
    cat(length(up), "Up deg\n")


    merged_df$DE = "no"
    merged_df[which(merged_df$gene_name%in% down), "DE"] = "down"
    merged_df[which(merged_df$gene_name%in% up), "DE"] = "up"

    merged_df$dar_logfc = dar[paste(merged_df$loc), "log2.fold_change."]

    down = dar[which( dar$log2.fold_change.< -.15 & dar$adjusted.p.value<0.1), ]
    down = down[1:min(nrow(down),3000),"feature.name"]
    cat(length(down), "Down dar\n")

    up = dar[which( dar$log2.fold_change.>.15& dar$adjusted.p.value<0.1), ]
    up = up[1:min(nrow(up),3000),"feature.name"]
    cat(length(up), "Up dar\n")

    merged_df$DAR = "no"
    merged_df[which(merged_df$loc%in% down), "DAR"] = "down"
    merged_df[which(merged_df$loc%in% up), "DAR"] = "up"

    if(oct == "STR D1 Gaba"){
        oct = "STR D12 Gaba"
    }
    merged_df$abc_diff = merged_df$abc_score_18mo-merged_df$abc_score_2mo
    merged_df$rpm_2mo = rpm[paste(merged_df$loc), paste(oct,":2mo", sep = "")]
    merged_df$rpm_9mo = rpm[paste(merged_df$loc), paste(oct,":9mo", sep = "")]
    merged_df$rpm_18mo = rpm[paste(merged_df$loc), paste(oct,":18mo", sep = "")]
    merged_df$rpm_diff = merged_df$rpm_18mo-merged_df$rpm_2mo 
    merged_df$log_rpm_diff = sign(merged_df$rpm_diff) * log2(abs(merged_df$rpm_diff)+1e-10)

    head(merged_df)
    table(merged_df$DAR)
    table(merged_df$DE)


    tab = merged_df[which(merged_df$DAR%in%c("up", "down", "no") & merged_df$DE%in%c("up", "down", "no")),]

    tab = tab[which((tab$abc_score_2mo | tab$abc_score_18mo )>.05),]


    pdf(paste(ct, "_DE_DAR_barplot.pdf", sep = ""))

    options(repr.plot.width=5, repr.plot.height=4)
    g1 = ggplot(tab, aes( x=DAR, fill=DE)) +
    geom_bar(position = "fill", alpha=0.5)+
    #geom_density(alpha=0.6)+ xlim(c(-50,50))+ 
    scale_color_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
    scale_fill_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
    labs(title="abc diff",x="DAR", y = "Density")+

    theme_classic()

    options(repr.plot.width=5, repr.plot.height=4)
    g2 = ggplot(tab, aes( x=DE, fill=DAR)) +
    geom_bar(position = "fill", alpha=0.5)+
    #geom_density(alpha=0.6)+ xlim(c(-50,50))+ 
    scale_color_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
    scale_fill_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
    labs(title="abc diff",x="DEG", y = "Density")+

    theme_classic()
    print(g1)
    print(g2)
    while (!is.null(dev.list())) dev.off()



    tab = tab[order(abs(tab$dar_logfc) , abs(tab$deg_logfc), decreasing = T),]

    genes = tab$gene_name[which(tab$DAR%in%c("up") & tab$DE%in%"up")]
    mat = RNA
    #pdf(paste(ct, "_DE_DAR.pdf", sep = ""))
    #print(plot_heatmap(genes,mat, TRUE,1, paste(ct, "Up")))
    generate_and_save_plots(genes,mat,1, paste(ct, "Up", sep = "_"))

    while (!is.null(dev.list())) dev.off()

    write.table(genes, file = paste(ct, "_up_DE_DAR.txt", sep = ""), sep = "\n", quote = F, row.names = F,col.names = F)

    genes = tab$gene_name[which(tab$DAR%in%c("down")& tab$DE%in%"down")]
    mat = RNA

    #plot_heatmap(genes,mat, TRUE,1, paste(ct, "Down"))
    #while (!is.null(dev.list())) dev.off()
    generate_and_save_plots(genes,mat, 1, paste(ct, "Down", sep = "_"))
    write.table(genes, file = paste(ct, "_down_DE_DAR.txt", sep = ""), sep = "\n", quote = F, row.names = F, col.names = F)
    
    
    downlocs = tab$loc[which(tab$DAR%in%c("down")& tab$DE%in%"down")]
    bed_df <- do.call(rbind, strsplit(downlocs, "[:-]"))
    bed_df <- as.data.frame(bed_df, stringsAsFactors = FALSE)
    head(bed_df)
    write.table(bed_df,  file = paste(ct, "_down_locs.bed", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

    uplocs = tab$loc[which(tab$DAR%in%c("up")& tab$DE%in%"up")]
    bed_df <- do.call(rbind, strsplit(uplocs, "[:-]"))
    bed_df <- as.data.frame(bed_df, stringsAsFactors = FALSE)
    head(bed_df)
    write.table(bed_df,  file = paste(ct, "_up_locs.bed", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

    
    
}

ct
oct


ct = 'Lamp5 Gaba'
RNA = av$RNA[,grep(ct, colnames(av$RNA))]
head(RNA)

head(mat)
genes
generate_and_save_plots(genes,mat,1, paste(ct, "Up", sep = "_"))



  good_groups = names(table(Idents(obj))[which(table(Idents(obj))>300)])
    good_groups = gsub("/", "-", good_groups)

good_groups
   # RNA= RNA[,which(colnames(RNA)%in%good_groups)]
 

'../h5ads_final/combined_diff/diff_csvs/diff_peaks_L6b-CT_ENT_Glut_2vs18_iter.csv

'../female_RNA/DEG_results_latent_rep_mito_together/Astro-NT_NN.csv

downlocs = tab$loc[which(tab$DAR%in%c("down")& tab$DE%in%"down")]
bed_df <- do.call(rbind, strsplit(downlocs, "[:-]"))
bed_df <- as.data.frame(bed_df, stringsAsFactors = FALSE)
head(bed_df)
write.table(bed_df, file = "output.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


'../h5ads_final/combined_diff/diff_csvs/diff_peaks_L2-3_IT_CTX_Glut_2vs18_iter.csv
