library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(edgeR)
library(enrichR)
meta = read.csv("~/projects/combined_all/female_RNA/meta_final.csv")
setwd("~/projects/fmouse_multiome/SoloTE_out/")

for (reg in c("HCA", "FC", "NAC","CP", "HCP")){
    setwd("~/projects/fmouse_multiome/SoloTE_out/")
    all_directories <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
    data_directories <- all_directories[grep("locus", all_directories)]
    data_directories = data_directories[grep(reg, data_directories)]
    if(reg == "CP"){
        data_directories = data_directories[-grep("HCP", data_directories)]
    }
    # Initialize an empty Seurat object to store the merged data
    directory = data_directories[1]
    setwd(directory)
    sample = strsplit(directory, "/", )[[1]][2]
    solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
    merged_seuratobj <- CreateSeuratObject(count = solote_matrix, min.features = 100, project = sample)
    merged_seuratobj <- CreateSeuratObject(count = solote_matrix, min.features = 100,project = sample)
    merged_seuratobj$sample = gsub("_SoloTE_output","",merged_seuratobj$orig.ident)
    merged_seuratobj$barcode = paste(merged_seuratobj$sample, sapply(strsplit(rownames(merged_seuratobj@meta.data), "-"), "[[", 1), sep  = ":")
    merged_seuratobj$barcode = gsub("02mo","2mo",merged_seuratobj$barcode)
    merged_seuratobj$barcode = gsub("09mo","9mo",merged_seuratobj$barcode)
    mat = match(merged_seuratobj$barcode, meta$barcode)
    merged_seuratobj$celltype = meta$celltype_final[mat]
    merged_seuratobj = merged_seuratobj[,which(!is.na(merged_seuratobj$celltype))]


    # Loop through each directory and load the data
    for (directory in data_directories[-1]) {
      setwd("~/projects/fmouse_multiome/SoloTE_out/")
      setwd(directory)
      sample = strsplit(directory, "/", )[[1]][2]
      cat(sample, "\n")
      solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
      if(dim(solote_matrix)[2]<25) {
          next
      }
      # Create a temporary Seurat object for the current directory
      temp_seuratobj <- CreateSeuratObject(count = solote_matrix, min.features = 100,project = sample)
      temp_seuratobj$sample = gsub("_SoloTE_output","",temp_seuratobj$orig.ident)
      temp_seuratobj$barcode = paste(temp_seuratobj$sample, sapply(strsplit(rownames(temp_seuratobj@meta.data), "-"), "[[", 1), sep  = ":")
      temp_seuratobj$barcode = gsub("02mo","2mo",temp_seuratobj$barcode)
      temp_seuratobj$barcode = gsub("09mo","9mo",temp_seuratobj$barcode)
      mat = match(temp_seuratobj$barcode, meta$barcode)
      temp_seuratobj$celltype = meta$celltype_final[mat]
      temp_seuratobj = temp_seuratobj[,which(!is.na(temp_seuratobj$celltype))]

      # Merge the current Seurat object with the merged_seuratobj
      merged_seuratobj <- merge(x=merged_seuratobj, y=temp_seuratobj)
        
    }
    merged_seuratobj$sample = gsub("_SoloTE_output","",merged_seuratobj$orig.ident)
    seuratobj = merged_seuratobj
    head(merged_seuratobj@meta.data)
    seuratobj=JoinLayers(seuratobj)
    seuratobj$barcode = paste(seuratobj$sample, sapply(strsplit(rownames(seuratobj@meta.data), "-"), "[[", 1), sep  = ":")
    seuratobj$barcode = gsub("02mo","2mo",seuratobj$barcode)
    seuratobj$barcode = gsub("09mo","9mo",seuratobj$barcode)
    mat = match(seuratobj$barcode, meta$barcode)
    seuratobj$celltype = meta$celltype_final[mat]
    seuratobj$age = sapply(strsplit(as.character(seuratobj$sample), "_"), "[[", 2)
    seuratobj = seuratobj[,which(!is.na(seuratobj$celltype))]
    rs <- rowSums(seuratobj@assays$RNA$counts)
    seuratobj = seuratobj[which(rs>=4),]
    seuratobj <- NormalizeData(seuratobj)
    seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(seuratobj)
    all.genes = all.genes[-grep("SoloTE", all.genes)]
    seuratobj <- ScaleData(seuratobj, features = all.genes)
    seuratobj <- RunPCA(seuratobj, features = VariableFeatures(object = seuratobj))
    seuratobj <- FindNeighbors(seuratobj, dims = 1:45)
    seuratobj <- FindClusters(seuratobj, resolution = 1)
    seuratobj <- RunUMAP(seuratobj, dims = 1:45)
    #DimPlot(seuratobj,reduction="umap")
    setwd("~/projects/combined_all/female_RNA/SoloTE/")

    pdf(paste(reg,"_UMAP.pdf", sep = ""), height = 5, width = 15)
    print(DimPlot(seuratobj,reduction="umap", group.by = "celltype", label = T))
    print(DimPlot(seuratobj,reduction="umap", group.by = "age", label = T))
    print(DimPlot(seuratobj,reduction="umap", group.by = "sample", label = T))
    dev.off()
    #saveRDS(seuratobj, paste('~/projects/combined_all/female_RNA/SoloTE/', reg, "_locus.RDS", sep = ""))
    ps <- AggregateExpression(seuratobj, assays = "RNA", return.seurat = T, group.by = c("celltype", "sample"))
    ps$age = sapply(strsplit(as.character(ps$sample), "-"), "[[", 2)
    ps$celltype_age = paste(ps$celltype,ps$age, sep = "_")
    saveRDS(ps, paste('~/projects/combined_all/female_RNA/SoloTE/', reg, "_locus_ps.RDS", sep = ""))
    counts = ps@assays$RNA$counts
    dir.create((paste(reg,"_locus", sep = "")))
    setwd((paste(reg,"_locus", sep = "")))
    for(ct in names(table(seuratobj$celltype)[which(table(seuratobj$celltype)>500)])){
        cat(ct, "\n")
        cl = gsub("/", "-", ct)
        cl = gsub(" ", "-", cl)
        cl = paste(reg, cl,sep = "_")
        colnames(counts)[grep(ct, colnames(counts))]
        counts_L23 = counts[,grep(ct, colnames(counts))]
        groups = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
        groups = sapply(strsplit(as.character(groups), "-"), "[[", 2)
        regs = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
        regs = sapply(strsplit(as.character(regs), "-"), "[[", 1)


        #hist(rowSums(counts_L23))
        min_reads = 50000
        colSums(counts_L23)
        
        bin.cov = log10(Matrix::rowSums(counts_L23)+1)
        #hist(bin.cov)
        bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.99)
        idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
        #length(idy)
        counts_L23=counts_L23[idy,]
        groups = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
        groups = sapply(strsplit(as.character(groups), "-"), "[[", 2)
        groups = factor(groups, levels = c("18mo","09mo", "02mo"))
        regs = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
        regs = sapply(strsplit(as.character(regs), "-"), "[[", 1)

        design <- model.matrix(~groups)
        head(design)
        # Create a DGEList object for the subset
        d <- DGEList(counts = counts_L23, group = groups)

        # Filter genes based on counts in the current subset
        cpm = cpm(d)


        keep = (rowSums(cpm>1)>=2)
        table(keep)
        d = d[keep,]
        # Normalize the filtered subset
        d_normalized <- calcNormFactors(d)
        d_disp <- estimateDisp(d_normalized)


        #plotMDS(d, gene.selection = "common", labels = groups)
        #plotMDS(d, gene.selection = "common", labels = regs)

        # Perform ANOVA analysis
        fit <- glmFit(d_disp, design)

        fit_contrasts <- glmLRT(fit, coef = 3)
        fit_contrasts$table$FDR <- p.adjust(fit_contrasts$table$PValue, method = "BH")  # Adjust p-values for FDR

        fit_contrasts$table = fit_contrasts$table[order(fit_contrasts$table$PValue),]
        head(fit_contrasts$table )
        sig = fit_contrasts$table[which(fit_contrasts$table$PValue<0.05),]
        nrow(sig)
        head(sig[grep("Solo", rownames(sig)),])
        if(nrow(sig)<2){
            next
        }
        write.table(sig, paste(cl, "_DEG_TE.txt",sep = ""), sep = "\t")

        dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "KEGG_2021_Human", "SynGO_2022")

        rn = rownames(sig[which(sig$logFC>0),])
        rn = sapply(strsplit(as.character(rn), ";"), tail, 1)

        enriched_down <- enrichr(rn, dbs)


        rn = rownames(sig[which(sig$logFC<0),])
        rn = sapply(strsplit(as.character(rn), ";"), tail, 1)

        enriched_up <- enrichr(rn, dbs)

        head(enriched_up[["GO_Biological_Process_2023"]],5)
        head(enriched_down[["GO_Biological_Process_2023"]],5)

        up = enriched_up[["GO_Biological_Process_2023"]]
        write.table(up[which(up$P.value<0.05),], paste(cl, "_GO_BP_up.txt",sep = ""), sep = "\t")


        down = enriched_down[["GO_Biological_Process_2023"]]
        write.table(down[which(down$P.value<0.05),], paste(cl, "_GO_BP_down.txt",sep = ""), sep = "\t")


        up = enriched_up[["SynGO_2022"]]
        write.table(up[which(up$P.value<0.05),], paste(cl, "_SynGO_up.txt",sep = ""), sep = "\t")


        down = enriched_down[["SynGO_2022"]]
        write.table(down[which(down$P.value<0.05),], paste(cl, "_SynGO_down.txt",sep = ""), sep = "\t")


        th  = cpm[rownames(sig),]
        type= rep("gene", nrow(th))
        type[grep("Solo", rownames(th))] = "TE"
        annotation_row = data.frame(
            type = type
          )
        rownames(annotation_row) = rownames(th)

        annotation_col= data.frame(
            age = groups,
            region = regs
          )
        rownames(annotation_col) = colnames(th)

        options(repr.plot.width=15, repr.plot.height=12)
        g=c()
        for(i in 2:length(groups)){
            g = c(g, groups[i]!=groups[i-1])
        }

        groups = factor(groups, levels = c("02mo","09mo", "18mo"))

        tt = t(scale(t(th[which(type == "gene"),order(groups)])))
        quantile_breaks <- function(xs, n = 10) {
          breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
          breaks[!duplicated(breaks)]
        }

        mat_breaks <- quantile_breaks(tt, n = 50)
        if(nrow(tt)>2){
        ph = pheatmap(tt, show_rownames = F,gaps_col = which(g),scale = "none", cluster_cols = F,annotation_row= annotation_row,annotation_col=annotation_col,color   = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
          breaks = mat_breaks)
        pdf(paste(cl, "_gene_heat.pdf",sep = ""), height =7, width =7)
        print(ph)
        dev.off()
        }
        options(repr.plot.width=12, repr.plot.height=12)

        tt = t(scale(t(th[which(type == "TE"),order(groups)])))
        quantile_breaks <- function(xs, n = 10) {
          breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
          breaks[!duplicated(breaks)]
        }

        mat_breaks <- quantile_breaks(tt, n = 50)
        sho=TRUE
        if(nrow(tt)>50){
            sho = FALSE
        }
        if(nrow(tt)>2){
        ph = pheatmap(tt, show_rownames = sho,gaps_col = which(g),scale = "none", cluster_cols = F,annotation_row= annotation_row,annotation_col=annotation_col,color   = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
          breaks= mat_breaks)
        pdf(paste(cl, "_SoloTE_heat.pdf",sep = ""), height =7, width =7)
        print(ph)
        dev.off()
        }
        options(repr.plot.width=8, repr.plot.height=7)

        markers = fit_contrasts$table
        markers$gene = rownames(markers)
        markers$type = ifelse(grepl("SoloTE", markers$gene), "SoloTE", "Gene")
        volcano_plot <- ggplot(markers, aes(x = -logFC, y = -log10(PValue+1e-1000))) +
          geom_point(aes(color = type), size = 2) +
          scale_color_manual(values = c("Gene" = "black", "SoloTE" = "red")) +
          labs(title = paste(ct, "Volcano Plot\n 2mo vs 18mo"), x = "Average Log2 Fold Change", y = "-log10(p-value)") +
          theme_minimal()+
          geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "blue")  # Add horizontal dotted line at y = -log10(0.01)


        top_significant_points <- markers[order(markers$PValue), ][1:25, ]

        # Add text labels to the top 15 significant points
        pdf(paste(cl, "_volcano.pdf",sep = ""), height =6, width =6)
        v1 = volcano_plot +
          geom_text_repel(data = top_significant_points, aes(x = -logFC, y = -log10(PValue), label = gene), 
                    vjust = -0.5, hjust = 0.5, size = 3, color = "black")+
          geom_point(data = subset(markers, grepl("SoloTE", gene)),
                     aes(x = -logFC, y = -log10(PValue + 1e-1000)), color = "red", size = 2)


        print(v1)
        dev.off()
    }

    
}

for (reg in c("HCA", "FC", "NAC","CP", "HCP", "ENT", "AMY", "RLP")){
    setwd("~/projects/fmouse_multiome/SoloTE_out/")
    all_directories <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
    data_directories <- all_directories[grep("subfamily", all_directories)]
    data_directories = data_directories[grep(reg, data_directories)]
    if(reg == "CP"){
        data_directories = data_directories[-grep("HCP", data_directories)]
    }
    # Initialize an empty Seurat object to store the merged data
    directory = data_directories[1]
    setwd(directory)
    sample = strsplit(directory, "/", )[[1]][2]
    solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
    merged_seuratobj <- CreateSeuratObject(count = solote_matrix, min.features = 100, project = sample)
    merged_seuratobj <- CreateSeuratObject(count = solote_matrix, min.features = 100,project = sample)
    merged_seuratobj$sample = gsub("_SoloTE_output","",merged_seuratobj$orig.ident)
    merged_seuratobj$barcode = paste(merged_seuratobj$sample, sapply(strsplit(rownames(merged_seuratobj@meta.data), "-"), "[[", 1), sep  = ":")
    merged_seuratobj$barcode = gsub("02mo","2mo",merged_seuratobj$barcode)
    merged_seuratobj$barcode = gsub("09mo","9mo",merged_seuratobj$barcode)
    mat = match(merged_seuratobj$barcode, meta$barcode)
    merged_seuratobj$celltype = meta$celltype_final[mat]
    merged_seuratobj = merged_seuratobj[,which(!is.na(merged_seuratobj$celltype))]


    # Loop through each directory and load the data
    for (directory in data_directories[-1]) {
      setwd("~/projects/fmouse_multiome/SoloTE_out/")
      setwd(directory)
      sample = strsplit(directory, "/", )[[1]][2]
      cat(sample, "\n")
      solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
      if(dim(solote_matrix)[2]<25) {
          next
      }
      # Create a temporary Seurat object for the current directory
      temp_seuratobj <- CreateSeuratObject(count = solote_matrix, min.features = 100,project = sample)
      temp_seuratobj$sample = gsub("_SoloTE_output","",temp_seuratobj$orig.ident)
      temp_seuratobj$barcode = paste(temp_seuratobj$sample, sapply(strsplit(rownames(temp_seuratobj@meta.data), "-"), "[[", 1), sep  = ":")
      temp_seuratobj$barcode = gsub("02mo","2mo",temp_seuratobj$barcode)
      temp_seuratobj$barcode = gsub("09mo","9mo",temp_seuratobj$barcode)
      mat = match(temp_seuratobj$barcode, meta$barcode)
      temp_seuratobj$celltype = meta$celltype_final[mat]
      temp_seuratobj = temp_seuratobj[,which(!is.na(temp_seuratobj$celltype))]

      # Merge the current Seurat object with the merged_seuratobj
      merged_seuratobj <- merge(x=merged_seuratobj, y=temp_seuratobj)
        
    }
    merged_seuratobj$sample = gsub("_SoloTE_output","",merged_seuratobj$orig.ident)
    seuratobj = merged_seuratobj
    head(merged_seuratobj@meta.data)
    seuratobj=JoinLayers(seuratobj)
    seuratobj$barcode = paste(seuratobj$sample, sapply(strsplit(rownames(seuratobj@meta.data), "-"), "[[", 1), sep  = ":")
    seuratobj$barcode = gsub("02mo","2mo",seuratobj$barcode)
    seuratobj$barcode = gsub("09mo","9mo",seuratobj$barcode)
    mat = match(seuratobj$barcode, meta$barcode)
    seuratobj$celltype = meta$celltype_final[mat]
    seuratobj$age = sapply(strsplit(as.character(seuratobj$sample), "_"), "[[", 2)
    seuratobj = seuratobj[,which(!is.na(seuratobj$celltype))]
    rs <- rowSums(seuratobj@assays$RNA$counts)
    seuratobj = seuratobj[which(rs>=4),]
    seuratobj <- NormalizeData(seuratobj)
    seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(seuratobj)
    all.genes = all.genes[-grep("SoloTE", all.genes)]
    seuratobj <- ScaleData(seuratobj, features = all.genes)
    seuratobj <- RunPCA(seuratobj, features = VariableFeatures(object = seuratobj))
    seuratobj <- FindNeighbors(seuratobj, dims = 1:45)
    seuratobj <- FindClusters(seuratobj, resolution = 1)
    seuratobj <- RunUMAP(seuratobj, dims = 1:45)
    #DimPlot(seuratobj,reduction="umap")
    setwd("~/projects/combined_all/female_RNA/SoloTE/")

    pdf(paste(reg,"_UMAP.pdf", sep = ""), height = 5, width = 15)
    print(DimPlot(seuratobj,reduction="umap", group.by = "celltype", label = T))
    print(DimPlot(seuratobj,reduction="umap", group.by = "age", label = T))
    print(DimPlot(seuratobj,reduction="umap", group.by = "sample", label = T))
    dev.off()
    #saveRDS(seuratobj, paste('~/projects/combined_all/female_RNA/SoloTE/', reg, "_locus.RDS", sep = ""))
    ps <- AggregateExpression(seuratobj, assays = "RNA", return.seurat = T, group.by = c("celltype", "sample"))
    ps$age = sapply(strsplit(as.character(ps$sample), "-"), "[[", 2)
    ps$celltype_age = paste(ps$celltype,ps$age, sep = "_")
    saveRDS(ps, paste('~/projects/combined_all/female_RNA/SoloTE/', reg, "_subfamily_ps.RDS", sep = ""))
    counts = ps@assays$RNA$counts
    dir.create((paste(reg,"_locus", sep = "")))
    setwd((paste(reg,"_locus", sep = "")))
    for(ct in names(table(seuratobj$celltype)[which(table(seuratobj$celltype)>500)])){
        cat(ct, "\n")
        cl = gsub("/", "-", ct)
        cl = gsub(" ", "-", cl)
        cl = paste(reg, cl,'subfamily',sep = "_")
        colnames(counts)[grep(ct, colnames(counts))]
        counts_L23 = counts[,grep(ct, colnames(counts))]
        groups = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
        groups = sapply(strsplit(as.character(groups), "-"), "[[", 2)
        regs = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
        regs = sapply(strsplit(as.character(regs), "-"), "[[", 1)


        #hist(rowSums(counts_L23))
        min_reads = 50000
        colSums(counts_L23)
        
        bin.cov = log10(Matrix::rowSums(counts_L23)+1)
        #hist(bin.cov)
        bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.99)
        idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
        #length(idy)
        counts_L23=counts_L23[idy,]
        groups = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
        groups = sapply(strsplit(as.character(groups), "-"), "[[", 2)
        groups = factor(groups, levels = c("18mo","09mo", "02mo"))
        regs = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
        regs = sapply(strsplit(as.character(regs), "-"), "[[", 1)

        design <- model.matrix(~groups)
        head(design)
        # Create a DGEList object for the subset
        d <- DGEList(counts = counts_L23, group = groups)

        # Filter genes based on counts in the current subset
        cpm = cpm(d)


        keep = (rowSums(cpm>1)>=2)
        table(keep)
        d = d[keep,]
        # Normalize the filtered subset
        d_normalized <- calcNormFactors(d)
        d_disp <- estimateDisp(d_normalized)


        #plotMDS(d, gene.selection = "common", labels = groups)
        #plotMDS(d, gene.selection = "common", labels = regs)

        # Perform ANOVA analysis
        fit <- glmFit(d_disp, design)

        fit_contrasts <- glmLRT(fit, coef = 3)
        fit_contrasts$table$FDR <- p.adjust(fit_contrasts$table$PValue, method = "BH")  # Adjust p-values for FDR

        fit_contrasts$table = fit_contrasts$table[order(fit_contrasts$table$PValue),]
        head(fit_contrasts$table )
        sig = fit_contrasts$table[which(fit_contrasts$table$PValue<0.05),]
        nrow(sig)
        head(sig[grep("Solo", rownames(sig)),])
        if(nrow(sig)<2){
            next
        }
        write.table(sig, paste(cl, "_DEG_TE.txt",sep = ""), sep = "\t")

        dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "KEGG_2021_Human", "SynGO_2022")

        rn = rownames(sig[which(sig$logFC>0),])
        rn = sapply(strsplit(as.character(rn), ";"), tail, 1)

        enriched_down <- enrichr(rn, dbs)


        rn = rownames(sig[which(sig$logFC<0),])
        rn = sapply(strsplit(as.character(rn), ";"), tail, 1)

        enriched_up <- enrichr(rn, dbs)

        head(enriched_up[["GO_Biological_Process_2023"]],5)
        head(enriched_down[["GO_Biological_Process_2023"]],5)

        up = enriched_up[["GO_Biological_Process_2023"]]
        write.table(up[which(up$P.value<0.05),], paste(cl, "_GO_BP_up.txt",sep = ""), sep = "\t")


        down = enriched_down[["GO_Biological_Process_2023"]]
        write.table(down[which(down$P.value<0.05),], paste(cl, "_GO_BP_down.txt",sep = ""), sep = "\t")


        up = enriched_up[["SynGO_2022"]]
        write.table(up[which(up$P.value<0.05),], paste(cl, "_SynGO_up.txt",sep = ""), sep = "\t")


        down = enriched_down[["SynGO_2022"]]
        write.table(down[which(down$P.value<0.05),], paste(cl, "_SynGO_down.txt",sep = ""), sep = "\t")


        th  = cpm[rownames(sig),]
        type= rep("gene", nrow(th))
        type[grep("Solo", rownames(th))] = "TE"
        annotation_row = data.frame(
            type = type
          )
        rownames(annotation_row) = rownames(th)

        annotation_col= data.frame(
            age = groups,
            region = regs
          )
        rownames(annotation_col) = colnames(th)

        options(repr.plot.width=15, repr.plot.height=12)
        g=c()
        for(i in 2:length(groups)){
            g = c(g, groups[i]!=groups[i-1])
        }

        groups = factor(groups, levels = c("02mo","09mo", "18mo"))

        tt = t(scale(t(th[which(type == "gene"),order(groups)])))
        quantile_breaks <- function(xs, n = 10) {
          breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
          breaks[!duplicated(breaks)]
        }

        mat_breaks <- quantile_breaks(tt, n = 50)
        if(nrow(tt)>2){
        ph = pheatmap(tt, show_rownames = F,gaps_col = which(g),scale = "none", cluster_cols = F,annotation_row= annotation_row,annotation_col=annotation_col,color   = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
          breaks = mat_breaks)
        pdf(paste(cl, "_gene_heat.pdf",sep = ""), height =7, width =7)
        print(ph)
        dev.off()
        }
        options(repr.plot.width=12, repr.plot.height=12)

        tt = t(scale(t(th[which(type == "TE"),order(groups)])))
        quantile_breaks <- function(xs, n = 10) {
          breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
          breaks[!duplicated(breaks)]
        }

        mat_breaks <- quantile_breaks(tt, n = 50)
        sho=TRUE
        if(nrow(tt)>50){
            sho = FALSE
        }
        if(nrow(tt)>2){
        ph = pheatmap(tt, show_rownames = sho,gaps_col = which(g),scale = "none", cluster_cols = F,annotation_row= annotation_row,annotation_col=annotation_col,color   = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
          breaks= mat_breaks)
        pdf(paste(cl, "_SoloTE_heat.pdf",sep = ""), height =7, width =7)
        print(ph)
        dev.off()
        }
        options(repr.plot.width=8, repr.plot.height=7)

        markers = fit_contrasts$table
        markers$gene = rownames(markers)
        markers$type = ifelse(grepl("SoloTE", markers$gene), "SoloTE", "Gene")
        volcano_plot <- ggplot(markers, aes(x = -logFC, y = -log10(PValue+1e-1000))) +
          geom_point(aes(color = type), size = 2) +
          scale_color_manual(values = c("Gene" = "black", "SoloTE" = "red")) +
          labs(title = paste(ct, "Volcano Plot\n 2mo vs 18mo"), x = "Average Log2 Fold Change", y = "-log10(p-value)") +
          theme_minimal()+
          geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "blue")  # Add horizontal dotted line at y = -log10(0.01)


        top_significant_points <- markers[order(markers$PValue), ][1:25, ]

        # Add text labels to the top 15 significant points
        pdf(paste(cl, "_volcano.pdf",sep = ""), height =6, width =6)
        v1 = volcano_plot +
          geom_text_repel(data = top_significant_points, aes(x = -logFC, y = -log10(PValue), label = gene), 
                    vjust = -0.5, hjust = 0.5, size = 3, color = "black")+
          geom_point(data = subset(markers, grepl("SoloTE", gene)),
                     aes(x = -logFC, y = -log10(PValue + 1e-1000)), color = "red", size = 2)


        print(v1)
        dev.off()
    }

    
}

tt[1:10,1:6]

table(sapply(strsplit(rownames(tt), "-"), "[[", 5))

type = sapply(strsplit(rownames(tt), "-"), "[[", 5)
type = sapply(strsplit(as.character(type), ":"), tail, 1)

type = rep("other",nrow(tt))
type[grep("LINE", rownames(tt))] = "LINE"
type[grep("SINE", rownames(tt))] = "SINE"



annotation_row = data.frame(
    type = type
)
rownames(annotation_row) = rownames(tt)



ph = pheatmap(tt, show_rownames = sho,gaps_col = which(g),scale = "none", cluster_cols = F,annotation_row= annotation_row,annotation_col=annotation_col,color   = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
          breaks= mat_breaks)



while(3>1){dev.off()}


