library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(edgeR)
library(enrichR)
meta = read.csv("~/projects/combined_all/female_RNA/meta_final.csv")
setwd("~/projects/combined_all/female_RNA/")

obj = readRDS("RNA_final.RDS")



for(reg in c("ENT", "HCA", "HCP", "NAC", "AMY", "CP", "RLP")){
    setwd("~/projects/combined_all/female_RNA/")
    seuratobj = subset(obj, cells = which(obj$region==reg))
    ps <- AggregateExpression(seuratobj, assays = "RNA", return.seurat = T, group.by = c("celltype", "orig.ident"))
    ps$age = sapply(strsplit(as.character(ps$orig.ident), "-"), "[[", 2)
    ps$celltype_age = paste(ps$celltype,ps$age, sep = "_")
    counts = ps@assays$RNA$counts
    dir.create((paste(reg,"_region_edgeR", sep = "")))
    setwd((paste(reg,"_region_edgeR", sep = "")))
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


    good_cts= names(table(seuratobj$celltype)[which(table(seuratobj$celltype)>500)])
    for(ct in good_cts){
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
        groups = factor(groups, levels = c("18mo","9mo", "2mo"))

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
        if(nrow(sig)<2){
            next
        }
        write.table(sig, paste(cl, "_DEG.txt",sep = ""), sep = "\t")

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

        

        up = enriched_up[["KEGG_2021_Human"]]
        write.table(up[which(up$P.value<0.05),], paste(cl, "_KEGG_2021_Human_up.txt",sep = ""), sep = "\t")


        down = enriched_down[["KEGG_2021_Human"]]
        write.table(down[which(down$P.value<0.05),], paste(cl, "_KEGG_2021_Human_down.txt",sep = ""), sep = "\t")

        
        th  = cpm[rownames(sig),]
        type= rep("gene", nrow(th))
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

        groups = factor(groups, levels = c("2mo","9mo", "18mo"))

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

        options(repr.plot.width=8, repr.plot.height=7)

        markers = fit_contrasts$table
        markers$gene = rownames(markers)
        markers$sig = "no"
        markers$sig[which(sig$PValue<0.01)] = "yes"
        volcano_plot <- ggplot(markers, aes(x = -logFC, y = -log10(PValue+1e-1000))) +
          geom_point(aes(color = sig), size = 2) + 
          scale_color_manual(values = c("no" = "black", "yes" = "red")) +
          labs(title = paste(ct, "Volcano Plot\n 2mo vs 18mo"), x = "Average Log2 Fold Change", y = "-log10(p-value)") +
          theme_minimal()+
          geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "blue")  # Add horizontal dotted line at y = -log10(0.01)


        top_significant_points <- markers[order(markers$PValue), ][1:25, ]

        # Add text labels to the top 15 significant points
        v1 = volcano_plot +
          geom_text_repel(data = top_significant_points, aes(x = -logFC, y = -log10(PValue), label = gene), 
                    vjust = -0.5, hjust = 0.5, size = 3, color = "red")

        pdf(paste(cl, "_volcano.pdf",sep = ""), height =6, width =6)
        print(v1)
        dev.off()
    }
}

enriched_down <- enrichr(rn, dbs)


while(3>1){dev.off()}

head(tab)

    head(cur)


reg = "NAC_"

setwd(paste("~/projects/combined_all/female_RNA/",reg,"region_edgeR/", sep = ""))
files = list.files(".", "DEG.txt", full.names = F)


curr = read.table(files[1])

curr$celltype = gsub(reg,"", files[1])
curr$celltype = gsub("_DEG.txt", "", curr$celltype)
curr$gene = rownames(curr)
curr$region = gsub("_", "", reg)
head(curr)

tab = curr 
for(i in 2:length(files)){
    cur = read.table(files[i])

    cur$celltype = gsub("_DEG.txt", "", files[i])
    cur$celltype = gsub(reg, "", cur$celltype)

    cur$gene = rownames(cur)
    cur$region = gsub("_", "", reg)

    tab = rbind(tab, cur)
        
}
print("OUT")
for(reg in c("FC_", "ENT_", "HCA_", "HCP_", "CP_", "AMY_", "RLP_")) {
setwd(paste("~/projects/combined_all/female_RNA/",reg,"region_edgeR/", sep = ""))
files = list.files(".", "DEG.txt", full.names = F)
print(reg)

for(i in 1:length(files)){
    cur = read.table(files[i])

    cur$celltype = gsub("_DEG.txt", "", files[i])
    cur$celltype = gsub(reg, "", cur$celltype)

    cur$gene = rownames(cur)
    cur$region = gsub("_", "", reg)
    head(cur)
    tab = rbind(tab, cur)
        
}
}
nrow(tab)
length(which(abs(tab$PValue)<0.01))

tab$direction = "Up in aging"
tab$direction[which(tab$logFC>0)] = "Down in aging"
tab$clade="Glut"
tab$clade[grep("NN", tab$celltype)]="NN"
tab$clade[grep("Gaba", tab$celltype)]="Gaba"
tab$clade[grep("Neur", tab$celltype)]="Gaba"
tab$clade[grep("Ex-IMN", tab$celltype)]="Glut"
tab$clade[grep("Inh-IMN", tab$celltype)]="Gaba"

ord = names(table(meta$celltype_final)[order(table(meta$celltype_final),decreasing = T)])
ord = gsub("/","-", ord)
ord = gsub(" ","-", ord)
reg = gsub("_","", reg)

options(repr.plot.width=17, repr.plot.height=4)


d = tab
d = d[which(d$PValue<0.01),]
d$celltype = factor(d$celltype, levels = rev(ord))


g = ggplot(d, aes(x = celltype, ,fill = direction)) +
  geom_bar(stat = "count",position = "dodge", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~region,scales = "free_y",ncol = 3)+ ggtitle(paste(reg, "P < 0.05 All DEGs"))+
  theme(text = element_text(size = 18))
    print(g)


options(repr.plot.width=25, repr.plot.height=15)

d = tab
d = d[which(d$FDR<0.05),]
d$celltype = factor(d$celltype, levels = rev(ord))


g = ggplot(d, aes(x = celltype, ,fill = direction)) +
  geom_bar(stat = "count",position = "dodge", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~region,ncol = 8)+ ggtitle(paste("FDR < 0.05 All DEGs"))+
  theme(text = element_text(size = 18))
    print(g)


tab[grep("Cr1", tab$gene),]

head(tab)


