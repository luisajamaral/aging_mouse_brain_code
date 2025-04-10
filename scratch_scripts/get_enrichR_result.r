library(data.table)
library(biomaRt)
library(enrichR)
library(tidyr)
library(ggplot2)
# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")



file = "Male.CGN.DMG.csv"
fh= read.table(file, sep = ",", header =T)
nrow(fh)
fh = fh[which(fh$corrected_pvalue<0.05 & fh$abs_change>0.05),]
nrow(fh)

data = as.data.frame(table(fh$celltype, fh$methylation_change>0))
pdf(gsub(".csv",".DMG_count.pdf", file),height =5 , width =8)
p = ggplot(data, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Category", y = "Count", fill = "Hyper in aging") +
  ggtitle(paste("# DMG per celltype (",gsub(".csv","",file),")\n corrected_pvalue<0.05, abs_change>0.05", sep = ""))
print(p)
dev.off()


table(fh$trend)


genes =  getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = gsub("\\..*","", fh$geneslop2k),
                      mart = ensembl)
if(length(which(duplicated(genes$ensembl_gene_id)))>0){
    genes = genes[-which(duplicated(genes$ensembl_gene_id)),]
}
rownames(genes) = genes$ensembl_gene_id

fh$geneid = genes[paste(gsub("\\..*","", fh$geneslop2k)),"external_gene_name"]

fh$celltype = gsub("/", "-", fh$celltype)
fh$celltype = gsub(" ", ".", fh$celltype)


for(ct in unique(fh$celltype)){
    gs  = fh[which(fh$celltype==ct & fh$methylation_change>0 ),"geneid"]
    if(length(gs[which(gs!="NA")])>1){
        enriched <- enrichr(gs, dbs)
        if(length(which(enriched$GO_Biological_Process_2023$Adjusted.P.value<0.05))>0){
            write.table(enriched$GO_Biological_Process_2023[which(enriched$GO_Biological_Process_2023$Adjusted.P.value<0.05),], 
                    file = paste(ct, gsub(".csv","_Hyper_GO_BP.txt", file)),sep = "_")
            pdf(paste(ct, gsub(".csv","_Hyper_GO_BP.pdf", file),sep = "_"))
            print(plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",title = paste(ct, "Hyper")))
            dev.off()
        }
        if(length(which(enriched$KEGG_2021_Human$Adjusted.P.value<0.05))>0){
            write.table(enriched$KEGG_2021_Human[which(enriched$KEGG_2021_Human$Adjusted.P.value<0.05),], 
                    file = paste(ct, gsub(".csv","_Hyper_KEGG.txt", file)),sep = "_")
            pdf(paste(ct, gsub(".csv","_Hyper_KEGG.pdf", file),sep = "_"))
            print(plotEnrich(enriched[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",title = paste(ct, "Hyper")))
            dev.off()
        }
    }
    
    if(length(gs[which(gs!="NA")])>1){
        gs  = fh[which(fh$celltype==ct & fh$methylation_change<0 ),"geneid"]
        enriched <- enrichr(gs, dbs)
        if(length(which(enriched$GO_Biological_Process_2023$Adjusted.P.value<0.05))>0){
            write.table(enriched$GO_Biological_Process_2023[which(enriched$GO_Biological_Process_2023$Adjusted.P.value<0.05),], 
                    file = paste(ct, gsub(".csv","_Hypo_GO_BP.txt", file)),sep = "_")
            pdf(paste(ct, gsub(".csv","_Hypo_GO_BP.pdf", file),sep = "_"))
            print(plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",title = paste(ct, "Hypo")))
            dev.off()
        }
        if(length(which(enriched$KEGG_2021_Human$Adjusted.P.value<0.05))>0){
            write.table(enriched$KEGG_2021_Human[which(enriched$KEGG_2021_Human$Adjusted.P.value<0.05),], 
                    file = paste(ct, gsub(".csv","_Hypo_KEGG.txt", file),sep = "_"))
            pdf(paste(ct, gsub(".csv","_Hypo_KEGG.pdf", file),sep = "_"))
            print(plotEnrich(enriched[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",title = paste(ct, "Hypo")))
            dev.off()
        }
    }
}




