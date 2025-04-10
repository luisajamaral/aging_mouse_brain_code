library(Seurat)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(ggrepel)

lkey=read.table("~/projects/combined_all/scTE/out.saf")
#lkey = lkey[-which(lkey$V5=="Retroposon"),]
ukey = lkey[-which(duplicated(lkey$V1)),]
rownames(ukey) = ukey$V1

bkey  = ukey[grep("_", rownames(ukey)), ]
rownames(bkey) = gsub("_", "-", rownames(bkey))
ukey = rbind(ukey,bkey)
ukey$V6 = gsub("[?]", "", ukey$V6)

tracksheet = read.csv("../../tracksheet_with_paths_oct23.csv")
rownames(tracksheet) = tracksheet$Sample_name


obj = readRDS("all_subfamily.RDS")

obj@meta.data$region= sapply(strsplit(as.character(obj@meta.data$sample), "_"), "[[", 1)
obj@meta.data$age= sapply(strsplit(as.character(obj@meta.data$sample), "_"), "[[", 2)
obj@meta.data$rep= sapply(strsplit(as.character(obj@meta.data$sample), "_"), "[[", 3)

rownames(tracksheet) = tracksheet$Sample_name
obj$track_Sample= tracksheet[gsub("0", "",paste(obj$sample)),"Sample"]
obj$track_Sample= tracksheet[gsub("0", "",paste(obj$sample)),"Sample"]
obj$is_Frozen = "no"
obj$is_Frozen[grep("Frozen", obj$track_Sample)] = "yes"
DimPlot(obj,reduction="umap",group.by = "is_Frozen",label = T)

# Subset genes
genes_obj <- subset(obj, features = rownames(obj)[-grep("SoloTE", rownames(obj))])

# Subset TEs
tes_obj <- subset(obj, features = rownames(obj)[grep("SoloTE", rownames(obj))])


genes_obj
tes_obj

head(tes_obj)

og = obj

genes_obj <- NormalizeData(genes_obj, normalization.method = "LogNormalize", scale.factor = 10000)
genes_obj <- FindVariableFeatures(genes_obj, selection.method = "vst")
genes_obj <- ScaleData(genes_obj)

#dir.create("../genes_only_MAST_latent_Frozen", showWarnings = FALSE)
setwd("../genes_only_MAST_latent_Frozen")

obj = genes_obj
obj$celltype_final=obj$celltype
obj$age = gsub("0", "", obj$age)
# Modify age_celltype_region to replace spaces with underscores
obj$age_celltype_region <- gsub(" ", "_", paste(obj$age, obj$celltype_final, obj$region, sep = "_"))

head(obj$age_celltype_region)
Idents(obj) <- "age_celltype_region"

for (r in c("CP", "RLP", "ENT", "HCA", "HCP", "NAC")) {
    for (celltype in unique(obj$celltype_final)) {
        # Define the file name for each cell type and region
        file_name <- paste(gsub("/", "-", gsub(" ", "_", celltype)), "--", r, ".csv", sep = "")

        # Skip if file already exists
        if (file.exists(file_name)) {
            message(file_name, " file exists")
            next
        }

        # Define the identities for the 2mo and 18mo groups
        ident1 <- gsub(" ", "_", paste("2mo", celltype, r, sep = "_"))
        ident2 <- gsub(" ", "_", paste("18mo", celltype, r, sep = "_"))

        # Check if both groups have more than 20 cells
         if (length(which(obj$age_celltype_region == ident1)) > 20 & 
            length(which(obj$age_celltype_region == ident2)) > 20 & 
           length(which(obj$age_celltype_region == ident1 & obj$rep ==1))>0 &
           length(which(obj$age_celltype_region == ident2 & obj$rep ==1))>0 &
           length(which(obj$age_celltype_region == ident1 & obj$rep ==2))>0&
           length(which(obj$age_celltype_region == ident2 & obj$rep ==2))>0) {
            message("Running ", file_name)   
            # Run differential expression analysis, excluding Malat1 and adding Gapdh as a latent variable
            de_results <- FindMarkers(
                obj,
                ident.1 = ident2,
                ident.2 = ident1,
                test.use = "MAST",
                logfc.threshold = 0.01, 
                min.pct = 0.01,
                latent.vars = c('percent.mt', "is_Frozen", "rep")
            )
            
            # Save the results to a CSV file for each cell type
            if (nrow(de_results) > 0) { # Save only if there are results
                write.csv(de_results, file = file_name, row.names = TRUE)
            } else {
                message("No significant markers found for ", file_name)
            }
        }
    }
}


for (r in c("AMY", "FC")) {
    for (celltype in unique(obj$celltype_final)) {
        # Define the file name for each cell type and region
        file_name <- paste(gsub("/", "-", gsub(" ", "_", celltype)), "--", r, ".csv", sep = "")

        # Skip if file already exists
        if (file.exists(file_name)) {
            message(file_name, " file exists")
            next
        }

        # Define the identities for the 2mo and 18mo groups
        ident1 <- gsub(" ", "_", paste("2mo", celltype, r, sep = "_"))
        ident2 <- gsub(" ", "_", paste("18mo", celltype, r, sep = "_"))

        # Check if both groups have more than 20 cells
         if (length(which(obj$age_celltype_region == ident1)) > 20 & 
            length(which(obj$age_celltype_region == ident2)) > 20 & 
           length(which(obj$age_celltype_region == ident1 & obj$rep ==1))>0 &
           length(which(obj$age_celltype_region == ident2 & obj$rep ==1))>0 &
           length(which(obj$age_celltype_region == ident1 & obj$rep ==2))>0&
           length(which(obj$age_celltype_region == ident2 & obj$rep ==2))>0) {
            message("Running ", file_name)   
            # Run differential expression analysis, excluding Malat1 and adding Gapdh as a latent variable
            de_results <- FindMarkers(
                obj,
                ident.1 = ident2,
                ident.2 = ident1,
                test.use = "MAST",
                logfc.threshold = 0.01, 
                min.pct = 0.01,
                latent.vars = c('percent.mt', "rep")
            )
            
            # Save the results to a CSV file for each cell type
            if (nrow(de_results) > 0) { # Save only if there are results
                write.csv(de_results, file = file_name, row.names = TRUE)
            } else {
                message("No significant markers found for ", file_name)
            }
        }
    }
}


tes_obj <- NormalizeData(tes_obj, normalization.method = "LogNormalize", scale.factor = 10000)
tes_obj <- FindVariableFeatures(tes_obj, selection.method = "vst")
tes_obj <- ScaleData(tes_obj)

cat("HI")

dir.create("../TEs_only_MAST_latent_Frozen_libnorm", showWarnings = FALSE)
setwd("../TEs_only_MAST_latent_Frozen_libnorm")
# Step 1: Calculate gene library size
gene_library_size <- colSums(GetAssayData(genes_obj, slot = "counts"))

# Step 2: Normalize TE counts using gene library size
te_counts <- GetAssayData(tes_obj, slot = "counts")
te_norm <- sweep(te_counts, 2, gene_library_size, FUN = "/") * 1e6  # CPM normalization
te_log_norm <- log1p(te_norm)

# Step 3: Update the Seurat object with normalized TE data
tes_obj <- SetAssayData(tes_obj, slot = "data", new.data = te_log_norm)

obj = tes_obj
obj$celltype_final=obj$celltype
obj$age = gsub("0", "", obj$age)
# Modify age_celltype_region to replace spaces with underscores
obj$age_celltype_region <- gsub(" ", "_", paste(obj$age, obj$celltype_final, obj$region, sep = "_"))

head(obj$age_celltype_region)
Idents(obj) <- "age_celltype_region"

for (r in c("CP", "RLP", "ENT", "HCA", "HCP", "NAC")) {
    for (celltype in unique(obj$celltype_final)) {
        # Define the file name for each cell type and region
        file_name <- paste(gsub("/", "-", gsub(" ", "_", celltype)), "--", r, ".csv", sep = "")

        # Skip if file already exists
        if (file.exists(file_name)) {
            message(file_name, " file exists")
            next
        }

        # Define the identities for the 2mo and 18mo groups
        ident1 <- gsub(" ", "_", paste("2mo", celltype, r, sep = "_"))
        ident2 <- gsub(" ", "_", paste("18mo", celltype, r, sep = "_"))

        # Check if both groups have more than 20 cells
        if (length(which(obj$age_celltype_region == ident1)) > 20 & 
            length(which(obj$age_celltype_region == ident2)) > 20 & 
           length(which(obj$age_celltype_region == ident1 & obj$rep ==1))>0 &
           length(which(obj$age_celltype_region == ident2 & obj$rep ==1))>0 &
           length(which(obj$age_celltype_region == ident1 & obj$rep ==2))>0&
           length(which(obj$age_celltype_region == ident2 & obj$rep ==2))>0) {
            message("Running ", file_name)
            
            # Run differential expression analysis, excluding Malat1 and adding Gapdh as a latent variable
            de_results <- FindMarkers(
                obj,
                ident.1 = ident2,
                ident.2 = ident1,
                test.use = "MAST",
                logfc.threshold = 0.01, 
                min.pct = 0.01,
                latent.vars = c('percent.mt', "is_Frozen", "rep")
            )
            
            # Save the results to a CSV file for each cell type
            if (nrow(de_results) > 0) { # Save only if there are results
                write.csv(de_results, file = file_name, row.names = TRUE)
            } else {
                message("No significant markers found for ", file_name)
            }
        }
    }
}


for (r in c("AMY", "FC")) {
    for (celltype in unique(obj$celltype_final)) {
        # Define the file name for each cell type and region
        file_name <- paste(gsub("/", "-", gsub(" ", "_", celltype)), "--", r, ".csv", sep = "")

        # Skip if file already exists
        if (file.exists(file_name)) {
            message(file_name, " file exists")
            next
        }

        # Define the identities for the 2mo and 18mo groups
        ident1 <- gsub(" ", "_", paste("2mo", celltype, r, sep = "_"))
        ident2 <- gsub(" ", "_", paste("18mo", celltype, r, sep = "_"))

        # Check if both groups have more than 20 cells
        if (length(which(obj$age_celltype_region == ident1)) > 20 & length(which(obj$age_celltype_region == ident2)) > 20) {
            message("Running ", file_name)
            
            # Run differential expression analysis, excluding Malat1 and adding Gapdh as a latent variable
            de_results <- FindMarkers(
                obj,
                ident.1 = ident2,
                ident.2 = ident1,
                test.use = "MAST",
                logfc.threshold = 0.01, 
                min.pct = 0.01,
                latent.vars = c('percent.mt', "rep")
            )
            
            # Save the results to a CSV file for each cell type
            if (nrow(de_results) > 0) { # Save only if there are results
                write.csv(de_results, file = file_name, row.names = TRUE)
            } else {
                message("No significant markers found for ", file_name)
            }
        }
    }
}


#dir.create("TEs_only_MAST_latent_Frozen", showWarnings = FALSE)
#setwd("TEs_only_MAST_latent_Frozen")

obj = tes_obj
obj$celltype_final=obj$celltype
obj$age = gsub("0", "", obj$age)
# Modify age_celltype_region to replace spaces with underscores
obj$age_celltype_region <- gsub(" ", "_", paste(obj$age, obj$celltype_final, obj$region, sep = "_"))

head(obj$age_celltype_region)
Idents(obj) <- "age_celltype_region"

for (r in c("CP", "RLP", "ENT",  "HCP", "NAC")) {
    for (celltype in unique(obj$celltype_final)) {
        # Define the file name for each cell type and region
        file_name <- paste(gsub("/", "-", gsub(" ", "_", celltype)), "9v18--", r, ".csv", sep = "")

        # Skip if file already exists
        if (file.exists(file_name)) {
            message(file_name, " file exists")
            next
        }

        # Define the identities for the 2mo and 18mo groups
        ident1 <- gsub(" ", "_", paste("9mo", celltype, r, sep = "_"))
        ident2 <- gsub(" ", "_", paste("18mo", celltype, r, sep = "_"))

        # Check if both groups have more than 20 cells
        if (length(which(obj$age_celltype_region == ident1)) > 20 & 
            length(which(obj$age_celltype_region == ident2)) > 20 & 
           length(which(obj$age_celltype_region == ident1 & obj$rep ==1))>0 &
           length(which(obj$age_celltype_region == ident2 & obj$rep ==1))>0 &
           length(which(obj$age_celltype_region == ident1 & obj$rep ==2))>0&
           length(which(obj$age_celltype_region == ident2 & obj$rep ==2))>0) {
            message("Running ", file_name)
            
            # Run differential expression analysis, excluding Malat1 and adding Gapdh as a latent variable
            de_results <- FindMarkers(
                obj,
                ident.1 = ident2,
                ident.2 = ident1,
                test.use = "MAST",
                logfc.threshold = 0.01, 
                min.pct = 0.01,
                latent.vars = c('percent.mt', "is_Frozen", "rep")
            )
            
            # Save the results to a CSV file for each cell type
            if (nrow(de_results) > 0) { # Save only if there are results
                write.csv(de_results, file = file_name, row.names = TRUE)
            } else {
                message("No significant markers found for ", file_name)
            }
        }
    }
}


for (r in c("HCA","AMY", "FC")) {
    for (celltype in unique(obj$celltype_final)) {
        # Define the file name for each cell type and region
        file_name <- paste(gsub("/", "-", gsub(" ", "_", celltype)), "9v18--", r, ".csv", sep = "")

        # Skip if file already exists
        if (file.exists(file_name)) {
            message(file_name, " file exists")
            next
        }

        # Define the identities for the 2mo and 18mo groups
        ident1 <- gsub(" ", "_", paste("9mo", celltype, r, sep = "_"))
        ident2 <- gsub(" ", "_", paste("18mo", celltype, r, sep = "_"))

        # Check if both groups have more than 20 cells
        if (length(which(obj$age_celltype_region == ident1)) > 20 & length(which(obj$age_celltype_region == ident2)) > 20) {
            message("Running ", file_name)
            
            # Run differential expression analysis, excluding Malat1 and adding Gapdh as a latent variable
            de_results <- FindMarkers(
                obj,
                ident.1 = ident2,
                ident.2 = ident1,
                test.use = "MAST",
                logfc.threshold = 0.01, 
                min.pct = 0.01,
                latent.vars = c('percent.mt', "rep")
            )
            
            # Save the results to a CSV file for each cell type
            if (nrow(de_results) > 0) { # Save only if there are results
                write.csv(de_results, file = file_name, row.names = TRUE)
            } else {
                message("No significant markers found for ", file_name)
            }
        }
    }
}


ident2
ident1
de_results <- FindMarkers(
                obj,
                ident.1 = ident2,
                ident.2 = ident1,
                test.use = "MAST",
                logfc.threshold = 0.01, 
                min.pct = 0.01,
                latent.vars = c('percent.mt', "rep")
            )

genes_obj <- NormalizeData(genes_obj, normalization.method = "LogNormalize", scale.factor = 10000)
genes_obj <- FindVariableFeatures(genes_obj, selection.method = "vst")
genes_obj <- ScaleData(genes_obj)



# Modify age_celltype_region to replace spaces with underscores
obj$age_celltype_region <- gsub(" ", "_", paste(obj$age, obj$celltype_final, obj$region, sep = "_"))

head(obj$age_celltype_region)
Idents(obj) <- "age_celltype_region"

for (r in c("CP", "RLP", "ENT", "HCA", "HCP", "NAC")) {
    for (celltype in unique(obj$celltype_final)) {
        # Define the file name for each cell type and region
        file_name <- paste(gsub("/", "-", gsub(" ", "_", celltype)), "2v9--", r, ".csv", sep = "")

        # Skip if file already exists
        if (file.exists(file_name)) {
            message(file_name, " file exists")
            next
        }

        # Define the identities for the 2mo and 18mo groups
        ident1 <- gsub(" ", "_", paste("2mo", celltype, r, sep = "_"))
        ident2 <- gsub(" ", "_", paste("9mo", celltype, r, sep = "_"))

        # Check if both groups have more than 20 cells
        if (length(which(obj$age_celltype_region == ident1)) > 20 & length(which(obj$age_celltype_region == ident2)) > 20) {
            message("Running ", file_name)
            
            # Run differential expression analysis, excluding Malat1 and adding Gapdh as a latent variable
            de_results <- FindMarkers(
                obj,
                ident.1 = ident2,
                ident.2 = ident1,
                test.use = "MAST",
                logfc.threshold = 0.01, 
                min.pct = 0.01,
                latent.vars = c('percent.mt', "is_Frozen", "rep")
            )
            
            # Save the results to a CSV file for each cell type
            if (nrow(de_results) > 0) { # Save only if there are results
                write.csv(de_results, file = file_name, row.names = TRUE)
            } else {
                message("No significant markers found for ", file_name)
            }
        }
    }
}


for (r in c("AMY", "FC")) {
    for (celltype in unique(obj$celltype_final)) {
        # Define the file name for each cell type and region
        file_name <- paste(gsub("/", "-", gsub(" ", "_", celltype)), "2v9--", r, ".csv", sep = "")

        # Skip if file already exists
        if (file.exists(file_name)) {
            message(file_name, " file exists")
            next
        }

        # Define the identities for the 2mo and 18mo groups
        ident1 <- gsub(" ", "_", paste("2mo", celltype, r, sep = "_"))
        ident2 <- gsub(" ", "_", paste("9mo", celltype, r, sep = "_"))

        # Check if both groups have more than 20 cells
        if (length(which(obj$age_celltype_region == ident1)) > 20 & length(which(obj$age_celltype_region == ident2)) > 20) {
            message("Running ", file_name)
            
            # Run differential expression analysis, excluding Malat1 and adding Gapdh as a latent variable
            de_results <- FindMarkers(
                obj,
                ident.1 = ident2,
                ident.2 = ident1,
                test.use = "MAST",
                logfc.threshold = 0.01, 
                min.pct = 0.01,
                latent.vars = c('percent.mt', "rep")
            )
            
            # Save the results to a CSV file for each cell type
            if (nrow(de_results) > 0) { # Save only if there are results
                write.csv(de_results, file = file_name, row.names = TRUE)
            } else {
                message("No significant markers found for ", file_name)
            }
        }
    }
}



# Modify age_celltype_region to replace spaces with underscores
obj$age_celltype_region <- gsub(" ", "_", paste(obj$age, obj$celltype_final, obj$region, sep = "_"))

head(obj$age_celltype_region)
Idents(obj) <- "age_celltype_region"

for (r in c("CP", "RLP", "ENT", "HCA", "HCP", "NAC")) {
    for (celltype in unique(obj$celltype_final)) {
        # Define the file name for each cell type and region
        file_name <- paste(gsub("/", "-", gsub(" ", "_", celltype)), "9v18--", r, ".csv", sep = "")

        # Skip if file already exists
        if (file.exists(file_name)) {
            message(file_name, " file exists")
            next
        }

        # Define the identities for the 2mo and 18mo groups
        ident1 <- gsub(" ", "_", paste("18mo", celltype, r, sep = "_"))
        ident2 <- gsub(" ", "_", paste("9mo", celltype, r, sep = "_"))

        # Check if both groups have more than 20 cells
        if (length(which(obj$age_celltype_region == ident1)) > 20 & length(which(obj$age_celltype_region == ident2)) > 20) {
            message("Running ", file_name)
            
            # Run differential expression analysis, excluding Malat1 and adding Gapdh as a latent variable
            de_results <- FindMarkers(
                obj,
                ident.1 = ident1,
                ident.2 = ident2,
                test.use = "MAST",
                logfc.threshold = 0.01, 
                min.pct = 0.01,
                latent.vars = c('percent.mt', "is_Frozen", "rep")
            )
            
            # Save the results to a CSV file for each cell type
            if (nrow(de_results) > 0) { # Save only if there are results
                write.csv(de_results, file = file_name, row.names = TRUE)
            } else {
                message("No significant markers found for ", file_name)
            }
        }
    }
}


for (r in c("AMY", "FC")) {
    for (celltype in unique(obj$celltype_final)) {
        # Define the file name for each cell type and region
        file_name <- paste(gsub("/", "-", gsub(" ", "_", celltype)), "9v18--", r, ".csv", sep = "")

        # Skip if file already exists
        if (file.exists(file_name)) {
            message(file_name, " file exists")
            next
        }

        # Define the identities for the 2mo and 18mo groups
        ident1 <- gsub(" ", "_", paste("18mo", celltype, r, sep = "_"))
        ident2 <- gsub(" ", "_", paste("9mo", celltype, r, sep = "_"))

        # Check if both groups have more than 20 cells
        if (length(which(obj$age_celltype_region == ident1)) > 20 & length(which(obj$age_celltype_region == ident2)) > 20) {
            message("Running ", file_name)
            
            # Run differential expression analysis, excluding Malat1 and adding Gapdh as a latent variable
            de_results <- FindMarkers(
                obj,
                ident.1 = ident1,
                ident.2 = ident2,
                test.use = "MAST",
                logfc.threshold = 0.01, 
                min.pct = 0.01,
                latent.vars = c('percent.mt', "rep")
            )
            
            # Save the results to a CSV file for each cell type
            if (nrow(de_results) > 0) { # Save only if there are results
                write.csv(de_results, file = file_name, row.names = TRUE)
            } else {
                message("No significant markers found for ", file_name)
            }
        }
    }
}

