library(Seurat)



tracksheet= read.csv("~/projects/fmouse_multiome/tracksheet_with_paths_oct23.csv")

tracksheet= read.csv("../../tracksheet_with_paths_oct23.csv")

head(tracksheet)

rownames(tracksheet) = tracksheet$Sample_name

obj = readRDS("../female_RNA/RNA_final_SCT.RDS")

FeaturePlot(obj, "Malat1")

obj_r = readRDS("../female_RNA/RNA_final.RDS")

obj_r@assays$RNA@counts[1:10,1:10]

DimPlot(obj)

DimPlot(obj, label = TRUE, group.by = "is_Frozen")

DimPlot(obj, label = TRUE, group.by = "celltype")+NoLegend()

obj$is_Frozen = factor(obj$is_Frozen)

head(obj$is_Frozen)

seuratobj = obj
# Normalize the data
seuratobj <- NormalizeData(seuratobj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify variable features
seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 2000)

# Scale data and regress out batch effect (is_Frozen) and library size
seuratobj <- ScaleData(seuratobj, vars.to.regress = c("is_Frozen", "nCount_RNA"), verbose = TRUE)

# Run PCA
seuratobj <- RunPCA(seuratobj, npcs = 30, verbose = FALSE)

# Elbow plot to determine dimensions
ElbowPlot(seuratobj)

# Clustering and UMAP
seuratobj <- FindNeighbors(seuratobj, dims = 1:50)  # Adjust dims based on ElbowPlot
seuratobj <- FindClusters(seuratobj, resolution = 0.5)  # Adjust resolution for granularity
seuratobj <- RunUMAP(seuratobj, dims = 1:50)

# Visualize UMAP
DimPlot(seuratobj, reduction = "umap", group.by = "seurat_clusters")
DimPlot(seuratobj, reduction = "umap", group.by = "is_Frozen")


DefaultAssay(obj) = "RNA"
if ("SCT" %in% Assays(obj)) {
  obj[["SCT"]] <- NULL
}

obj <- SCTransform(
    obj, 
    assay = "RNA", 
    vars.to.regress = c("is_Frozen", "percent.mt"), 
    n_genes = 5000,     # Increase the number of variable genes to consider
    n_cells = 50000,    # Increase the number of cells to use in the model
    verbose = TRUE
)

obj@reductions

obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:50)  # Adjust dims based on ElbowPlot
obj <- FindClusters(obj, resolution = 0.5)  # Adjust resolution as needed
obj <- RunUMAP(obj, dims = 1:50, reduction = "pca")
DimPlot(obj, reduction = "umap", group.by = "seurat_clusters")


DimPlot(obj, reduction = "umap", group.by = "is_Frozen")


DimPlot(obj, group.by = "celltype",label=T) + NoLegend()


saveRDS(file = "SCT_is_Frozen_ptmt.RDS", obj)

DimPlot(obj, label = TRUE)

obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:50, verbose = FALSE)

obj <- FindNeighbors(obj, dims = 1:50, verbose = FALSE)
obj <- FindClusters(obj, verbose = FALSE)
DimPlot(obj, label = TRUE)

DimPlot(obj, label = TRUE, group.by = "is_Frozen")

DimPlot(obj, label = TRUE, group.by = "is_Frozen")

DimPlot(obj, label = TRUE, group.by = "celltype_final")+NoLegend()
DimPlot(obj, label = TRUE, group.by = "region")
DimPlot(obj, label = TRUE, group.by = "age")
DimPlot(obj, label = TRUE, group.by = "is_Frozen")
DimPlot(obj, label = TRUE, group.by = "orig.ident")

DimPlot(obj_r, label = TRUE, group.by = "age")+NoLegend()

DimPlot(obj_r, label = TRUE, group.by = "orig.ident")

DimPlot(obj_r, label = TRUE, group.by = "orig.ident")

saveRDS(obj_r, file = "RNA_SCT_latent_orig_ident.RDS") 

obj$track_Sample= tracksheet[gsub("0", "",paste(obj$orig.ident)),"Sample"]

obj_r$track_Sample= tracksheet[gsub("0", "",paste(obj_r$orig.ident)),"Sample"]
obj_r$is_Frozen = "no"
obj_r$is_Frozen[grep("Frozen", obj_r$track_Sample)] = "yes"

DimPlot(obj_r,reduction="umap",group.by = "is_Frozen",label = T)


FeaturePlot(obj_r,reduction="umap","Xist",label = T)
FeaturePlot(obj,reduction="umap","Xist",label = T)
FeaturePlot(obj,reduction="umap","Tsix",label = T)


obj_r


obj

head(obj@meta.data)

obj$is_Frozen = "no"
obj$is_Frozen[grep("Frozen", obj$track_Sample)] = "yes"


DimPlot(obj,reduction="umap",group.by = "is_Frozen",label = T)


Idents(obj) = "is_Frozen"
de_froz_genes = FindMarkers(obj, `ident.1` = "no", `ident.2` = "yes")

head(de_froz_genes,50
    )


rm(obj)

mkdir("DEG_results_latent_isFrozen_percentmt")

setwd("DEG_results_latent_isFrozen_percentmt")

setwd("DEG_results_latent_isfrozen_sample_removed_lowq_samples")

table(obj$orig.ident) 

obj = obj[,-which(obj$orig.ident %in% c("AMY_2mo_1", "AMY_18mo_1", "CP_9mo_2", "ENT_18mo_1", "ENT_9mo_2",
                            "HCP_9mo_1"))]


table(obj$orig.ident) 

head(obj@meta.data)

obj$Gapdh_expr <- FetchData(obj, vars = "Gapdh")


obj <- obj[!rownames(obj) %in% "Malat1", ]


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MAST",force = TRUE)



remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)


library(Seurat)
SessionInfo()

celltype = "Oligo NN"
r = "HCP"


ident1 <- gsub(" ", "_", paste("2mo", celltype, r, sep = "_"))
ident2 <- gsub(" ", "_", paste("18mo", celltype, r, sep = "_"))


ident1

unique(obj$orig.ident)

fm = FindMarkers(
                obj_r,
                ident.1 = ident1,
                ident.2 = ident2,
                test.use = "MAST",
                logfc.threshold = 0.25, 
                min.pct = 0.05,
                latent.vars = c('percent.mt', "orig.ident", "Gapdh_expr", "Xist_expr")
            )

setwd("../DEG_results_MAST_latent_Frozen_pc_mt/")


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
                latent.vars = c('percent.mt', "is_Frozen")
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



