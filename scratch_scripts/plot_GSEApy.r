library(pheatmap)
library(data.table)
library(matrixStats)
library(IRdisplay)
library(ggplot2)

setwd("~/projects/combined_all/female_RNA/DEG_results_latent_rep_mito/EnrichR/")


options(repr.plot.width=12, repr.plot.height=8)

for (ct in c("--AMY","--ENT", "--HCP", "--NAC", "--RLP", "--HCA", "--FC", "--CP")) {
    database <- "Go_Bp"
    dir <- "Up"
    files <- list.files(pattern = ".csv", path = ".")
    files <- files[grep(database, files)]
    files <- files[which(file.exists(files))]
    files <- files[grep(ct, files)]
    files <- files[grep(dir, files)]
    pval <- 0.05

    gkey <- c()
    for (f in files) {
        if (!file.exists(f)) {
            cat(f, " NO FILE\n")
            next()
        }
        curr <- fread(f, sep = ",", fill = TRUE)
        curr <- as.data.frame(curr)
        curr <- curr[which(curr$`Adjusted P-value` < pval), ]
        gkey <- c(gkey, unique(curr$Term))
    }
    gkey = unique(gkey)
    if (length(gkey) < 2) {
        next()
    }
    
    mat <- matrix(0, nrow = length(gkey), ncol = length(files))
    gene_mat <- matrix(0, nrow = length(gkey), ncol = length(files))
    
    for (i in 1:length(gkey)) {
        for (j in 1:length(files)) {
            f <- files[j]
            curr <- fread(f, sep = ",", fill = TRUE)
            curr <- as.data.frame(curr)
            names <- curr$Term
            curr <- curr[!duplicated(names), ]
            curr <- curr[!is.na(curr$Term), ]
            if (gkey[i] %in% curr$Term) {
                mat[i, j] <- -log10(curr[curr$Term == gkey[i], "Adjusted P-value"] + 1e-100)
                genes <- unlist(strsplit(as.character(curr[curr$Term == gkey[i], "Genes"]), ";"))
                gene_mat[i, j] <- length(genes)
            } else {
                mat[i, j] <- NA
                gene_mat[i, j] <- NA
            }
        }
    }
 
    cn <- gsub(database, "", files)
    cn <- gsub(".txt", "", cn)
    cn <- gsub("/", "", cn)
    cn <- gsub("Up_", "", cn)
    cn <- gsub("Down_", "", cn)
    cn <- gsub("_.csv", "", cn)

    colnames(mat) <- cn
    colnames(gene_mat) <- cn
    rownames(mat) <- gkey
    rownames(gene_mat) <- gkey
    mat[which(is.na(mat), arr.ind = TRUE)] <- 0
    gene_mat[which(is.na(gene_mat), arr.ind = TRUE)] <- 0

    nmat <- mat
    gene_nmat <- gene_mat
    gene_nmat <- gene_nmat[, which(colMaxs(nmat) > 1.5)]
    nmat <- nmat[, which(colMaxs(nmat) > 1.5)]
    gene_nmat <- gene_nmat[order(rowMaxs(nmat), decreasing = TRUE), ]
    nmat <- nmat[order(rowMaxs(nmat), decreasing = TRUE), ]
    gene_nmat <- gene_nmat[,order(colMaxs(nmat), decreasing = TRUE)]
    nmat <- nmat[,order(colMaxs(nmat), decreasing = TRUE)]

    if (nrow(nmat) > 30) {
        nmat <- nmat[1:30, ]
        gene_nmat <- gene_nmat[1:30, ]
    }
    nmat[which(nmat > 30, arr.ind = TRUE)] <- 30

    df <- expand.grid(Celltype = colnames(nmat), Term = rownames(nmat))
    df$Size <- as.vector(t(gene_nmat))
    df$Adjusted_P_value <- as.vector(t(nmat))
    
    p <- ggplot(df, aes(x = Celltype, y = Term, size = Size, color = Adjusted_P_value)) +
        geom_point() +
        scale_size_continuous(range = c(0, 10)) +
        scale_color_gradient(low = "blue", high = "red") +
        theme_minimal() +
        ggtitle(paste(gsub("--", "", ct), dir, database)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    pdf(paste(dir, ct, database,"dotplot","pdf", sep = "."),height = 9, width = 14)
    print(p)
    dev.off()
    print(p)
}


library(ggplot2)
options(repr.plot.width=12, repr.plot.height=8)

for (ct in c("OPC_NN", "Oligo_NN", "Microglia_NN", "Astro", "DG_Glut", "STR_D12_Gaba", "CA3_Glut")) {
    database <- "KEGG"
    dir <- "Down"
    files <- list.files(pattern = ".csv", path = ".")
    files <- files[grep(database, files)]
    files <- files[which(file.exists(files))]
    files <- files[grep(ct, files)]
    files <- files[grep(dir, files)]
    pval <- 0.05

    gkey <- c()
    for (f in files) {
        if (!file.exists(f)) {
            cat(f, " NO FILE\n")
            next()
        }
        curr <- fread(f, sep = ",", fill = TRUE)
        curr <- as.data.frame(curr)
        curr <- curr[which(curr$`Adjusted P-value` < pval), ]
        gkey <- c(gkey, unique(curr$Term))
    }
    gkey = unique(gkey)
    if (length(gkey) < 2) {
        next()
    }
    
    mat <- matrix(0, nrow = length(gkey), ncol = length(files))
    gene_mat <- matrix(0, nrow = length(gkey), ncol = length(files))
    
    for (i in 1:length(gkey)) {
        for (j in 1:length(files)) {
            f <- files[j]
            curr <- fread(f, sep = ",", fill = TRUE)
            curr <- as.data.frame(curr)
            names <- curr$Term
            curr <- curr[!duplicated(names), ]
            curr <- curr[!is.na(curr$Term), ]
            if (gkey[i] %in% curr$Term) {
                mat[i, j] <- -log10(curr[curr$Term == gkey[i], "Adjusted P-value"] + 1e-100)
                genes <- unlist(strsplit(as.character(curr[curr$Term == gkey[i], "Genes"]), ";"))
                gene_mat[i, j] <- length(genes)
            } else {
                mat[i, j] <- NA
                gene_mat[i, j] <- NA
            }
        }
    }
 
    cn <- gsub(database, "", files)
    cn <- gsub(".txt", "", cn)
    cn <- gsub("/", "", cn)
    cn <- gsub("Up_", "", cn)
    cn <- gsub("Down_", "", cn)
    cn <- gsub("_.csv", "", cn)

    colnames(mat) <- cn
    colnames(gene_mat) <- cn
    rownames(mat) <- gkey
    rownames(gene_mat) <- gkey
    mat[which(is.na(mat), arr.ind = TRUE)] <- 0
    gene_mat[which(is.na(gene_mat), arr.ind = TRUE)] <- 0

    nmat <- mat
    gene_nmat <- gene_mat
    gene_nmat <- gene_nmat[, which(colMaxs(nmat) > 1.5)]
    nmat <- nmat[, which(colMaxs(nmat) > 1.5)]
    gene_nmat <- gene_nmat[order(rowMaxs(nmat), decreasing = TRUE), ]
    nmat <- nmat[order(rowMaxs(nmat), decreasing = TRUE), ]
    gene_nmat <- gene_nmat[,order(colMaxs(nmat), decreasing = TRUE)]
    nmat <- nmat[,order(colMaxs(nmat), decreasing = TRUE)]

    if (nrow(nmat) > 30) {
        nmat <- nmat[1:30, ]
        gene_nmat <- gene_nmat[1:30, ]
    }
    nmat[which(nmat > 30, arr.ind = TRUE)] <- 30

    df <- expand.grid(Celltype = colnames(nmat), Term = rownames(nmat))
    df$Size <- as.vector(t(gene_nmat))
    df$Adjusted_P_value <- as.vector(t(nmat))
    
    p <- ggplot(df, aes(x = Celltype, y = Term, size = Size, color = Adjusted_P_value)) +
        geom_point() +
        scale_size_continuous(range = c(0, 10)) +
        scale_color_gradient(low = "blue", high = "red") +
        theme_minimal() +
        ggtitle(paste(gsub("--", "", ct), dir, database)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    pdf(paste(dir, ct, database,"dotplot","pdf", sep = "."),height = 9, width = 14)
    print(p)
    dev.off()
    print(p)
}


#setwd("~/ps-renlab2/projects/combined_all/female_RNA/DEG_results_region//")
for(ct in c("--ENT", "--HCP", "--AMY", "--NAC", "--RLP", "--HCA", "--FC", "--CP")) {
    database = "Go_Bp"
    dir = "Down"
    #ct = "--HCA"
    #setwd(paste( dir,"/GO",sep = ""))
    files = list.files(pattern = ".csv", path = ".")
    files = files[grep(database,files)]
    files = files[which(file.exists(files))]
    files=files[grep(ct, files)]
    files=files[grep(dir, files)]
    pval = .01
    f = files[1]
    curr = fread(f, sep = ",",fill = T)
    curr = as.data.frame(curr)
    curr = curr[which(curr$`Adjusted P-value` < pval),]
    gkey = paste(curr[,c(3)])
    for(f in files[2:length(files)]) {
      if(!file.exists(f)){
        cat(f, " NO FILE\n")
        next()
      }
      curr = fread(f, sep = ",",fill = T)
      curr = as.data.frame(curr)
      #  show(head(curr))
      curr = curr[which(curr$`Adjusted P-value` < pval),]
      gkey = c(gkey ,paste(curr[,3]))
    }

    gkey = gkey[which(!duplicated(gkey))]
    cat(length(gkey))
    if(length(gkey)<2){
        next
    }

   # files=files[-length(files)]


    # getting significance for each term + file
    mat = matrix(0, nrow=length(gkey),ncol=length(files))
    smat = matrix('', nrow=length(gkey),ncol=length(files))

    x = 1
    for(i in 1:length(gkey)) {
      y = 1
      for(f in files) {
        curr = fread(f, sep = ",", fill=T)
        curr = as.data.frame(curr)
        names = curr$Term
        which(duplicated(names))
        curr = curr[which(!duplicated(names)),]
        curr = curr[which(!is.na(curr$Term)),]

        rownames(curr) = curr$Term
        if(paste(gkey[i]) %in% paste(rownames(curr))) {
           mat[x,y]  = -log10(curr[paste(gkey[i]), "Adjusted P-value"]+1e-100)
          #mat[x,y]  = -curr[paste(gkey[i]),"logP"]
          #  cat(mat[x,y], " ")
            if(curr[paste(gkey[i]), "Adjusted P-value"]<0.05){
                smat[x,y]= '*'
            }
            if(curr[paste(gkey[i]), "Adjusted P-value"]<0.01){
                smat[x,y] = '**'
            }
            if(curr[paste(gkey[i]), "Adjusted P-value"]<0.001){
                smat[x,y] = '***'
            }
        } else {
          mat[x,y]  = NA
        }
        y= y+1
      }
      x=x+1
    }

    options(repr.plot.width=15, repr.plot.height=10)
    rownames(mat) = gkey
    cn = gsub(database, "", files)
    cn = gsub(".txt", "", cn)
    cn = gsub("/", "", cn)
    cn = gsub("Up_", "",cn)
    cn = gsub("Down_", "",cn)
    cn = gsub("_.csv", "",cn)


    colnames(mat) = cn
    mat[which(is.na(mat), arr.ind = T)] = 0
    nmat = mat
    #nmat = nmat[which(rowMaxs(nmat)>2), ]
    smat=smat[,which(colMaxs(nmat)>1.5)]
    nmat=nmat[,which(colMaxs(nmat)>1.5)]

    smat = smat[order(rowMaxs(nmat), decreasing = T),]
    nmat = nmat[order(rowMaxs(nmat), decreasing = T),]
    if(nrow(nmat)>30) {
        nmat= nmat[1:30,]
        smat= smat[1:30,]

    }
    #nmat = nmat[,which(colMaxs(nmat)>7)]
    nmat[which(nmat>20, arr.ind = T)] = 20

    clade = colnames(nmat)
    clade[grep("Neur", clade)] = "Gaba"
    clade[grep("Gaba", clade)] = "Gaba"
    clade[grep("Glut", clade)] = "Glut"
    clade[grep("NN", clade)] = "NN"
    clade[grep("IMN", clade)] = "IMN"

    nannotation_col = data.frame(
        Clade = clade
      )
    rownames(nannotation_col) = colnames(nmat) 

    heatmap <- pheatmap(nmat, main = paste(gsub("--", "",ct), dir, database), annotation_col = nannotation_col,display_numbers = smat)
    print(heatmap)
    pdf(paste(dir, ct, database,"pdf", sep = "."),height = 9, width = 14)
    print(heatmap)
    dev.off()
}

    dev.off()


#setwd("~/ps-renlab2/projects/combined_all/female_RNA/DEG_results_region//")
for(ct in c("OPC_NN", "Oligo_NN", "Microglia_NN", "Astro", "DG_Glut", "STR_D12_Gaba", "CA3_Glut")) {
    database = "Go_Bp"
    dir = "Down"
    #ct = "--HCA"
    #setwd(paste( dir,"/GO",sep = ""))
    files = list.files(pattern = ".csv", path = ".")
    files = files[grep(database,files)]
    files = files[which(file.exists(files))]
    files=files[grep(ct, files)]
    files=files[grep(dir, files)]
    pval = .05
    f = files[1]
    curr = fread(f, sep = ",",fill = T)
    curr = as.data.frame(curr)
    curr = curr[which(curr$`Adjusted P-value` < pval),]
    gkey = paste(curr[,c(3)])
    for(f in files[2:length(files)]) {
      if(!file.exists(f)){
        cat(f, " NO FILE\n")
        next()
      }
      curr = fread(f, sep = ",",fill = T)
      curr = as.data.frame(curr)
      #  show(head(curr))
      curr = curr[which(curr$`Adjusted P-value` < pval),]
      gkey = c(gkey ,paste(curr[,3]))
    }

    gkey = gkey[which(!duplicated(gkey))]
    cat(length(gkey))
    if(length(gkey)<2){
        next
    }

   # files=files[-length(files)]


    # getting significance for each term + file
    mat = matrix(0, nrow=length(gkey),ncol=length(files))
    smat = matrix('', nrow=length(gkey),ncol=length(files))

    x = 1
    for(i in 1:length(gkey)) {
      y = 1
      for(f in files) {
        curr = fread(f, sep = ",", fill=T)
        curr = as.data.frame(curr)
        names = curr$Term
        which(duplicated(names))
        curr = curr[which(!duplicated(names)),]
        curr = curr[which(!is.na(curr$Term)),]

        rownames(curr) = curr$Term
        if(paste(gkey[i]) %in% paste(rownames(curr))) {
           mat[x,y]  = -log10(curr[paste(gkey[i]), "Adjusted P-value"]+1e-100)
          #mat[x,y]  = -curr[paste(gkey[i]),"logP"]
          #  cat(mat[x,y], " ")
            if(curr[paste(gkey[i]), "Adjusted P-value"]<0.05){
                smat[x,y]= '*'
            }
            if(curr[paste(gkey[i]), "Adjusted P-value"]<0.01){
                smat[x,y] = '**'
            }
            if(curr[paste(gkey[i]), "Adjusted P-value"]<0.001){
                smat[x,y] = '***'
            }
        } else {
          mat[x,y]  = NA
        }
        y= y+1
      }
      x=x+1
    }

    options(repr.plot.width=15, repr.plot.height=10)
    rownames(mat) = gkey
    cn = gsub(database, "", files)
    cn = gsub(".txt", "", cn)
    cn = gsub("/", "", cn)
    cn = gsub("Up_", "",cn)
    cn = gsub("Down_", "",cn)
    cn = gsub("_.csv", "",cn)


    colnames(mat) = cn
    mat[which(is.na(mat), arr.ind = T)] = 0
    nmat = mat
    #nmat = nmat[which(rowMaxs(nmat)>2), ]
   # smat=smat[,which(colMaxs(nmat)>1.5)]
    #nmat=nmat[,which(colMaxs(nmat)>1.5)]

    smat = smat[order(rowMaxs(nmat), decreasing = T),]
    nmat = nmat[order(rowMaxs(nmat), decreasing = T),]
    if(nrow(nmat)>30) {
        nmat= nmat[1:30,]
        smat= smat[1:30,]

    }
    #nmat = nmat[,which(colMaxs(nmat)>7)]
    nmat[which(nmat>20, arr.ind = T)] = 20

    clade = colnames(nmat)
    clade[grep("Neur", clade)] = "Gaba"
    clade[grep("Gaba", clade)] = "Gaba"
    clade[grep("Glut", clade)] = "Glut"
    clade[grep("NN", clade)] = "NN"
    clade[grep("IMN", clade)] = "IMN"

    nannotation_col = data.frame(
        Clade = clade
      )
    rownames(nannotation_col) = colnames(nmat) 

    heatmap <- pheatmap(nmat, main = paste(gsub("--", "",ct), dir, database), annotation_col = nannotation_col,display_numbers = smat)
    print(heatmap)
    pdf(paste(ct, dir,database,"pdf", sep = "."),height =7,width=10)
    print(heatmap)
    dev.off()
}


