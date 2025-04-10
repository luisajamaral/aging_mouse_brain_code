library(dplyr)
library(smoother)

setwd("../../region_DARs_redone/")
chr.len = read.table("/tscc/projects/ps-renlab/share/bwa_indices/mm10.fa.fai")

pval <- 0.01
win_exp <- 5
win_size <- 10^win_exp

files <- list.files(path = ".", pattern = "diff_peaks", full.names = TRUE)
chr.len <- rbind(chr.len, chr.len)
chr.len$logFC <- rep(c(-1, 1), each = nrow(chr.len) / 2)
colnames(chr.len)[1] <- "chr"

for (direction in c("Up", "Down")) {
  smthscore_dict <- list()
  diff_cluster_dict <- list()
  direction_bool <- ifelse(direction == "Up", TRUE, FALSE)
  
  for (file in files) {
    print(file)
    a <- read.csv(file)
    a$chr <- sub("(.*):(.*)-(.*)", "\\1", a$feature.name)
    if (direction_bool){
        a <- a[a$adjusted.p.value < pval & a$chr %in% chr.len$chr & (a$log2.fold_change.) > 0.25, ]
    } else {
        a <- a[a$adjusted.p.value < pval & a$chr %in% chr.len$chr & (a$log2.fold_change.) < -0.25, ]

    }
    if (nrow(a) > 0) {
      if (nrow(a) > 5000) {
        a <- a[1:5000, ]
      }
      a$sample <- file
      a$pos <- floor(as.integer(sub("(.*):(.*)-(.*)", "\\2", a$feature.name)) / win_size)
      #b <- a %>% count(chr, pos)
      b <- a %>% count(chr, pos, abs(log2.fold_change.) > 0)
      colnames(b)[3] <- "logFC"
      b <- data.frame(b)
      print(direction)
      print(head(b))
      clust <- gsub("diff_peaks_", "", file)
      clust <- gsub("_2vs18.csv", "", clust)
      clust <- gsub("./", "", clust)
      print(clust)
      
      pdf(paste("smooth/",clust, "_", direction, "_", pval, "_smoothed_diffplots.pdf", sep = ""))
      
      locations <- c()
      counts <- c()
      smth20_score <- c()
      
      for (chr in unique(chr.len$chr)) {
        bo <- b[b$chr == chr, ]
        if (nrow(bo) == 0) next
        
        rownames(bo) <- paste0(bo$chr, ":", bo$pos)
        min_boa <- max(0, min(bo$pos) - 11)
        max_boa <- min(round(chr.len[chr.len$chr == chr, 2] / win_size), max(bo$pos) + 11)
        
        boa <- data.frame(chr = chr, pos = min_boa:max_boa, count = 0)
        rownames(boa) <- paste0(boa$chr, ":", boa$pos)
        boa[rownames(bo), "count"] <- bo$n
        print(length(boa$count))
        print(sum(boa$count))
        # Check for sufficient data before smoothing
        if (length(boa$count) <= 20 || all(boa$count == 0)) {
          cat("Skipping smoothing: insufficient data for chr", chr, "\n")
          next
        }
        
        # Apply smoothing
        if (length(boa$count) >= 5 && any(boa$count > 0)) {
          smth5 <- smth.gaussian(boa$count, window = 5, tails = TRUE)
        } else {
          smth5 <- rep(0, length(boa$count))
        }

        if (length(boa$count) >= 10 && any(boa$count > 0)) {
          smth10 <- smth.gaussian(boa$count, window = 10, tails = TRUE)
        } else {
          smth10 <- rep(0, length(boa$count))
        }

        if (length(boa$count) > 20 && any(boa$count > 0)) {
          smth20 <- smth.gaussian(boa$count, window = 20, tails = TRUE)
        } else {
          smth20 <- rep(0, length(boa$count))
        }

        # Plot smoothed scores
        col <- ifelse(smth20 > 0.2, "red", "black")
        plot(boa$pos, smth20, main = chr, col = col, ylim = c(0, 0.75))
        
        boa <- cbind(boa, smth20, clust)
        smthscore_dict[[paste(chr, clust, direction, sep = "_")]] <- boa
        boa <- cbind(boa, smth10, smth20, smth5)
        boa <- boa[order(boa$smth20, decreasing = TRUE), ]
        
        ordered_pos <- rownames(boa)[boa$smth20 > 0.2]
        boa <- boa[order(boa$pos), ]
        
        while (length(ordered_pos) > 0) {
          mid <- ordered_pos[1]
          mid_pos <- boa[mid, "pos"]
          
          # Determine cluster boundaries
          d <- 10
          end_loc <- mid_pos
          c <- 1
          while (d > 0) {
            d <- d - 1
            curr_pos <- mid_pos + c
            if (curr_pos > max(boa$pos)) break
            if (boa[boa$pos == curr_pos, "count"] > 0) {
              end_loc <- curr_pos
              d <- 10
            }
            c <- c + 1
          }
          
          d <- 10
          start_loc <- mid_pos
          c <- 1
          while (d > 0) {
            d <- d - 1
            curr_pos <- mid_pos - c
            if (curr_pos < min(boa$pos)) break
            if (boa[boa$pos == curr_pos, "count"] > 0) {
              start_loc <- curr_pos
              d <- 10
            }
            c <- c + 1
          }
          
          locations <- c(locations, paste(chr, ":", start_loc, paste(rep("0", win_exp), collapse = ""), "-", end_loc + 1, paste(rep("0", win_exp), collapse = ""), sep = ""))
          counts <- c(counts, sum(boa[boa$pos %in% start_loc:end_loc, "count"]))
          smth20_score <- c(smth20_score, boa[mid, "smth20"])
          
          # Remove processed positions
          ordered_pos <- ordered_pos[!ordered_pos %in% rownames(boa[boa$pos %in% start_loc:end_loc, ])]
        }
      }
      dev.off()
      
      if (length(locations) > 0) {
        ret <- data.frame(locations, counts, smth20_score, clust)
        ret <- ret[ret$counts > 2, ]
        diff_cluster_dict[[paste(clust, direction, sep = "_")]] <- ret
      }
    }
  }
  
  # Save smoothed scores and clusters
  smthscore_table <- do.call(rbind, smthscore_dict)
  if(nrow(smthscore_table)>0){
  colnames(smthscore_table) <- c("chr", "window", "count", "smth20_score", "cell type")
  write.table(smthscore_table, row.names = FALSE, quote = FALSE, sep = "\t",
              file = paste("smooth/","smth20_scores_table_top5k", direction, ".txt", sep = ""))
  }
  diff_cluster_table <- do.call(rbind, diff_cluster_dict)
  diff_cluster_table <- diff_cluster_table[order(diff_cluster_table$counts, diff_cluster_table$locations, decreasing = TRUE), ]
  write.table(diff_cluster_table, row.names = FALSE, quote = FALSE, sep = "\t",
              file = paste("smooth/","clustered_diff_peaks_top5k", direction, "_", pval, ".txt", sep = ""))
}


smth.gaussian(boa$count, window = 20, tails = TRUE)

chr

    head(boa)



