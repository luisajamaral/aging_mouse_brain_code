# Detect Aging-Associated Genomic Hotspots from Differential Peaks
# Dependencies: smoother, dplyr

library(dplyr)
library(smoother)

# Set working directory and input parameters
setwd("combined_DARs_redone/")
chr.len <- read.table("/tscc/projects/ps-renlab/share/bwa_indices/mm10.fa.fai")
pval <- 0.001
win_exp <- 5
win_size <- 10^win_exp

# Input differential peak files
files <- list.files(path = ".", pattern = "diff_peaks", full.names = TRUE)
chr.len <- rbind(chr.len, chr.len)
chr.len$logFC <- rep(c(-1, 1), each = nrow(chr.len) / 2)
colnames(chr.len)[1] <- "chr"

# Analyze Up and Down directions
for (direction in c("Up", "Down")) {
  smthscore_dict <- list()
  diff_cluster_dict <- list()
  direction_bool <- direction == "Up"

  for (file in files) {
    a <- read.csv(file)
    a$chr <- sub("(.*):(.*)-(.*)", "\\1", a$feature.name)

    # Filter based on direction and thresholds
    if (direction_bool) {
      a <- a[a$adjusted.p.value < pval & a$chr %in% chr.len$chr & a$log2.fold_change. > 0.25, ]
    } else {
      a <- a[a$adjusted.p.value < pval & a$chr %in% chr.len$chr & a$log2.fold_change. < -0.25, ]
    }

    if (nrow(a) == 0) next
    a <- a[1:min(nrow(a), 2000), ]

    a$sample <- file
    a$pos <- floor(as.integer(sub("(.*):(.*)-(.*)", "\\2", a$feature.name)) / win_size)
    b <- a %>% count(chr, pos)
    colnames(b)[3] <- "logFC"

    clust <- gsub("diff_peaks_|_2vs18.csv|\\./", "", file)
    pdf(paste0("smooth/", clust, "_", direction, "_", pval, "_smoothed_diffplots.pdf"))

    locations <- c()
    counts <- c()
    smth20_score <- c()

    for (chr in unique(chr.len$chr)) {
      bo <- b[b$chr == chr, ]
      if (nrow(bo) == 0) next

      rownames(bo) <- paste0(bo$chr, ":", bo$pos)
      min_pos <- max(0, min(bo$pos) - 11)
      max_pos <- min(round(chr.len[chr.len$chr == chr, 2] / win_size), max(bo$pos) + 11)
      boa <- data.frame(chr = chr, pos = min_pos:max_pos, count = 0)
      rownames(boa) <- paste0(boa$chr, ":", boa$pos)
      boa[rownames(bo), "count"] <- bo$n

      # Apply Gaussian smoothing
      smth20 <- if (length(boa$count) >= 20 && any(boa$count > 0)) {
        smth.gaussian(boa$count, window = 20, tails = TRUE)
      } else rep(0, length(boa$count))

      plot(boa$pos, smth20, main = chr,
           col = ifelse(smth20 > 0.2, "red", "black"), ylim = c(0, 0.75))

      boa <- cbind(boa, smth20, clust)
      smthscore_dict[[paste(chr, clust, direction, sep = "_")]] <- boa

      # Cluster regions with smth20 > 0.2
      ordered_pos <- rownames(boa)[boa$smth20 > 0.2]
      boa <- boa[order(boa$pos), ]

      while (length(ordered_pos) > 0) {
        mid <- ordered_pos[1]
        mid_pos <- boa[mid, "pos"]

        get_bound <- function(pos, step, limit_fun, dir = 1) {
          d <- 10
          bound <- pos
          c <- 1
          while (d > 0) {
            curr_pos <- pos + dir * c
            if (limit_fun(curr_pos)) break
            if (boa[boa$pos == curr_pos, "count"] > 0) {
              bound <- curr_pos
              d <- 10
            }
            c <- c + 1
            d <- d - 1
          }
          bound
        }

        start_loc <- get_bound(mid_pos, -1, function(x) x < min(boa$pos), dir = -1)
        end_loc <- get_bound(mid_pos, 1, function(x) x > max(boa$pos), dir = 1)

        locations <- c(locations, paste0(chr, ":", start_loc, paste(rep("0", win_exp), collapse = ""), "-", end_loc + 1, paste(rep("0", win_exp), collapse = "")))
        counts <- c(counts, sum(boa[boa$pos %in% start_loc:end_loc, "count"]))
        smth20_score <- c(smth20_score, boa[mid, "smth20"])

        ordered_pos <- ordered_pos[!ordered_pos %in% rownames(boa[boa$pos %in% start_loc:end_loc, ])]
      }
    }

    dev.off()

    if (length(locations) > 0) {
      ret <- data.frame(locations, counts, smth20_score, clust)
      diff_cluster_dict[[paste(clust, direction, sep = "_")]] <- ret[ret$counts > 2, ]
    }
  }

  # Save results
  smthscore_table <- do.call(rbind, smthscore_dict)
  if (nrow(smthscore_table) > 0) {
    colnames(smthscore_table) <- c("chr", "window", "count", "smth20_score", "celltype")
    write.table(smthscore_table, row.names = FALSE, quote = FALSE, sep = "\t",
                file = paste0("smooth/smth20_scores_table_", direction, ".txt"))
  }

  diff_cluster_table <- do.call(rbind, diff_cluster_dict)
  diff_cluster_table <- diff_cluster_table[order(diff_cluster_table$counts, decreasing = TRUE), ]
  write.table(diff_cluster_table, row.names = FALSE, quote = FALSE, sep = "\t",
              file = paste0("smooth/clustered_diff_peaks_", direction, "_", pval, ".txt"))
}
