library(dplyr)
library(smoother)

setwd("../../h5ads_final/region_sex_DARs/diff_csvs/")


smthscore_table

pval = 0.05
direction = "Down" 
locations = c()

win_exp = 5
win_size = 10^win_exp
 # setwd(tissue)
  files = list.files(path=".",pattern="iter.csv",full.names=T)
  chr.len = read.table("/projects/ps-renlab/share/bwa_indices/mm10.fa.fai")
  chr.len = rbind(chr.len,chr.len)
  chr.len$logFC = rep(c(-1,1),each=nrow(chr.len)/2)
  colnames(chr.len)[1] = "chr"
  smthscore_dict = list()
  diff_cluster_dict = list()
  for (file in files) {
    print(file)
    a = read.csv(file) 
    a$chr = sub("(.*):(.*)-(.*)", "\\1", a$feature.name)

    #if (nrow(a) > 0) {
    a = a[which(a$adjusted.p.value < pval & a$chr %in% chr.len$chr & abs(a$log2.fold_change.)>.2), ]
    if (nrow(a) > 0) {
      if(nrow(a)>500){
       a = a[1:500,] 
      }
      a$sample = file
      a$pos = floor(as.integer(sub("(.*):(.*)-(.*)", "\\2", a$feature.name)) / win_size)
      b = a %>% count(chr, pos, log2.fold_change. < 0)
      colnames(b)[3] = "logFC"
      b = data.frame(b)
      clust = gsub("diff_csvs/", "", file)
      clust = gsub(".edger.txt", "", clust)
      clust = gsub("_2vs18_iter.csv", "", clust)
      clust = gsub("diff_peaks_", "", clust)

      #pdf(paste(clust,"_",direction,"_",pval,"_smoothed_diffplots.pdf", sep = ""))
      if(direction == "Up") {
        direction_bool = TRUE
      } else {
        direction_bool = FALSE
      }
      locations = c()
      counts = c()
      smth20_score = c()
      for (chr in unique(chr.len$chr)) {
        bo = b[which(b$logFC == direction_bool & b$chr==chr),]
        if(nrow(bo) == 0) {
          next
        }
        rownames(bo) = paste0(bo$chr, ":",bo$pos)
        min_boa = max(0,min(bo$pos)-11)
        max_boa = min(round(chr.len[which(chr.len$chr==chr)[1],2]/win_size),max(bo$pos)+11)
        cat(chr,length(min_boa:max_boa), "\n")
        if (length(min_boa:max_boa)<=20) {
          max_boa = max_boa + 20
          min_boa = min_boa - 20
        }
        boa = data.frame(chr = chr, pos = min_boa:max_boa, count = 0 )
        rownames(boa) = paste0(boa$chr, ":",boa$pos)
        boa[rownames(bo),"count"] = bo$n  
        
        smth5 = smth.gaussian(boa$count, window = 5, tails = T)
        smth10 = smth.gaussian(boa$count, window = 10, tails = T)
        smth20 = smth.gaussian(boa$count, window = 20, tails = T)

        col = rep("black",length(smth20))
        col[which(smth20>.2)] = "red"
        #print(plot(boa$pos, smth20, main = chr, col = col, ylim = c(0,.75)) )
        boa = cbind(boa,smth20, clust)
        smthscore_dict[[paste(chr,clust, sep = "_")]] = boa
        boa = cbind(boa, smth10, smth20, smth5)
        boa = boa[order(boa$smth20,decreasing = T),]
        # .2 is the arbitrary score I chose based on the plots when window_size is 10^5
        ordered_pos = rownames(boa)[which(boa$smth20>.2)]
        boa = boa[order(boa$pos),]
        while (length(ordered_pos)>0) { 
          mid = ordered_pos[1]
          mid_pos = boa[mid,"pos"]
          # d is the maximum distance (in window_size) between two differential peaks for them to 
          # be considered in the same domain
          d = 10
          end_loc = mid_pos
          c = 1
          # gets end position
          while (d > 0) {
            d = d-1
            curr_pos = mid_pos+c
            if (curr_pos>max(boa$pos)){
              curr_count = 0
              if (end_loc == -1) {
                end_loc = curr_pos-1
              }
              break
            } else{
              curr_count = boa[which(boa$pos == curr_pos),"count"]
            }
            if (curr_count > 0) {
              cat(curr_pos, curr_count, d, end_loc, "\n")
              end_loc = curr_pos
              d = 10
            }
            c = c+1
          }
          
          d = 10
          start_loc = mid_pos
          c = 1
          # gets start position
          while (d > 0) {
            d = d-1
            curr_pos = mid_pos-c
            if (curr_pos<min(boa$pos)){
              curr_count = 0
              if (end_loc == -1) {
                end_loc = curr_pos+1
              }
              break
            } else{
            curr_count = boa[which(boa$pos == curr_pos),"count"]
            }
            if (curr_count > 0) {
              cat(curr_pos, curr_count, d, start_loc, "\n")
              start_loc = curr_pos
              d = 10
            }
            c = c+1
          }
          locations = c(locations, paste(chr,":",start_loc,paste(rep("0",win_exp),collapse=""),"-", end_loc+1,paste(rep("0",win_exp),collapse=""), sep = ""))
          sum_count = sum(boa[which(boa$pos %in% start_loc:end_loc),"count"])
          counts = c(counts, sum_count)
          smth20_score = c(smth20_score,boa[mid,"smth20"])
          
          # remove all regions near mid
          idx_mid = which(boa$pos == mid_pos)
          idx_0s = which(boa$smth20<0.05)
          if (length(which(idx_0s>idx_mid))>0) {
            end_region = idx_0s[which(idx_0s>idx_mid)[1]] 
          } else {
            end_region = idx_mid+1
          }
          if (length(which(idx_0s<idx_mid))>0) {
            start_region = idx_0s[which(idx_0s<idx_mid)[length(which(idx_0s<idx_mid))]]
          } else {
            start_region = idx_mid-1
          }
          remove_from_pos = rownames(boa[start_region:end_region,])
          ordered_pos = ordered_pos[-which(ordered_pos %in% remove_from_pos)]
        }
      }
    #dev.off()
    }
    if (length(locations)>0) {
      ret = data.frame(locations, counts,smth20_score, clust)
      ret = data.frame(ret)
      if(length(which(ret$counts <=2))>0){
        ret = ret[-which(ret$counts <=2),]
      }
      show(ret)
      diff_cluster_dict[[clust]] = ret
    } 
  }
  smthscore_table = do.call(rbind, smthscore_dict)
  colnames(smthscore_table) = c("chr", "window", "count","smth20_score", "cell type")
  write.table(smthscore_table,row.names = F, quote = F, sep = "\t", 
              file = paste("smth20_scores_table_500_", direction,".txt",sep = ""))
#                                  direction,"_",pval, ".txt",sep = ""))
  diff_cluster_table = do.call(rbind, diff_cluster_dict)
  diff_cluster_table = diff_cluster_table[order(diff_cluster_table$counts,diff_cluster_table$locations,decreasing = T),]
  write.table(diff_cluster_table, row.names = F, quote = F, sep = "\t", file = paste("clustered_diff_peaks",
                                                                                     direction,"_500_",pval, ".txt",sep = ""))
  return(diff_cluster_dict)

head(a)


