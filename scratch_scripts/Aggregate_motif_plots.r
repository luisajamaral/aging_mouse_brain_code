library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table) 
setwd("~/projects/combined_all/aggregate_motif_plots/")


region_color = read.csv("../Figures/color_scheme/MajorRegion.csv")
rownames(region_color) = region_color$Region
region_color$Region_short = ""
region_color$Region_short[which(region_color$Region=="Anterior_Hippocampus")] = "HCA"
region_color$Region_short[which(region_color$Region=="Posterior_Hippocampus")] = "HCP"
region_color$Region_short[which(region_color$Region=="Entorhinal_Cortex")] = "ENT"
region_color$Region_short[which(region_color$Region=="Amygdala")] = "AMY"
region_color$Region_short[which(region_color$Region=="Frontal_Cortex")] = "FC"
region_color$Region_short[which(region_color$Region=="Nucleus_accumbens")] = "NAC"
region_color$Region_short[which(region_color$Region=="Caudate_Putamen")] = "CP"
region_color$Region_short[which(region_color$Region=="PAG-PCG")] = "RLP"


color_vector <- setNames(region_color$Color, region_color$Region_short)


region_color

files <- list.files(pattern = "*_all_motifs_ann_hist.txt")

cts  =sapply(strsplit(as.character(files), ":"), `[`, 1)
cts  =cts[grep("background_peaks_", cts)]
cts = gsub("background_peaks_","",cts)

cts

for(ct in cts){ 
    for(motif in motifs) {
        curr_files = files[grep(ct ,files)]
        # Load the data into a single data frame
        all_data <- lapply(curr_files, function(file) {
          df <- fread(file, header = TRUE)
          df$peak_type <- ifelse(grepl("up_peaks", file), "up", 
                                 ifelse(grepl("down_peaks", file), "down", "background"))
          df$region <-  sapply(strsplit(as.character(file), ":"), `[`, 2)
          df$region <-  sapply(strsplit(as.character(df$region), "_"), `[`, 1)

          df$file <- file
          colnames(df)[1] = c("Distance" ) 

          return(df)
        }) %>% bind_rows()
        
        all_data = as.data.frame(all_data)
        if(nrow(all_data)<100){
            print(all_data)
            next
        }
        options(repr.plot.width=6.5, repr.plot.height=2.5)

        curr = all_data[,c(1,grep(motif, colnames(all_data),fixed = T),(ncol(all_data)-2):ncol(all_data))] 
        colnames(curr)[2:4] = c("combined", "+", "-" )
        p = ggplot(curr, aes(x = as.numeric(Distance), y = combined, color = region, group = region)) +
          geom_line() + facet_wrap(~peak_type) + 
          xlim(c(-400,400))+  scale_color_manual(values = color_vector) +
          labs(x = "Distance from center (bp)", y = "Sites per bp per peak", 
               title = paste("Aggregate Plot of ",motif," Motif in ",ct ,sep = "")) +
          theme_minimal()

        ggsave(p, file = paste(ct, "_",motif, ".svg", sep = ""), height = 2.5, width = 6.5)
        }
}


nrow(all_data)

ct = "L6_CT"
motif= "CTCF"



curr_files = files[grep(ct ,files)]
# Load the data into a single data frame
all_data <- lapply(curr_files, function(file) {
  df <- fread(file, header = TRUE)
  df$peak_type <- ifelse(grepl("up_peaks", file), "up", 
                         ifelse(grepl("down_peaks", file), "down", "background"))
  df$region <-  sapply(strsplit(as.character(file), ":"), `[`, 2)
  df$region <-  sapply(strsplit(as.character(df$region), "_"), `[`, 1)

  df$file <- file
  colnames(df)[1] = c("Distance" ) 

  return(df)
}) %>% bind_rows()

all_data = as.data.frame(all_data)

options(repr.plot.width=6.5, repr.plot.height=2.5)

curr = all_data[,c(1,grep(motif, colnames(all_data)),(ncol(all_data)-2):ncol(all_data))] 
colnames(curr)[2:4] = c("combined", "+", "-" )
p = ggplot(curr, aes(x = as.numeric(Distance), y = combined, color = region, group = region)) +
  geom_line() + facet_wrap(~peak_type) + 
  xlim(c(-400,400))+  scale_color_manual(values = color_vector) +
  labs(x = "Distance from center (bp)", y = "Sites per bp per peak", 
       title = paste("Aggregate Plot of ",motif," Motif in ",ct ,sep = "")) +
  theme_minimal()
p
#ggsave(p, file = paste(ct, "_",motif, ".svg", sep = ""), height = 2.5, width = 6.5)

ct = "Glut"
motif= "CTCF(Zf)"



curr_files = files[c(grep("L6_CT" ,files), grep("DG_Glut" ,files), grep("L2-3_IT_ENT_Glut" ,files))]
# Load the data into a single data frame
all_data <- lapply(curr_files, function(file) {
  df <- fread(file, header = TRUE)
  df$peak_type <- ifelse(grepl("up_peaks", file), "up", 
                         ifelse(grepl("down_peaks", file), "down", "background"))
  df$region <-  sapply(strsplit(as.character(file), ":"), `[`, 2)
  df$region <-  sapply(strsplit(as.character(df$region), "_"), `[`, 1)

  df$file <- file
  colnames(df)[1] = c("Distance" ) 

  return(df)
}) %>% bind_rows()

all_data = as.data.frame(all_data)

options(repr.plot.width=6.5, repr.plot.height=2.5)

curr = all_data[,c(1,grep(motif, colnames(all_data), fixed = T),(ncol(all_data)-2):ncol(all_data))] 
colnames(curr)[2:4] = c("combined", "+", "-" )
p = ggplot(curr, aes(x = as.numeric(Distance), y = combined, color = region, group = region)) +
  geom_line() + facet_wrap(~peak_type) + 
  xlim(c(-400,400))+  scale_color_manual(values = color_vector) +
  labs(x = "Distance from center (bp)", y = "Sites per bp per peak", 
       title = paste("Aggregate Plot of ",motif," Motif in ",ct ,sep = "")) +
  theme_minimal()
p
ggsave(p, file = paste(ct, "_",motif, ".svg", sep = ""), height = 2.5, width = 6.5)

p

options(repr.plot.width=6.5, repr.plot.height=2.5)

curr = all_data[,c(1,grep(motif, colnames(all_data)),(ncol(all_data)-2):ncol(all_data))] 
colnames(curr)[2:4] = c("combined", "+", "-" )
p = ggplot(curr, aes(x = as.numeric(Distance), y = combined, color = region, group = region)) +
  geom_line() + facet_wrap(~peak_type) + 
  xlim(c(-400,400))+  scale_color_manual(values = color_vector) +
  labs(x = "Distance from center (bp)", y = "Sites per bp per peak", 
       title = paste("Aggregate Plot of ",motif," Motif in ",ct ,sep = "")) +
  theme_minimal()
p
ggsave(p, file = paste(ct, motif, ".svg", sep = ""), height = 2.5, width = 6.5)

p 

motifs = sapply(strsplit(as.character(colnames(all_data)), "/"), `[`, 1)
motifs = unique(motifs)
motifs = motifs[grep("[(]", motifs)]
motifs

head(all_data)


