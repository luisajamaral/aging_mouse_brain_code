library(data.table)
library(ggplot2)
library(dplyr)

library(enrichR)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("Enrichr") # Human genes   
}

setwd("../combined_DARs_redone/")
files = list.files(".", "diff_peaks")


ct = "Oligo_NN"

f =paste( "background_peaks_",ct,"_2vs18_ann_motifs_mscore.txt", sep = "")


anno  = fread(f)
anno$ann = sapply(strsplit(as.character(anno$Annotation), " "), `[`, 1)
anno=as.data.frame(anno)
rownames(anno) = paste(anno$Chr,":",anno$Start-1,"-",anno$End, sep = "")


colnames(anno) = sapply(strsplit(colnames(anno), "[(]"), "[[", 1)

f = paste("diff_peaks_",ct,"_2vs18.csv", sep = "")
cd = fread(f)
cd$`log2(fold_change)` = cd$`log2(fold_change)`


back_pt  = fread(paste("background_peaks_",ct,"_2vs18.csv", sep = ""))
back_pt = as.data.frame(back_pt)
rownames(back_pt) = back_pt[,1]
#cd$pct_2mo = back_pt[paste(cd$`feature name`), 'pct_cells_2mo']
#cd$pct_18mo = back_pt[paste(cd$`feature name`), 'pct_cells_18mo']
#cd$annotation = anno[paste(cd$`feature name`),"ann"]


#cd = cbind(cd, anno[paste(cd$`feature name`),22:35])

head(cd)

anno$pct_2mo = back_pt[rownames(anno), 'pct_cells_2mo']
anno$pct_18mo = back_pt[rownames(anno), 'pct_cells_18mo']


head(anno)

anno$`feature name` = rownames(anno)

# Perform left join to merge cd into anno
merged_data <- anno %>%
  left_join(cd, by = "feature name")

# Rename columns appropriately
merged_data <- merged_data %>%
  rename(adj_pvalue = `adjusted p-value`, logFC = `log2(fold_change)`)

# Ensure missing values are set to default (if needed)
merged_data$adj_pvalue[is.na(merged_data$adj_pvalue)] <- 1
merged_data$logFC[is.na(merged_data$logFC)] <- 0

anno <- merged_data 


# Inspect the resulting anno data frame
head(anno)

max(anno$logFC)

diff = anno 
diff$dir = "ns"
diff$dir[which(diff$adj_pvalue<0.05 & diff$logFC>0.25)] = "up"
diff$dir[which(diff$adj_pvalue<0.05 & diff$logFC<0.25)] = "down"
table(diff$dir)

diff$type = "non-TSS"
diff$type[which(diff$ann=="promoter-TSS")]  = "promoter-TSS"

diff$AP1_binary = "no AP1"
diff$AP1_binary[which(diff$`AP-1`>8)] = "AP1"
table(diff$AP1_binary)

head(diff$ann)



diff$pct = (diff$pct_2mo+diff$pct_18mo)/2
diff$pct_binary = "low accessibility"
diff$pct_binary[which(diff$pct>2)] = "high accessibility"
table(diff$pct_binary)


options(repr.plot.width=8.5, repr.plot.height=4.5)

ggplot(diff, aes(x = ann,y = pct, fill = dir
                )) +ggtitle(paste(ct, "Peaks"))+ 
  geom_boxplot() + ylab("% all cells")+

  theme_minimal()



options(repr.plot.width=6, repr.plot.height=4.5)

ggplot(diff[which(diff$dir != "ns"),], aes(x = ann, fill = dir)) +
  geom_bar(stat = "count", position = "dodge", color = "black") +
  labs(title = paste(ct,"DARs"),
       x = "Peak Type",
       y = "Count",
       fill = "Direction") +
  theme_minimal() + coord_flip() 

head(diff$Sox6)

library(dplyr)
library(ggplot2)

process_and_plot <- function(data, column_name, limit = 2) {
  # Ensure the column exists
  if (!column_name %in% names(data)) {
    stop("Column not found in the data frame.")
  }
  
  # Calculate z-scores
  data <- data %>%
    mutate(z_score = (get(column_name) - mean(get(column_name), na.rm = TRUE)) / sd(get(column_name), na.rm = TRUE))

  # Create binary column based on z-score limit
  binary_col_name <- paste0(column_name, "_binary")
  data[[binary_col_name]] <- paste("no", column_name)
  data[[binary_col_name]][which(data$z_score > limit)] <- column_name

  # Plot histogram with vertical line
  hist(data$z_score, main = paste("Histogram of z-scores for", column_name), xlab = "z-score", col = "lightblue", border = "black")
  abline(v = limit, col = "red", lwd = 2, lty = 2)

  # Display table of binary classification
  print(table(data[[binary_col_name]]))

  # Create bar plot
  bar_plot <- ggplot(data, aes(x = type, fill = dir)) +
    geom_bar(position = "fill") +
    labs(title = paste("Age change in peaks containing", column_name),
         x = "Type",
         y = "Count",
         fill = "Category") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) + facet_wrap(as.formula(paste("~ `", binary_col_name, "`", sep = "")))

  # Set plot size and display the plot
  options(repr.plot.width = 6.5, repr.plot.height = 4)
  print(bar_plot)
}

# Example usage
options(repr.plot.width=5, repr.plot.height=5)
head(diff$`Mef2a`)
process_and_plot(diff, "Mef2a")


process_and_plot(diff, "Sox6")


process_and_plot(diff, "KLF1")


process_and_plot(diff, "AP-1")


head(diff$`Gene Name`)

get_enrich <-function(limit,dir,column_name, diff, max_dist_tss) {
    dbs <- c("GO_Molecular_Function_2023", "GO_Biological_Process_2023", "Reactome_2022","WikiPathway_2023_Human")
    data <- diff %>%
    mutate(z_score = (get(column_name) - mean(get(column_name), na.rm = TRUE)) / sd(get(column_name), na.rm = TRUE))

  # Create binary column based on z-score limit
  binary_col_name <- paste0(column_name, "_binary")
  data[[binary_col_name]] <- paste("no", column_name)
  data[[binary_col_name]][which(data$z_score > limit)] <- column_name

  # Plot histogram with vertical line
  

    limit = min(limit, quantile(data$z_score, .95 ))
    hist(data$z_score, main = paste("Histogram of z-scores for", column_name), xlab = "z-score", col = "lightblue", border = "black")
    abline(v = limit, col = "red", lwd = 2, lty = 2)

    head(data$`Gene Name`[which(data$z_score>limit & diff$dir == dir)])

    #cat(data$`Gene Name`[which(data$z_score>2 & data$dir == "up"  & data$`Distance to TSS`< 2000)], sep = "\n")
    gene_list = data$`Gene Name`[which(data$z_score>limit& data$dir == dir  & data$`Distance to TSS`< max_dist_tss)]

   # cat(length(unique(diff$`Gene Name`)))
    cat(length(unique(gene_list)))

    enriched <- enrichr(gene_list, dbs)
    head(enriched$GO_Molecular_Function_2023)
    options(repr.plot.width=9, repr.plot.height=4.5)

    print(plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",title = paste(ct, dir , column_name, "\nnearby gene GO_Molecular_Function_2023")))
    print(plotEnrich(enriched[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",title = paste(ct, dir , column_name,"\nnearby gene GO_Biological_Process_2023")))
    print(plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",title = paste(ct, dir , column_name,"\nnearby gene Reactome_2022")))
    print(plotEnrich(enriched[[4]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",title = paste(ct, dir , column_name,"\nnearby gene WikiPathway_2023_Human")))

}

get_enrich(2,"up", "CTCF", diff, 5000)

get_enrich(1,"down", "Sox6", diff, 2500)

get_enrich(1.5,"up", "AP-1", diff, 5000)

get_enrich(2,"up", "CTCF", diff, 10000)

get_enrich(2,"down", "Mef2a", diff, 10000)

get_enrich(2,"up", "KLF1", diff, 10000)

dbs <- c("GO_Molecular_Function_2023", "GO_Biological_Process_2023", "Reactome_2022","WikiPathway_2023_Human")

enriched <- enrichr(gene_list, dbs)
head(enriched$GO_Molecular_Function_2023)

 plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

head(enriched$GO_Biological_Process_2023)

head(enriched$Reactome_2022)

head(enriched$WikiPathway_2023_Human)

process_and_plot(diff, "AP-1", limit =1.5)


process_and_plot(diff, "JunB")


head(diff$

library(tidyr)

# Define a function to find modes in a vector and return a maximum of two
find_modes <- function(x) {
  density_est <- density(x, na.rm = TRUE)
  peaks <- which(diff(sign(diff(density_est$y))) == -2)
  mode_values <- density_est$x[peaks]
  return(head(mode_values, 3))  # Return only the first two modes if they exist
}
# Calculate modes for each group
modes <- diff %>%
  group_by(type, dir) %>%
  summarise(modes = list(find_modes(`CTCF`))) %>%
  unnest(cols = c(modes))

# Create the violin plot with mode points
ggplot(diff, aes(x = type, y = `CTCF`, fill = dir)) +
  geom_violin(position = position_dodge(0.9)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", color = "black", position = position_dodge(0.9)) +
  geom_point(data = modes, aes(x = type, y = modes, color = dir), position = position_dodge(0.9), size = 2, shape = 21, fill = "white") +
  labs(title = paste(ct, "CTCF by Type and Dir"),
       x = "Type",
       y = "CTCF",
       fill = "Direction",
       color = "Direction") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) 

# Define a function to find modes in a vector and return a maximum of two
find_modes <- function(x) {
  density_est <- density(x, na.rm = TRUE)
  peaks <- which(diff(sign(diff(density_est$y))) == -2)
  mode_values <- density_est$x[peaks]
  return(head(mode_values, 3))  # Return only the first two modes if they exist
}
# Calculate modes for each group
modes <- diff %>%
  group_by(type, dir) %>%
  summarise(modes = list(find_modes(`Sox6`))) %>%
  unnest(cols = c(modes))

# Create the violin plot with mode points
ggplot(diff, aes(x = type, y = `Sox6`, fill = dir)) +
  geom_violin(position = position_dodge(0.9)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", color = "black", position = position_dodge(0.9)) +
  geom_point(data = modes, aes(x = type, y = modes, color = dir), position = position_dodge(0.9), size = 2, shape = 21, fill = "white") +
  labs(title = paste(ct, "Sox6 by Type and Dir"),
       x = "Type",
       y = "Sox6",
       fill = "Direction",
       color = "Direction") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) 

# Define a function to find modes in a vector and return a maximum of two
find_modes <- function(x) {
  density_est <- density(x, na.rm = TRUE)
  peaks <- which(diff(sign(diff(density_est$y))) == -2)
  mode_values <- density_est$x[peaks]
  return(head(mode_values, 3))  # Return only the first two modes if they exist
}
# Calculate modes for each group
modes <- diff %>%
  group_by(type, dir) %>%
  summarise(modes = list(find_modes(`JunB`))) %>%
  unnest(cols = c(modes))

# Create the violin plot with mode points
ggplot(diff, aes(x = type, y = `JunB`, fill = dir)) +
  geom_violin(position = position_dodge(0.9)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", color = "black", position = position_dodge(0.9)) +
  geom_point(data = modes, aes(x = type, y = modes, color = dir), position = position_dodge(0.9), size = 2, shape = 21, fill = "white") +
  labs(title = paste(ct, "JunB by Type and Dir"),
       x = "Type",
       y = "JunB",
       fill = "Direction",
       color = "Direction") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) 

# Define a function to find modes in a vector and return a maximum of two
find_modes <- function(x) {
  density_est <- density(x, na.rm = TRUE)
  peaks <- which(diff(sign(diff(density_est$y))) == -2)
  mode_values <- density_est$x[peaks]
  return(head(mode_values, 3))  # Return only the first two modes if they exist
}
# Calculate modes for each group
modes <- diff %>%
  group_by(type, dir) %>%
  summarise(modes = list(find_modes(`AP-1`))) %>%
  unnest(cols = c(modes))

# Create the violin plot with mode points
ggplot(diff, aes(x = type, y = `AP-1`, fill = dir)) +
  geom_violin(position = position_dodge(0.9)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", color = "black", position = position_dodge(0.9)) +
  geom_point(data = modes, aes(x = type, y = modes, color = dir), position = position_dodge(0.9), size = 2, shape = 21, fill = "white") +
  labs(title = paste(ct, "AP-1 by Type and Dir"),
       x = "Type",
       y = "AP-1",
       fill = "Direction",
       color = "Direction") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) 

# Define a function to find modes in a vector and return a maximum of two
find_modes <- function(x) {
  density_est <- density(x, na.rm = TRUE)
  peaks <- which(diff(sign(diff(density_est$y))) == -2)
  mode_values <- density_est$x[peaks]
  return(head(mode_values, 3))  # Return only the first two modes if they exist
}
# Calculate modes for each group
modes <- diff %>%
  group_by(type, dir) %>%
  summarise(modes = list(find_modes(`Mef2a`))) %>%
  unnest(cols = c(modes))

# Create the violin plot with mode points
ggplot(diff, aes(x = type, y = `Mef2a`, fill = dir)) +
  geom_violin(position = position_dodge(0.9)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", color = "black", position = position_dodge(0.9)) +
  geom_point(data = modes, aes(x = type, y = modes, color = dir), position = position_dodge(0.9), size = 2, shape = 21, fill = "white") +
  labs(title = paste(ct, "Mef2a by Type and Dir"),
       x = "Type",
       y = "Mef2a",
       fill = "Direction",
       color = "Direction") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) 

# Define a function to find modes in a vector and return a maximum of two
find_modes <- function(x) {
  density_est <- density(x, na.rm = TRUE)
  peaks <- which(diff(sign(diff(density_est$y))) == -2)
  mode_values <- density_est$x[peaks]
  return(head(mode_values, 3))  # Return only the first two modes if they exist
}
# Calculate modes for each group
modes <- diff %>%
  group_by(type, dir) %>%
  summarise(modes = list(find_modes(`KLF1`))) %>%
  unnest(cols = c(modes))

# Create the violin plot with mode points
ggplot(diff, aes(x = type, y = `KLF1`, fill = dir)) +
  geom_violin(position = position_dodge(0.9)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", color = "black", position = position_dodge(0.9)) +
  geom_point(data = modes, aes(x = type, y = modes, color = dir), position = position_dodge(0.9), size = 2, shape = 21, fill = "white") +
  labs(title = paste(ct, "KLF1 by Type and Dir"),
       x = "Type",
       y = "KLF1",
       fill = "Direction",
       color = "Direction") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) 

ggplot(diff, aes(x = type, y = `JunB`, fill = dir)) +
  geom_violin(position = position_dodge(0.9)) +
  #geom_jitter(aes(color = dir), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), alpha = 0.5) +  # Add jittered points within corresponding violins
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", color = "black", position = position_dodge(0.9)) +  # Add summary statistics
  labs(title = paste(ct, "AP-1 by Type and Dir"),
       x = "Type",
       y = "AP-1(bZIP) Odds ratio",
       fill = "Direction",
       color = "Direction") +  # Add labels and title
  theme_minimal() +  # Apply minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center title
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggplot(diff, aes(x = type, y = `CTCF`, fill = dir)) +
  geom_violin(position = position_dodge(0.9)) +
  #geom_jitter(aes(color = dir), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), alpha = 0.5) +  # Add jittered points within corresponding violins
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", color = "black", position = position_dodge(0.9)) +  # Add summary statistics
  labs(title = "Violin Plot of AP-1 odds ratio by Type and Dir",
       x = "Type",
       y = "AP-1(bZIP) Odds ratio",
       fill = "Direction",
       color = "Direction") +  # Add labels and title
  theme_minimal() +  # Apply minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center title
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + facet_wrap(~pct_binary)

options(repr.plot.width=9, repr.plot.height=4.5)

ggplot(diff, aes(x = type,y = `NeuroD1(bHLH)`, fill = dir
                )) +
  geom_boxplot() +

  theme_minimal()


setwd("../combined_DARs_redone/")
files = list.files(".", "diff_peaks")

all_df  = list()
for(f in files){
    cd = fread(f)
    cd$`log2(fold_change)` = cd$`log2(fold_change)`
    ct = gsub("diff_peaks_", "", f)
    ct = gsub("_2vs18.csv", "", ct)
    anno  = fread(paste("background_peaks_",ct,"_2vs18_ann.txt", sep = ""))
    anno$ann = sapply(strsplit(as.character(anno$Annotation), " "), `[`, 1)
    anno=as.data.frame(anno)
    rownames(anno) = paste(anno$Chr,":",anno$Start-1,"-",anno$End, sep = "")
    back_pt  = fread(paste("background_peaks_",ct,"_2vs18.csv", sep = ""))
    back_pt = as.data.frame(back_pt)
    rownames(back_pt) = back_pt[,1]
    cd$pct_2mo = back_pt[paste(cd$`feature name`), 'pct_cells_2mo']
    cd$pct_18mo = back_pt[paste(cd$`feature name`), 'pct_cells_18mo']
    cd$annotation = anno[paste(cd$`feature name`),"ann"]
    cd$annotation = anno[paste(cd$`feature name`),"AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer Best Motif log-odds Score"]
    cat(ct)

    #cat(cd$`feature name`[which(! paste(cd$`feature name`) %in% rownames(anno))])
    cd$celltype = ct
    #head(cd)
    all_df[[ct]] = cd
}
    

head(diff)

write.table(diff, file  = "~/projects/combined_all/Figures/Figure3-NN/Oligo_annotated_peaks.txt", sep = "\t")


