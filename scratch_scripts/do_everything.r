library(data.table)

setwd("../histone_overlap/")

path = "/projects/ps-renlab2/szu/projects/amb_pairedtag/04.peakcalling/out/reproPeak_nolambda"


cts = c("DG_Glut", "L2_3_IT_CTX_Glut"  , "L2_3_IT_PPP_Glut", "L5_IT_CTX_Glut", 
        "Microglia_NN" , "Astro_TE_NN","Astro_NT_NN", "CA1_ProS_Glut", "MOL_NN","COP_NN","NFOL_NN", "MFOL_NN", 
        "STR_D1_Sema5a_Gaba", "STR_D1_Gaba", "STR_PAL_Chst9_Gaba", "LA_BLA_BMA_PA_Glut",
       "OPC_NN" ,"L6b_CT_ENT_Glut", "Vip_Gaba")

for(ct in cts) {
    files = list.files(path, ct, full.names = T) 
    for(f in files){
        system(paste("cp", f, "all_histone/"))
    }
}


files

files

for(ct in cts) {
    files = list.files(path, ct, full.names = T) 
    for(f in files){
        t = fread(f, header = T)
        head(t)
        t = t[which(t$ScorePerMillion>3),]
       # c = gsub('/projects/ps-renlab2/szu/projects/amb_pairedtag/04.peakcalling/out/reproPeak_nolambda/', '',f)
        c= gsub('.reproPeak', '',f)
       # c= gsub('_5-', '-',c)
       # c= gsub('_4-', '-',c)
       # c= gsub('_2-', '-',c)
       # c= gsub('_3-', '-',c)
       # c= gsub('_1-', '-',c)
        #c= gsub('MOL', 'Oligo',c)
       # c= gsub('D1', 'D12',c)


        write.table(t[,c(1,2,3)], file = paste(c,".bed", sep = ""),row.names = F, sep = "\t", quote = F, col.names = F)
    }
}


path = "../../combined_DARs_redone"

cts = c("DG_Glut", "L2-3_IT_CTX_Glut"  , "L2-3_IT_PPP_Glut", "L5_IT_CTX_Glut", "Microglia_NN" , "Astro-TE_NN",  
        "Astro-NT_NN","CA3_Glut", "CA1-ProS_Glut","Oligo_NN", "STR_D12_Gaba", "STR_D1_Sema5a_Gaba","STR-PAL_Chst9_Gaba",
       "LA-BLA-BMA-PA_Glut", "OPC_NN","L6b-CT_ENT_Glut", "Vip_Gaba")



for(ct in cts) {
    files = list.files(path, ct, full.names = T) 
    files = files[grep(".bed", files)]
    if(length(grep("ann", files))>0){
        files = files[-grep("ann", files)]
    }
    if(length(grep("csv", files))>0){
        files = files[-grep("csv", files)]
    }

    for(f in files){
        t = fread(f, header = T)
        head(t)
        c = gsub('../combined_DARs_redone/', '',f)
        c = gsub('.bed', '',c)
        c = gsub('_2vs18', '',c)
        c = gsub('-', '_',c)
        c = gsub('__', '_',c)

        cat(c)
        t = t[which(
        if(length(grep("background", c))==0){
            write.table(t[1:min(5000,nrow(t)),c(1,2,3)], file = paste(c,"_5k.bed", sep = ""),row.names = F, sep = "\t", quote = F, col.names = F)
        } else{
            write.table(t[,c(1,2,3)], file = paste(c,".bed", sep = ""),row.names = F, sep = "\t", quote = F, col.names = F)

            }
    }
}




setwd("all_histone/")

# Load required libraries
library(ggplot2)
library(dplyr)

# Read the data and set column names explicitly
col_names <- c("celltype", "histone", "type", "overlap", "total_peaks", "total_merged")
overlap_results <- read.csv("overlap_results_5k.csv", header = FALSE, col.names = col_names)


head(overlap_results)

options(repr.plot.width=7.5 , repr.plot.height=5
       )

plist = list()

# Iterate over each histone mark
for (hist in c("H3K9me3", "H3K4me1", "H3K27me3", "H3K27ac")) {
  
  # Filter for specific histone
  hist_results <- overlap_results %>%
    filter(histone == hist)
  
  # Separate up, down, and background data
  hist_up <- hist_results %>% filter(type == "up")
  hist_down <- hist_results %>% filter(type == "down")
  hist_background <- hist_results %>% filter(type == "background")
  
  # Fisher's exact test function (keeping this to calculate significance)
  run_fishers_test <- function(overlap, total_peaks, overlap_background, total_background) {
    contingency_table <- matrix(c(overlap, total_peaks - overlap, 
                                  overlap_background, total_background - overlap_background), 
                                nrow = 2)
    fisher_result <- fisher.test(contingency_table)
    return(c(fisher_result$p.value))
  }
  
  # Calculate fold change and Fisher's test for up vs background and down vs background
  hist_up$fold_change <- (hist_up$overlap / hist_up$total_peaks) / (hist_background$overlap / hist_background$total_peaks)
  hist_down$fold_change <- (hist_down$overlap / hist_down$total_peaks) / (hist_background$overlap / hist_background$total_peaks)
  
  # Log2 transform the fold change
  hist_up$log2_fold_change <- log2(hist_up$fold_change)
  hist_down$log2_fold_change <- log2(hist_down$fold_change)
  
  # Run Fisher's exact test to get p-values
  up_p_values <- mapply(run_fishers_test, 
                        hist_up$overlap, 
                        hist_up$total_peaks, 
                        hist_background$overlap, 
                        hist_background$total_peaks)
  
  down_p_values <- mapply(run_fishers_test, 
                          hist_down$overlap, 
                          hist_down$total_peaks, 
                          hist_background$overlap, 
                          hist_background$total_peaks)
  
  # Assign p-values and significance stars
  hist_up$p_value <- up_p_values
  hist_down$p_value <- down_p_values

  # Combine up and down results for plotting
  hist_combined <- bind_rows(hist_up, hist_down)
  hist_combined$adj_p_value <- p.adjust(hist_combined$p_value, method = "BH")

  hist_combined$stars = cut(hist_combined$adj_p_value, 
                         breaks = c(-Inf, 1e-30, 1e-20, 1e-10, Inf), 
                         labels = c("***", "**", "*", ""))
    

  # Add cell type category (Gaba, Glut, NN)
  hist_combined$category <- ifelse(grepl("Gaba", hist_combined$celltype), "Gaba",
                                   ifelse(grepl("Glut", hist_combined$celltype), "Glut", "NN"))
  hist_combined$log2_fold_change = as.numeric(hist_combined$log2_fold_change)
  # Traditional log fold change plot using geom_point
  p <- ggplot(hist_combined, aes(x = celltype, y = fold_change, color = type)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    theme_minimal() +
    labs(title = paste("Fold Change of", hist, "Overlap vs background"),
         x = "Cell Type", y = "Fold Change") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # Reference line at log2(fold change) = 0
    theme(axis.text.x = element_text(angle = 65, hjust = 1),text = element_text(size = 16),plot.title = element_text(size = 14) ) +
    facet_grid(category~. , scales = "free_y") +  # Facet by cell type and "up" or "down"
    ylim(c(0,5))+
    scale_color_manual(values=c( "blue", "red")) + coord_flip()+ 

    geom_text(aes(label = stars), 
              position = position_dodge(width = 0.5), 
              vjust = -0.5, size = 3)  # Add significance stars
  
  print(p)
  plist[[hist]] = p
}


head(hist_combined)

options(repr.plot.width=7.5 , repr.plot.height=5
       )

plist = list()

# Iterate over each histone mark
for (hist in c("H3K9me3", "H3K4me1", "H3K27me3", "H3K27ac")) {
  
  # Filter for specific histone
  hist_results <- overlap_results %>%
    filter(histone == hist)
  
  # Separate up, down, and background data
  hist_up <- hist_results %>% filter(type == "up")
  hist_down <- hist_results %>% filter(type == "down")
  hist_background <- hist_results %>% filter(type == "background")
  
  # Fisher's exact test function (keeping this to calculate significance)
  run_fishers_test <- function(overlap, total_peaks, overlap_background, total_background) {
    contingency_table <- matrix(c(overlap, total_peaks - overlap, 
                                  overlap_background, total_background - overlap_background), 
                                nrow = 2)
    fisher_result <- fisher.test(contingency_table)
    return(c(fisher_result$p.value))
  }
  
  # Calculate fold change and Fisher's test for up vs background and down vs background
  hist_up$fold_change <- (hist_up$overlap / hist_up$total_peaks) / (hist_background$overlap / hist_background$total_peaks)
  hist_down$fold_change <- (hist_down$overlap / hist_down$total_peaks) / (hist_background$overlap / hist_background$total_peaks)
  
  # Log2 transform the fold change
  hist_up$log2_fold_change <- log2(hist_up$fold_change)
  hist_down$log2_fold_change <- log2(hist_down$fold_change)
  
  # Run Fisher's exact test to get p-values
  up_p_values <- mapply(run_fishers_test, 
                        hist_up$overlap, 
                        hist_up$total_peaks, 
                        hist_background$overlap, 
                        hist_background$total_peaks)
  
  down_p_values <- mapply(run_fishers_test, 
                          hist_down$overlap, 
                          hist_down$total_peaks, 
                          hist_background$overlap, 
                          hist_background$total_peaks)
  
  # Assign p-values and significance stars
  hist_up$p_value <- up_p_values
  hist_down$p_value <- down_p_values

  # Combine up and down results for plotting
  hist_combined <- bind_rows(hist_up, hist_down)
  hist_combined$adj_p_value <- p.adjust(hist_combined$p_value, method = "BH")

  hist_combined$stars = cut(hist_combined$adj_p_value, 
                         breaks = c(-Inf, 1e-30, 1e-20, 1e-10, Inf), 
                         labels = c("***", "**", "*", ""))
    

  # Add cell type category (Gaba, Glut, NN)
  hist_combined$category <- ifelse(grepl("Gaba", hist_combined$celltype), "Gaba",
                                   ifelse(grepl("Glut", hist_combined$celltype), "Glut", "NN"))
  hist_combined$log2_fold_change = as.numeric(hist_combined$log2_fold_change)
  # Traditional log fold change plot using geom_point
  p <- ggplot(hist_combined, aes(x = celltype, y = fold_change, fill = type, color=type)) +
    geom_bar(position = position_dodge(width = 1), stat = "identity" ) +
    theme_minimal() +
    labs(title = paste( hist, "Overlap"),
         x = "Cell Type", y = "Fold Change (vs all bg peaks)") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # Reference line at log2(fold change) = 0
    theme(axis.text.x = element_text(angle = 65, hjust = 1),text = element_text(size = 16),plot.title = element_text(size = 14) ) +
    facet_grid(category~. , scales = "free_y", space = "free") +  # Facet by cell type and "up" or "down"
    ylim(c(0,25))+
    scale_color_manual(values=c( "blue", "red"),name = "Age-diff (18mo/2mo)") +
    scale_fill_manual(values=c( "blue", "red"),name = "Age-diff (18mo/2mo)") + coord_flip()+ 

    geom_text(aes(label = stars), 
              position = position_dodge(width = 0.5), 
              vjust = 0.5, hjust = -.25,size = 3)  # Add significance stars
  
  print(p)
  plist[[hist]] = p
}


options(repr.plot.width=25 , repr.plot.height=7
       )
library(ggplot2)
library(patchwork)

plot1 <- plist[[1]]
plot2 <- plist[[2]]

# Combine the plots side by side
combined_plot <- plot1+ theme(legend.position = "none", strip.text = element_blank()) + plot2+ theme(legend.position = "none", strip.text = element_blank(),axis.title.y = element_blank(), 
                            axis.text.y = element_blank())  + plist[[3]]+ theme(legend.position = "none", strip.text = element_blank(),axis.title.y = element_blank(), 
                            axis.text.y = element_blank()) +plist[[4]]+ theme(axis.title.y = element_blank(), 
                            axis.text.y = element_blank())

# Display the combined plot
combined_plot+ plot_layout(nrow = 1)

pdf("~/projects/combined_all/Figures/Figure6-Het-TEs/histone_overlap_bar.pdf", width = 21, height = 6)
combined_plot+ plot_layout(nrow = 1)

dev.off()

options(repr.plot.width=25 , repr.plot.height=7
       )
library(ggplot2)
library(patchwork)

plot1 <- plist[[1]]
plot2 <- plist[[2]]

# Combine the plots side by side
combined_plot <- plot1+ theme(legend.position = "none", strip.text = element_blank()) + plot2+ theme(legend.position = "none", strip.text = element_blank(),axis.title.y = element_blank(), 
                            axis.text.y = element_blank())  + plist[[3]]+ theme(legend.position = "none", strip.text = element_blank(),axis.title.y = element_blank(), 
                            axis.text.y = element_blank()) +plist[[4]]+ theme(axis.title.y = element_blank(), 
                            axis.text.y = element_blank())

# Display the combined plot
combined_plot+ plot_layout(nrow = 1)


combined_plot+ plot_layout(nrow = 1)


# Plot the overlap results for H3K9me3
ggplot(h3k9me3_results, aes(x = celltype, y = fraction, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Overlap of Peaks with H3K9me3 by Cell Type", 
       x = "Cell Type", y = "Number of Overlaps") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

h3k9me3_combined


