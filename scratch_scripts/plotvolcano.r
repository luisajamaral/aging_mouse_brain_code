library(ggplot2)
library(ggrepel)

setwd("TEs_only_MAST_latent_Frozen/")

options(repr.plot.width=6, repr.plot.height=6)

files = list.files(".", ".csv")
files = files[grep("L5", files)]

for (f in files){ 
    
    cur = gsub(".csv", "", f)
    cat(cur)
    data=read.csv(f)
    head(data)
    data$gene = data$X
    # Calculate -log10 of the p-values for plotting
    data$log_p_val <- -log10(data$p_val)

    # Define significance thresholds
    p_value_threshold <- -log10(1e-10)
    logFC_threshold <- 0.2

    # Add a column to indicate significance and direction
    data$significance <- with(data, ifelse(log_p_val > p_value_threshold & avg_log2FC > logFC_threshold, "Upregulated",
                                  ifelse(log_p_val > p_value_threshold & avg_log2FC < -logFC_threshold, "Downregulated", "Not Significant")))

    # Create the volcano plot
    g = ggplot(data, aes(x = avg_log2FC, y = log_p_val)) +
      geom_point(aes(color = significance), alpha = 0.7) +
      # Define colors for points
      scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
      # Add labels to selected points
      geom_text_repel(data = subset(data, significance != "Not Significant"), 
                      aes(label = gene), 
                      size = 3, 
                      max.overlaps = 10) +
      # Custom aesthetics
      theme_minimal() +
      labs(title = cur,
           x = "Log2 Fold Change",
           y = "-Log10 P-value",
           color = "Significance") +
      theme(plot.title = element_text(hjust = 0.5))
    print(g)
    
    }


