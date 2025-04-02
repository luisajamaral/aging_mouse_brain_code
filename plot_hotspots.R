# Chromosomal Hotspot Plotting Script
# Description: Plots aging-associated chromatin hotspots across the genome.

library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)

# ==== 1. Load hotspot data ====
hotspots <- read.csv("hotspot_summary.csv")  # contains chr, start, end, region, clade, score, etc.
h3k9me3 <- read.csv("h3k9me3_overlap.csv")    # contains chr, start, end of regions overlapping H3K9me3

# Add logical column for H3K9me3 overlap
hotspots$h3k9me3_overlap <- with(hotspots, paste(chr, start, end) %in% paste(h3k9me3$chr, h3k9me3$start, h3k9me3$end))

# ==== 2. Plot 1: Split by Region, Color by Clade ====
p1 <- ggplot(hotspots, aes(x = start / 1e6, y = score, color = clade)) +
  geom_point(alpha = 0.7, size = 1.2) +
  facet_wrap(~region, scales = "free_y") +
  theme_classic() +
  labs(x = "Genomic Position (Mb)", y = "Hotspot Score", title = "Aging Hotspots by Region and Clade") +
  scale_color_brewer(palette = "Dark2")

ggsave("hotspot_region_clade_plot.pdf", p1, width = 10, height = 6)

# ==== 3. Plot 2: Highlight H3K9me3 Overlap, Split by Clade ====
p2 <- ggplot(hotspots, aes(x = start / 1e6, y = score)) +
  geom_point(data = hotspots[!hotspots$h3k9me3_overlap, ], color = "gray", alpha = 0.5, size = 1) +
  geom_point(data = hotspots[hotspots$h3k9me3_overlap, ], color = "red", size = 1.2) +
  facet_wrap(~clade, scales = "free_y") +
  theme_minimal() +
  labs(x = "Genomic Position (Mb)", y = "Hotspot Score", title = "Hotspots Overlapping H3K9me3 by Clade")

ggsave("hotspot_h3k9me3_by_clade.pdf", p2, width = 10, height = 6)

# ==== 4. Optional: Save hotspot data with overlap ====
write.csv(hotspots, "hotspots_annotated.csv", row.names = FALSE)
