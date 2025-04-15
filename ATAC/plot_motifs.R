# === Motif Enrichment Plotting Script (Cleaned for Submission) ===
# Description: Plot HOMER motif enrichment results (up- and down-regulated DARs)
# Output: PDF dot plots, colored by -log10(p-value), sized by % target sequences

library(tidyverse)
library(ggplot2)
library(stringr)
library(ggrepel)
library(gridExtra)

# === Load Helper Function ===
load_motif_file <- function(file_path, direction) {
  df <- read_tsv(file_path, col_names = TRUE)
  df$Direction <- direction
  df <- df %>% 
    rename(
      Motif = `Motif Name`,
      Consensus = `Consensus`,
      PValue = `P-value`,
      QValue = `q-value (Benjamini)`,
      TargetPerc = `% of Target Sequences with Motif`
    ) %>%
    mutate(
      Motif = str_remove(Motif, "\t.*"),
      CellType = str_extract(file_path, "(?<=motifs/).*(?=_up|_down)"),
      LogP = -log10(PValue)
    )
  return(df)
}

# === Define Motif Lists and Input Files ===
motif_dir <- "motifs/"
up_files <- list.files(motif_dir, pattern = "_up\.motif", full.names = TRUE)
down_files <- list.files(motif_dir, pattern = "_down\.motif", full.names = TRUE)

# === Combine All Files ===
motif_up <- map_dfr(up_files, load_motif_file, direction = "Up")
motif_down <- map_dfr(down_files, load_motif_file, direction = "Down")
motif_all <- bind_rows(motif_up, motif_down)

# === Filter for P-value Threshold ===
motif_all <- motif_all %>% 
  filter(PValue < 0.01) %>%
  mutate(TargetPerc = as.numeric(str_replace(TargetPerc, "%", "")))

# === Select Key Motifs to Plot ===
selected_motifs <- c("Sox2", "Sox9", "Olig2", "Lhx2", "Lhx6", "NF1", "Jun", "Fos", "Atf3", "Egr1", "Mef2c", "Mef2d", "Hes5")
motif_all <- motif_all %>% filter(Consensus %in% selected_motifs)

# === Plotting Function ===
plot_motif_dotplot <- function(df, direction, file_out) {
  ggplot(df %>% filter(Direction == direction), 
         aes(x = CellType, y = Consensus, size = TargetPerc, color = LogP)) +
    geom_point() +
    facet_wrap(~ Direction) +
    scale_color_gradient(low = "gray90", high = "red3") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    labs(
      title = paste("Motif Enrichment:", direction),
      x = "Cell Type",
      y = "Motif",
      size = "% Target Sequences",
      color = "-log10(p-value)"
    )
  ggsave(file_out, width = 10, height = 6)
}

# === Output Plots ===
plot_motif_dotplot(motif_all, direction = "Up", file_out = "motif_up_final.pdf")
plot_motif_dotplot(motif_all, direction = "Down", file_out = "motif_down_final.pdf")
