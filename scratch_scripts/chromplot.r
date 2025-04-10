library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(org.Mm.eg.db)
library(viridis)
library(ggrepel)
library(AnnotationDbi)
library(biomaRt)
library(GenomicRanges)

setwd("../../region_DARs_redone/smooth/")


h3k9me3 = read.table("../../histone_volcano/H3K9me3_reprocessed_P0_forebrain.bed", header = F)

colnames(h3k9me3) = c("GeneID", "chr", "start", "end", "strand", "score")

nrow(h3k9me3)

h3k9me3 = h3k9me3[which(h3k9me3$score>5000),]

h3k9_gr <- GRanges(seqnames = h3k9me3$chr,
                   ranges = IRanges(start = h3k9me3$start, end = h3k9me3$end))
h3k9_gr

setwd("../../../../h5ads_final/region_sex_DARs/diff_csvs/smooth/")




setwd("../../../region_DARs_redone/smooth/")


head(smoothed_all)

head(smoothed_up)

library(ggplot2)
library(dplyr)
library(data.table)

# Load full smoothed scores
smoothed_up <- fread("./smth20_scores_table_top1k_1logfcUp.txt")
smoothed_down <- fread("./smth20_scores_table_top1k_1logfcDown.txt")

# Add direction labels
smoothed_up$direction <- "Down"
smoothed_down$direction <- "Up"
smoothed_up$smth20 <- -smoothed_up$smth20

# Combine both datasets
smoothed_all <- rbind(smoothed_up, smoothed_down)
head(smoothed_all)
# Identify the first and last row for each chromosome
first_last <- smoothed_all %>%
  group_by(chr) %>%
  filter(row_number() == 1 | row_number() == n()) %>%
  ungroup()

# Filter rows: Keep first/last row for each chromosome, otherwise remove rows where smth20 == 0
smoothed_all <- smoothed_all %>%
  filter(abs(smth20) > 0 | chr %in% first_last$chr & pos %in% first_last$pos)

# Assign cell-type clades
smoothed_all$clade <- "NN"
smoothed_all$clade[grep("Glut", smoothed_all$cluster_id)] <- "Glut"
smoothed_all$clade[grep("Gaba", smoothed_all$cluster_id)] <- "Gaba"
smoothed_all$clade[grep("Neur", smoothed_all$cluster_id)] <- "Gaba"

# Make low values gray
smoothed_all$clade[which(abs(smoothed_all$smth20) < 0.2)] <- "ns"

# Set chromosome ordering
smoothed_all$chr <- factor(smoothed_all$chr, levels=paste0("chr", c(1:19, "X", "Y")))

smoothed_all$sex = "Male"
smoothed_all$sex[grep("Female" , smoothed_all$cluster_id) ] = "Female"

# Plot
g <- ggplot(smoothed_all, aes(x = pos, y = smth20, color = clade)) +
  geom_point(alpha = .8, size = .8) +
  
  scale_color_manual(values = c("ns" = "gray70", 
                                "Glut" = "firebrick1", 
                                "Gaba" = "blue3", 
                                "NN" = "green2")) + 
  
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + 
  
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour="grey")) +
  
  guides(colour = guide_legend(override.aes = list(size=2))) +
  facet_grid( direction ~ chr, scales = "free_x", space = "fixed") + 
  ggtitle("Chromosomal Hotspots")

print(g)


options(repr.plot.width=5, repr.plot.height=3)


g <- ggplot(smoothed_all[which(smoothed_all$chr%in%c("chrY", "chrX")),], aes(x = pos, y = smth20, color = clade)) +
  geom_point(alpha = .8, size = .8) +
  
  scale_color_manual(values = c("ns" = "gray70", 
                                "Glut" = "firebrick1", 
                                "Gaba" = "blue3", 
                                "NN" = "green2")) + 
  
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + 
  
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour="grey")) +
  
  guides(colour = guide_legend(override.aes = list(size=2))) +
  facet_grid(sex ~ chr, scales = "free_x", space = "free") + 
  ggtitle("Sex Chromosomes Hotspots")

print(g)

#ggsave(g , file = "../../../../Figures/Figure6-Het-TEs/sex_chrom.pdf", height = 3, width = 5)

head(smoothed_all)

smoothed_all$region = sapply(strsplit(as.character(smoothed_all$cluster_id), ":"), `[`, 3)
smoothed_all$region = sapply(strsplit(as.character(smoothed_all$region), "_"), `[`, 1)
smoothed_all$region = factor(smoothed_all$region, levels = c("FC", "ENT", "HCA", "HCP", "AMY", "NAC","RLP",  "CP"))

options(repr.plot.width=7, repr.plot.height=5)


g <- ggplot(smoothed_all[which(smoothed_all$chr%in%c("chrY", "chrX")),], aes(x = pos, y = smth20, color = clade)) +
  geom_point(alpha = .25, size = 1) +
  
  scale_color_manual(values = c("ns" = "gray70", 
                                "Glut" = "firebrick1", 
                                "Gaba" = "blue3", 
                                "NN" = "green2")) + 
  
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + 
  
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour="grey")) +
  
  guides(colour = guide_legend(override.aes = list(size=2))) +
  facet_grid(region ~ chr + sex, scales = "free_x", space = "free") + 
  ggtitle("Sex Chromosomes Hotspots")

print(g)

#ggsave(g , file = "../../../../Figures/Figure6-Het-TEs/sex_chrom.pdf", height = 3, width = 5)


x[which(x$cluster_id == "Oligo_NN:Female:CP"),]

x = smoothed_all[order(smoothed_all$smth20, decreasing = T),]
head(x[which(x$chr == "chrX"),], 5)

x$cluster_id = gsub("_2vs18_iter.csv" , "",x$cluster_id) 
table(x[which(x$chr == "chrX" & x$smth20>0.5),'cluster_id'])

options(repr.plot.width=10, repr.plot.height=8)

g <- ggplot(smoothed_all[,], aes(x = pos, y = smth20, color = clade)) +
  geom_point(alpha = .8, size = .8) +
  
  scale_color_manual(values = c("ns" = "gray70", 
                                "Glut" = "firebrick1", 
                                "Gaba" = "blue3", 
                                "NN" = "green2")) + 
  
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + 
  
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour="grey")) +
  
  guides(colour = guide_legend(override.aes = list(size=2))) +
  facet_grid(sex ~ chr, scales = "free_x", space = "fixed") + 
  ggtitle("Chromosomal Hotspots")

print(g)

head(smoothed_all)


smoothed_all$clade[which(abs(smoothed_all$smth20)<0.2)] = "ns"

options(repr.plot.width=10, repr.plot.height=8)

g <- ggplot() +
  # Background: Plot full chromosome smoothed scores (light gray)
  geom_line(data = smoothed_all , 
            aes(x = midpoint, y = smth20, group = chr), 
            color = "gray70", size = 0.5, alpha = 0.5) +
  
  # Overlay: Highlight hotspot points
  geom_point(data = smoothed_all, 
             aes(x = midpoint, y = smth20_score, color = clade), 
             alpha = .8, size = .8) + 
  facet_grid(sex ~ chr, scales = "free_x", space = "fixed")
  
  scale_color_manual(values = c("green2", "firebrick1", "blue3")) + 
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + ylim(c(-1.5,1.5)) +
  
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour="grey")) +
  
  guides(colour = guide_legend(override.aes = list(size=2))) +
  facet_grid(sex ~ chr, scales = "free_x", space = "fixed") + 
  ggtitle("Chromosomal Hotspots with Full Representation")

print(g)

files <- list.files(path=".", pattern="smth20_scores_table_top1k_", full.names=TRUE)

down = fread(files[1])
up = fread(files[2])
down$direction = "Down"
up$direction = "Up"
out = rbind(up,down)

out <- out %>%
  mutate(
    chr = sub("(.*):(.*)-(.*)", "\\1", locations),
    start = as.numeric(sub("(.*):(.*)-(.*)", "\\2", locations)),
    end = as.numeric(sub("(.*):(.*)-(.*)", "\\3", locations)),
    midpoint = (start + end) / 2
)
out$region = sapply(strsplit(as.character(out$clust), ":"), `[`, 2)
out$region = factor(out$region, levels = c("FC", "ENT", "HCA", "HCP", "AMY", "NAC","RLP",  "CP"))

out$clade = "NN"
out$clade[grep("Glut", out$clust)] = "Glut"
out$clade[grep("Gaba", out$clust)] = "Gaba"
out$clade[grep("Neur", out$clust)] = "Gaba"

out$direction = factor(out$direction  , levels = c("Up", "Down"))

out$smth20_score[which(out$direction=="Down")] = -out$smth20_score[which(out$direction=="Down")]
out$chr = factor(out$chr, levels=paste0("chr",c(1:19,"X","Y")))
h3k9me3$chr = factor(h3k9me3$chr, levels=paste0("chr",c(1:19,"X","Y")))
options(repr.plot.width=12, repr.plot.height=8)
options(repr.plot.width=10, repr.plot.height=8)
out$smth20_score = -out$smth20_score

#out$clade = factor(out$clade, levels = c("NN", "Gaba", "Glut"))
g = ggplot(out[which(! out$chr %in% c("chrX", "chrY")),]) +
  geom_point(aes(x=midpoint,y=smth20_score,color=clade),alpha=.8,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
#scale_color_manual(values=c( brewer.pal(n = 3, name = "Dark2")[1],brewer.pal(n = 3, name = "Dark2")[3],brewer.pal(n = 3, name = "Dark2")[2]
                   # ))+ 
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + ylim(c(-2.5,2.5))+
 # Add annotation bar for regions
  # geom_rect(data = h3k9me3[! which(h3k9me3$chr %in% c("chrX", "chrY")),], aes(
  #   xmin = start-100000, xmax = end+100000,fill = (score),
  #   ymin = -4, ymax = -3  # Position the bar beneath the plot
  # ), alpha = 0.5) +
  # scale_fill_gradient(low = "white", high = "black")+
  
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("chromosomal hotspots")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 

g


files <- list.files(path=".", pattern="500_0.05", full.names=TRUE)

down = fread(files[1])
up = fread(files[2])
down$direction = "Down"
up$direction = "Up"
out = rbind(up,down)

out <- out %>%
  mutate(
    chr = sub("(.*):(.*)-(.*)", "\\1", locations),
    start = as.numeric(sub("(.*):(.*)-(.*)", "\\2", locations)),
    end = as.numeric(sub("(.*):(.*)-(.*)", "\\3", locations)),
    midpoint = (start + end) / 2
)
out$region = sapply(strsplit(as.character(out$clust), ":"), `[`, 2)
out$region = factor(out$region, levels = c("FC", "ENT", "HCA", "HCP", "AMY", "NAC","RLP",  "CP"))

out$clade = "NN"
out$clade[grep("Glut", out$clust)] = "Glut"
out$clade[grep("Gaba", out$clust)] = "Gaba"
out$clade[grep("Neur", out$clust)] = "Gaba"

out$direction = factor(out$direction  , levels = c("Up", "Down"))

out$smth20_score[which(out$direction=="Down")] = -out$smth20_score[which(out$direction=="Down")]
out$chr = factor(out$chr, levels=paste0("chr",c(1:19,"X","Y")))
h3k9me3$chr = factor(h3k9me3$chr, levels=paste0("chr",c(1:19,"X","Y")))
options(repr.plot.width=12, repr.plot.height=8)
options(repr.plot.width=10, repr.plot.height=2)
out$smth20_score = -out$smth20_score

#out$clade = factor(out$clade, levels = c("NN", "Gaba", "Glut"))
g = ggplot(out[,]) +
  geom_point(aes(x=midpoint,y=smth20_score,color=clade),alpha=.8,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
#scale_color_manual(values=c( brewer.pal(n = 3, name = "Dark2")[1],brewer.pal(n = 3, name = "Dark2")[3],brewer.pal(n = 3, name = "Dark2")[2]
                   # ))+ 
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + ylim(c(-2.5,2.5))+
 # Add annotation bar for regions
  # geom_rect(data = h3k9me3[! which(h3k9me3$chr %in% c("chrX", "chrY")),], aes(
  #   xmin = start-100000, xmax = end+100000,fill = (score),
  #   ymin = -4, ymax = -3  # Position the bar beneath the plot
  # ), alpha = 0.5) +
  # scale_fill_gradient(low = "white", high = "black")+
  
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("chromosomal hotspots")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 

g


options(repr.plot.width=10, repr.plot.height=3)

out$sex = "Male"
out$sex[grep("Female", out$clust)] = "Female"
g = ggplot(out[which(out$chr %in% c("chrX", "chrY")),]) +
  geom_point(aes(x=midpoint,y=smth20_score,color=clade),alpha=.8,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
#scale_color_manual(values=c( brewer.pal(n = 3, name = "Dark2")[1],brewer.pal(n = 3, name = "Dark2")[3],brewer.pal(n = 3, name = "Dark2")[2]
                   # ))+ 
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + ylim(c(-2,2))+
   geom_rect(data = h3k9me3[ which(h3k9me3$chr %in% c("chrX", "chrY")),], aes(
     xmin = start-100000, xmax = end+100000,fill = (score),
     ymin = -4, ymax = -3  # Position the bar beneath the plot
   ), alpha = 0.5) +
   scale_fill_gradient(low = "black", high = "black")+
  
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(sex~chr,scales="free_x",space="fixed") + ggtitle("chromosomal hotspots")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 

g


out$sex = "female"
fem = out

files <- list.files(path=".", pattern="top1k_1logfc", full.names=TRUE)

down = fread(files[1])
up = fread(files[2])
down$direction = "Down"
up$direction = "Up"
out = rbind(up,down)

out <- out %>%
  mutate(
    chr = sub("(.*):(.*)-(.*)", "\\1", locations),
    start = as.numeric(sub("(.*):(.*)-(.*)", "\\2", locations)),
    end = as.numeric(sub("(.*):(.*)-(.*)", "\\3", locations)),
    midpoint = (start + end) / 2
)
out$region = sapply(strsplit(as.character(out$clust), ":"), `[`, 2)
out$region = factor(out$region, levels = c("FC", "ENT", "HCA", "HCP", "AMY", "NAC","RLP",  "CP"))

out$clade = "NN"
out$clade[grep("Glut", out$clust)] = "Glut"
out$clade[grep("Gaba", out$clust)] = "Gaba"
out$clade[grep("Neur", out$clust)] = "Gaba"

out$direction = factor(out$direction  , levels = c("Up", "Down"))

out$smth20_score[which(out$direction=="Down")] = -out$smth20_score[which(out$direction=="Down")]
out$chr = factor(out$chr, levels=paste0("chr",c(1:19,"X","Y")))
h3k9me3$chr = factor(h3k9me3$chr, levels=paste0("chr",c(1:19,"X","Y")))
options(repr.plot.width=12, repr.plot.height=8)
options(repr.plot.width=10, repr.plot.height=2)
out$smth20_score = -out$smth20_score

#out$clade = factor(out$clade, levels = c("NN", "Gaba", "Glut"))
g = ggplot(out[which(! out$chr %in% c("chrY")),]) +
  geom_point(aes(x=midpoint,y=smth20_score,color=clade),alpha=.8,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
#scale_color_manual(values=c( brewer.pal(n = 3, name = "Dark2")[1],brewer.pal(n = 3, name = "Dark2")[3],brewer.pal(n = 3, name = "Dark2")[2]
                   # ))+ 
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + ylim(c(-1.5,1.5))+
 # Add annotation bar for regions
  # geom_rect(data = h3k9me3[! which(h3k9me3$chr %in% c("chrX", "chrY")),], aes(
  #   xmin = start-100000, xmax = end+100000,fill = (score),
  #   ymin = -4, ymax = -3  # Position the bar beneath the plot
  # ), alpha = 0.5) +
  # scale_fill_gradient(low = "white", high = "black")+
  
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("chromosomal hotspots")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 

g


out$sex = "male"


out = rbind(out,fem)

g = ggplot(out[which(! out$chr %in% c("chrY")),]) +
  geom_point(aes(x=midpoint,y=smth20_score,color=clade),alpha=.8,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
#scale_color_manual(values=c( brewer.pal(n = 3, name = "Dark2")[1],brewer.pal(n = 3, name = "Dark2")[3],brewer.pal(n = 3, name = "Dark2")[2]
                   # ))+ 
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + ylim(c(-2,2))+
 # Add annotation bar for regions
   geom_rect(data = h3k9me3[which(h3k9me3$chr != "chrY"),], aes(
     xmin = start-100000, xmax = end+100000,fill = (score),
     ymin = -4, ymax = -3  # Position the bar beneath the plot
   ), alpha = 0.5) +
   scale_fill_gradient(low = "white", high = "black")+
  
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(sex~chr,scales="free_x",space="fixed") + ggtitle("chromosomal hotspots")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 

g


out$smth20_score = -out$smth20_score



options(repr.plot.width=10, repr.plot.height=8)


out$region = sapply(strsplit(as.character(out$clust), ":"), `[`, 3)
head(out)
out$region = gsub("_2vs18_iter.csv", "", out$region)
out$region = factor(out$region, levels = c("FC", "ENT", "HCA", "HCP", "AMY", "NAC","RLP",  "CP"))

out$sex = "Male"
out$sex[grep("Female",out$clust)] = "Female"

out$clade = "NN"
out$clade[grep("Glut", out$clust)] = "Glut"
out$clade[grep("Gaba", out$clust)] = "Gaba"
out$clade[grep("Neur", out$clust)] = "Gaba"
head(out)
g = ggplot(out[,]) +
  geom_point(aes(x=midpoint,y=smth20_score,color=clade),alpha=.8,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
#scale_color_manual(values=c( brewer.pal(n = 3, name = "Dark2")[1],brewer.pal(n = 3, name = "Dark2")[3],brewer.pal(n = 3, name = "Dark2")[2]
                   # ))+ 
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + ylim(c(-2.5,2.5))+
 # Add annotation bar for regions
  # geom_rect(data = h3k9me3[which(h3k9me3$chr %in% c("chrX")),], aes(
  #   xmin = start-100000, xmax = end+100000,fill = (score),
  #   ymin = -4, ymax = -2  # Position the bar beneath the plot
  # ), alpha = 0.5) +
  scale_fill_gradient(low = "black", high = "black")+
  
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(sex~chr,scales="free_x",space="fixed") + ggtitle("chromosomal hotspots")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 

g

unique(out[which(out$chr == "chrY"),"locations"])

unique(out[which(out$chr == "chrY"),"locations"])

g = ggplot(out[which(out$chr %in% c("chrX", "chrY")),]) +
  geom_point(aes(x=midpoint,y=smth20_score,color=clade),alpha=.8,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
#scale_color_manual(values=c( brewer.pal(n = 3, name = "Dark2")[1],brewer.pal(n = 3, name = "Dark2")[3],brewer.pal(n = 3, name = "Dark2")[2]
                   # ))+ 
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + #ylim(c(-2.5,2.5))+
 # Add annotation bar for regions
  geom_rect(data = h3k9me3[which(h3k9me3$chr %in% c("chrX", "chrY")),], aes(
    xmin = start-100000, xmax = end+100000,fill = (score),
    ymin = -4, ymax = -3  # Position the bar beneath the plot
  ), alpha = 0.5) +
  scale_fill_gradient(low = "white", high = "black")+
  
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(~chr,scales="free_x",space="fixed") + ggtitle("chromosomal hotspots")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 

g

g = ggplot(out[which(out$sex=="Male" & out$chr != "chrY"),]) +
  geom_point(aes(x=midpoint,y=smth20_score,color=clade),alpha=.8,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
#scale_color_manual(values=c( brewer.pal(n = 3, name = "Dark2")[1],brewer.pal(n = 3, name = "Dark2")[3],brewer.pal(n = 3, name = "Dark2")[2]
                   # ))+ 
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + ylim(c(-2.5,2.5))+
 # Add annotation bar for regions
  geom_rect(data = h3k9me3[! which(h3k9me3$chr %in% c("chrX", "chrY")),], aes(
    xmin = start-100000, xmax = end+100000,fill = (score),
    ymin = -4, ymax = -3  # Position the bar beneath the plot
  ), alpha = 0.5) +
  scale_fill_gradient(low = "white", high = "black")+
  
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("chromosomal hotspots")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 

g

ggsave(g, file = "../../Figures/Figure6-Het-TEs/chrom_hotspot_region_og_colors.pdf", height =9 , width = 9 )

up_file

setwd("../../../../combined_DARs_redone/smooth/")
files <- list.files(path=".", pattern="clustered_diff_peaks_top1k_1logfc", full.names=TRUE)

down = fread(files[1])
up = fread(files[2])
down$direction = "Down"
up$direction = "Up"
out = rbind(up,down)

out <- out %>%
  mutate(
    chr = sub("(.*):(.*)-(.*)", "\\1", locations),
    start = as.numeric(sub("(.*):(.*)-(.*)", "\\2", locations)),
    end = as.numeric(sub("(.*):(.*)-(.*)", "\\3", locations)),
    midpoint = (start + end) / 2
)
out$region = sapply(strsplit(as.character(out$clust), ":"), `[`, 2)
out$region = factor(out$region, levels = c("FC", "ENT", "HCA", "HCP", "AMY", "NAC","RLP",  "CP"))

out$clade = "NN"
out$clade[grep("Glut", out$clust)] = "Glut"
out$clade[grep("Gaba", out$clust)] = "Gaba"
out$clade[grep("Neur", out$clust)] = "Gaba"

out$direction = factor(out$direction  , levels = c("Up", "Down"))

out$smth20_score[which(out$direction=="Down")] = -out$smth20_score[which(out$direction=="Down")]
out$chr = factor(out$chr, levels=paste0("chr",c(1:19,"X","Y")))
h3k9me3$chr = factor(h3k9me3$chr, levels=paste0("chr",c(1:19,"X","Y")))
options(repr.plot.width=12, repr.plot.height=8)
options(repr.plot.width=10, repr.plot.height=15)

g = ggplot(out) +
  geom_point(aes(x=midpoint,y=smth20_score,color=clade),alpha=.8,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
  geom_hline(yintercept = 0, linetype = "solid") +
  geom_hline(yintercept = 0.2, linetype = "dashed") +

  theme_bw() + ylim(c(-4,4))+
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(clade~chr,scales="free_x",space="fixed") + ggtitle("DEG in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 
options(repr.plot.width=9, repr.plot.height=3)

g


hist(abs(out$smth20_score))

quantile(abs(out$smth20_score), 0.01)

g = ggplot(out) +
  geom_point(aes(x=midpoint,y=smth20_score,color=clade),alpha=.8,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + ylim(c(-5.5,5.5))+
 # Add annotation bar for regions
  geom_rect(data = h3k9me3, aes(
    xmin = start, xmax = end,fill = (score),
    ymin = -5.5, ymax = -4.5  # Position the bar beneath the plot
  ), alpha = 0.5) +
  scale_fill_gradient(low = "white", high = "black")+
  
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("DEG in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 
options(repr.plot.width=9, repr.plot.height=5)

g


peaks_gr = GRanges(seqnames = out$chr,
                   ranges = IRanges(start = out$start, end = out$end))


overlaps <- findOverlaps(peaks_gr, h3k9_gr)


h3k9_gr
table(h3k9_gr@seqnames)

peaks_gr

table(out[which(out$smth20_score<.2),"mark"])

50/3031

1943/2183


1943/(2183  +  1943 )

out$mark = ""
out$mark[queryHits(overlaps)] <- "H3K9me3"
overlaps


h3k9me3$score[which(h3k9me3$score>10000)] = 10000

out$clade = factor(out$clade, levels = c("NN", "Gaba", "Glut"))

g1 = ggplot(out[which(!out$chr %in% c("chrY", "chrX")),]) +
  geom_point(aes(x=midpoint,y=smth20_score,color=mark),alpha=.8,size =.8) + 
  scale_color_manual(values=c("gray","firebrick1", "blue3")) + 
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + ylim(c(-5,5))+
 # Add annotation bar for regions
  # geom_rect(data = h3k9me3[which(!h3k9me3$chr %in% c("chrY", "chrX")),], aes(
  #   xmin = start-100000, xmax = end+100000,
  #   ymin = -5.5, ymax = -4.5  # Position the bar beneath the plot
  # ), alpha = 1) +
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(clade~chr,scales="free_x",space="fixed") + ggtitle("Chromosomal Hotspots, overlap with H3K9me3 domains ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 
options(repr.plot.width=9, repr.plot.height=3)
g1
g = ggplot(out[which(!out$chr %in% c("chrY", "chrX")),]) +
  geom_point(aes(x=midpoint,y=smth20_score,color=mark),alpha=.8,size =.8) + 
  scale_color_manual(values=c("gray","firebrick1", "blue3")) + 
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + ylim(c(-5,5))+
 # Add annotation bar for regions
  # geom_rect(data = h3k9me3[which(!h3k9me3$chr %in% c("chrY", "chrX")),], aes(
  #   xmin = start-100000, xmax = end+100000,
  #   ymin = -5.5, ymax = -4.5  # Position the bar beneath the plot
  # ), alpha = 1) +
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(~chr,scales="free_x",space="fixed") + ggtitle("Chromosomal Hotspots, overlap with H3K9me3 domains ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 
options(repr.plot.width=9, repr.plot.height=2)


g


ggsave(g1, file = "../../Figures/Figure6-Het-TEs/chrom_hotspot_h3k9_clade.pdf", height =3 , width = 9 )

options(repr.plot.width=9, repr.plot.height=2.5)
#out$mark[which(abs(out$smth20_score)<1)] = ""
g = ggplot(out[which(!out$chr %in% c("chrY")),]) +
  geom_point(aes(x=midpoint,y=smth20_score,color=mark),alpha=.8,size =.8) + 
  scale_color_manual(values=c("grey","firebrick1"
                             )) + 
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme_bw() +#ylim(c(-5.5,5.5))+
 # Add annotation bar for regions
  #geom_rect(data = h3k9me3[which(!h3k9me3$chr %in% c("chrY", "chrX")),], aes(
  #  xmin = start, xmax = end,fill = (score),
  #  ymin = -5.5, ymax = -4.5  # Position the bar beneath the plot
  #), alpha = 0.5) +
  scale_fill_gradient(low = "white", high = "black")+
  
  theme( axis.ticks.x= element_blank(),legend.position = "top",
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(~chr,scales="free_x",space="fixed") + ggtitle("Grouped DARs in Aging (Hotspots)")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 

g

brewer.pal(n = 3, name = "Set1")


pdf("~/projects/combined_all/Figures/Figure6-Het-TEs/H3K9me3_chromplot.pdf", height =2.5, width = 10)
print(g)
dev.off()

options(repr.plot.width=7, repr.plot.height=6)

#out$clade[which(abs(out$smth20_score)<1)] = "ns"

g = ggplot(out[which(out$chr == "chr18" ),]) +
  geom_point(aes(x=midpoint,y=smth20_score,color=clade),alpha=.8,size =.8) + 
  scale_color_manual(values=c( brewer.pal(n = 3, name = "Dark2")[1],brewer.pal(n = 3, name = "Dark2")[3],brewer.pal(n = 3, name = "Dark2")[2]
                    )) + 
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_bw() + ylim(c(-3,3))+
 # Add annotation bar for regions
  #geom_rect(data = h3k9me3[which(h3k9me3$chr == "chr13"),], aes(
  #  xmin = start, xmax = end,fill = (score),
  #  ymin = -4.5, ymax = -3.5  # Position the bar beneath the plot
  #), alpha = 0.5) +
  scale_fill_gradient(low = "white", high = "black")+
  
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("DAR in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 

g

ggsave(g, file = "../../Figures/Figure6-Het-TEs/chrom_hotspot_region_chr13.pdf", height =6 , width = 6 )




locs = unique(out[which(!out$chr %in% c("chrY", "chrX")& out$smth20_score>1 ),"locations"])


setwd("../../female_RNA/DEG_results_rm_rpl_rps/")


files <- list.files(path=".", pattern=".csv", full.names=TRUE)
chr.len <- read.table("/projects/ps-renlab2/ps-renlab/share/bwa_indices/mm10.fa.fai")
chr.len <- rbind(chr.len, chr.len)
chr.len$logFC <- rep(c(-1,1), each=nrow(chr.len)/2)
colnames(chr.len)[1] <- "chr"
chr.len$color = c(rep(c("black","grey"),10),"black")
pval <- 0.1
locations <- c()
win_exp <- 5
win_size <- 10^win_exp
chr.window = NULL
size = 10^5
for (i in 1:nrow(chr.len)){
    chr.window = c(chr.window, paste0(chr.len[i,1],":",0:ceiling(chr.len[i,2]/size) ))
    } 

alldf = list()
pval = 1.1
for(f in files){
    a = read.csv(f)
    a <- a[which(a$p_val_adj < pval & abs(a$avg_log2FC) >= 0), ]
    #if(nrow(a)>15000){
    #    a = a[1:15000,]
    #}
    if(nrow(a)<1){
        next
    }
   # a$pos <- floor(as.integer(sub("(.*):(.*)-(.*)", "\\2", a$feature.name)) / win_size)
   # a$sig = log(a$adjusted.p.value+1e-200)*sign(a$log2.fold_change.)
    a$avg_log2FC = -a$avg_log2FC 
    a$dir = "Up"
    a$dir[which(a$avg_log2FC<0)] = "Down"
    a$file = gsub("./", "", f)
    a$file = gsub(".csv", "", a$file)
    alldf[[f]]=a
}

out = do.call(rbind,alldf)
head(out)

library(AnnotationDbi)
library(biomaRt)

all_genes = read.table("~/projects/combined_all/Figures/all_genes.txt")
all_genes = all_genes$V1
# Connect to Ensembl
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")  # Change to "hsapiens_gene_ensembl" for human
# Retrieve gene types for all genes in your dataset
  gene_annotations <- getBM(attributes = c("external_gene_name", "gene_biotype","transcript_length","transcript_gencode_basic", 
                                          "chromosome_name", "start_position", "end_position" , "percentage_gene_gc_content"
                                          ),
                          filters = "external_gene_name", values = all_genes,
                          mart = mart)
    

  # Filter for lncRNAs
  lncRNA_genes <- unique(gene_annotations$external_gene_name[gene_annotations$gene_biotype == "lncRNA"])
  protein_coding_genes <- unique(gene_annotations$external_gene_name[gene_annotations$gene_biotype == "protein_coding"])
  


head(gene_annotations)
gene_annotations_small = gene_annotations[-which(duplicated(gene_annotations$external_gene_name)),]
rownames(gene_annotations_small) = gene_annotations_small$external_gene_name

out$annotation <- "lncRNA + pseudogenes"
#deg_results$annotation[deg_results$X %in% lncRNA_genes] <- "lncRNA"
out$annotation[out$X %in% protein_coding_genes] <- "protein_coding"
head(out)

out$chr = paste("chr", gene_annotations_small[paste(out$X), "chromosome_name"], sep = "")
out$start = gene_annotations_small[paste(out$X), "start_position"]
out$end = gene_annotations_small[paste(out$X), "end_position"]
head(out)

library(tidyr)
locs_df <- as.data.frame(locs)

# Separate "chr" from "start-end"
locs_df <- locs_df %>%
  separate(col = "locations", into = c("chr", "positions"), sep = ":") %>%
  separate(col = "positions", into = c("start", "end"), sep = "-")

# Convert start/end to numeric
locs_df$start <- as.integer(locs_df$start)
locs_df$end   <- as.integer(locs_df$end)
nrow(locs_df)

nrow(out[which(is.na(out$start)),])
nrow(out)



# Convert h3k9 to a GRanges object
locs_gr <- GRanges(seqnames = locs_df$chr,
                   ranges = IRanges(start = locs_df$start, end = locs_df$end))

# Convert gene_annotations to a GRanges object
gene_annotations_gr <- GRanges(seqnames = paste0("chr", gene_annotations$chromosome_name),
                               ranges = IRanges(start = gene_annotations$start_position, 
                                                end = gene_annotations$end_position))

# Find overlaps
overlaps <- findOverlaps(gene_annotations_gr, h3k9_gr)

# Add mark column, default to NA
gene_annotations$mark <- NA

# Assign "h3k9me3" to genes that overlap with any h3k9 region
gene_annotations$mark[queryHits(overlaps)] <- "hotspot"

# View updated gene_annotations
head(gene_annotations)


gene_annotations$differential = "no"
gene_annotations$differential[which(gene_annotations$external_gene_name %in%
                                    out$X[which(out$dir=="Up")])] = "Up"
gene_annotations$differential[which(gene_annotations$external_gene_name %in%
                                    out$X[which(out$dir=="Down")])] = "Down"


out$hotspot = ""
out$hotspot[which(out$X%in%gene_annotations$external_gene_name[
    which(gene_annotations$mark == "hotspot")])] = "hotspot"

out = out[order(out$p_val),]
out[which(out$dir == "Up" & out$hotspot=="hotspot"),]

locs_df

gene_annotations[which(gene_annotations$mark == "hotspot")[1:5],]

h3k9_genes = unique(gene_annotations$external_gene_name[which(gene_annotations$mark == "hotspot")])

h3k9_genes

options(repr.plot.width=7, repr.plot.height=8)

library(ggplot2)
library(dplyr)
library(tidyr)
# If available, consider readr::read_csv for faster file reading.

# Define directory containing DEG results
deg_dir <- "./../DEG_results_rm_rpl_rps//"

# Get all CSV files in the directory
deg_files <- list.files(deg_dir, full.names = TRUE, pattern = "*.csv")

# Placeholder for results
deg_summary <- data.frame()

# Loop through each DEG file
for (file in deg_files) {
  
  # Extract cell type name from file name
  celltype <- gsub(".csv", "", basename(file))
  
  # Read data
  deg_results <- read.csv(file)
  
  # Reverse the sign of avg_log2FC if needed
  deg_results$avg_log2FC <- -deg_results$avg_log2FC
    
  
  # Annotate gene types
  deg_results$annotation <- "lncRNA + pseudogenes"
  #deg_results$annotation[deg_results$X %in% lncRNA_genes] <- "lncRNA"
  deg_results$annotation[deg_results$X %in% protein_coding_genes] <- "protein_coding"
  deg_results$mark <- "Other"
  deg_results$mark[deg_results$X %in% h3k9_genes] <- "H3K9me3 domain"
  
  n = table(deg_results$annotation)[[1]]/(table(deg_results$annotation)[[2]]+table(deg_results$annotation)[[1]])
    
  # Filter significant genes (p_val_adj < 0.05) to decide if there is enough signal
  sig_degs <- deg_results %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) >0.25)
  
  # Only proceed if there are at least 100 significant DEGs
 # if (nrow(sig_degs) < 100) {
 #   next
 # }
  
  # Select top 100 upregulated and top 100 downregulated genes from the significant ones
  top_up <- sig_degs[which(sig_degs$p_val_adj< 0.05),] %>% 
    arrange(-avg_log2FC) %>% 
    head(1000) %>% 
    mutate(category = "Upregulated")
  
  top_down <- sig_degs[which(sig_degs$p_val_adj< 0.05),] %>% 
    arrange(avg_log2FC) %>% 
    head(1000) %>% 
    mutate(category = "Downregulated")
  
  # Instead of filtering only for non-significant genes, select all genes that are not in the top lists.
  # This will include genes that are significant (but not in the top 100) and non-significant genes.
  all_others <- deg_results %>% 
    filter(!(X %in% c(top_up$X, top_down$X))) %>% 
    mutate(category = "All Genes")
  
  # Combine the three groups: top up, top down, and all other genes.
  combined_genes <- bind_rows(top_up, top_down, all_others)
  
  # Count genes in each category per annotation type
  count_summary <- combined_genes %>%
    group_by(celltype = celltype, category, annotation, mark) %>%
    summarise(count = n(), .groups = "drop")
  
  # Append to overall summary
  deg_summary <- bind_rows(deg_summary, count_summary)
}

# Assign clade labels based on celltype patterns
deg_summary$clade <- "Glut"
deg_summary$clade[grep("NN", deg_summary$celltype)] <- "NN"
deg_summary$clade[grep("IMN", deg_summary$celltype)] <- "NN"
deg_summary$clade[grep("Gaba", deg_summary$celltype)] <- "Gaba"
deg_summary$clade[grep("Neur", deg_summary$celltype)] <- "Gaba"

deg_summary$clade= factor(deg_summary$clade, levels = c("NN", "Gaba", "Glut"))

deg_summary = deg_summary[-which(deg_summary$category=="All Genes"),]
#deg_summary = deg_summary[which(deg_summary$celltype %in% cts), ] 

deg_summary$celltype = factor(as.character(deg_summary$celltype))
deg_summary$celltype = factor(deg_summary$celltype , levels = rev(levels(deg_summary$celltype)))

deg_summary$annotation = factor(deg_summary$annotation , levels = c("protein_coding", "lncRNA + pseudogenes"))

# Plot gene counts per cell type
g=ggplot(deg_summary, aes(x = celltype, y = count, fill = annotation)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(clade ~ category, scales = "free", space = "free_y") +
  theme_bw() +  
  scale_fill_manual(values = c("H3K9me3 domain" = "tomato", "Other"= "grey", "lncRNA + pseudogenes" = "tomato", "protein_coding" = "lightgrey" ))+
  labs(title = "Top 100 Age-DEGs",
       x = "Cell Type",
       y = "Gene Proportion",
       fill = "Gene Annotation") +
  coord_flip() +
  geom_hline(yintercept = n
             , linetype = "dashed") +  # Horizontal line at 0.2

theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        strip.text.y = element_blank(),
        axis.text.y = element_text(size = 14), text = element_text(size = 14),
        legend.position = "top") 
g

g=ggplot(deg_summary, aes(x = celltype, y = count, fill = mark)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(clade ~ category, scales = "free", space = "free_y") +
  theme_bw() +  
  scale_fill_manual(values = c("H3K9me3 domain" = "tomato", "Other"= "grey", "lncRNA + pseudogenes" = "tomato", "protein_coding" = "lightgrey" ))+
  labs(title = "Top 100 Age-DEGs",
       x = "Cell Type",
       y = "Gene Proportion",
       fill = "Gene Annotation") +
  coord_flip() +
  geom_hline(yintercept = n
             , linetype = "dashed") +  # Horizontal line at 0.2

theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        strip.text.y = element_blank(),
        axis.text.y = element_text(size = 14), text = element_text(size = 14),
        legend.position = "top") 
g

length(unique(gene_annotations$external_gene_name[
    which(gene_annotations$gene_biotype == "lncRNA")]))

length(unique(gene_annotations$external_gene_name[
    which(gene_annotations$gene_biotype == "protein_coding")]))

length(unique(gene_annotations$external_gene_name[
    which(gene_annotations$gene_biotype == "protein_coding" &
         gene_annotations$mark == "hotspot")]))

length(unique(gene_annotations$external_gene_name[
    which(gene_annotations$gene_biotype == "lncRNA" &
         gene_annotations$mark == "hotspot")]))


length(unique(gene_annotations$external_gene_name))

17132/24480

111/7308
354/17132


gene_annotations_small = gene_annotations[-which(duplicated(gene_annotations$external_gene_name)),]

table(gene_annotations_small[,"gene_biotype"])

out$hotspot = NA
out$hotspot[which(out$X%in%
                  gene_annotations[which(gene_annotations$mark=="hotspot"),"external_gene_name"])] = "hotspot"

out = out[order(out$p_val),]
(out[which(out$hotspot == "hotspot"),])


