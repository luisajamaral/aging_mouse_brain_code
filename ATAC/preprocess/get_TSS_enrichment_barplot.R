library(ggplot2)
# make large tss enrichment plot 
# ls *18*/*.txt *09

setwd("~/projects/brain_aging_mouse/analysis/by_sample.snapATAC/")
files = system("ls */snapFiles/*7.meta.txt", intern = T)
#files = system("ls *18*/snapFiles/*7.meta.txt *9mo*/snapFiles/*7.meta.txt", intern = T)

alltss = list()
for(f in files) {
   cur = read.csv(file = f, sep = "\t", header = T)
   alltss[[paste(cur$x.sp.sample[1])]] = cur
 }

t = do.call("rbind", alltss)
#t = t[-which(t$V5>=50), ]
#t = t[-which(t$V5<10), ]
t$sample = gsub("03", "8wk", t$sample)
pdf("all_TSS_QC_plot.pdf", height = 7, width = 15)

p <- ggplot(t, aes(x=x.sp.sample, y=TSS_enrich)) + 
  geom_boxplot(fill='deepskyblue1', outlier.shape = NA,notch=TRUE) + 
  scale_y_continuous(limits = c(10,40)) + ylab(label = "TSS Enrichment") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p #+ stat_summary(fun.y=mean, geom="point", shape=23, size=4)


p <- ggplot(t, aes(x=x.sp.sample, y=log10(UQ))) + 
  geom_boxplot(fill='pink', outlier.shape = NA,notch=TRUE) + 
   ylab(label = "log10UMI") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p


g <- ggplot(t, aes(x.sp.sample)) +
  ylab(label = "#cells") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
g + geom_bar()
dev.off()


 setwd("~/projects/brain_aging_mouse/data/snATAC/qc.TSSvsFrag.out/")
 files = system("ls */*.txt", intern = T)
# 
# 
 alltss = list()
 for(f in files) {
   cur = read.csv(file = f, sep = "\t", header = F)
   cur$sample = gsub("/stat.txt", "", f)
   alltss[[gsub("/stat.txt", "", f)]] = cur
 }

t = do.call("rbind", alltss)
#t = t[-which(t$V5>=50), ]
#t = t[-which(t$V5<10), ]

p <- ggplot(t, aes(x=sample, y=V5)) + 
  geom_boxplot(fill='deepskyblue1', outlier.shape = NA,notch=TRUE) + 
  scale_y_continuous(limits = c(0,40)) + ylab(label = "TSS Enrich before filtering") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p #+ stat_summary(fun.y=mean, geom="point", shape=23, size=4)
dev.off()




p <- ggplot(sumF, aes(x=sample, y=doublet_scores)) + 
  geom_boxplot(fill='deepskyblue1', outlier.shape = NA,notch=TRUE) + ylab(label = "Doublet scores") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p #
