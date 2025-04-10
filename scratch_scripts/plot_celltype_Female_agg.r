library(ggplot2)

meta= read.csv("../final_meta.csv")
head(meta)
meta$ct = gsub(" ", "_" , meta$celltype_final)
good_cts = names(table(meta$ct)[which(table(meta$ct)>5000)])
good_cts = gsub("/", "_", good_cts)

good_cts


a=read.table("Male_L1MA5A_clade_age.matrix.mat.gz",skip=1)
b=readLines("Male_L1MA5A_clade_age.matrix.mat.gz",n=1)


head(b)

samples = unlist(strsplit(b,","))[grep("merge", unlist(strsplit(b,",")))]
samples = gsub("\"", "", samples)
samples = gsub(".all", "", samples)
samples = gsub("]", "", samples)
samples = gsub("sample_labels:[", "", samples,fixed = T)
samples = samples[grep("merge",samples)]

samples

a2= a[,-c(1:6)]
mean = colSums(a2)/nrow(a2)

len = length(a2)/length(samples)

dat = data.frame(sample=rep(samples,each=len),
                 pos = rep(1:len,length(samples)),
                 mean =mean
)


WIN = 10

dat2 = dat[which(dat$pos>WIN & dat$pos < 1000-WIN),]

for ( sample in unique(dat2$sample)) {
  for ( pos in unique(dat2$pos)) {
  #  print(as.character(sample))
   # print(pos)
    dat2$mean[which(dat2$sample==sample & dat2$pos==pos)] = mean(dat$mean[which(dat$sample==sample & dat$pos %in% (pos-WIN):(pos+WIN))]) 
  }
}

head(dat2)


#dat2 = dat2[which( dat2$pos < 115), ]
#dat2 = dat2[which( dat2$pos >35), ]

dat2$age = rep("2mo" , nrow(dat2))
dat2$age[grep("18mo", dat2$sample)] = "18mo"
dat2$age[grep("9mo", dat2$sample)] = "9mo"
dat2$age = factor(dat2$age , levels = c("2mo", "9mo", "18mo"))

dat2$clade = rep("Gaba" , nrow(dat2))
dat2$clade[grep("Glut", dat2$sample)] = "Glut"
dat2$clade[grep("NN", dat2$sample)] = "NN"
dat2$clade[grep("IMN", dat2$sample)] = "IMN"
dat2$clade = factor(dat2$clade, levels = c("Glut", "Gaba", "NN", "IMN"))

dat2$celltype = dat2$sample
dat2$celltype = gsub("Female_", "", dat2$celltype)
dat2$celltype = gsub("_merged", "", dat2$celltype)
dat2$celltype = gsub("_18mo", "", dat2$celltype)
dat2$celltype = gsub("_2mo", "", dat2$celltype)
dat2$celltype = gsub("_9mo", "", dat2$celltype)


head(dat2)


head(dat2)

#dat2 = dat2[which(dat2$celltype%in% good_cts),]
dat2$celltype = factor(dat2$celltype)

head(subset(dat2,pos%%100==0))


unique(dat2$sample)

dat2 = dat2[-grep("Oligo", dat2$sample),]

dat2$region = sapply(strsplit(as.character(dat2$sample), "_"), `[`, 1)


sub = subset(dat2[which(dat2$clade!= "IMN"),],pos%%73==0)
levs = unique(sub[order(sub$mean, decreasing = T), "region"])
dat2$region = factor(dat2$region, levels = levs)

dat2$region=factor(dat2$region,levels = c('FC','EC','AMY','HCA','HCP','RLP','NAC','CP'))

unique(dat2$region)

options(repr.plot.width=7, repr.plot.height=4)

g = ggplot(dat2[which(dat2$clade!= "IMN"),]) + geom_line(aes(x=pos,y=mean,color=clade,group=sample, linetype = region)) +
# geom_hline(yintercept = 0, linetype = "dotted") + # Add dotted line at y = 0
    xlim(c(20,120))+

  #geom_point(data=subset(dat2[which(dat2$clade!= "IMN"),],pos%%73==0), aes(x=pos,y=mean,color=clade),size=5) + 
  facet_grid(~age) + theme_classic() + ggtitle("L1MA5A Accessibility")

options(repr.plot.width=6.5, repr.plot.height=3.5)


g

pdf("~/projects/combined_all/Figures/Figure6-Het-TEs/L1MA5A_small.pdf", height = 3.5 , 7)
print(g)
dev.off()

options(repr.plot.width=6.5, repr.plot.height=6)

g = ggplot(dat2[which(dat2$clade!= "IMN"),]) + geom_line(aes(x=pos,y=mean,color=region,group=sample)) +
 geom_hline(yintercept = 0, linetype = "dotted") + # Add dotted line at y = 0

  #geom_point(data=subset(dat2[which(dat2$clade!= "IMN"),],pos%%73==0), aes(x=pos,y=mean,color=clade),size=5) + 
  facet_grid(clade~age) + theme_classic() + ggtitle("L1MA5A Accessibility - Female")

g



options(repr.plot.width=6.5, repr.plot.height=10)

g = ggplot(dat2[which(dat2$clade!= "IMN"),]) + geom_line(aes(x=pos,y=mean,color=clade,group=sample)) +
 geom_hline(yintercept = 0, linetype = "dotted") + # Add dotted line at y = 0

  #geom_point(data=subset(dat2[which(dat2$clade!= "IMN"),],pos%%73==0), aes(x=pos,y=mean,color=clade),size=5) + 
  facet_grid(region~age) + theme_classic() + ggtitle("L1MA5A Accessibility - Female")

g

ggsave(g, width =6.5, height = 10 , file = "Female_clade_L1MA5A.pdf")



dats = dat2[grep("L2_3_IT_CTX", dat2$sample),]

options(repr.plot.width=25, repr.plot.height=6)

ggplot(dat2) + geom_line(aes(x=pos,y=mean, color = clade, shape = celltype)) +
  geom_point(data=subset(dat2,pos%%222==0), aes(x=pos,y=mean),size=5) + 
  facet_wrap(~age) + 
  theme_classic()

ggplot(dat2) + geom_line(aes(x=pos,y=mean)) +
  geom_point(data=subset(dat2,pos%%222==0), aes(x=pos,y=mean),size=5) + 
  facet_wrap(~sample) + 
  theme_classic()

library(ggplot2)


