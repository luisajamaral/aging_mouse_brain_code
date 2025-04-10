library(dplyr)
library(data.table)
library(ggplot2)

setwd("../../h5ads_final/region_sex_DARs/diff_csvs/")

files <- list.files(path=".", pattern="iter.csv", full.names=TRUE)
chr.len <- read.table("/projects/ps-renlab2/ps-renlab/share/bwa_indices/mm10.fa.fai")
chr.len <- rbind(chr.len, chr.len)
chr.len$logFC <- rep(c(-1,1), each=nrow(chr.len)/2)
colnames(chr.len)[1] <- "chr"
chr.len$color = c(rep(c("black","grey"),10),"black")

pval <- 0.05
locations <- c()
win_exp <- 5
win_size <- 10^win_exp

alldf = list()

for(f in files){
    a = read.csv(f)
    a$chr <- sub("(.*):(.*)-(.*)", "\\1", a$feature.name)
    a <- a[which(a$adjusted.p.value < pval & a$chr %in% chr.len$chr & abs(a$log2.fold_change.) > 0.2), ]
    if(nrow(a)>10000){
        a = a[1:10000,]
    }
    if(nrow(a)<1){
        next
    }
    a$pos <- floor(as.integer(sub("(.*):(.*)-(.*)", "\\2", a$feature.name)) / win_size)
    a$sig = log(a$adjusted.p.value+1e-200)*sign(a$log2.fold_change.)
    a$dir = "Up"
    a$dir[which(a$sig<0)] = "Down"
    a$file = gsub("./diff_peaks_", "", f)
    a$file = gsub("_2vs18_iter.csv", "", a$file)
    alldf[[f]]=a
}

out = do.call(rbind,alldf)

out$clade = "NN"
out$clade[grep("Glut", out$file)] = "Glut"
out$clade[grep("Gaba", out$file)] = "Gaba"
out$clade[grep("Neur", out$file)] = "Gaba"

out$sex = "Male"
out$sex[grep("Female", out$file)] = "Female"

out$region = sapply(strsplit(as.character(out$file), ":"), `[`, 3)

head(out)


out$region = factor(out$region, levels = c("FC", "ENT", "AMY", "RLP", "HCA","HCP","NAC","CP"))
out$chr = factor(out$chr, levels=paste0("chr",c(1:19,"X","Y")))
out = out[order(out$adjusted.p.value, decreasing = F),]

options(repr.plot.width=12, repr.plot.height=6)

ggplot(out[which(out$adjusted.p.value<0.001&  abs(out$log2.fold_change.)>1 ),]) +
  geom_point(aes(x=pos,y=-log2.fold_change.,color=clade),alpha=.8,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("logFC in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 


head(out)


ggplot(out[which(out$adjusted.p.value<0.001&  abs(out$log2.fold_change.)>1 & out$cel),]) +
  geom_point(aes(x=pos,y=-log2.fold_change.,color=clade),alpha=.8,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("logFC in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 


options(repr.plot.width=12, repr.plot.height=3)

g1 = ggplot(out[which(  out$sig>5 & abs(out$log2.fold_change.)>1.5 ),]) +
  geom_point(aes(x=pos,y=-log2.fold_change.,color=clade),alpha=1,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
    #  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
 # geom_jitter(data=sig_up,aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +
  #  geom_text(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.35,label=pos),angle=45,color="red") +
  #geom_hline(yintercept=2,linetype="dashed") +
  #ylim(-.5,2) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(sex~chr,scales="free_x",space="fixed") + ggtitle("Up in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 

ggsave(g1, file = "~/projects/combined_all/Figures/sex_logFC_chromplot_up.pdf", height = 3 , width = 12)

g2 = ggplot(out[which(  out$sig<(-5) & abs(out$log2.fold_change.)>1.5 ),]) +
  geom_point(aes(x=pos,y=log2.fold_change.,color=clade),alpha=1,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
    #  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
 # geom_jitter(data=sig_up,aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +
  #  geom_text(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.35,label=pos),angle=45,color="red") +
  #geom_hline(yintercept=2,linetype="dashed") +
  #ylim(-.5,2) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(sex~chr,scales="free_x",space="fixed") + ggtitle("Down in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 
ggsave(g2, file = "~/projects/combined_all/Figures/sex_logFC_chromplot_down.pdf", height = 3 , width = 12)


unique((out[order(out$log2.fold_change.),"feature.name"]))[1:200]

#### options(repr.plot.width=12, repr.plot.height=4)

ggplot(out[which(  out$sig>5 & abs(out$log2.fold_change.)>1.5 ),]) +
  geom_point(aes(x=pos,y=-log2.fold_change.,color=clade),alpha=1,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
    #  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
 # geom_jitter(data=sig_up,aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +
  #  geom_text(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.35,label=pos),angle=45,color="red") +
  #geom_hline(yintercept=2,linetype="dashed") +
  #ylim(-.5,2) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(sex~chr,scales="free_x",space="fixed") + ggtitle("Up in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 



ggplot(out[which(  out$sig<(-5) & abs(out$log2.fold_change.)>1.5 ),]) +
  geom_point(aes(x=pos,y=log2.fold_change.,color=clade),alpha=1,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
    #  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
 # geom_jitter(data=sig_up,aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +
  #  geom_text(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.35,label=pos),angle=45,color="red") +
  #geom_hline(yintercept=2,linetype="dashed") +
  #ylim(-.5,2) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(sex~chr,scales="free_x",space="fixed") + ggtitle("Down in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 


options(repr.plot.width=12, repr.plot.height=8)

ggplot(out[which(out$sex =="Female" & out$adjusted.p.value<0.001 & abs(out$log2.fold_change.)>1 ),]) +
  geom_point(aes(x=pos,y=-log2.fold_change.,color=clade),alpha=1,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
    #  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
 # geom_jitter(data=sig_up,aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +
  #  geom_text(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.35,label=pos),angle=45,color="red") +
  #geom_hline(yintercept=2,linetype="dashed") +
  #ylim(-.5,2) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("LogFC in Aging - Female")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 


ggplot(out[which(out$sex =="Male" & out$adjusted.p.value<0.001 & abs(out$log2.fold_change.)>1 ),]) +
  geom_point(aes(x=pos,y=-log2.fold_change.,color=clade),alpha=1,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
    #  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
 # geom_jitter(data=sig_up,aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +
  #  geom_text(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.35,label=pos),angle=45,color="red") +
  #geom_hline(yintercept=2,linetype="dashed") +
  #ylim(-,2) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("LogFC in Aging - Male")#+  geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 





ggplot(out[which(out$sex =="Male" & out$sig>5 & abs(out$log2.fold_change.)>1 ),]) +
  geom_point(aes(x=pos,y=-log2.fold_change.,color=clade),alpha=1,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
    #  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
 # geom_jitter(data=sig_up,aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +
  #  geom_text(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.35,label=pos),angle=45,color="red") +
  #geom_hline(yintercept=2,linetype="dashed") +
  #ylim(-.5,2) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("Up in Aging - Male")#+  geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 


ggplot(out[which(out$sex =="Male" & out$sig<(-5) & abs(out$log2.fold_change.)>1 ),]) +
  geom_point(aes(x=pos,y=-log2.fold_change.,color=clade),alpha=1,size =.8) + 
  scale_color_manual(values=c("green2","firebrick1", "blue3")) + 
    #  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
 # geom_jitter(data=sig_up,aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +
  #  geom_text(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.35,label=pos),angle=45,color="red") +
  #geom_hline(yintercept=2,linetype="dashed") +
  #ylim(-.5,2) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("Down in Aging - Male")#+  geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 





