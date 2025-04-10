library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)

p=system("pwd", intern = T)
p

setwd("DEG_results_rm_rpl_rps/")
#setwd("../DEG_results_.01_.01//")

gtf = read.table("/mnt/tscc/yanxiao/annotations/mm10/gencode.vM10.annotation.gene.tss1k.bed")

gtf = read.table("../gencode.vM10.annotation.gene.tss1k.bed")

head(gtf)

length(gtf$V4[which(duplicated(gtf$V4))])

gtf[which(gtf$V4=="Flg"),]
gtf=gtf[-which(duplicated(gtf$V4)),]
rownames(gtf)=gtf$V4

      
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

a= read.csv(files[1])
a <- a[which(a$p_val_adj < pval & abs(a$avg_log2FC) > 0.1), ]

head(a)

files[1]

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

length(unique(out$X[which(out$p_val_adj<0.000001 & abs(out$avg_log2FC)>0.5)]))

write(unique(out$X[which(out$p_val_adj<0.000001 & abs(out$avg_log2FC)>0.5)]), file = "DEGs_0.5LFC_fdr_1e6.txt", sep = "\n")

table(out$file[which(out$p_val_adj<0.000001 & abs(out$avg_log2FC)>0.5)])

out$position = (gtf[paste(out$X),2]+gtf[paste(out$X),3])/2

out$chr = gtf[paste(out$X),1]

head(out)


out = out[order(out$p_val_adj),]

unannot = out[which(is.na(out$chr)),]

out[which(out$X == "Gm47283"), "chr"] = "chrY"
out[which(out$X == "Gm47283"), "position"] = (90796007+90827734)/2

out$pos = floor(out$position/10^5)

head(out[which(is.na(out$chr)),])

out$clade = "NN"
out$clade[grep("Glut", out$file)] = "Glut"
out$clade[grep("Gaba", out$file)] = "Gaba"
out$clade[grep("Neur", out$file)] = "Gaba"


out$region = sapply(strsplit(as.character(out$file), "--"), `[`, 2)
out$celltype = sapply(strsplit(as.character(out$file), "--"), `[`, 1)



head(out)


out$region = factor(out$region, levels = c("FC", "ENT", "AMY", "RLP", "HCA","HCP","NAC","CP"))
out$chr = factor(out$chr, levels=paste0("chr",c(1:19,"X","Y")))


options(repr.plot.width=12, repr.plot.height=6)

ggplot(out[which(out$p_val_adj<0.001&  abs(out$avg_log2FC)>1 ),]) +
  geom_point(aes(x=position,y=(log(p_val_adj+1e-200)*sign(avg_log2FC)),color=clade),alpha=.8,size =.8) + 
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
  facet_grid(~chr,scales="free_x",space="fixed") + ggtitle("DEG in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 


options(repr.plot.width=12, repr.plot.height=6)

ggplot(out[which(out$p_val_adj<0.001&  abs(out$avg_log2FC)>1 ),]) +
  geom_point(aes(x=position,y=(log(p_val_adj+1e-200)*sign(avg_log2FC)),color=clade),alpha=.8,size =.8) + 
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
  facet_grid(~chr,scales="free_x",space="fixed") + ggtitle("DEG in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 


options(repr.plot.width=12, repr.plot.height=6)

ggplot(out[which(out$p_val_adj<0.001&  abs(out$avg_log2FC)>1 ),]) +
  geom_point(aes(x=position,y=(log(p_val_adj+1e-200)*sign(avg_log2FC)),color=clade),alpha=.8,size =.8) + 
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
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("DEG in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 


ggplot(out[which(out$p_val_adj<0.001&  abs(out$avg_log2FC)>.5 & !is.na(out$chr)),]) +
  geom_point(aes(x=position,y=(log(p_val_adj+1e-600)*sign(avg_log2FC)),color=clade),alpha=.8,size =.8) + 
  scale_color_manual(values=c( "green2","blue3", "firebrick1")) + 
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("DEG in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 


ggplot(out[which(out$p_val_adj<0.001&  abs(out$avg_log2FC)>.5 & !is.na(out$chr)),]) +
  geom_point(aes(x=position,y=(log(p_val_adj+1e-600)*sign(avg_log2FC)),color=dir),alpha=.8,size =.8) + 
  scale_color_manual(values=c( "blue3", "firebrick1")) + 
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(clade~chr,scales="free_x",space="fixed") + ggtitle("DEG in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 


head(out[which(out$p_val_adj<0.01&  abs(out$avg_log2FC)>.75 & !is.na(out$chr) & out$celltype=="L2-3_IT_CTX_Glut"),])

unique(out$celltype)


ggplot(out[which(out$p_val_adj<0.01&  abs(out$avg_log2FC)>.75 & !is.na(out$chr) & out$celltype=="STR_D12_Gaba"),]) +
  geom_point(aes(x=pos,y=(log(p_val_adj+1e-500)*sign(avg_log2FC)),color=dir),alpha=.8,size =.8) + 
  scale_color_manual(values=c( "blue3", "firebrick1")) + 
  theme_bw() + geom_hline(yintercept=0,linetype="dashed") +
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("DEG in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 




out$family = ""
out$family[grep("Rp", out$X)]="Rp"
out$family[grep("Atp", out$X)]="Atp"
out$family[grep("Nduf", out$X)]="Nduf"
out$family[grep("Mrp", out$X)]="Mrp"
out$family[grep("Cox", out$X)]="Cox"



#out$family[which(out$X %in% fam)]=fam
#out$family[which(out$X %in% unclustered_protos)]="unclustered_protocadherins"
#out$family[which(out$X %in% major_cadherins)]="major_cadherins"

#out$family[grep("Atp", out$X)]="Atp"
#out$family[grep("Cdh", out$X)]="Cdh"

out$family = factor(out$family, levels = c( "Rp", "Atp", "Nduf", "Mrp","Cox",""))
out = out[order(out$family,decreasing = T),]


ggplot(out[which(out$p_val_adj<0.01&  abs(out$avg_log2FC)>.25 & !is.na(out$chr)),]) +
  geom_point(aes(x=position,y=(log(p_val_adj+1e-600)*sign(avg_log2FC)),color=family),alpha=.8,size =.8) + 
  scale_color_manual(values=c( "green2","blue3", "firebrick1","magenta", "orange","grey")) + 
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+ geom_hline(yintercept=0,linetype="dashed") +
  facet_grid(clade~chr,scales="free_x",space="fixed") + ggtitle("DEG in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 


options(repr.plot.width=12, repr.plot.height=10)


out$family = ""
out$family[grep("Rp", out$X)]="Rp"
out$family[grep("Atp", out$X)]="Atp"
out$family[grep("Nduf", out$X)]="Nduf"
out$family[grep("Mrp", out$X)]="Mrp"
out$family[grep("Cox", out$X)]="Cox"



#out$family[which(out$X %in% fam)]=fam
#out$family[which(out$X %in% unclustered_protos)]="unclustered_protocadherins"
#out$family[which(out$X %in% major_cadherins)]="major_cadherins"

#out$family[grep("Atp", out$X)]="Atp"
#out$family[grep("Cdh", out$X)]="Cdh"

out$family = factor(out$family, levels = c( "Rp", "Atp", "Nduf", "Mrp","Cox",""))
out = out[order(out$family,decreasing = T),]


ggplot(out[which(out$p_val_adj<0.01&  abs(out$avg_log2FC)>.25 & !is.na(out$chr)),]) +
  geom_point(aes(x=position,y=(log(p_val_adj+1e-600)*sign(avg_log2FC)),color=family),alpha=.8,size =.8) + 
  scale_color_manual(values=c( "green2","blue3", "firebrick1","magenta", "orange","grey")) + 
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+ geom_hline(yintercept=0,linetype="dashed") +
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("DEG in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 


#unique(out[which(out$family!=""),"X"]) 
#write(unique(out[which(out$family!=""),"X"])  , sep = "\n" ,file = "oxphos_ribo_genes.txt") 

rps = out[grep("Atp", out$X),]


ggplot(rps[which(rps$p_val_adj<0.05 & !is.na(rps$chr) ),]) +
  geom_point(aes(x=pos,y=(log(p_val_adj+1e-500)*sign(avg_log2FC)),color=dir),alpha=.8,size =.8) + 
  scale_color_manual(values=c( "blue3", "firebrick1")) + 
  theme_bw() + geom_hline(yintercept=0,linetype="dashed") +
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("DEG in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 

rps = out[grep("Rp", out$X),]


ggplot(rps[which(rps$p_val_adj<0.05 & !is.na(rps$chr) ),]) +
  geom_point(aes(x=pos,y=(log(p_val_adj+1e-500)*sign(avg_log2FC)),color=dir),alpha=.8,size =.8) + 
  scale_color_manual(values=c( "blue3", "firebrick1")) + 
  theme_bw() + geom_hline(yintercept=0,linetype="dashed") +
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("DEG in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 


rps = out[grep("Nduf", out$X),]


ggplot(rps[which(rps$p_val_adj<0.05 & !is.na(rps$chr) ),]) +
  geom_point(aes(x=pos,y=(log(p_val_adj+1e-500)*sign(avg_log2FC)),color=dir),alpha=.8,size =.8) + 
  scale_color_manual(values=c( "blue3", "firebrick1")) + 
  theme_bw() + geom_hline(yintercept=0,linetype="dashed") +
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+guides(colour = guide_legend(override.aes = list(size=2)))+
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("DEG in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 



rps[which(is.na(rps$region)),]

options(repr.plot.width=12, repr.plot.height=6)

ggplot(out[which(out$p_val_adj<0.001&  abs(out$avg_log2FC)>1 ),]) +
  geom_point(aes(x=pos,y=(log(p_val_adj+1e-200)*sign(avg_log2FC)),color=clade),alpha=.8,size =.8) + 
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
  facet_grid(region~chr,scales="free_x",space="fixed") + ggtitle("DEG in Aging ")#+ # geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 


rps = out[grep("Rps", out$X),]
nrow(rps)

cadherins= c("Cdh16", "Cdh17", "Cdhr1", "Fat3", "Fat4", "Clstn1", "Clstn2", "Clstn3", "Pcdh15", "Ret", "Cdhr2", "Cdh23", "Cdhr3", "Cdhr4", "Cdhr5", "Dchs1", "Dchs2", "Fat1", "Fat2", "Celsr1", "Celsr2", "Celsr3", "Pcdha1", "Pcdha10", "Pcdha11", "Pcdha12", "Pcdha13", "Pcdha2", "Pcdha3", "Pcdha4", "Pcdha5", "Pcdha6", "Pcdha7", "Pcdha8", "Pcdha9", "Pcdhac1", "Pcdhac2", "Pcdhact", "Pcdhb1", "Pcdhb10", "Pcdhb11", "Pcdhb12", "Pcdhb13", "Pcdhb14", "Pcdhb15", "Pcdhb16", "Pcdhb17", "Pcdhb18", "Pcdhb19", "Pcdhb2", "Pcdhb3", "Pcdhb4", "Pcdhb5", "Pcdhb6", "Pcdhb7", "Pcdhb8", "Pcdhb9", "Pcdhg1", "Pcdhg10", "Pcdhg11", "Pcdhg12", "Pcdhg2", "Pcdhg3", "Pcdhg4", "Pcdhg5", "Pcdhg6", "Pcdhg7", "Pcdhg8", "Pcdhg9", "Pcdhga1", "Pcdhga10", "Pcdhga11", "Pcdhga12", "Pcdhga2", "Pcdhga3", "Pcdhga4", "Pcdhga5", "Pcdhga6", "Pcdhga7", "Pcdhga8", "Pcdhga9", "Pcdhgb1", "Pcdhgb2", "Pcdhgb3", "Pcdhgb4", "Pcdhgb5", "Pcdhgb6", "Pcdhgb7", "Pcdhgc3", "Pcdhgc4", "Pcdhgc5", "Pcdhgc")

clustered_protocadherins = c("Pcdha1", "Pcdha2", "Pcdha3", "Pcdha4", "Pcdha5", "Pcdha6", "Pcdha7", "Pcdha8", "Pcdha9", "Pcdha10", "Pcdha11", "Pcdha12", "Pcdha13", "Pcdha14", "Pcdhb1", "Pcdhb2", "Pcdhb3", "Pcdhb4", "Pcdhb5", "Pcdhb6", "Pcdhb7", "Pcdhb8", "Pcdhb9", "Pcdhb10", "Pcdhb11", "Pcdhb12", "Pcdhb13", "Pcdhb14", "Pcdhb15", "Pcdhb16", "Pcdhb17", "Pcdhb18", "Pcdhb19", "Pcdhg2", "Pcdhga1", "Pcdhga2", "Pcdhga3", "Pcdhga4", "Pcdhga5", "Pcdhga6", "Pcdhga7", "Pcdhga8", "Pcdhga9", "Pcdhga10", "Pcdhga11", "Pcdhga12", "Pcdhgb1", "Pcdhgb2", "Pcdhgb3", "Pcdhgb4", "Pcdhgb5", "Pcdhgb6", "Pcdhgb7", "Pcdhgb8", "Pcdhgb9", "Pcdhgc3", "Pcdhgc4", "Pcdhgc5")

unclustered_protos = c("Pcdh1", "Pcdh7", "Pcdh8", "Pcdh9", "Pcdh10", "Pcdh11x", "Pcdh11y", "Pcdh12", "Pcdh17", "Pcdh18", "Pcdh19", "Pcdh20")

major_cadherins = c("Cdh16", "Cdh17", "Celsr1", "Celsr2", "Celsr3", "Dsc1", "Dsc2", "Dsc3", "Dsg1", "Dsg2", "Dsg3", "Dsg4", "Cdh13", "Cdh26", "Cdh1", "Cdh15", "Cdh2", "Cdh3", "Cdh4", "Cdh10", "Cdh11", "Cdh12", "Cdh18", "Cdh19", "Cdh20", "Cdh22", "Cdh24", "Cdh5", "Cdh6", "Cdh7", "Cdh8", "Cdh9")

hsp = c("Bbs10","Bbs12","Tcp1","Cct2","Cct3","Cct4","Cct5","Cct6a","Cct6b","Cct7","Cct8","Clpb","Hspd1","Hspe1","Mkks","Dnaja1","Dnaja2","Dnaja3","Dnaja4","Dnajb1","Dnajb11","Dnajb12","Dnajb13","Dnajb14","Dnajb2","Dnajb3","Dnajb4","Dnajb5","Dnajb6","Dnajb7","Dnajb8","Dnajb9","Dnajc1","Dnajc10","Dnajc11","Dnajc12","Dnajc13","Dnajc14","Dnajc15","Dnajc16","Dnajc17","Dnajc18","Dnajc19","Dnajc2","Hscb","Dnajc21","Dnajc22","Sec63","Dnajc24","Dnajc25","Gak","Dnajc27","Dnajc28","Sacs","Dnajc3","Dnajc30","Dnajc4","Dnajc5","Dnajc5b","Dnajc5g","Dnajc6","Dnajc7","Dnajc8","Dnajc9","Hspa12a","Hspa12b","Hspa13","Hspa14","Hspa1a","Hspa1b","Hspa1l","Hspa2","Hspa4","Hspa4l","Hspa5","Hspa6","Hspa7","Hspa8","Hspa9","Hsph1","Hyou1","Hsp90aa1","Hsp90aa3p","Hsp90ab1","Hsp90b1","Trap1","Hspb1","Odf1","Ift25","Hspb2","Hspb3","Cryaa"
)

library(RColorBrewer)
library(ggrepel)

head(out$region[grep("mt-", out$X)])

length(rownames(x = obj))

rbps = c(rownames(x = obj)[grep("Rps", rownames(x = obj))] , rownames(x = obj)[grep("Rpl", rownames(x = obj))])

mts = c(rownames(x = obj)[grep("mt-", rownames(x = obj))] , rownames(x = obj)[grep("mt-", rownames(x = obj))])

mts = unique(mts)
mts

obj = AddModuleScore(obj, mts, assay = "RNA", name = "mts")

options(repr.plot.width=10, repr.plot.height=4.5)

Idents(obj)="region"
VlnPlot(obj, features = "mt-Co3", assay= "RNA", split.by="sample", pt.size = 0)

library(RColorBrewer)
library(ggrepel)

options(repr.plot.width=10, repr.plot.height=5)

fam = "Hdac1"

out$family = ""
out$family[grep("Hdac1$", out$X)]="Hdac1"
out$family[grep("Hdac11", out$X)]=""

#out$family[which(out$X %in% fam)]=fam
#out$family[which(out$X %in% unclustered_protos)]="unclustered_protocadherins"
#out$family[which(out$X %in% major_cadherins)]="major_cadherins"

#out$family[grep("Atp", out$X)]="Atp"
#out$family[grep("Cdh", out$X)]="Cdh"

out$family = factor(out$family, levels = c(  fam,""))
out = out[order(out$family,decreasing = T),]
#out = out[order(out$p_val,decreasing = F),]

sig = out[which(out$p_val_adj<1e-2 & out$family!= "" & abs(out$avg_log2FC)>.1),]
#sig = out[which(out$X==" "),]




ggplot(out, aes(x = (out$`pct.1`+out$`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +  
#geom_text_repel(data = sig, aes(label = sig$celltype, x = (sig$`pct.1`+sig$`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3) +
facet_wrap(~clade,nrow = 1)

options(repr.plot.width=16, repr.plot.height=5)

ggplot(out, aes(x = (out$`pct.1`+out$`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +   
geom_text_repel(data = sig, aes(label = sig$X, x = (sig$`pct.1`+sig$`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3)+
facet_wrap(~region, nrow = 1)

options(repr.plot.width=16, repr.plot.height=5)

ggplot(out, aes(x = avg_log2FC, y = -log10(p_val_adj+1e-100), color = family)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c(brewer.pal(3, "Set1")[1], "lightgrey")) +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Log Fold Change", y = "-log10(Adjusted p-value)", color = "Gene Family") +
  geom_text_repel(data = sig, aes(label = sig$X, x = avg_log2FC, y = -log10(p_val_adj)), 
                  color = "black", size = 3) +
  facet_wrap(~region, nrow = 1)

options(repr.plot.width=10, repr.plot.height=4.5)

fam = "Ets"

out$family = ""
out$family[grep("Elk", out$X)]="Ets"
out$family[grep("Ets", out$X)]="Ets"
out$family[grep("Elf", out$X)]="Ets"
out$family[grep("Fli1", out$X)]="Ets"
out$family[grep("Erg$", out$X)]="Ets"
out$family[grep("Elfn", out$X)]=""

#out$family[which(out$X %in% fam)]=fam
#out$family[which(out$X %in% unclustered_protos)]="unclustered_protocadherins"
#out$family[which(out$X %in% major_cadherins)]="major_cadherins"

#out$family[grep("Atp", out$X)]="Atp"
#out$family[grep("Cdh", out$X)]="Cdh"

out$family = factor(out$family, levels = c(  fam,""))
out = out[order(out$family,decreasing = T),]
#out = out[order(out$p_val,decreasing = F),]

sig = out[which(out$p_val_adj<1e-2 & out$family!= "" & abs(out$avg_log2FC)>.1),]
#sig = out[which(out$X==" "),]




ggplot(out, aes(x = (out$`pct.1`+out$`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +  
geom_text_repel(data = sig, aes(label = sig$X, x = (sig$`pct.1`+sig$`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3) +
facet_wrap(~clade,nrow = 1)

options(repr.plot.width=16, repr.plot.height=5)

ggplot(out, aes(x = (out$`pct.1`+out$`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +   
geom_text_repel(data = sig, aes(label = sig$X, x = (sig$`pct.1`+sig$`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3)+
facet_wrap(~region,nrow = 1)

library(ggrepel)
library(RColorBrewer)

brewer.pal(4, "Set1")

ggplot(out, aes(x = (out$`pct.1`+out$`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1"),"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +  
geom_text_repel(data = sig, aes(label = sig$X, x = (sig$`pct.1`+sig$`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3) +
facet_wrap(~clade,nrow = 1)


head(out)

myelin_sheath = fread("~/projects/combined_all/Figures/GO_term_summary_myelin_sheath.txt")

myelin_sheath = fread("~/projects/combined_all/Figures/GO_term_summary_ensheathment_of_neurons.txt")

myelin_sheath = unique(myelin_sheath$`MGI Gene/Marker ID`)
myelin_sheath

options(repr.plot.width=10, repr.plot.height=4.5)


out$family = "other or NS"



out$family[which(out$X %in% myelin_sheath)]="myelin_sheath"
#out$family[which(out$X %in% unclustered_protos)]="unclustered_protocadherins"


#out$family[grep("Cdh", out$X)]="Cdh"
#out$family[grep("Zkscan", out$X)]="Zkscan"
out$family[which(out$p_val_adj>0.05)] = "other or NS"

out$family = factor(out$family, levels = c(  "myelin_sheath","other or NS"))
out = out[order(out$family,decreasing = T),]
#out = out[order(out$p_val,decreasing = F),]

sig = out[which(out$p_val_adj<1e-10 & out$family!= "other or NS" & abs(out$avg_log2FC)>1),]
head(sig)
#sig = out[which(out$X==" "),]




ggplot(out, aes(x = (`pct.1`+`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(4, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +  
geom_text_repel(data = sig,max.overlaps = 100, aes(label =X, x = (`pct.1`+`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3) +
facet_wrap(~clade,nrow = 1)

options(repr.plot.width=16, repr.plot.height=5)
nn = out[which(out$clade == "NN"),]
sign = sig[which(sig$clade=="NN"), ] 
ggplot(nn, aes(x = (`pct.1`+`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(4, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +   
geom_text_repel(data = sign,max.overlaps = 100, aes(label = X, x = (`pct.1`+`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3)+
facet_wrap(~region,nrow = 1)


glut = out[which(out$clade == "Glut"),]
sign = sig[which(sig$clade=="Glut"), ] 

ggplot(glut, aes(x = (`pct.1`+`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(4, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +   
geom_text_repel(data = sign,max.overlaps = 100, aes(label = X, x = (`pct.1`+`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3)+
facet_wrap(~region,nrow = 1)

gaba = out[which(out$clade == "Gaba"),]
sign = sig[which(sig$clade=="Gaba"), ] 

ggplot(gaba, aes(x = (`pct.1`+`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(4, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +   
geom_text_repel(data = sign, max.overlaps = 100, aes(label = X, x = (`pct.1`+`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3)+
facet_wrap(~region,nrow = 1)

nrow(nn)
nrow(nn_combined)
nrow(nn_filtered)

head(nn)
nn[which(nn$X== "Ntrk3"),]

options(repr.plot.width=10
        , repr.plot.height=5)

# Extract all points for Oligo_NN
nn = out[which(out$celltype == "Oligo_NN"),]

# Remove overlapping non-significant points with p-value > 0.5 or abs(logFC) < 0.1
nn_filtered = nn %>%
  filter(!(p_val > 0.05 | abs(avg_log2FC) < 0.25))  # Removing points with high p-values or low log fold change

significant_myelin_sheath = nn_filtered[which(nn_filtered$family!="other or NS"),]
significant_myelin_sheath = significant_myelin_sheath[which(abs(significant_myelin_sheath$avg_log2FC)>1),]

# Define top 5 significant myelin_sheath genes for labeling per region
significant_myelin_sheath = significant_myelin_sheath %>%
  group_by(region) %>%
  arrange(log10(p_val_adj)*abs(avg_log2FC)) %>%  
  slice_head(n = 4) %>%    # Select top 5 points per region
  ungroup()

#significant_myelin_sheath = significant_myelin_sheath[which(significant_myelin_sheath$family!="other or NS"),]

#significant_myelin_sheath = nn_filtered[grep("Sox10", nn_filtered$X),]
nn_filtered$family = as.character(nn_filtered$family)
nn_filtered$family[grep("Sox10", nn_filtered$X)] = "Sox10"
nn_filtered$family[grep("Neat1", nn_filtered$X)] = "Neat1"
nn_filtered$family[grep("Cdh8", nn_filtered$X)] = "Cdh8"
nn_filtered$family[grep("Hdac", nn_filtered$X)] = "Hdac"


nn_filtered$family[grep("myelin", nn_filtered$family)] = "CNS myelination genes"
nn_filtered$family  = factor(nn_filtered$family , levels = c( "Sox10","Snca","Neat1","Cdh8","CNS myelination genes","other or NS"))

# Plot the filtered dataset with significant myelin_sheath genes labeled
myelin_oligo = ggplot(nn_filtered, aes(x = (`pct.1` + `pct.2`) / 2, y = -avg_log2FC, color = family)) +
  geom_point(alpha = 0.7) + 
  scale_color_manual(values = c(brewer.pal(4, "Set1")[1:5], "lightgrey")) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed") + ylim(-6.5, 6.5) +
  labs(x = "pct.cells", y = "Log Fold Change (18mo/2mo)", color = "") +   
 # geom_point(data = nn_filtered[which(nn_filtered$family=="Sox10"),], aes(x = (`pct.1` + `pct.2`) / 2, y = -avg_log2FC), color = brewer.pal(4, "Set1")[2]) +

  geom_text_repel(nudge_y = -0.05,data = significant_myelin_sheath, max.overlaps = 25, aes(label = X, x = (`pct.1` + `pct.2`) / 2, y = -avg_log2FC), color = "black", size = 3) +
  facet_wrap(~region, ncol = 4) + 
  ggtitle("Oligo DEGs MA") +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1))

# Display plot
print(myelin_oligo)

options(repr.plot.width=10
        , repr.plot.height=5)

# Extract all points for Oligo_NN
nn = out[which(out$celltype == "Oligo_NN"),]

# Remove overlapping non-significant points with p-value > 0.5 or abs(logFC) < 0.1
nn_filtered = nn %>%
  filter(!(p_val > 0.05 | abs(avg_log2FC) < 0.25))  # Removing points with high p-values or low log fold change

significant_myelin_sheath = nn_filtered[which(nn_filtered$family!="other or NS"),]
significant_myelin_sheath = significant_myelin_sheath[which(abs(significant_myelin_sheath$avg_log2FC)>1),]

# Define top 5 significant myelin_sheath genes for labeling per region
significant_myelin_sheath = significant_myelin_sheath %>%
  group_by(region) %>%
  arrange(log10(p_val_adj)*abs(avg_log2FC)) %>%  
  slice_head(n = 4) %>%    # Select top 5 points per region
  ungroup()

#significant_myelin_sheath = significant_myelin_sheath[which(significant_myelin_sheath$family!="other or NS"),]

#significant_myelin_sheath = nn_filtered[grep("Sox10", nn_filtered$X),]
nn_filtered$family = as.character(nn_filtered$family)
#nn_filtered$family[grep("Sox10", nn_filtered$X)] = "Sox10"
#nn_filtered$family[grep("Neat1", nn_filtered$X)] = "Neat1"
#nn_filtered$family[grep("Cdh8", nn_filtered$X)] = "Cdh8"
nn_filtered$family[grep("Gtf3", nn_filtered$X)] = "Klf"


nn_filtered$family[grep("myelin", nn_filtered$family)] = "CNS myelination genes"
nn_filtered$family  = factor(nn_filtered$family , levels = (c( "Sox10","Klf","Neat1","Cdh8","CNS myelination genes","other or NS")))

# Plot the filtered dataset with significant myelin_sheath genes labeled
myelin_oligo = ggplot(nn_filtered, aes(x = (`pct.1` + `pct.2`) / 2, y = -avg_log2FC, color = family)) +
  geom_point(alpha = 0.7) + 
  scale_color_manual(values = c(brewer.pal(4, "Set1")[1:2], "lightgrey")) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed") + ylim(-6.5, 6.5) +
  labs(x = "pct.cells", y = "Log Fold Change (18mo/2mo)", color = "") +   
  geom_point(data = nn_filtered[which(nn_filtered$family=="Klf"),], aes(x = (`pct.1` + `pct.2`) / 2, y = -avg_log2FC), color = brewer.pal(4, "Set1")[1]) +

  geom_text_repel(nudge_y = -0.05,data = significant_myelin_sheath, max.overlaps = 25, aes(label = X, x = (`pct.1` + `pct.2`) / 2, y = -avg_log2FC), color = "black", size = 3) +
  facet_wrap(~region, ncol = 4) + 
  ggtitle("Oligo DEGs MA") +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1))

# Display plot
print(myelin_oligo)

Klf = nn[grep("Atf", nn$X),]
Klf[order(Klf$p_val_adj),]

myelin_oligo

ggsave(myelin_oligo, height = 5, width = 7, file = "~/projects/combined_all/Figures/Figure3-NN/Oligo_MA.pdf")


myelin_oligo

ggsave(myelin_oligo, height = 5, width = 7, file = "~/projects/combined_all/Figures/Figure3-NN/Oligo_myelin.pdf")


options(repr.plot.width=10, repr.plot.height=4.5)


out$family = "other or NS"



#out$family[which(out$X %in% fam)]=fam
#out$family[which(out$X %in% unclustered_protos)]="unclustered_protocadherins"
#out$family[which(out$X %in% major_cadherins)]="major_cadherins"
out$family[grep("Rps", out$X)]="Rps"
out$family[grep("Atp", out$X)]="Atp"
out$family[grep("Nduf", out$X)]="Nduf"
out$family[grep("Cox", out$X)]="Cox"
out$family[grep("Mettl", out$X)]="Mettl"
out$family[grep("Ctbp", out$X)]="Ctbp"
out$family[grep("Kmt2", out$X)]="Kmt2"
out$family[grep("Yy1", out$X)]="Yy1"



#out$family[grep("Cdh", out$X)]="Cdh"
#out$family[grep("Zkscan", out$X)]="Zkscan"
out$family[which(out$p_val_adj>0.2)] = "other or NS"

out$family = factor(out$family, levels = c("Mettl", "Ctbp", "Kmt2","Yy1","Rps","Atp", "Nduf", "Cox","other or NS"))
out = out[order(out$family,decreasing = T),]
#out = out[order(out$p_val,decreasing = F),]

sig = out[which(out$p_val_adj<1e-5 & out$family!= "" & abs(out$avg_log2FC)>4),]
#sig = out[which(out$X==" "),]




ggplot(out, aes(x = (`pct.1`+`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(8, "Set1"),"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +  
#geom_text_repel(data = sig, aes(label = sig$X, x = (sig$`pct.1`+sig$`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3) +
facet_wrap(~clade,nrow = 1)

options(repr.plot.width=16, repr.plot.height=5)
nn = out[which(out$clade == "NN"),]
sign = sig[which(sig$clade=="NN"), ] 
ggplot(nn, aes(x = (`pct.1`+`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(8, "Set1"),"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +   
#geom_text_repel(data = sign, aes(label = X, x = (`pct.1`+`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3)+
facet_wrap(~region,nrow = 1)


glut = out[which(out$clade == "Glut"),]
ggplot(glut, aes(x = (`pct.1`+`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(8, "Set1"),"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +   
#geom_text_repel(data = sign, aes(label = X, x = (`pct.1`+`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3)+
facet_wrap(~region,nrow = 1)

gaba = out[which(out$clade == "Gaba"),]
ggplot(gaba, aes(x = (`pct.1`+`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(8, "Set1"),"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +   
#geom_text_repel(data = sign, aes(label = X, x = (`pct.1`+`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3)+
facet_wrap(~region,nrow = 1)



options(repr.plot.width=10, repr.plot.height=4.5)

fam = "Rps"

out$family = ""
out$family[grep("Rps", out$X)]="Rps"


#out$family[which(out$X %in% fam)]=fam
#out$family[which(out$X %in% unclustered_protos)]="unclustered_protocadherins"
#out$family[which(out$X %in% major_cadherins)]="major_cadherins"

#out$family[grep("Atp", out$X)]="Atp"
#out$family[grep("Cdh", out$X)]="Cdh"

out$family = factor(out$family, levels = c(  fam,""))
out = out[order(out$family,decreasing = T),]
#out = out[order(out$p_val,decreasing = F),]

sig = out[which(out$p_val_adj<1e-2 & out$family!= "" & abs(out$avg_log2FC)>.1),]
#sig = out[which(out$X==" "),]




ggplot(out, aes(x = (out$`pct.1`+out$`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +  
geom_text_repel(data = sig, aes(label = sig$X, x = (sig$`pct.1`+sig$`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3) +
facet_wrap(~clade,nrow = 1)

options(repr.plot.width=16, repr.plot.height=5)

ggplot(out, aes(x = (out$`pct.1`+out$`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +   
geom_text_repel(data = sig, aes(label = sig$X, x = (sig$`pct.1`+sig$`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3)+
facet_wrap(~region,nrow = 1)

table(out$family, out$dir)

nadh = c('Sdha', 'Sdhb', 'Sdhc', 'Sdhd', 'Uqcrc1', 'Uqcr10', 'Uqcr11', 'Uqcrc2', 'Mt-cyb', 'Cyc1', 'Uqcrfs1', 'Uqcrh', 'Uqcrb', 'Uqcrq', 'Cox4i1', 'Cox4i2', 'Cox5a', 'Cox5b', 'Cox6a1', 'Cox6a2', 'Cox6b1', 'Cox6b2', 'Cox6c', 'Cox7a1', 'Cox7a2', 'Cox7b', 'Cox7b2', 'Cox7c', 'Cox8a', 'Cox8c', 'Mt-co1', 'Mt-co2', 'Mt-co3', 'Atp5f1a', 'Atp5f1b', 'Atp5f1c', 'Atp5f1d', 'Atp5f1e', 'Atp5mc1', 'Atp5mc2', 'Atp5mc3', 'Atp5me', 'Atp5mf', 'Atp5mg', 'Atp5mj', 'Atp5mk', 'Mt-atp6', 'Mt-atp8', 'Atp5pb', 'Atp5pd', 'Atp5pf', 'Atp5po', 'Atp5if1', 'Mt-nd1', 'Mt-nd2', 'Mt-nd3', 'Mt-nd4', 'Mt-nd4l', 'Mt-nd5', 'Mt-nd6', 'Ndufs1', 'Ndufs2', 'Ndufs3', 'Ndufs7', 'Ndufs8', 'Ndufv1', 'Ndufv2', 'Ndufa1', 'Ndufa10', 'Ndufa11', 'Ndufa12', 'Ndufa13', 'Ndufa2', 'Ndufa3', 'Ndufa5', 'Ndufa6', 'Ndufa7', 'Ndufa8', 'Ndufa9', 'Ndufab1', 'Ndufb1', 'Ndufb10', 'Ndufb11', 'Ndufb2')

hdac = c('Dnttip1', 'Hdac1', 'Hdac2', 'Mideas', 'Gps2', 'Hdac3', 'Ncor1', 'Ncor2', 'Tbl1x', 'Tbl1xr1', 'Chd3', 'Chd4', 'Chd5', 'Hdac1', 'Hdac2', 'Mbd2', 'Mbd3', 'Mta1', 'Mta2', 'Mta3', 'Rbbp4', 'Rbbp7', 'Hdac1', 'Hdac2', 'Rbbp4', 'Rbbp7', 'Sap18', 'Sap30', 'Sin3a', 'Sin3b', 'Suds3')

lyt = c('Suv39h1', 'Suv39h2', 'Ehmt2', 'Ehmt1', 'Setdb1', 'Setdb2', 'Kmt2a', 'Kmt2b', 'Kmt2c', 'Kmt2d', 'Kmt2e', 'Setd1a', 'Setd1b', 'Ash1l', 'Setd2', 'Nsd1', 'Smyd2', 'Smyd1', 'Smyd3', 'Nsd3', 'Nsd2', 'Dot1l', 'Kmt5a', 'Kmt5b', 'Kmt5c', 'Ezh2', 'Ezh1', 'Setd7', 'Prdm2', 'Prdm9', 'Prdm6', 'Prdm8', 'Mecom', 'Prdm16')

ribo = c('Rpl10', 'Rpl10a', 'Rpl10l', 'Rpl11', 'Rpl12', 'Rpl13', 'Rpl13a', 'Rpl14', 'Rpl15', 'Rpl17', 'Rpl18', 'Rpl18a', 'Rpl19', 'Rpl21', 'Rpl22', 'Rpl22l1', 'Rpl23', 'Rpl23a', 'Rpl24', 'Rpl26', 'Rpl26l1', 'Rpl27', 'Rpl27a', 'Rpl28', 'Rpl29', 'Rpl3', 'Rpl30', 'Rpl31', 'Rpl32', 'Rpl34', 'Rpl35', 'Rpl35a', 'Rpl36', 'Rpl36a', 'Rpl36al', 'Rpl37', 'Rpl37a', 'Rpl38', 'Rpl39', 'Rpl39l', 'Rpl3l', 'Rpl4', 'Uba52', 'Rpl41', 'Rpl5', 'Rpl6', 'Rpl7', 'Rpl7a', 'Rpl7l1', 'Rpl8', 'Rpl9', 'Rplp0', 'Rplp1', 'Rplp2', 'Mrpl1', 'Mrpl10', 'Mrpl11', 'Mrpl12', 'Mrpl13', 'Mrpl14', 'Mrpl15', 'Mrpl16', 'Mrpl17', 'Mrpl18', 'Mrpl19', 'Mrpl2', 'Mrpl20', 'Mrpl21', 'Mrpl22', 'Mrpl23', 'Mrpl24', 'Mrpl27', 'Mrpl28', 'Mrpl3', 'Mrpl30', 'Mrpl32', 'Mrpl33', 'Mrpl34', 'Mrpl35', 'Mrpl36', 'Mrpl37', 'Mrpl38', 'Mrpl39', 'Mrpl4', 'Mrpl40', 'Mrpl41', 'Mrpl42', 'Mrpl43', 'Mrpl44', 'Mrpl45', 'Mrpl46', 'Mrpl47', 'Mrpl48', 'Mrpl49', 'Mrpl50', 'Mrpl51', 'Mrpl52', 'Mrpl53', 'Mrpl54', 'Mrpl55', 'Mrpl58', 'Mrpl57', 'Gadd45gip1', 'Mrps30', 'Mrps18a', 'Mrpl1', 'Mrpl10', 'Mrpl11', 'Mrpl12', 'Mrpl13', 'Mrpl14', 'Mrpl15', 'Mrpl16', 'Mrpl17', 'Mrpl18', 'Mrpl19', 'Mrpl2', 'Mrpl20', 'Mrpl21', 'Mrpl22', 'Mrpl23', 'Mrpl24', 'Mrpl27', 'Mrpl28', 'Mrpl3', 'Mrpl30', 'Mrpl32', 'Mrpl33', 'Mrpl34', 'Mrpl35', 'Mrpl36', 'Mrpl37', 'Mrpl38', 'Mrpl39', 'Mrpl4', 'Mrpl40', 'Mrpl41', 'Mrpl42', 'Mrpl43', 'Mrpl44', 'Mrpl45', 'Mrpl46', 'Mrpl47', 'Mrpl48', 'Mrpl49', 'Mrpl50', 'Mrpl51', 'Mrpl52', 'Mrpl53', 'Mrpl54', 'Mrpl55', 'Mrpl58', 'Mrpl57', 'Gadd45gip1', 'Mrps30', 'Mrps18a', 'Mrpl9', 'Mrps10', 'Mrps11', 'Mrps12', 'Mrps14', 'Mrps15', 'Mrps16', 'Mrps17', 'Mrps18c', 'Mrps2', 'Mrps21', 'Mrps22', 'Mrps23', 'Mrps24', 'Mrps25', 'Mrps26', 'Mrps27', 'Mrps28', 'Dap3', 'Mrps31', 'Mrps33', 'Mrps34', 'Mrps35', 'Chchd1', 'Aurkaip1', 'Ptcd3', 'Mrps18b', 'Mrps5', 'Mrps6', 'Mrps7', 'Mrps9', 'Rps10', 'Rps11', 'Rps12', 'Rps13', 'Rps14', 'Rps15', 'Rps15a', 'Rps16', 'Rps17', 'Rps18', 'Rps19', 'Rps2', 'Rps20', 'Rps21', 'Rps23', 'Rps24', 'Rps25', 'Rps26', 'Rps27', 'Rps27a', 'Rps27l', 'Rps28', 'Rps29', 'Rps3', 'Fau', 'Rps3a', 'Rps4x', 'Rps4y1', 'Rps4y2', 'Rps5', 'Rps6', 'Rps7', 'Rps8', 'Rps9', 'Rpsa', 'Mrps10', 'Mrps11', 'Mrps12', 'Mrps14', 'Mrps15', 'Mrps16', 'Mrps17', 'Mrps18c', 'Mrps2', 'Mrps21', 'Mrps22', 'Mrps23', 'Mrps24', 'Mrps25', 'Mrps26', 'Mrps27', 'Mrps28', 'Mrps31', 'Mrps33', 'Mrps34', 'Mrps35', 'Chchd1', 'Aurkaip1', 'Ptcd3', 'Mrps18b'
)

mito = c("Flvcr1", "Macrod1", "Wwox", "Oxr1", "Hivep1", "Pnkp", "Pgam5", "Mtch2", "Adcy10", "Prodh", "Uqcrq", "Arglu1", "Mpv17l", "Mrps10", "Plpbp", "Immp1l", "Dnajc19", "Ndufa10", "Lrrk2", "Crebzf", "Mrpl2", "Miga2", "Ckmt2", "Tmem186", "Pstk", "Ndrg4", "Mrpl41", "Dnm3", "P2ry12", "Dact2", "Tmem126a", "Emc2", "Taco1", "Cox7b", "Rmdn1", "Cox6a1", "Timmdc1", "Slc25a27", "Sqstm1", "Aifm3", "Cwc15", "Srebf2")

sen = c("Cdkn2a", )

obj=readRDS("../RNA_final_SCT.RDS")
obj

head(out)

library(Seurat)

obj

Idents(obj) = "celltype_final"
ae = AverageExpression(obj) 

colnames(ae$RNA) = gsub("/", "-", colnames(ae$RNA))

unique(out$celltype)[which(!unique(out$celltype)%in%colnames(ae$RNA))]

out = out[-which(out$celltype %in% unique(out$celltype)[which(!unique(out$celltype)%in%colnames(ae$RNA))]),]

colnames(ae$RNA)

out$average_expression <- ae$RNA[cbind(out$X, out$celltype)]
head(out)

out$logSum=geneLogSums[paste(out$X)]

display.brewer.all()

options(repr.plot.width=8, repr.plot.height=3.5)

fam = "nadh"

out$family = ""
out$family[which(out$X %in% nadh)]=fam
out$family[which(out$X %in% ribo)]="ribo"

#out$family[grep( "mt-", out$X)]=fam

#out$family[which(out$X %in% unclustered_protos)]="unclustered_protocadherins"
#out$family[which(out$X %in% major_cadherins)]="major_cadherins"

#out$family[grep("Atp", out$X)]="Atp"
#out$family[grep("Cdh", out$X)]="Cdh"

out$family = factor(out$family, levels = c(  fam,"ribo",""))
out = out[order(out$family,decreasing = T),]
#out = out[order(out$p_val,decreasing = F),]

sig = out[which(out$p_val_adj<1e-5 & out$family!= "" & abs(out$avg_log2FC)>2.5),]
sig = out[which(out$X==" "),]

out$clade = factor(out$clade , levels = c("Glut", "Gaba", "NN"))


ggplot(out, aes(x = log(average_expression+1), y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.5) + scale_color_manual(values=c(brewer.pal(3, "Set1")[c(1,3)],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ ylim(-10,10)+
  labs(x = "average expression", y = "Log Fold Change", color = "Gene Family") +  
geom_text_repel(data = sig, aes(label = sig$X, x = log(average_expression+1), y = -avg_log2FC), color = "black", size = 3,max.overlaps = 50) +
facet_wrap(~clade,nrow = 1)

options(repr.plot.width=15, repr.plot.height=4)

ggplot(out, aes(x = log(average_expression+1), y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.5) + scale_color_manual(values=c(brewer.pal(3, "Set1")[c(1,3)],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ylim(-10,10)+
  labs(x = "average expression", y = "Log Fold Change", color = "Gene Family") +   
geom_text_repel(data = sig, aes(label = sig$X, x = log(average_expression+1), y = -avg_log2FC), color = "black", size = 3,max.overlaps = 50)+
facet_wrap(~region,nrow = 1)


library(RColorBrewer)
options(repr.plot.width=15, repr.plot.height=14)
colors <- (colorRampPalette(brewer.pal(9, "Oranges")[4:9])(100))
#breaks <- c(0, 5,10, 50,100,200)


cdh8 = out[which(out$X=="Sox4" & out$p_val_adj<0.1),]
ggplot(cdh8, aes(x = celltype, y = -avg_log2FC, fill = -log10(p_val_adj))) +
  geom_bar(stat = "identity") +  
  scale_fill_gradientn(colors = colors) + theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = "Oligo Cdh8 Age LogFC (18mo/2mo)", x = "Region", y = "Avg LogFC (18mo/2mo)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+facet_wrap(~region, scales = "free")  # Rotate x-axis labels for better readability

library(RColorBrewer)
options(repr.plot.width=8, repr.plot.height=5)
colors <- (colorRampPalette(brewer.pal(9, "Oranges")[4:9])(100))
#breaks <- c(0, 5,10, 50,100,200)


cdh8 = out[which(out$X=="E2f3" & out$p_val_adj<0.05),]
ggplot(cdh8, aes(x = celltype, y = -avg_log2FC, fill = -log10(p_val_adj))) +
  geom_bar(stat = "identity") +  
  scale_fill_gradientn(colors = colors) + theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = "E2f3 Age LogFC (18mo/2mo)", x = "Region", y = "Avg LogFC (18mo/2mo)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+facet_wrap(~region,nrow = 1, scales = "free_x"
                                                                            )  # Rotate x-axis labels for better readability

max(out$p_val_adj
   )

library(RColorBrewer)
options(repr.plot.width=4, repr.plot.height=4)
colors <- (colorRampPalette(brewer.pal(9, "Oranges")[4:9])(100))
#breaks <- c(0, 5,10, 50,100,200)


cdh8 = out[which(out$X=="Eif5a" & out$clade == "NN"),]
cdh8 = cdh8[-grep("IMN", cdh8$file),]

cdh8 = cdh8[which(cdh8$p_val<.25),]

g = ggplot(cdh8, aes(x = file, y = -avg_log2FC, fill = -log10(p_val_adj))) +
  geom_bar(stat = "identity") +  
  scale_fill_gradientn(colors = colors) + theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = "RNA Eif5a Age LogFC (18mo/2mo)", x = "Region", y = "Avg LogFC (18mo/2mo)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability
pdf("~/projects/combined_all/Figures/Figure3/eif5a_rna.pdf" , height=4, width = 4)
print(g)
dev.off()

g

options(repr.plot.width=5, repr.plot.height=4)

gene = "Jund"
cdh8 = out[which(out$X==gene ),]
ggplot(cdh8, aes(x = file, y = -avg_log2FC, fill = -log10(p_val_adj+1e-100))) +
  geom_bar(stat = "identity") +  
  scale_fill_gradientn(colors = colors) + theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = paste(gene,"Age LogFC (18mo/2mo)"), x = "Region", y = "Avg LogFC (18mo/2mo)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability


options(repr.plot.width=15, repr.plot.height=7)
colors <- (colorRampPalette(brewer.pal(9, "Oranges")[4:9])(100))
#breaks <- c(0, 5,10, 50,100,200)


cdh8 = out[which(out$X=="AC149090.1" ),]
ggplot(cdh8, aes(x = file, y = -avg_log2FC, fill = -log10(p_val_adj+1e-100))) +
  geom_bar(stat = "identity") +  
  scale_fill_gradientn(colors = colors) + theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = "Oligo Cdh8 Age LogFC (18mo/2mo)", x = "Region", y = "Avg LogFC (18mo/2mo)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability


gene = "Xist"
cdh8 = out[which(out$X==gene ),]
ggplot(cdh8, aes(x = file, y = -avg_log2FC, fill = -log10(p_val_adj+1e-100))) +
  geom_bar(stat = "identity") +  
  scale_fill_gradientn(colors = colors) + theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = paste(gene,"Age LogFC (18mo/2mo)"), x = "Region", y = "Avg LogFC (18mo/2mo)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability


head(cdh8)
ggplot(cdh8, aes(x = file, y = -avg_log2FC, fill = -log10(p_val_adj+1e-100))) +
  geom_bar(stat = "identity") +   theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = paste(gene,"Age LogFC (18mo/2mo)"), x = "Region", y = "Avg LogFC (18mo/2mo)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability


library(RColorBrewer)

colors <- (colorRampPalette(brewer.pal(9, "PuBu")[4:9])(100))

gene = "Il11"
cdh8 = out[which(out$X==gene ),]
ggplot(cdh8, aes(x = file, y = -avg_log2FC, fill = -log10(p_val_adj+1e-100))) +
  geom_bar(stat = "identity") +  
  scale_fill_gradientn(colors = colors) + theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = paste(gene,"Age LogFC (18mo/2mo)"), x = "Region", y = "Avg LogFC (18mo/2mo)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability


out$X[grep("Il11", out$X)]

options(repr.plot.width=5, repr.plot.height=5)


gene = "Crlf2"
cdh8 = out[which(out$X==gene ),]
ggplot(cdh8, aes(x = file, y = -avg_log2FC, fill = -log10(p_val_adj+1e-100))) +
  geom_bar(stat = "identity") +  
  scale_fill_gradientn(colors = colors) + theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = paste(gene,"Age LogFC (18mo/2mo)"), x = "Region", y = "Avg LogFC (18mo/2mo)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability


options(repr.plot.width=9, repr.plot.height=5)

colors <- (colorRampPalette(brewer.pal(9, "PuBu")[4:9])(100))

gene = "Prkn"
cdh8 = out[which(out$X==gene  ),]
ggplot(cdh8, aes(x = file, y = -avg_log2FC, fill = -log10(p_val_adj+1e-100))) +
  geom_bar(stat = "identity") +  
  scale_fill_gradientn(colors = colors) + theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = paste(gene,"Age LogFC (18mo/2mo)"), x = "Region", y = "Avg LogFC (18mo/2mo)", fill = "-log10(adj p-value)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability



gene = "Gm47283"
cdh8 = out[which(out$X==gene ),]
ggplot(cdh8, aes(x = file, y = -avg_log2FC, fill = -log10(p_val_adj+1e-100))) +
  geom_bar(stat = "identity") +  
  scale_fill_gradientn(colors = colors) + theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = paste(gene,"Age LogFC (18mo/2mo)"), x = "Region", y = "Avg LogFC (18mo/2mo)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability


#obj$celltype_final[which(obj$final_clusters_final=="IOL")] = "IOL NN"

obj@assays$RNA$counts[1:10,1:10]
obj@assays$SCT$counts[1:10,1:10]

Idents(obj) = "celltype_final"
ae = AggregateExpression(obj, assays = "RNA",normalization.method = "RC", scale.factor = 1e6,return.seurat = TRUE)

ae

hist(ae@assays$RNA$data[,"Oligo NN"], breaks = 50000, xlim = c(0,.5))

length(which(ae@assays$SCT$data[,"Oligo NN"]>0))

ae@assays$SCT$scale.data[1:5,1:5]

options(repr.plot.width=10, repr.plot.height=5)


gene = "Rgn"
cdh8 = out[which(out$X==gene ),]
ggplot(cdh8, aes(x = file, y = -avg_log2FC, fill = -log10(p_val_adj+1e-100))) +
  geom_bar(stat = "identity") +  
  scale_fill_gradientn(colors = colors) + theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = paste(gene,"Age LogFC (18mo/2mo)"), x = "Region", y = "Avg LogFC (18mo/2mo)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability



head(mtab)



updiff = out[which(out$avg_log2FC>(.15) & out$p_val_adj<0.05),]
tab = (table(updiff$X,updiff$file)[order(rowSums(table(updiff$X,updiff$file)), decreasing =T),])
mtab = melt(tab[1:40,])
mtab = mtab[which(mtab$value>0),]
mtab$region = sapply(strsplit(as.character(mtab$Var2), "--"), `[`, 2)
mtab$celltype = sapply(strsplit(as.character(mtab$Var2), "--"), `[`, 1)
mtab$clade = sapply(strsplit(as.character(mtab$celltype), "_"), tail, n = 1)
mtab$clade[which(mtab$clade %in%c("Dopa", "Glut-Sero"))] = "Glut"
mtab$clade[which(mtab$clade %in%c("Neur"))] = "Gaba"


options(repr.plot.width=5, repr.plot.height=7)

ggplot(mtab, aes(x = Var1, fill=region)) +
  geom_bar() + theme_minimal()+coord_flip()

ggplot(mtab, aes(x = Var1, fill=clade)) +
  geom_bar() + theme_minimal()+coord_flip()+labs(x = "Gene", y = "Occurances DE down in aging", color = "Clade") 

updiff = out[which(out$avg_log2FC<(-.1)& out$p_val_adj<0.05),]
head(table(updiff$X)[order(table(updiff$X) , decreasing = T)],50)
up_tab = table(updiff$X)[order(table(updiff$X) , decreasing = T)]
hist(table(updiff$X)[order(table(updiff$X) , decreasing = T)])

updiff


gene = "Dohh"
cdh8 = out[which(out$X==gene ),]
ggplot(cdh8, aes(x = file, y = -avg_log2FC, fill = -log10(p_val_adj+1e-100))) +
  geom_bar(stat = "identity") +  
  scale_fill_gradientn(colors = colors) + theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = paste(gene,"Age LogFC (18mo/2mo)"), x = "Region", y = "Avg LogFC (18mo/2mo)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability


options(repr.plot.width=10, repr.plot.height=5)

fam = "ribo"

out$family = ""
out$family[which(out$X %in% ribo)]=fam
#out$family[grep( "mt-", out$X)]=fam

#out$family[which(out$X %in% unclustered_protos)]="unclustered_protocadherins"
#out$family[which(out$X %in% major_cadherins)]="major_cadherins"

#out$family[grep("Atp", out$X)]="Atp"
#out$family[grep("Cdh", out$X)]="Cdh"

out$family = factor(out$family, levels = c(  fam,""))
out = out[order(out$family,decreasing = T),]
#out = out[order(out$p_val,decreasing = F),]

sig = out[which(out$p_val_adj<1e-3 & out$family!= "" & abs(out$avg_log2FC)>1),]
#sig = out[which(out$X==" "),]




ggplot(out, aes(x = log(average_expression+1), y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ ylim(-10,10)+
  labs(x = "average expression", y = "Log Fold Change", color = "Gene Family") +  
geom_text_repel(data = sig, aes(label = sig$X, x = log(average_expression+1), y = -avg_log2FC), color = "black", size = 3) +
facet_wrap(~clade,nrow = 1)

options(repr.plot.width=16, repr.plot.height=6)

ggplot(out, aes(x = log(average_expression+1), y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ylim(-10,10)+
  labs(x = "average expression", y = "Log Fold Change", color = "Gene Family") +   
geom_text_repel(data = sig, aes(label = sig$X, x = log(average_expression+1), y = -avg_log2FC), color = "black", size = 3)+
facet_wrap(~region,nrow = 1)

options(repr.plot.width=26, repr.plot.height=26)

ggplot(out, aes(x = log(average_expression+1), y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +  
geom_text_repel(data = sig, aes(label = sig$X, x = log(average_expression+1), y = -avg_log2FC), color = "black", size = 3) +
facet_wrap(~celltype,nrow = 6)

options(repr.plot.width=10, repr.plot.height=5)

fam = "ribo"

out$family = ""
out$family[which(out$X %in% ribo)]=fam
#out$family[grep( "mt-", out$X)]=fam

#out$family[which(out$X %in% unclustered_protos)]="unclustered_protocadherins"
#out$family[which(out$X %in% major_cadherins)]="major_cadherins"

#out$family[grep("Atp", out$X)]="Atp"
#out$family[grep("Cdh", out$X)]="Cdh"

out$family = factor(out$family, levels = c(  fam,""))
out = out[order(out$family,decreasing = T),]
#out = out[order(out$p_val,decreasing = F),]

sig = out[which(out$p_val_adj<1e-3 & out$family!= "" & abs(out$avg_log2FC)>1),]
#sig = out[which(out$X==" "),]




ggplot(out, aes(x = log(average_expression+1), y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ ylim(-10,10)+
  labs(x = "average expression", y = "Log Fold Change", color = "Gene Family") +  
geom_text_repel(data = sig, aes(label = sig$X, x = log(average_expression+1), y = -avg_log2FC), color = "black", size = 3) +
facet_wrap(~clade,nrow = 1)

options(repr.plot.width=16, repr.plot.height=6)

ggplot(out, aes(x = log(average_expression+1), y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ylim(-10,10)+
  labs(x = "average expression", y = "Log Fold Change", color = "Gene Family") +   
geom_text_repel(data = sig, aes(label = sig$X, x = log(average_expression+1), y = -avg_log2FC), color = "black", size = 3)+
facet_wrap(~region,nrow = 1)

options(repr.plot.width=26, repr.plot.height=26)

ggplot(out, aes(x = log(average_expression+1), y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +  
geom_text_repel(data = sig, aes(label = sig$X, x = log(average_expression+1), y = -avg_log2FC), color = "black", size = 3) +
facet_wrap(~celltype,nrow = 6)

options(repr.plot.width=16, repr.plot.height=16)

diff = out 
diff= diff[order(diff$family,decreasing = T),]
ggplot(diff, 
  aes(x = celltype, y = -avg_log2FC, color = family))+
  #color = ifelse(p_val_adj < 0.01, ifelse(avg_log2FC > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) + ylim(c(-10,10))+
  labs(title = "Age Differential Expression",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_minimal() + coord_flip() + facet_wrap(~region, scales = "free", ncol = 3)

options(repr.plot.width=26, repr.plot.height=26)

ggplot(out, aes(x = (out$`pct.1`+out$`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +  
geom_text_repel(data = sig, aes(label = sig$X, x = (sig$`pct.1`+sig$`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3) +
facet_wrap(~celltype,nrow = 6)

options(repr.plot.width=10, repr.plot.height=5)

out$family = ""
out$family[which(out$X %in% clustered_protocadherins)]="clustered_protocadherins"
out$family[which(out$X %in% unclustered_protos)]="unclustered_protocadherins"
out$family[which(out$X %in% major_cadherins)]="major_cadherins"

#out$family[grep("Atp", out$X)]="Atp"
#out$family[grep("Cdh", out$X)]="Cdh"

out$family = factor(out$family, levels = c(  "clustered_protocadherins","major_cadherins","unclustered_protocadherins",""))
out = out[order(out$family,decreasing = T),]
#out = out[order(out$p_val,decreasing = F),]

#sig = out[which(out$p_val_adj<1e-80 & out$family!= "" & abs(out$avg_log2FC)>1),]
sig = out[which(out$X==" "&out$celltype=="Oligo_NN"),]




ggplot(out, aes(x = (out$`pct.1`+out$`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1:3],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +  
geom_text_repel(data = sig, aes(label = sig$X, x = (sig$`pct.1`+sig$`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3) +
facet_wrap(~clade,nrow = 1)

options(repr.plot.width=16, repr.plot.height=6)

ggplot(out, aes(x = (out$`pct.1`+out$`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1:3],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +   
geom_text_repel(data = sig, aes(label = sig$X, x = (sig$`pct.1`+sig$`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3)+
facet_wrap(~region,nrow = 1)

out$family = ""
out$family[grep("Rps", out$X)]="Rps"
out$family[grep("Rpl", out$X)]="Rpl"
out$family[grep("Mrpl", out$X)]="Mrpl"

out$family[grep("Nduf", out$X)]="Nduf"
#out$family[grep("Pcdhg", out$X)]="Pcdhg"
#out$family[grep("Olf", out$X)]="Olf"
out$family[grep("Atp", out$X)]="Atp"
#out$family[grep("Cdh", out$X)]="Cdh"
#out$family[grep("Prdm", out$X)]="Prdm"

#out$family[grep("Hist1h", out$X)]="Hist1h"


unique(out$X[grep("Mir", out$X)])

head(out[which(out$X=="Rnasek"),],20)

out$family = factor(out$family, levels = c(  "Mrpl", "Rpl", "Rps", "Nduf", "Atp", ""))

out$family = factor(out$family, levels = c( "Cdh", "Pcdhg", "Rpl", "Rps", "Nduf", "Atp", ""))

table(out$family,out$dir)

out = out[order(out$family,decreasing = T),]

library(RColorBrewer)


options(repr.plot.width=12, repr.plot.height=6)


ggplot(out, aes(x = (out$`pct.1`+out$`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.3) + scale_color_manual(values=c(brewer.pal(5, "Set1"),"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") + facet_wrap(~clade,nrow = 1)

head(out)


