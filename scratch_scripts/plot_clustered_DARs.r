library(data.table)
library(ggplot2)

setwd(".")

files = list.files(path=".",pattern="smth20_scores",recursive=F,full.names=T)
files = files[grepl("300",files)]
files

#files = files[grepl("up",files)]

file = files[2]
chr.len = read.table("/projects/ps-renlab/share/bwa_indices/mm10.fa.fai")
chr.len$color = c(rep(c("black","grey"),10),"black")
colnames(chr.len)[1] = "chr"

chr.window = NULL
size = 10^5
for (i in 1:nrow(chr.len)){
    chr.window = c(chr.window, paste0(chr.len[i,1],":",0:ceiling(chr.len[i,2]/size) ))
    }

diff_dict = list()


diff_dict = list()

for(file in files){
print(file)
a=data.frame(fread(file))
a$pos = factor(paste0(a$chr,":",a$window),levels=chr.window)
a$ct = strsplit(file, ":")[[1]][1]
a$direct = strsplit(file, "_")[[1]][4]
#a$color = chr.len$color[match(a$chr,chr.len$chr)]
a$color = ifelse(grepl("Up",file),"red","blue")
name = paste(strsplit(file, "/")[[1]][2], strsplit(file, "_")[[1]][4])

diff_dict[[name]] = a
}


out = do.call(rbind,diff_dict)
out = out[!is.na(out$pos),]
out$posN = as.numeric(out$pos)
out$chr = factor(out$chr, levels=paste0("chr",c(1:19,"X","Y")))


sig_up = read.delim("clustered_diff_peaksUp_0.001.txt")
sig_up$chr = sub("(.*):(.*)-(.*)","\\1",sig_up$locations)
sig_up$pos1 = as.numeric(as.character(sub("(.*):(.*)-(.*)","\\2",sig_up$locations)))
sig_up$pos2 = as.numeric(sub("(.*):(.*)-(.*)","\\3",sig_up$locations))
#sig_up = sig_up[-which(is.na(sig_up$pos1)),]

sig_up$pos = factor(paste0(sig_up$chr, ":",as.numeric(ceiling( (sig_up$pos1+sig_up$pos2)/2/size)), sep = ""), levels=chr.window)
sig_up$posN = factor(as.numeric(ceiling( (sig_up$pos1+sig_up$pos2)/2/size)))

sig_up$chr = factor(sig_up$chr, levels=paste0("chr",c(1:19,"X","Y")))
## plot H3K9me3 domain. 
k9 = read.table("/mnt/tscc/yanxiao/projects/mouse_aging/data/encode/histone/peaks/FB.P0.H3K9me3-W1000-G3000-FDR0.01-island.bed")
colnames(k9) = c("chr","start","end","score")
k9$start_pos = factor(paste0(k9$chr,":",ceiling(k9$start/size)),levels=chr.window)
k9$end_pos = factor(paste0(k9$chr,":",ceiling(k9$end/size)),levels=chr.window)
k9$size = k9$end-k9$start
k9$start_posN = as.numeric(k9$start_pos)
k9$end_posN = as.numeric(k9$end_pos)
k9 = k9[which(k9$size>=100e3),]
k9$chr = factor(k9$chr, levels=paste0("chr",c(1:19,"X","Y")))


dim(out)
out$posN=as.numeric(out$posN)
out$smth20_score=as.numeric(out$smth20_score)

sig_up$posN=as.numeric(sig_up$posN)
sig_up$smth20_score=as.numeric(sig_up$smth20_score)
out$clade = "NN"
out$clade[grep("Glut", out$cell.type)] = "Glut"
out$clade[grep("Gaba", out$cell.type)] = "Gaba"
out$clade[grep("Neur", out$cell.type)] = "Gaba"
head(out)


table(out$chr)

options(repr.plot.width=15, repr.plot.height=12)

ggplot(out[which(out$smth20_score>0.2),]) +
  geom_line(aes(x=posN,y=smth20_score,color=color),alpha=0.5) + 
  scale_color_manual(values=c("blue","red")) + 
#  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
#      geom_jitter(data=sig_up,aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +

    geom_hline(yintercept=0.4,linetype="dashed") +
  ylim(-0.5,2) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
  axis.text.x = element_blank(),
  panel.spacing = unit(0, "lines"),
  panel.grid = element_blank(),
  panel.border = element_rect(colour="grey")
  )+
  facet_grid(clade~chr,scales="free_x",space="free_x")  +  geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 





options(repr.plot.width=10
        , repr.plot.height=5)

for(ct in unique(out$cell.type)){
if(nrow(out[which(out$smth20_score>0.5 & out$cell.type == ct),])<1){
    next
    }
p=ggplot(out[which(out$smth20_score>0.2 & out$cell.type == ct),]) +
  geom_line(aes(x=posN,y=smth20_score,color=color),alpha=0.5) + 
  scale_color_manual(values=c("blue","red")) + 
#  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
#      geom_jitter(data=sig_up,aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +

    geom_hline(yintercept=0.4,linetype="dashed") +
  ylim(-0.5,2) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
  axis.text.x = element_blank(),
  panel.spacing = unit(0, "lines"),
  panel.grid = element_blank(),
  panel.border = element_rect(colour="grey")
  )+ggtitle(ct)+
  facet_grid(~chr,scales="free_x",space="free_x")  +  geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5)
    
    print(p)
}




options(repr.plot.width=15, repr.plot.height=12)

ggplot(out[which(out$smth20_score>0.2),]) +
  geom_line(aes(x=posN,y=smth20_score,color=cell.type),alpha=0.5) + 
 # scale_color_manual(values=c("blue","red")) + 
#  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
    geom_hline(yintercept=0.4,linetype="dashed") +
  ylim(-0.5,2) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
  axis.text.x = element_blank(),
  panel.spacing = unit(0, "lines"),
  panel.grid = element_blank(),
  panel.border = element_rect(colour="grey")
  )+
  facet_grid(clade~chr,scales="free_x",space="free_x") 

options(repr.plot.width=15, repr.plot.height=12)

ggplot(out[which(out$cell.type=="Oligo_NN:Male"),]) +
  geom_line(aes(x=posN,y=smth20_score,color=cell.type),alpha=0.5) + 
 # scale_color_manual(values=c("blue","red")) + 
#  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
    geom_hline(yintercept=0.4,linetype="dashed") +
  ylim(-0.5,2) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
  axis.text.x = element_blank(),
  panel.spacing = unit(0, "lines"),
  panel.grid = element_blank(),
  panel.border = element_rect(colour="grey")
  )+
  facet_grid(~chr,scales="free_x",space="free_x") 

# write the large k9 domains. 
write.table(k9[,c("chr","start","end","score","size")],"FB.P0.H3K9me3.gt100k_domain.bed",row.names=F,col.names=F,quote=F,sep="\t")



#png("DH_FC.diff_peaks.chr_location.png",units="in",width=16,height=6,res=300)
pdf("all_tissue_sep.diff_peaks.chr_location.pdf",height=10,width=16)
ggplot(out2) +
  geom_line(aes(x=posN,y=smth20_score,color=color),alpha=0.5) + 
  scale_color_manual(values=c("blue","red")) + 
#  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
  geom_jitter(data=sig_up,aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +
#  geom_text(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.35,label=pos),angle=45,color="red") +
  geom_hline(yintercept=0.4,linetype="dashed") +
  ylim(-0.5,2) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
  axis.text.x = element_blank(),
  panel.spacing = unit(0, "lines"),
  panel.grid = element_blank(),
  panel.border = element_rect(colour="grey")
  )+
  facet_grid(tissue~chr,scales="free_x",space="free_x") +  geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 
dev.off()

out2$clade = ckey[paste(out2$cell.type), "L1Annot"]
out2 = out2[-which(is.na(out2$clade)),]
out2$clade[which(out2$clade == "SERO;GABA")] = "GABA"

pdf("all_tissue_sep.diff_peaks.chr_location_ct_up_down_nh3_c2.pdf",height=7.5,width=12)
ggplot(out2[which(out2$direct=="Up"),]) +
  geom_line(aes(x=posN,y=smth20_score,color=clade),alpha=0.9) + 
  scale_color_manual(values=c("darkorange1","firebrick1", "blue3")) + 
    #  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
 # geom_jitter(data=sig_up,aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +
  #  geom_text(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.35,label=pos),angle=45,color="red") +
  geom_hline(yintercept=0.4,linetype="dashed") +
  ylim(0,3) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+
facet_grid(tissue~chr,scales="free_x",space="free_x") + ggtitle("Up in Aging")#+  geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 

ggplot(out2[which(out2$direct=="Up"),]) +
  geom_line(aes(x=posN,y=smth20_score,color=clade),alpha=0.9) + 
  scale_color_manual(values=c("darkorange1","firebrick1", "blue3")) + 
  #  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
  # geom_jitter(data=sig_up,aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +
  #  geom_text(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.35,label=pos),angle=45,color="red") +
  geom_hline(yintercept=0.4,linetype="dashed") +
  ylim(-0.5,3) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+
  facet_grid(tissue~chr,scales="free_x",space="free_x") + ggtitle("Up in Aging")+  geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 

ggplot(out2[which(out2$direct=="Down"),]) +
  geom_line(aes(x=posN,y=smth20_score,color=clade),alpha=0.9) + 
  scale_color_manual(values=c("darkorange1","firebrick1", "blue3")) + 
  #  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
  #  geom_jitter(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +
  #  geom_text(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.35,label=pos),angle=45,color="red") +
  geom_hline(yintercept=0.4,linetype="dashed") +
  ylim(0,3) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+
  facet_grid(tissue~chr,scales="free_x",space="free_x") +  ggtitle("Down in Aging") #+geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 

ggplot(out2[which(out2$direct=="Down"),]) +
  geom_line(aes(x=posN,y=smth20_score,color=clade),alpha=0.9) + 
  scale_color_manual(values=c("darkorange1","firebrick1", "blue3")) + 
  #  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
  #  geom_jitter(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +
  #  geom_text(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.35,label=pos),angle=45,color="red") +
  geom_hline(yintercept=0.4,linetype="dashed") +
  ylim(-0.7,3) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
         axis.text.x = element_blank(),
         panel.spacing = unit(0, "lines"),
         panel.grid = element_blank(),
         panel.border = element_rect(colour="grey")
  )+
facet_grid(tissue~chr,scales="free_x",space="free_x") +  ggtitle("Down in Aging") +geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=7) 
dev.off()

#out = out[which(out$chr=="chr9"),]
#sig = sig[which(sig$chr=="chr9"),]
#k9 = k9[which(k9$chr=="chr9"),]
#k9$start_pos = factor(paste0(k9$chr,":",ceiling(k9$start/size)),levels=chr.window)
#k9$end_pos = factor(paste0(k9$chr,":",ceiling(k9$end/size)),levels=chr.window)


#ggplot(out) +geom_line(aes(x=pos,y=smth20_score,color=color)) +
#  geom_point(data=subset(sig,direct=="Down"),aes(x=pos,y=-0.25),shape=25,color="blue",fill="blue") +
#  geom_point(data=subset(sig,direct=="Up"),aes(x=pos,y=-0.25),shape=24,color="red",fill="red") +
#  geom_segment(data=k9,aes(x=start_pos,y=-0.5,xend=end_pos,yend=-0.5),size=5)





#ggplot(sig) +geom_point(aes(x=as.numeric(pos),y=-1),shape=2)
#ggplot(sig) +geom_point(aes(x=as.numeric(pos),y=paste0(tissue,clust),color=tissue))



ggplot(out) +
  geom_line(aes(x=posN,y=smth20_score,color=color),alpha=0.5) + 
  scale_color_manual(values=c("blue","red"))

nrow(out)

while(3>1){dev.off()}

ggplot(out) +
  geom_line(aes(x=posN,y=smth20_score,color=color),alpha=0.5) + 
  scale_color_manual(values=c("blue","red")) + 
#  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
  geom_jitter(data=sig_up,aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +
#  geom_text(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.35,label=pos),angle=45,color="red") +
  geom_hline(yintercept=0.4,linetype="dashed") +
  ylim(-0.5,2) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
  axis.text.x = element_blank(),
  panel.spacing = unit(0, "lines"),
  panel.grid = element_blank(),
  panel.border = element_rect(colour="grey")
  )+
  facet_grid(~chr,scales="free_x",space="free_x") +  geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) 
dev.off()


