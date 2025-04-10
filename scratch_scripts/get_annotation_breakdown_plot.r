options(scipen=200)
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("ggplot2"))
#suppressPackageStartupMessages(library("ggpubr"))
#suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("reshape2"))
#suppressPackageStartupMessages(library("circlize"))
#suppressPackageStartupMessages(library("ggforce"))
suppressPackageStartupMessages(library("ggrepel"))
library(repr)
options(repr.plot.width=5, repr.plot.height=4)


cl = "D12MSN"

# * load peak annotation
up_male <- read.table(paste("~/projects/combined/diff_peaks_iter/",cl,":Male_up_.peaks.anno.txt", sep = ""), sep="\t", quote="", fill=TRUE, header=F)
down_male <- read.table(paste("~/projects/combined/diff_peaks_iter/",cl,":Male_down_.peaks.anno.txt", sep = ""), sep="\t", quote="", fill=TRUE, header=F)
all_male <- read.table(paste("~/projects/combined/diff_peaks_iter/",cl,":Male_peaks_selected_.peaks.anno.txt", sep = ""), sep="\t", quote="", fill=TRUE, header=F)

up_female <- read.table(paste("~/projects/combined/diff_peaks_iter/",cl,":Female_up_.peaks.anno.txt", sep = ""), sep="\t", quote="", fill=TRUE, header=F)
down_female <- read.table(paste("~/projects/combined/diff_peaks_iter/",cl,":Female_down_.peaks.anno.txt", sep = ""), sep="\t", quote="", fill=TRUE, header=F)
all_female <- read.table(paste("~/projects/combined/diff_peaks_iter/",cl,":Female_peaks_selected_.peaks.anno.txt", sep = ""), sep="\t", quote="", fill=TRUE, header=F)

up_male$dir ="Up"
down_male$dir = "Down"
all_male$dir = "All"

up_male$sex ="Male"
down_male$sex = "Male"
all_male$sex = "Male"

up_female$dir ="Up"
down_female$dir = "Down"
all_female$dir = "All"

up_female$sex ="Female"
down_female$sex = "Female"
all_female$sex = "Female"

up_male$V2 = as.numeric(up_male$V2)
down_male$V2 = as.numeric(down_male$V2)
all_male$V2 = as.numeric(all_male$V2)

up_female$V2 = as.numeric(up_female$V2)
down_female$V2 = as.numeric(down_female$V2)
all_female$V2 = as.numeric(all_female$V2)

peak.anno = rbind(up_male, down_male, all_male,up_female,down_female, all_female )

peak.anno$V1 = trimws(paste(peak.anno$V1))
peak.anno$Repeat = "Not repeat"
peak.anno$Repeat[which(peak.anno$V1%in%c("LINE", "SINE", "LTR", "Other Repeats"))] = "Repeat"

ggplot(peak.anno, aes(fill=Repeat, y=V2, x=dir)) + 
    geom_bar(position="fill", stat="identity") + scale_fill_brewer(palette = "Set2") + facet_wrap(~sex)

ggplot(peak.anno, aes(fill=V1, y=V2, x=dir)) + 
    geom_bar(position="fill", stat="identity") + scale_fill_brewer(palette = "Set3") + facet_wrap(~sex)



#annotatePeaks.pl ITRSPL23GL2:Female_up.bed mm10 -annStats ITRSPL23GL2:Female_up.stats -p 10
for i in *homer_annotation.txt; do 

cut -f 9 $i | sed 's/ (.*)//; s/.*|\([^|]*\)|.*/\1/g; s/-[0-9][0-9]*//; s/\..*//; s/?$//; s/'\ /'-/; s/CpG/CpG island/; s/exon/Exon/; s/intron/Intron/; s/non-coding/Non-coding/; s/promoter-TSS/Promoter-TSS/' > peak.anno.clean.txt

cat peak.anno.clean.txt | sed 's/^Simple_repeat$/Other Repeats/; s/^Low_complexity$/Other Repeats/; s/^DNA$/Other Repeats/; s/^snRNA$/Other Repeats/; s/^Other$/Other Repeats/; s/^Satellite$/Other Repeats/; s/^srpRNA$/Other Repeats/; s/^tRNA$/Other Repeats/; s/^RNA$/Other Repeats/; s/^scRNA$/Other Repeats/; s/^RC$/Other Repeats/; s/^rRNA$/Other Repeats/; s/^Retroposon$/Other Repeats/; s/^Unknown$/Other Repeats/' > peak.anno.merge.txt

cat peak.anno.merge.txt | sort | uniq -c | grep -v Annotation | awk 'BEGIN{OFS="\t"}{print $2" "$3, $1}' > ${i%homer_annotation.txt}.peaks.anno.txt
done
