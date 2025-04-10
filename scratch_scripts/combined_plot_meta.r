library(ggplot2)
setwd("/home/lamaral/ps-renlab2/projects/combined_all/")

#meta = read.csv("before_filter_clusters/meta.csv")

#meta = read.csv("combined_meta.csv")

meta = read.csv("final_meta.csv")

head(meta$replicate)

meta = meta[-which(meta$doub == "TRUE"),]

head(meta$batch)

table(meta[which(meta$batch == "Female"), "keep"])

nrow(meta)

meta$age = ""
meta$age[which(meta$batch == "Female")] = sapply(strsplit(meta$sample[which(meta$batch == "Female")],'_'),'[[',2) 
meta$age[which(meta$batch == "Male")] = sapply(strsplit(meta$sample[which(meta$batch == "Male")],'_'),'[[',1) 
meta$age = gsub("Male:","",meta$age)
meta$age[which(meta$age == "8wk")] = "2mo" 
meta$region = ""
meta$region[which(meta$batch == "Female")] = sapply(strsplit(meta$sample[which(meta$batch == "Female")],'_'),'[[',1) 
meta$region = gsub("Female:", "", meta$region)
meta$region[which(meta$batch == "Male")] = sapply(strsplit(meta$sample[which(meta$batch == "Male")],'_'),'[[',2) 
meta$region[which(meta$region %in% c("11E", "11E-11F-12E", "11F", "12E"))] = "HCP"
meta$region[which(meta$region %in% c("12D", "12D-13B", "13B"))] = "ENT"
meta$region[which(meta$region %in% c("13D", "13D-14C", "14C"))] = "RLP"
meta$region[which(meta$region %in% c("2A", "2A-3A", "3A"))] = "FC"
meta$region[which(meta$region %in% c("3F", "3F-4E", "4E"))] = "NAC"
meta$region[which(meta$region %in% c("7H", "8H", "7H-8H-9G", "9G"))] = "AMY"
meta$region[which(meta$region %in% c("5E", "5E-6E", "6E"))] = "CP"
meta$region[which(meta$region %in% c("8E", "9H","8E-9H-8J-9J", "8J", "9J"))] = "HCA"
meta$age = factor(meta$age, levels = c("2mo", "9mo", "18mo"))
meta$sample_name = sapply(strsplit(meta$sample,':'),'[[',2) 
asamp = unique(meta$sample_name)
asamp = asamp[c(grep("_1$", asamp),grep("_rep1$", asamp),grep("_2$", asamp), grep("_rep2$", asamp),grep("_3$", asamp))]
asamp = asamp[c(grep("8wk", asamp),grep("2mo", asamp), grep("9mo", asamp),grep("18mo", asamp))]
meta$sample_name = factor(meta$sample_name, levels = asamp)

unique_samples <- meta$sample
last_characters <- substr(unique_samples, nchar(unique_samples), nchar(unique_samples))
unique(last_characters)
meta$replicate = last_characters

head(meta)

write.csv(meta, file = "meta.csv", row.names=F)

meta$age = factor(meta$age, levels = c("2mo", "9mo", "18mo"))

meta$diss_rep = gsub("18mo_", "", meta$sample_name)
meta$diss_rep = gsub("09mo_", "", meta$diss_rep)
meta$diss_rep = gsub("02mo_", "", meta$diss_rep)


library(repr)
options(repr.plot.width=12, repr.plot.height=6)

ggplot(meta, aes(fill = age , x = sample_name
                )) +
  geom_bar(position = "stack", stat = "count") + theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~batch+region, nrow = 2,scales = "free_x")




table(meta$keep)

ggplot(meta, aes(fill = keep , x = sample_name)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~batch+region, nrow = 2,scales = "free_x") + geom_text(stat='count', aes(label=..count..), vjust=-1)


meta$overTsse7 = "No"
meta$overTsse7[which(meta$tsse>=7)] = "Yes"
ggplot(meta, aes(fill = overTsse7 , x = sample_name)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~batch+region, nrow = 2,scales = "free_x")

ggplot(meta, aes(fill = batch , x = factor(leiden_3))) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(meta[which(meta$leiden_3==68),], aes(x = sample)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

meta$keep = "yes"
meta[which(meta$leiden_3 %in% c(34,45,54)),"keep"] = "no"

ggplot(meta[which(meta$keep=="yes"),], aes(x = sample)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(meta, aes( x = factor(leiden_3), y = tsse)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + facet_wrap(~batch, nrow = 2)

ggplot(meta, aes(fill = overTsse7 , x = factor(leiden_3))) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(meta, aes(fill = over2500 , x = factor(leiden_2))) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


ggplot(meta, aes(x=sample, y=n_fragment, color = batch)) + 
  geom_boxplot(notch=TRUE)


head(metaf)

new_df

library(dplyr)
#predictions
metaf = meta
metaf = metaf[which(metaf$male_ct!=''),]
predictions <- table(metaf$leiden_3,metaf$male_ct)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

new_df <- predictions %>%
  group_by(Var1) %>%
  filter(Freq == max(Freq)) %>%
  select(Var1, Var2, Freq)

new_df = as.data.frame(new_df)
rownames(new_df) = new_df$Var1
meta$best_celltype = new_df[paste(meta$leiden_3), "Var2"]
meta$best_celltype = as.character(meta$best_celltype)
#metab$best_celltype[which(metab$`leiden_2.5`%in%c(44,69,39))] = "LQ"


head(meta)

#meta$best_celltype[which(meta$`leiden_3`%in%c(41,49,71,72))] = "remove"
#meta$best_celltype[which(meta$`leiden_3`%in%c(62))] = "LQ"
meta$best_celltype[which(meta$`best_celltype`%in%c("OGC1", "OGC2"))] = "OGC"
meta$best_celltype[which(meta$`best_celltype`%in%c("ASCTENT", "ASCTE1"))] = "ASCTE"

#meta$best_celltype[which(meta$`leiden_2`== 38 & meta$`leiden_3`== 50)] = "RGL"

write.csv(meta, "meta_bestcelltype.csv", row.names = F)

ggplot(meta, aes(fill = batch , x = factor(celltype_final))) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + facet_wrap(~age, nrow=3)

ggplot(meta, aes(fill = batch , x = factor(leiden_3))) +
  geom_bar(position = "fill", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

ggplot(meta, aes(fill = region, color = age , x = factor(leiden_1))) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  + facet_wrap(~batch, nrow=2)

ggplot(meta, aes(fill = over2500, color = age , x = factor(leiden_1))) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  + facet_wrap(~batch, nrow=2)

nrow(meta)


ggplot(meta, aes(fill = overTsse7, x = factor(best_celltype))) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  + facet_wrap(~batch, nrow=2)

ggplot(meta, aes(fill = batch , x = factor(best_celltype))) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

length(unique(meta$best_celltype))

write.csv(meta, "meta_bestcelltype.csv")

ggplot(meta, aes(fill = batch, color = age ,y = n_fragment, x = region)) +
  geom_boxplot() 

library(repr)
options(repr.plot.width=16, repr.plot.height=8)

metaf = meta[-which(meta$best_celltype %in% c("remove", "LQ")),] 

ggplot(metaf, aes(fill = age , x = sample_name)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~batch+region, nrow = 2,scales = "free_x") + geom_text_repel(stat='count', aes(label=..count..), vjust=-1)




ggplot(metaf, aes(fill = age , x = sample_name)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~batch+region, nrow = 2,scales = "free_x") 


ggplot(meta, aes(fill = age , x = sample_name)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~batch+region, nrow = 2,scales = "free_x") + geom_text_repel(stat='count', aes(label=..count..), vjust=-1)



meta = read.csv("meta_bestcelltype.csv")

nrow(meta)

asamp = unique(meta$sample_name)
asamp = asamp[c(grep("_1$", asamp),grep("_rep1$", asamp),grep("_2$", asamp), grep("_rep2$", asamp))]
asamp = asamp[c(grep("8wk", asamp),grep("2mo", asamp), grep("9mo", asamp),grep("18mo", asamp))]
meta$sample_name = factor(meta$sample_name, levels = asamp)
meta$age = factor(meta$age, levels = c("2mo", "9mo", "18mo"))

ggplot(meta, aes(fill = age , x = sample_name)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~batch+region, nrow = 2,scales = "free_x") 


metaf = meta[-which(meta$best_celltype %in% c("remove", "LQ")),] 

ggplot(metaf, aes(fill = age , x = sample_name)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~batch+region, nrow = 2,scales = "free_x") + geom_text(stat='count', aes(label=..count..), vjust=-1)



ggplot(metaf, aes(fill = age , x = sample_name)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~batch+region, nrow = 2,scales = "free_x") 


ggplot(metaf, aes(fill = age, y = tsse, x = sample_name)) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~batch+region, nrow = 2,scales = "free_x") 


nrow(metaf)

library(ggrepel)

ggplot(meta, aes(fill = age , x = sample_name)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~batch+region, nrow = 2,scales = "free_x") + geom_text_repel(stat='count', aes(label=..count..), vjust=-1)


ggplot(metaf, aes(fill = region, color = age , x = best_celltype)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  + facet_wrap(~batch, nrow=2)

pdf("age_comp_stack.pdf", height =6 , width = 30)
ggplot(metaf, aes(fill = age ,x = best_celltype)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))   + facet_wrap(~batch+region, nrow=2)
dev.off()

ggplot(metaf, aes(fill = age ,x = best_celltype)) +
  geom_bar(position = "fill", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))   + facet_wrap(~batch+region, nrow=2)



