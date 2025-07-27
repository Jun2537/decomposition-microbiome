setwd("E:/abs/Git/Figure2")

#Figure2
#select bacterial key decomposer
library(WGCNA)
library(ggplot2)  
library(vegan)
library(psych)

rm(list=ls())
otu<-read.csv('bac_genera_MS_QMP.csv',row.names=1)
otu<-t(otu)
group<-read.csv('Group.csv',row.names=1)
MS_group<-group[group$Bodysite=='MS',]
otu_group<-cbind(MS_group$Day,otu)

cor<- corAndPvalue(otu_group,use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')

r[p > 0.001 | r < 0.60] = 0
r<-data.frame(r)[-1,]
r<-r[r$V1!=0,]
rownames(r) #obtained names of bacterial key decomposer

#select fungal key decomposer
rm(list=ls())
otu<-read.csv('fun_genera_MS_QMP.csv',row.names=1)
otu<-t(otu)
group<-read.csv('Group.csv',row.names=1)
MS_group<-group[group$Bodysite=='MS',]
otu_group<-cbind(MS_group$Day,otu)

cor<- corAndPvalue(otu_group,use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')

r[p > 0.001 | r < 0.60] = 0
r<-data.frame(r)[-1,]
r<-r[r$V1!=0,]
rownames(r) #obtained names of fungal key decomposer

#Draw a heatmap
library(pheatmap)
ann_col=data.frame(Sample=factor(c(rep("Day1",6),rep("Day3",6),rep("Day7",6),rep("Day14",6),rep("Day21",6),rep("Day28",6),rep("Day35",6))))
row.names(ann_col)=colnames(data)

#bacteria
data<-read.csv('bac_genera_key.csv',row.names=1)
ann_color = list(Sample=c(Day1="#c0a6ff",Day3="#ff9289",Day7="#ff81f2",Day14="#70cf07",Day21="#00d1ff",Day28="#00dcaf",Day35="#e0b400"))
pheatmap(data, scale = "row",          
         cluster_rows = T,cluster_cols =F, 
         border_color = "grey60", 
         fontsize_number = 6,
         number_color = "grey30",
         fontsize =10,
         show_rownames = T,
         show_colnames = F,
         color = colorRampPalette(c("#4387b5","white","#e64b35"))(100),
         annotation_col = ann_col,
         annotation = NA, annotation_colors = ann_color,
         annotation_legend = TRUE,
         annotation_names_row = TRUE, annotation_names_col = TRUE) 


#fungi
data<-read.csv('fun_genera_key.csv',row.names=1)
ann_color = list(Sample=c(Day1="#c0a6ff",Day3="#ff9289",Day7="#ff81f2",Day14="#70cf07",Day21="#00d1ff",Day28="#00dcaf",Day35="#e0b400"))
pheatmap(data, scale = "row",          
         cluster_rows = T,cluster_cols =F, 
         border_color = "grey60", 
         fontsize_number = 6,
         number_color = "grey30",
         fontsize =10,
         show_rownames = T,
         show_colnames = F,
         color = colorRampPalette(c("#4387b5","white","#e64b35"))(100),
         annotation_col = ann_col,
         annotation = NA, annotation_colors = ann_color,
         annotation_legend = TRUE,
         annotation_names_row = TRUE, annotation_names_col = TRUE) 