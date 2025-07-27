setwd('E:/abs/Git/Figure4')

#normalizing the data of metabolites
rm(list=ls())
library(MetaboAnalystR)

mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "metabolites.csv", "rowu", "disc")
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "SumNorm", "NULL", "ParetoNorm", ratio=FALSE)
write.csv(mSet$dataSet$norm,'metabolites_normal.csv',row.names=T)

#Figure4A_PLS-DA
library(ropls)
data<-read.csv('metabolites_normal.csv',row.names=1)
group<-read.csv('metabolites.csv',row.names=1)

df1_plsda <- opls(data, group$Group, orthoI = 0)
df <- as.data.frame(df1_plsda@scoreMN)

df$Group = group$Group
df$samples = rownames(df)

x_lab <- df1_plsda@modelDF[1, "R2X"] * 100
y_lab <- df1_plsda@modelDF[2, "R2X"] * 100

library(ggplot2)
library(RColorBrewer)
mycolor<- brewer.pal(9,"Paired")
p <- ggplot(df, aes(p1, p2,colour =Group))+
  geom_point(size = 4,alpha=0.7) +
scale_fill_manual(values = mycolor) +
scale_colour_manual(values = mycolor) +
  xlab(paste("Comp1 (", x_lab, " %)", sep = "")) +
  ylab(paste("PCo2 (", y_lab, " %)", sep = "")) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw()

#Figure4B_boxplot
library(RColorBrewer)
library(ggsci)
library(ggpubr)

p1_group<-cbind(df$p1,group)
colnames(p1_group)[1] <- "p1"

mycolor<- brewer.pal(9,"Paired")
ggplot(p1_group, aes(x = Group, y = p1, color=Group)) + 
stat_boxplot(geom="errorbar")+
geom_jitter(size = 1.5,alpha = 0.4,width = 0.2)+ 
scale_colour_manual(values = mycolor) +
geom_boxplot() +
labs(y = "Comp1", x = "Postmoterm interval (days)")+
theme_bw()+
coord_flip() 


#Figure4C
#Sort the metabolites according to their VIP values
data_VIP <- df1_plsda@vipVn
data_VIP_select <- data_VIP[data_VIP > 1.5]
name<-read.csv('metabolites_name.csv',row.names=1)
select_name<-name[rownames(data.frame(data_VIP_select)),]
df<-cbind(data.frame(data_VIP_select),data.frame(select_name))

df = df[order(df[,1],decreasing =F),]
df$select_name=factor(df$select_name,levels =df$select_name,ordered=TRUE)

ggplot(data = df, aes(x=select_name,y=data_VIP_select)) + 
  geom_bar(stat="identity",fill='steelblue4')+coord_flip()+ theme_bw()+theme(panel.grid=element_line(colour=NA))+
  scale_y_continuous(name="Variable importance in projection") 


#Figure_4D

library(WGCNA)
library(ggplot2)  
library(vegan)
library(psych)

#select mebolites which significant related to PMI
data<-read.csv('metabolites_normal.csv',row.names=1)
group<-read.csv('Group_MS.csv',row.names=1)
data_group<-cbind(group,data)

cor<- corAndPvalue(data_group,use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')

r[p > 0.001 |abs(r) < 0.60] = 0
r<-data.frame(r)[-1,]
r<-r[r$Day!=0,]
rownames(r) 
meta_select<-data[,rownames(r)]


#Draw a heatmap
library(pheatmap)
ann_col=data.frame(Sample=factor(c(rep("Day1",6),rep("Day3",6),rep("Day7",6),rep("Day14",6),rep("Day21",6),rep("Day28",6),rep("Day35",6))))
row.names(ann_col)=rownames(meta_select)

ann_color = list(Sample=c(Day1="#c0a6ff",Day3="#ff9289",Day7="#ff81f2",Day14="#70cf07",Day21="#00d1ff",Day28="#00dcaf",Day35="#e0b400"))
pheatmap(t(meta_select), scale = "row",          
         cluster_rows = T,cluster_cols =F, 
         cutree_rows=2,
         border_color = "grey60", 
         fontsize_number = 6,
         fontsize =10,
         show_rownames = F,
         show_colnames = F,
         color = colorRampPalette(c("#4387b5","white","#e64b35"))(100),
         annotation_col = ann_col,
         annotation = NA, annotation_colors = ann_color,
         annotation_legend = TRUE,
         annotation_names_row = TRUE, annotation_names_col = TRUE) 

#Figure4E
library(RColorBrewer)
library(ggplot2)

data<-read.csv('stack.csv',row.names=1)
mycolor<- brewer.pal(8,"Paired")

ggplot(data, aes( x = Cluster,y=100 * Abundance,fill = Superclass))+
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  labs(x="",y="Percentage (%)")+
  scale_fill_manual(values=mycolor)+
  theme_bw() 

