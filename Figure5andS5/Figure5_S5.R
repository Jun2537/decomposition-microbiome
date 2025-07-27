setwd('E:/abs/Git/Figure5andS5')
library(vegan)


#Figure5A_bacteria
metab<-read.csv('metabolites_normal.csv',row.names=1)
metab_pca<-rda(metab,scale=F)

otu<-read.csv('bac_QMP.csv',row.names=1)
otu<-t(otu)
otu_hel<-decostand(otu,method='hellinger')
otu_pca<-rda(otu_hel,scale=F)

biplot(metab_pca,choices=c(1,2),scaling=1,main='metab PCA',col=c('red','blue'))
biplot(otu_pca,choices=c(1,2),scaling=1,main='otu PCA',col=c('red','blue'))

site_metab<-summary(metab_pca,scaling=1)$site
site_otu<-summary(otu_pca,scaling=1)$site

proc<-procrustes(X=otu_pca,Y=metab_pca,symmetric=TRUE)


names(proc)
head(proc$Yrot)#Procrustes y_lab
head(proc$X)#Procrustes x_lab
proc$ss #M2
proc$rotation 

#Fig.S5_bacteria
plot(proc,kind=2)
residuals(proc)

set.seed(123)
prot<-protest(X=metab_pca,Y=otu_pca,permutations=how(nperm=999))
prot

names(prot)
prot$signif
prot$ss

library(ggplot2)
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)

group <- read.csv('Group.csv',row.names=1)
Y <- cbind(group,Y)
library(RColorBrewer)
mycolor<- brewer.pal(9,"Paired")

#Fig.5A_Bacteria
p <- ggplot(Y) +
geom_point(aes(X1, X2, color = Group), size = 4, shape = 16) +
geom_point(aes(PC1, PC2, color = Group), size = 4, shape = 17) +
scale_color_manual(values =mycolor) +
geom_segment(aes(x = X1, y = X2, xend = PC1, yend = PC2,color=Group), arrow = arrow(length = unit(0.1, 'cm')),
    , size = 0.3) +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
    legend.key = element_rect(fill = 'transparent')) +
labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) 

#Figure5A_fungi
metab<-read.csv('metabolites_normal.csv',row.names=1)
metab_pca<-rda(metab,scale=F)

otu<-read.csv('fun_QMP.csv',row.names=1)
otu<-t(otu)
otu_hel<-decostand(otu,method='hellinger')
otu_pca<-rda(otu_hel,scale=F)


biplot(metab_pca,choices=c(1,2),scaling=1,main='metab PCA',col=c('red','blue'))
biplot(otu_pca,choices=c(1,2),scaling=1,main='otu PCA',col=c('red','blue'))

site_metab<-summary(metab_pca,scaling=1)$site
site_otu<-summary(otu_pca,scaling=1)$site

proc<-procrustes(X=otu_pca,Y=metab_pca,symmetric=TRUE)


names(proc)
head(proc$Yrot)#Procrustes y_lab
head(proc$X)#Procrustes x_lab
proc$ss #M2
proc$rotation 

#Fig.S5_fungi
plot(proc,kind=2)
residuals(proc)

set.seed(123)
prot<-protest(X=metab_pca,Y=otu_pca,permutations=how(nperm=999))
prot

names(prot)
prot$signif
prot$ss

library(ggplot2)
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)

group <- read.csv('Group.csv',row.names=1)
Y <- cbind(group,Y)
library(RColorBrewer)
mycolor<- brewer.pal(9,"Paired")

#Fig.5A_fungi
p <- ggplot(Y) +
geom_point(aes(X1, X2, color = Group), size = 4, shape = 16) +
geom_point(aes(PC1, PC2, color = Group), size = 4, shape = 17) +
scale_color_manual(values =mycolor) +
geom_segment(aes(x = X1, y = X2, xend = PC1, yend = PC2,color=Group), arrow = arrow(length = unit(0.1, 'cm')),
    , size = 0.3) +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
    legend.key = element_rect(fill = 'transparent')) +
labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) 

#Figure5B for mantel test
bac<-read.csv('bac_QMP.csv',row.names=1)
bac<-t(bac)
fun<-read.csv('fun_QMP.csv',row.names=1)
fun<-t(fun)
bac_dist<-vegdist(bac,'bray')
fun_dist<-vegdist(fun,'bray')

metab<-read.csv('metabolites_normal.csv',row.names=1)
metab_dist<-dist(metab)

#mantel test
mantel(metab_dist,bac_dist,method='spearman',permutations=999)
mantel(metab_dist,fun_dist,method='spearman',permutations=999)

#partial mantel test
mantel.partial(metab_dist,bac_dist,fun_dist,method='spearman',permutations=999)
mantel.partial(metab_dist,fun_dist,bac_dist,method='spearman',permutations=999)

summary(lm(as.numeric(metab_dist)~as.numeric(bac_dist)))
summary(lm(as.numeric(metab_dist)~as.numeric(fun_dist)))

df1<-data.frame(cbind(as.numeric(metab_dist),as.numeric(bac_dist)))
df2<-data.frame(cbind(as.numeric(metab_dist),as.numeric(fun_dist)))

#Figure5B_bacteria
ggplot(df1, aes(x=X2, y=X1)) +
    geom_point(size = 1.5,alpha=0.2) +    
    geom_smooth(method=lm,color='#E64B35')+
    xlab(paste("Bray_curtis distance of bacterial community", sep = "")) +
    ylab(paste("Euclidean distance of metabolites", sep = "")) +
    theme_bw()

#Figure5B_fungi
ggplot(df2, aes(x=X2, y=X1)) +
    geom_point(size = 1.5,alpha=0.2) +     
    geom_smooth(method=lm,color='#3c5488')+
    scale_x_continuous(limits = c(0.37, 1), breaks = seq(0.4, 1, 0.2))+
    xlab(paste("Bray_curtis distance of fungal community", sep = "")) +
    ylab(paste("Euclidean distance of metabolites", sep = "")) +
    theme_bw()


#Figure5C
#Redundancy analysis
library(vegan)
#select microbial factors significantly correlated to metabolites by RDA
metabo<-read.csv('metabolites_normal.csv',row.names=1)
micro<-read.csv('micro1.csv',row.names=1)
micro<-data.frame(t(micro))

rda<-rda(metabo~.,micro,scale=F)
result=summary(rda)
envsig<-envfit(rda,micro)
envsig

#reperform RDA using selected microbial factors
metabo<-read.csv('metabolites_normal.csv',row.names=1)
micro<-read.csv('micro_select.csv',row.names=1)
micro<-data.frame(t(micro))

rda<-rda(metabo~.,micro,scale=F)
result=summary(rda)

eig_percent<-round(result$eig/sum(result$eig)*100,3)
rda1<-round(rda$CCA$eig[1]/sum(rda$CCA$eig)*100,2)
rda2<-round(rda$CCA$eig[2]/sum(rda$CCA$eig)*100,2)

spe=as.data.frame(result$species[,1:2])*4
head(spe)

sit=as.data.frame(result$sites[,1:2])
head(sit)

env2=as.data.frame(result$biplot[,1:2])*5
head(env2)

envsig<-envfit(rda,micro)
envsig

grp=as.data.frame(c(rep("MS01",6),rep("MS03",6),rep("MS07",6),rep("MS14",6),rep("MS21",6),rep("MS28",6),
rep("MS35",6)))

colnames(grp)="group"

library(RColorBrewer)
library(ggplot2)
library(ggrepel)
mycolor<-brewer.pal(9,"Paired")

p1<-ggplot()+
geom_point(data=sit,aes(RDA1,RDA2,color=grp$group),size=4,alpha=0.6)+
scale_colour_manual(values=mycolor)+
geom_segment(data=env2[1:5,],aes(x=0,y=0,xend=RDA1,yend=RDA2),
arrow=arrow(angle=22.5,length=unit(0.35,"cm"),
type="closed"),linetype=1,size=0.6,colour="#E64B35")+
geom_segment(data=env2[6:10,],aes(x=0,y=0,xend=RDA1,yend=RDA2),
arrow=arrow(angle=22.5,length=unit(0.35,"cm"),
type="closed"),linetype=1,size=0.6,colour="#3C5488")+
geom_text_repel(data=env2,aes(RDA1,RDA2,label=row.names(env2)))+
labs(x="RDA1 (66.22%)",y="RDA2 (12.52%)")+
geom_hline(yintercept=0,linetype=3,size=1)+
geom_vline(xintercept=0,linetype=3,size=1)+
guides(shape=guide_legend(title=NULL,color="black"),
fill=guide_legend(title=NULL))+
theme_bw()+theme(panel.grid=element_blank())


#Figure5D for randomforest indicator
setwd('E:/abs/Git/Figure5/Figure5D')
library(vegan)

rm(list=ls())

#select metabolites significantly correlated to bacterial community 
data<-read.csv('bac_QMP.csv',row.names=1)
data<-t(data)
data<-decostand(data,'hellinger')

pca_result<-prcomp(data,center=T, scale = F)
pc<-pca_result$x

metab<-read.csv('class_normal.csv',row.names=1)

library(randomForest)
set.seed(1)
imp = replicate(100,importance(randomForest(pc[,1]~.,data=metab, ntree=500, proximity=TRUE, importance=TRUE))[,1],simplify=F)

merged_result<-imp[[1]]
for (i in 2:100){
current_list<-imp[[i]]
merged_result<-cbind(merged_result,current_list)
}

VIP_mean<-data.frame(rowMeans(merged_result))
colnames(VIP_mean)[1]<-'VIP'
VIP_mean<-cbind(rownames(VIP_mean),VIP_mean)
bac_VIP_mean<-VIP_mean[order(-VIP_mean$VIP),]
write.csv(bac_VIP_mean[1:20,],'bac_VIP_top20.csv',row.names=T)


#select metabolites significantly correlated to fungal community 
data<-read.csv('fun_QMP.csv',row.names=1)
data<-t(data)
data<-decostand(data,'hellinger')

pca_result<-prcomp(data, center=T,scale = F)
pc<-pca_result$x

metab<-read.csv('class_normal.csv',row.names=1)

library(randomForest)
set.seed(1)
imp = replicate(100,importance(randomForest(pc[,1]~.,data=metab, ntree=500, proximity=TRUE, importance=TRUE))[,1],simplify=F)

merged_result<-imp[[1]]
for (i in 2:100){
current_list<-imp[[i]]
merged_result<-cbind(merged_result,current_list)
}

VIP_mean<-data.frame(rowMeans(merged_result))
colnames(VIP_mean)[1]<-'VIP'
VIP_mean<-cbind(rownames(VIP_mean),VIP_mean)
fun_VIP_mean<-VIP_mean[order(-VIP_mean$VIP),]
write.csv(fun_VIP_mean[1:20,],'fun_VIP_top20.csv',row.names=T)

#select the intersection of metabolites most relevant (top 20) to bacteria and fungi
inter<-intersect(x=bac_VIP_mean[1:20,1],y=fun_VIP_mean[1:20,1])
write.csv(inter,'inter.csv',row.names=T)


#Manually convert 'bac_VIP_top20.csv', 'fun_VIP_top20.csv','inter.csv' to facilitate drawing the graph

library(ggplot2)
library(tidyr)
library(dplyr)

hd1<-read.csv('bac_VIP_top20_trans.csv',row.names=1)
hd1$ID<-factor(hd1$ID,levels=rev(unique(hd1$ID)),ordered=T)

p1<-ggplot(hd1,aes(-IncMSE,ID))+
geom_col(fill='#7EB4C6',width=0.7)+
labs(x="%IncMSE",y="")+
scale_x_continuous(position="bottom",
expand=expansion(add=c(0.2,0.2)),
limits=c(-10,0),
breaks=c(-10,-5,0),
label=c("10","5","0"))+
scale_y_discrete(position="right")+
theme(axis.text=element_text(size=13))+
theme(panel.grid.minor.x=element_blank(),
plot.title=element_text(hjust=0.5),
axis.title.x=element_text(size=13),
panel.background=element_blank(),
legend.position="none")


hd2<-read.csv('fun_VIP_top20_trans.csv',row.names=1)
hd2$ID<-factor(hd2$ID,levels=rev(unique(hd2$ID)),ordered=T)

p2<-ggplot(hd2,aes(IncMSE,ID))+
geom_col(fill='#AD93B8',width=0.7)+
labs(x="%IncMSE",y="")+
scale_x_continuous(position="bottom",
expand=expansion(add=c(0.2,0.2)),
limits=c(0,11),
breaks=c(0,5,10),
label=c("0","5","10"))+
scale_y_discrete(position="left")+
theme(axis.text=element_text(size=13))+
theme(panel.grid.minor.x=element_blank(),
plot.title=element_text(hjust=0.5),
axis.title.x=element_text(size=13),
panel.background=element_blank(),
legend.position="none")

net<-read.csv('inter_trans.csv',row.names=1)

h1<-hd1$ID
sou<-net$Source
y1<-NULL
for(i in sou){
ys<-which(h1==i)
y1<-c(y1,ys)
}
net_y1<-length(h1)+1-y1

h2<-hd2$ID
tar<-net$Target
y2<-NULL
for(i in tar){
yt<-which(h2==i)
y2<-c(y2,yt)
}
net_y2<-length(h2)+1-y2


x1<-rep(2,length(tar))
x2<-rep(3,length(sou))

df<-data.frame(sou,x1,net_y1,tar,x2,net_y2)
df

n=length(hd1$ID)

p_mid<-ggplot(df)+geom_segment(aes(x1,net_y1,xend=x2,yend=net_y2),
size=0.5,color="red")+
geom_point(aes(x=x1,y=net_y1),size=3.5,
fill="#7EB4C6",
color="#7EB4C6",
stroke=1,
shape=21)+
geom_point(aes(x=x2,y=net_y2),size=3.5,
fill="#AD93B8",
color="#AD93B8",
stroke=1,
shape=21)+
scale_y_continuous(limits=c(1,n),expand=expansion(add=c(0.5,0.7)))+
scale_x_continuous(expand=expansion(0,0.1))+
theme_void()
p_mid

library(patchwork)
p1+p_mid+p2




