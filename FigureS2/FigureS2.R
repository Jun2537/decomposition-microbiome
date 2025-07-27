setwd('E:/abs/Git/FigureS2')

#Figure S2
#upset for bacterial communities in gravesoil 
#changes in QMP data of bacterial communities 
data<-read.csv("QMP_bac_genera_GS_top100.csv",row.names=1)
data<-t(data)

group<-read.csv('Group_GS.csv',row.names=1)

data_group<-cbind(group,data)

Fresh<-data_group[data_group$Stage=='Fresh',-1:-2]
Aged<-data_group[data_group$Stage=='Aged',-1:-2]

Result=data.frame()
for(i in 1:100){
  wil <- wilcox.test(Fresh[,i],Aged[,i],alternative='two.sided',correct=T)
  dif<-mean(Aged[,i])-mean(Fresh[,i])
  Result=rbind(Result,
                 cbind(id=colnames(data)[i],
                 dif,
                 pvalue=wil$p.value)
                  )
}

library(dplyr)
Result_QMP<-mutate(Result,rd = cut(as.numeric(Result$dif), breaks = c(-Inf, 0, Inf),
                  labels = c("Decreased", "Increased")),
         pd = cut(as.numeric(Result$pvalue), breaks = c(-Inf, 0.05, Inf),
                  labels = c("< 0.05", "ns")))

#changes in RMP data of bacterial communities 
data<-read.csv("RMP_bac_genera_GS_top100.csv",row.names=1)
data<-t(data)

group<-read.csv('Group_GS.csv',row.names=1)

data_group<-cbind(group,data)

Fresh<-data_group[data_group$Stage=='Fresh',-1:-2]
Aged<-data_group[data_group$Stage=='Aged',-1:-2]

Result=data.frame()
for(i in 1:100){
  wil <- wilcox.test(Fresh[,i],Aged[,i],alternative='two.sided',correct=T)
  dif<-mean(Aged[,i])-mean(Fresh[,i])
  Result=rbind(Result,
                 cbind(id=colnames(data)[i],
                 dif,
                 pvalue=wil$p.value)
                  )
}

library(dplyr)
Result_RMP<-mutate(Result,rd = cut(as.numeric(Result$dif), breaks = c(-Inf, 0, Inf),
                  labels = c("Decreased", "Increased")),
         pd = cut(as.numeric(Result$pvalue), breaks = c(-Inf, 0.05, Inf),
                  labels = c("< 0.05", "ns")))


upset<-cbind(data.frame(Result_QMP[,c(1,4,5)]),data.frame(Result_RMP[,4:5]))
colnames(upset)<-c('ID','QMP_rd','QMP_pd','RMP_rd','RMP_pd')

QMP_increased<-upset[upset$QMP_rd=='Increased'&upset$QMP_pd=='< 0.05',1]
QMP_decreased<-upset[upset$QMP_rd=='Decreased'&upset$QMP_pd=='< 0.05',1]
QMP_ns<-upset[upset$QMP_pd=='ns',1]
RMP_increased<-upset[upset$RMP_rd=='Increased'&upset$RMP_pd=='< 0.05',1]
RMP_decreased<-upset[upset$RMP_rd=='Decreased'&upset$RMP_pd=='< 0.05',1]
RMP_ns<-upset[upset$RMP_pd=='ns',1]

input<-list(QMP_increased=QMP_increased,QMP_decreased=QMP_decreased,QMP_ns=QMP_ns,RMP_increased=RMP_increased,RMP_decreased=RMP_decreased,RMP_ns=RMP_ns)

library(UpSetR)
upset(fromList(input),
sets = c("RMP_ns", "QMP_ns", "RMP_increased", "QMP_increased",  "RMP_decreased","QMP_decreased"),
nsets=6,mb.ratio = c(0.6, 0.4),
keep.order = TRUE)







#upset for bacterial communities in Tissue 
#changes in QMP data of bacterial communities 
data<-read.csv("QMP_bac_genera_MS_top100.csv",row.names=1)
data<-t(data)

group<-read.csv('Group_MS.csv',row.names=1)

data_group<-cbind(group,data)

Fresh<-data_group[data_group$Stage=='Fresh',-1:-2]
Aged<-data_group[data_group$Stage=='Aged',-1:-2]

Result=data.frame()
for(i in 1:100){
  wil <- wilcox.test(Fresh[,i],Aged[,i],alternative='two.sided',correct=T)
  dif<-mean(Aged[,i])-mean(Fresh[,i])
  Result=rbind(Result,
                 cbind(id=colnames(data)[i],
                 dif,
                 pvalue=wil$p.value)
                  )
}

library(dplyr)
Result_QMP<-mutate(Result,rd = cut(as.numeric(Result$dif), breaks = c(-Inf, 0, Inf),
                  labels = c("Decreased", "Increased")),
         pd = cut(as.numeric(Result$pvalue), breaks = c(-Inf, 0.05, Inf),
                  labels = c("< 0.05", "ns")))

#changes in RMP data of bacterial communities 
data<-read.csv("RMP_bac_genera_MS_top100.csv",row.names=1)
data<-t(data)

group<-read.csv('Group_MS.csv',row.names=1)

data_group<-cbind(group,data)

Fresh<-data_group[data_group$Stage=='Fresh',-1:-2]
Aged<-data_group[data_group$Stage=='Aged',-1:-2]

Result=data.frame()
for(i in 1:100){
  wil <- wilcox.test(Fresh[,i],Aged[,i],alternative='two.sided',correct=T)
  dif<-mean(Aged[,i])-mean(Fresh[,i])
  Result=rbind(Result,
                 cbind(id=colnames(data)[i],
                 dif,
                 pvalue=wil$p.value)
                  )
}

library(dplyr)
Result_RMP<-mutate(Result,rd = cut(as.numeric(Result$dif), breaks = c(-Inf, 0, Inf),
                  labels = c("Decreased", "Increased")),
         pd = cut(as.numeric(Result$pvalue), breaks = c(-Inf, 0.05, Inf),
                  labels = c("< 0.05", "ns")))


upset<-cbind(data.frame(Result_QMP[,c(1,4,5)]),data.frame(Result_RMP[,4:5]))
colnames(upset)<-c('ID','QMP_rd','QMP_pd','RMP_rd','RMP_pd')

QMP_increased<-upset[upset$QMP_rd=='Increased'&upset$QMP_pd=='< 0.05',1]
QMP_decreased<-upset[upset$QMP_rd=='Decreased'&upset$QMP_pd=='< 0.05',1]
QMP_ns<-upset[upset$QMP_pd=='ns',1]
RMP_increased<-upset[upset$RMP_rd=='Increased'&upset$RMP_pd=='< 0.05',1]
RMP_decreased<-upset[upset$RMP_rd=='Decreased'&upset$RMP_pd=='< 0.05',1]
RMP_ns<-upset[upset$RMP_pd=='ns',1]

input<-list(QMP_increased=QMP_increased,QMP_decreased=QMP_decreased,QMP_ns=QMP_ns,RMP_increased=RMP_increased,RMP_decreased=RMP_decreased,RMP_ns=RMP_ns)

library(UpSetR)
upset(fromList(input),
sets = c("RMP_ns", "QMP_ns", "RMP_increased", "QMP_increased",  "RMP_decreased","QMP_decreased"),
nsets=6,mb.ratio = c(0.6, 0.4),
keep.order = TRUE)






#upset for fungal communities in gravesoil 
#changes in QMP data of fungal communities 
data<-read.csv("QMP_fungi_genera_GS_top100.csv",row.names=1)
data<-t(data)

group<-read.csv('Group_GS.csv',row.names=1)

data_group<-cbind(group,data)

Fresh<-data_group[data_group$Stage=='Fresh',-1:-2]
Aged<-data_group[data_group$Stage=='Aged',-1:-2]

Result=data.frame()
for(i in 1:100){
  wil <- wilcox.test(Fresh[,i],Aged[,i],alternative='two.sided',correct=T)
  dif<-mean(Aged[,i])-mean(Fresh[,i])
  Result=rbind(Result,
                 cbind(id=colnames(data)[i],
                 dif,
                 pvalue=wil$p.value)
                  )
}

library(dplyr)
Result_QMP<-mutate(Result,rd = cut(as.numeric(Result$dif), breaks = c(-Inf, 0, Inf),
                  labels = c("Decreased", "Increased")),
         pd = cut(as.numeric(Result$pvalue), breaks = c(-Inf, 0.05, Inf),
                  labels = c("< 0.05", "ns")))

#changes in RMP data of fungal communities 
data<-read.csv("RMP_fungi_genera_GS_top100.csv",row.names=1)
data<-t(data)

group<-read.csv('Group_GS.csv',row.names=1)

data_group<-cbind(group,data)

Fresh<-data_group[data_group$Stage=='Fresh',-1:-2]
Aged<-data_group[data_group$Stage=='Aged',-1:-2]

Result=data.frame()
for(i in 1:100){
  wil <- wilcox.test(Fresh[,i],Aged[,i],alternative='two.sided',correct=T)
  dif<-mean(Aged[,i])-mean(Fresh[,i])
  Result=rbind(Result,
                 cbind(id=colnames(data)[i],
                 dif,
                 pvalue=wil$p.value)
                  )
}

library(dplyr)
Result_RMP<-mutate(Result,rd = cut(as.numeric(Result$dif), breaks = c(-Inf, 0, Inf),
                  labels = c("Decreased", "Increased")),
         pd = cut(as.numeric(Result$pvalue), breaks = c(-Inf, 0.05, Inf),
                  labels = c("< 0.05", "ns")))


upset<-cbind(data.frame(Result_QMP[,c(1,4,5)]),data.frame(Result_RMP[,4:5]))
colnames(upset)<-c('ID','QMP_rd','QMP_pd','RMP_rd','RMP_pd')

QMP_increased<-upset[upset$QMP_rd=='Increased'&upset$QMP_pd=='< 0.05',1]
QMP_decreased<-upset[upset$QMP_rd=='Decreased'&upset$QMP_pd=='< 0.05',1]
QMP_ns<-upset[upset$QMP_pd=='ns',1]
RMP_increased<-upset[upset$RMP_rd=='Increased'&upset$RMP_pd=='< 0.05',1]
RMP_decreased<-upset[upset$RMP_rd=='Decreased'&upset$RMP_pd=='< 0.05',1]
RMP_ns<-upset[upset$RMP_pd=='ns',1]

input<-list(QMP_increased=QMP_increased,QMP_decreased=QMP_decreased,QMP_ns=QMP_ns,RMP_increased=RMP_increased,RMP_decreased=RMP_decreased,RMP_ns=RMP_ns)

library(UpSetR)
upset(fromList(input),
sets = c("RMP_ns", "QMP_ns", "RMP_increased", "QMP_increased",  "RMP_decreased","QMP_decreased"),
nsets=6,mb.ratio = c(0.6, 0.4),
keep.order = TRUE)







#upset for fungal communities in Tissue 
#changes in QMP data of fungal communities 
data<-read.csv("QMP_fungi_genera_MS_top100.csv",row.names=1)
data<-t(data)

group<-read.csv('Group_MS.csv',row.names=1)

data_group<-cbind(group,data)

Fresh<-data_group[data_group$Stage=='Fresh',-1:-2]
Aged<-data_group[data_group$Stage=='Aged',-1:-2]

Result=data.frame()
for(i in 1:100){
  wil <- wilcox.test(Fresh[,i],Aged[,i],alternative='two.sided',correct=T)
  dif<-mean(Aged[,i])-mean(Fresh[,i])
  Result=rbind(Result,
                 cbind(id=colnames(data)[i],
                 dif,
                 pvalue=wil$p.value)
                  )
}

library(dplyr)
Result_QMP<-mutate(Result,rd = cut(as.numeric(Result$dif), breaks = c(-Inf, 0, Inf),
                  labels = c("Decreased", "Increased")),
         pd = cut(as.numeric(Result$pvalue), breaks = c(-Inf, 0.05, Inf),
                  labels = c("< 0.05", "ns")))

#changes in RMP data of fungal communities 
data<-read.csv("RMP_fungi_genera_MS_top100.csv",row.names=1)
data<-t(data)

group<-read.csv('Group_MS.csv',row.names=1)

data_group<-cbind(group,data)

Fresh<-data_group[data_group$Stage=='Fresh',-1:-2]
Aged<-data_group[data_group$Stage=='Aged',-1:-2]

Result=data.frame()
for(i in 1:100){
  wil <- wilcox.test(Fresh[,i],Aged[,i],alternative='two.sided',correct=T)
  dif<-mean(Aged[,i])-mean(Fresh[,i])
  Result=rbind(Result,
                 cbind(id=colnames(data)[i],
                 dif,
                 pvalue=wil$p.value)
                  )
}

library(dplyr)
Result_RMP<-mutate(Result,rd = cut(as.numeric(Result$dif), breaks = c(-Inf, 0, Inf),
                  labels = c("Decreased", "Increased")),
         pd = cut(as.numeric(Result$pvalue), breaks = c(-Inf, 0.05, Inf),
                  labels = c("< 0.05", "ns")))


upset<-cbind(data.frame(Result_QMP[,c(1,4,5)]),data.frame(Result_RMP[,4:5]))
colnames(upset)<-c('ID','QMP_rd','QMP_pd','RMP_rd','RMP_pd')

QMP_increased<-upset[upset$QMP_rd=='Increased'&upset$QMP_pd=='< 0.05',1]
QMP_decreased<-upset[upset$QMP_rd=='Decreased'&upset$QMP_pd=='< 0.05',1]
QMP_ns<-upset[upset$QMP_pd=='ns',1]
RMP_increased<-upset[upset$RMP_rd=='Increased'&upset$RMP_pd=='< 0.05',1]
RMP_decreased<-upset[upset$RMP_rd=='Decreased'&upset$RMP_pd=='< 0.05',1]
RMP_ns<-upset[upset$RMP_pd=='ns',1]

input<-list(QMP_increased=QMP_increased,QMP_decreased=QMP_decreased,QMP_ns=QMP_ns,RMP_increased=RMP_increased,RMP_decreased=RMP_decreased,RMP_ns=RMP_ns)

library(UpSetR)
upset(fromList(input),
sets = c("RMP_ns", "QMP_ns", "RMP_increased", "QMP_increased",  "RMP_decreased","QMP_decreased"),
nsets=6,mb.ratio = c(0.6, 0.4),
keep.order = TRUE)



