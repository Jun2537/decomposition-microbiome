setwd('E:/abs/Git/FigureS3')

#time-decay for QMP bacteria
data<-read.csv('QMP_bac_ASV.csv',row.names=1)
data<-t(data)

GS_QMP<-data[49:90,]
MS_QMP<-data[91:132,]

group<-read.csv("Group.csv",row.names=1)
day<-group[49:90,2]

GS_sim_QMP<-1-vegdist(GS_QMP,"bray")
MS_sim_QMP<-1-vegdist(MS_QMP,"bray")

day_dist<-dist(day)

summary(lm(as.vector(log10(GS_sim_QMP+0.001))~as.vector(day_dist+0.1)))
summary(lm(as.vector(log10(MS_sim_QMP+0.001))~as.vector(day_dist+0.1)))

#time-decay for RMP bacteria
data<-read.csv('RMP_bac_ASV.csv',row.names=1)
data<-t(data)

GS_RMP<-data[49:90,]
MS_RMP<-data[91:132,]

group<-read.csv("Group.csv",row.names=1)
day<-group[49:90,2]

GS_sim_RMP<-1-vegdist(GS_RMP,"bray")
MS_sim_RMP<-1-vegdist(MS_RMP,"bray")

day_dist<-dist(day)

summary(lm(as.vector(log10(GS_sim_RMP+0.001))~as.vector(day_dist+0.1)))
summary(lm(as.vector(log10(MS_sim_RMP+0.001))~as.vector(day_dist+0.1)))

df1<-data.frame(cbind(as.vector(log10(GS_sim_QMP+0.001)),as.vector(log10(GS_sim_RMP+0.001)),as.vector(day_dist+0.1)))
colnames(df1)<-c('QMP','RMP','Day')
df2<-data.frame(cbind(as.vector(log10(MS_sim_QMP+0.001)),as.vector(log10(MS_sim_RMP+0.001)),as.vector(day_dist+0.1)))
colnames(df2)<-c('QMP','RMP','Day')

library(tidyr)
df1<-data.frame(pivot_longer(df1,cols=1:2,names_to='type',values_to='Similarity'))
df2<-data.frame(pivot_longer(df2,cols=1:2,names_to='type',values_to='Similarity'))

mycolor<-c('#E64B35','#3C5488')

#time-decay for bacterial communities in gravesoil
ggplot(df1, aes(x=Day, y=Similarity,color=type)) +
    geom_point(size = 1.5,alpha=0.2) + 
   scale_color_manual(values = mycolor)+      
    geom_smooth(method=lm)+
  xlab(paste("Day", sep = "")) +
  ylab(paste("Log10(Similarity)", sep = "")) +
    theme_bw()

#time-decay for bacterial communities in tissue
ggplot(df2, aes(x=Day, y=Similarity,color=type)) +
    geom_point(size = 1.5,alpha=0.2) + 
   scale_color_manual(values = mycolor)+      
    geom_smooth(method=lm)+
  xlab(paste("Day", sep = "")) +
  ylab(paste("Log10(Similarity)", sep = "")) +
    theme_bw()


#time-decay for QMP fungi
data<-read.csv('QMP_fun_ASV.csv',row.names=1)
data<-t(data)

GS_QMP<-data[49:90,]
MS_QMP<-data[91:132,]

group<-read.csv("Group.csv",row.names=1)
day<-group[49:90,2]

GS_sim_QMP<-1-vegdist(GS_QMP,"bray")
MS_sim_QMP<-1-vegdist(MS_QMP,"bray")

day_dist<-dist(day)

summary(lm(as.vector(log10(GS_sim_QMP+0.001))~as.vector(day_dist+0.1)))
summary(lm(as.vector(log10(MS_sim_QMP+0.001))~as.vector(day_dist+0.1)))

#time-decay for RMP fungi
data<-read.csv('RMP_fun_ASV.csv',row.names=1)
data<-t(data)

GS_RMP<-data[49:90,]
MS_RMP<-data[91:132,]

group<-read.csv("Group.csv",row.names=1)
day<-group[49:90,2]

GS_sim_RMP<-1-vegdist(GS_RMP,"bray")
MS_sim_RMP<-1-vegdist(MS_RMP,"bray")

day_dist<-dist(day)

summary(lm(as.vector(log10(GS_sim_RMP+0.001))~as.vector(day_dist+0.1)))
summary(lm(as.vector(log10(MS_sim_RMP+0.001))~as.vector(day_dist+0.1)))

df1<-data.frame(cbind(as.vector(log10(GS_sim_QMP+0.001)),as.vector(log10(GS_sim_RMP+0.001)),as.vector(day_dist+0.1)))
colnames(df1)<-c('QMP','RMP','Day')
df2<-data.frame(cbind(as.vector(log10(MS_sim_QMP+0.001)),as.vector(log10(MS_sim_RMP+0.001)),as.vector(day_dist+0.1)))
colnames(df2)<-c('QMP','RMP','Day')

library(tidyr)
df1<-data.frame(pivot_longer(df1,cols=1:2,names_to='type',values_to='Similarity'))
df2<-data.frame(pivot_longer(df2,cols=1:2,names_to='type',values_to='Similarity'))

mycolor<-c('#E64B35','#3C5488')

#time-decay for fungal communities in gravesoil
ggplot(df1, aes(x=Day, y=Similarity,color=type)) +
    geom_point(size = 1.5,alpha=0.2) + 
   scale_color_manual(values = mycolor)+      
    geom_smooth(method=lm)+
  xlab(paste("Day", sep = "")) +
  ylab(paste("Log10(Similarity)", sep = "")) +
    theme_bw()

#time-decay for fungal communities in tissue
ggplot(df2, aes(x=Day, y=Similarity,color=type)) +
    geom_point(size = 1.5,alpha=0.2) + 
   scale_color_manual(values = mycolor)+      
    geom_smooth(method=lm)+
  xlab(paste("Day", sep = "")) +
  ylab(paste("Log10(Similarity)", sep = "")) +
    theme_bw()


