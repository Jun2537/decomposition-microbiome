setwd('E:/abs/Git/Figure1')
library(ggplot2)

#Figure1A_bacteria
data<-read.csv('copies_shannon_bacteria.csv',row.names=1)
data$Day<-as.numeric(as.character(data$Day))

ggplot(data, aes(x=Day, y=log10(copies),color=Bodysite,group=Bodysite))+
  geom_smooth(size = 1) +
  geom_jitter(width = 0.1,alpha = .5,size = 1) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#3C5488')) +
  theme_bw() +
  ylab('Log10 (copies per gram)') +
  xlab("Postmortem interval (day) ") +
  scale_x_continuous(breaks = c(0,1,3,7,21,28,35)) +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA))

#Figure1A_fungi
data<-read.csv('copies_shannon_fungi.csv',row.names=1)
data$Day<-as.numeric(as.character(data$Day))

ggplot(data, aes(x=Day, y=log10(copies),color=Bodysite,group=Bodysite))+
  geom_smooth(size = 1) +
  geom_jitter(width = 0.1,alpha = .5,size = 1) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#3C5488')) +
  theme_bw() +
  ylab('Log10 (copies per gram)') +
  xlab("Postmortem interval (day) ") +
  scale_x_continuous(breaks = c(0,1,3,7,21,28,35)) +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA))


#Figure1B_bacteria
data<-read.csv('copies_shannon_bacteria.csv',row.names=1)
data$Day<-as.numeric(as.character(data$Day))
ggplot(data, aes(x=Day, y=Shannon,color=Bodysite,group=Bodysite))+
  geom_smooth(size = 1) +
  geom_jitter(width = 0.1,alpha = .5,size = 1) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#3C5488')) +
  theme_bw() +
  ylab('Shannon index') +
  xlab("Postmortem interval (day) ") +
  scale_x_continuous(breaks = c(0,1,3,7,21,28,35)) +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) 



#Figure1B_fungi
data<-read.csv('copies_shannon_fungi.csv',row.names=1)
data$Day<-as.numeric(as.character(data$Day))
ggplot(data, aes(x=Day, y=Shannon,color=Bodysite,group=Bodysite))+
  geom_smooth(size = 1) +
  geom_jitter(width = 0.1,alpha = .5,size = 1) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#3C5488')) +
  theme_bw() +
  ylab('Shannon index') +
  xlab("Postmortem interval (day) ") +
  scale_x_continuous(breaks = c(0,1,3,7,21,28,35)) +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) 


#Figure1C_bacteria
library(ape)
library(RColorBrewer)
library(ggplot2)
library(vegan)

data<-read.csv('bac_ASV_QMP.csv',row.names=1)
data<-t(data)
group<-read.csv('Group.csv',row.names=1)
data_group<-cbind(data,group)
dis<-vegdist(data,'bray')

pcoa<-pcoa(dis)
axes<-as.data.frame(pcoa$vectors)
axes<-cbind(axes,group)

eigval <- round(pcoa$values$Relative_eig * 100, digits = 2)
eigval <- data.frame( PC = 1:length(eigval), Eigval = eigval)

mycolor<- brewer.pal(9,"Paired")
p <- ggplot(axes, aes(Axis.1, Axis.2,colour =Day2,fill=Day2,shape=Bodysite))+
  geom_point(size = 4,alpha=0.7) +
scale_fill_manual(values = mycolor) +
scale_colour_manual(values = mycolor) +
 scale_shape_manual(values=c(8,16,17))+
   xlab(paste("PCo1 (", eigval$Eigval[1], " %)", sep = "")) +
  ylab(paste("PCo2 (", eigval$Eigval[2], " %)", sep = "")) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw()

#Figure1C_fungi
data<-read.csv('Fun_ASV_QMP.csv',row.names=1)
data<-t(data)
group<-read.csv('Group.csv',row.names=1)
data_group<-cbind(data,group)
dis<-vegdist(data,'bray')

pcoa<-pcoa(dis)
axes<-as.data.frame(pcoa$vectors)
axes<-cbind(axes,group)

eigval <- round(pcoa$values$Relative_eig * 100, digits = 2)
eigval <- data.frame( PC = 1:length(eigval), Eigval = eigval)

mycolor<- brewer.pal(9,"Paired")
p <- ggplot(axes, aes(Axis.1, Axis.2,colour =Day2,fill=Day2,shape=Bodysite))+
  geom_point(size = 4,alpha=0.7) +
scale_fill_manual(values = mycolor) +
scale_colour_manual(values = mycolor) +
 scale_shape_manual(values=c(8,16,17))+
   xlab(paste("PCo1 (", eigval$Eigval[1], " %)", sep = "")) +
  ylab(paste("PCo2 (", eigval$Eigval[2], " %)", sep = "")) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw()




