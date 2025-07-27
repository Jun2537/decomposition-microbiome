setwd('E:/abs/Git/Figure3')
library(ggplot2)

#Figure1A_Gravesoil_bac
data<-read.csv('phylum_bacteria.csv',row.names=1)
data<-data[data$Type=="GS",]
data$Day<-as.numeric(as.character(data$Day))

#Pseudomonadota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Pseudomonadota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Pseudomonadota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 

#Bacillota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Bacillota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Bacillota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 

#Actinobacteriota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Actinobacteriota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Actinobacteriota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 


#Bacteroidota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Bacteroidota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Bacteroidota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 

#Acidobacteriota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Acidobacteriota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Acidobacteriota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 

#Figure1A_Tissue_bac
data<-read.csv('phylum_bacteria.csv',row.names=1)
data<-data[data$Type=="MS",]
data$Day<-as.numeric(as.character(data$Day))

#Pseudomonadota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Pseudomonadota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Pseudomonadota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 

#Bacillota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Bacillota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Bacillota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 

#Actinobacteriota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Actinobacteriota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Actinobacteriota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 


#Bacteroidota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Bacteroidota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Bacteroidota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 

#Acidobacteriota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Acidobacteriota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Acidobacteriota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 





#Figure1B_Gravesoil_Fungi
data<-read.csv('phylum_fungi.csv',row.names=1)
data<-data[data$Type=="GS",]
data$Day<-as.numeric(as.character(data$Day))

#Ascomycota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Ascomycota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Ascomycota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 

#Mortierellomycota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Mortierellomycota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Mortierellomycota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 

#Basidiomycota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Basidiomycota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Basidiomycota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 


#Mucoromycota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Mucoromycota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Mucoromycota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 

#Chytridiomycota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Chytridiomycota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Chytridiomycota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 


#Figure1B_Tissue_Fungi
data<-read.csv('phylum_fungi.csv',row.names=1)
data<-data[data$Type=="MS",]
data$Day<-as.numeric(as.character(data$Day))

#Ascomycota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Ascomycota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Ascomycota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 

#Mortierellomycota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Mortierellomycota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Mortierellomycota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 

#Basidiomycota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Basidiomycota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Basidiomycota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 


#Mucoromycota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Mucoromycota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Mucoromycota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 

#Chytridiomycota
ggplot(data, aes(x=Day)) +
  geom_smooth(aes(y=RA_Chytridiomycota*100),color='#3C5488')+
  geom_smooth(aes(y = log10(AA_Chytridiomycota+1)*10), color = '#E64B35') + 
  scale_y_continuous(name = "Relative abundance (%)",sec.axis = sec_axis( trans=~./10, name="Log10 (copies per gram)")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = NA)) +
  scale_x_continuous(breaks = c(1,3,7,14,21,28,35)) 




