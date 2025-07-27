setwd('E:/abs/Git/FigureS4')

library(WGCNA)
library(ggplot2)  
library(vegan)
library(psych)
library(igraph)
library(reshape2)

rm(list=ls())

#generate QMP bacterial network in gravesoil
otu<-read.csv('QMP_bac_genera_GS_top100.csv',row.names=1)

cor<- corAndPvalue(t(otu),use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')
r[p > 0.01 | abs(r) < 0.60] = 0

g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) 
g = delete.vertices(g, names(degree(g)[degree(g) == 0]))   

E(g)$corr = E(g)$weight            
E(g)$weight = abs(E(g)$weight)  

tax_all = read.csv('bac_genera_tax.csv',row.names=1)  
tax<-tax_all[rownames(otu),]
tax<-data.frame(cbind(rownames(otu),tax))
colnames(tax)<-c('Genera','Phylum')
rownames(tax)<-tax[,1]

tax = tax[as.character(V(g)$name), ]                  

V(g)$Phylum = tax$Phylum                              


node_list = data.frame(
  label = names(V(g)),
  phylum = V(g)$Phylum
  )  

head(node_list)
write.csv(node_list, 'network.node_QMP_bac_GS.csv')   


edge = data.frame(as_edgelist(g))                      
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

edge_list$PN<- cut(edge_list$correlation, breaks = c(-Inf,0, Inf), labels = c("neg","pos"), right=FALSE)

df_QMP_bac<-data.frame(cbind(paste(edge_list$source,edge_list$target,sep='_'),edge_list$PN))


head(edge_list)
write.csv(edge_list, 'network.edge_list_QMP_bac_GS.csv')  
write_graph(g, 'network_QMP_bac_GS.graphml', format = 'graphml') 


#generate RMP bacterial network in gravesoil
otu<-read.csv('RMP_bac_genera_GS_top100.csv',row.names=1)

cor<- corAndPvalue(t(otu),use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')
r[p > 0.01 | abs(r) < 0.60] = 0

g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) 
g = delete.vertices(g, names(degree(g)[degree(g) == 0]))   

E(g)$corr = E(g)$weight            
E(g)$weight = abs(E(g)$weight)  

tax_all = read.csv('bac_genera_tax.csv',row.names=1)  
tax<-tax_all[rownames(otu),]
tax<-data.frame(cbind(rownames(otu),tax))
colnames(tax)<-c('Genera','Phylum')
rownames(tax)<-tax[,1]

tax = tax[as.character(V(g)$name), ]                  

V(g)$Phylum = tax$Phylum                              


node_list = data.frame(
  label = names(V(g)),
  phylum = V(g)$Phylum
  )  

head(node_list)
write.csv(node_list, 'network.node_RMP_bac_GS.csv')   


edge = data.frame(as_edgelist(g))                      
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

edge_list$PN<- cut(edge_list$correlation, breaks = c(-Inf,0, Inf), labels = c("neg","pos"), right=FALSE)
df_RMP_bac<-data.frame(cbind(paste(edge_list$source,edge_list$target,sep='_'),edge_list$PN))


head(edge_list)
write.csv(edge_list, 'network.edge_list_RMP_bac_GS.csv')  
write_graph(g, 'network_RMP_bac_GS.graphml', format = 'graphml') 

#generate confusion matrix
df1<-data.frame(df_QMP_bac[,2])
df2<-data.frame(df_RMP_bac[,2])

rownames(df1)<-df_QMP_bac[,]$X1
rownames(df2)<-df_RMP_bac[,]$X1

colnames(df1)<-'QMP'
colnames(df2)<-'RMP'

merge<-merge(df1,df2,by=0,all=T)
merge<-replace(merge,merge=='2','1Pos')
merge<-replace(merge,merge=='1','2Neg')
merge[is.na(merge)]<-'3NS'

cm<-table(merge$QMP,merge$RMP)

confusion_matrix_df <- as.data.frame.matrix(cm)
confusion_matrix_df $Response <- rownames(confusion_matrix_df )
draw_data <- melt(confusion_matrix_df )

colnames(draw_data)<-c("QMP","RMP","num")
library(ggplot2)
ggplot(draw_data,aes(QMP,RMP,fill = num))+ 
geom_tile() +
geom_text(aes(label = num),col='red',size=6)+
scale_fill_gradientn(colors = rev(hcl.colors(10,"Blues")))+ 
coord_fixed() + 
theme_minimal()



rm(list=ls())
#generate QMP bacterial network in tissue
otu<-read.csv('QMP_bac_genera_MS_top100.csv',row.names=1)

cor<- corAndPvalue(t(otu),use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')
r[p > 0.01 | abs(r) < 0.60] = 0

g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) 
g = delete.vertices(g, names(degree(g)[degree(g) == 0]))   

E(g)$corr = E(g)$weight            
E(g)$weight = abs(E(g)$weight)  

tax_all = read.csv('bac_genera_tax.csv',row.names=1)  
tax<-tax_all[rownames(otu),]
tax<-data.frame(cbind(rownames(otu),tax))
colnames(tax)<-c('Genera','Phylum')
rownames(tax)<-tax[,1]

tax = tax[as.character(V(g)$name), ]                  

V(g)$Phylum = tax$Phylum                              


node_list = data.frame(
  label = names(V(g)),
  phylum = V(g)$Phylum
  )  

head(node_list)
write.csv(node_list, 'network.node_QMP_bac_MS.csv')   


edge = data.frame(as_edgelist(g))                      
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

edge_list$PN<- cut(edge_list$correlation, breaks = c(-Inf,0, Inf), labels = c("neg","pos"), right=FALSE)

df_QMP_bac<-data.frame(cbind(paste(edge_list$source,edge_list$target,sep='_'),edge_list$PN))


head(edge_list)
write.csv(edge_list, 'network.edge_list_QMP_bac_MS.csv')  
write_graph(g, 'network_QMP_bac_MS.graphml', format = 'graphml') 


#generate RMP bacterial network in tissue
otu<-read.csv('RMP_bac_genera_MS_top100.csv',row.names=1)

cor<- corAndPvalue(t(otu),use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')
r[p > 0.01 | abs(r) < 0.60] = 0

g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) 
g = delete.vertices(g, names(degree(g)[degree(g) == 0]))   

E(g)$corr = E(g)$weight            
E(g)$weight = abs(E(g)$weight)  

tax_all = read.csv('bac_genera_tax.csv',row.names=1)  
tax<-tax_all[rownames(otu),]
tax<-data.frame(cbind(rownames(otu),tax))
colnames(tax)<-c('Genera','Phylum')
rownames(tax)<-tax[,1]

tax = tax[as.character(V(g)$name), ]                  

V(g)$Phylum = tax$Phylum                              


node_list = data.frame(
  label = names(V(g)),
  phylum = V(g)$Phylum
  )  

head(node_list)
write.csv(node_list, 'network.node_RMP_bac_MS.csv')   


edge = data.frame(as_edgelist(g))                      
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

edge_list$PN<- cut(edge_list$correlation, breaks = c(-Inf,0, Inf), labels = c("neg","pos"), right=FALSE)
df_RMP_bac<-data.frame(cbind(paste(edge_list$source,edge_list$target,sep='_'),edge_list$PN))


head(edge_list)
write.csv(edge_list, 'network.edge_list_RMP_bac_MS.csv')  
write_graph(g, 'network_RMP_bac_MS.graphml', format = 'graphml') 

#generate confusion matrix
df1<-data.frame(df_QMP_bac[,2])
df2<-data.frame(df_RMP_bac[,2])

rownames(df1)<-df_QMP_bac[,]$X1
rownames(df2)<-df_RMP_bac[,]$X1

colnames(df1)<-'QMP'
colnames(df2)<-'RMP'

merge<-merge(df1,df2,by=0,all=T)
merge<-replace(merge,merge=='2','1Pos')
merge<-replace(merge,merge=='1','2Neg')
merge[is.na(merge)]<-'3NS'

cm<-table(merge$QMP,merge$RMP)

confusion_matrix_df <- as.data.frame.matrix(cm)
confusion_matrix_df $Response <- rownames(confusion_matrix_df )
draw_data <- melt(confusion_matrix_df )

colnames(draw_data)<-c("QMP","RMP","num")
library(ggplot2)
ggplot(draw_data,aes(QMP,RMP,fill = num))+ 
geom_tile() +
geom_text(aes(label = num),col='red',size=6)+
scale_fill_gradientn(colors = rev(hcl.colors(10,"Blues")))+ 
coord_fixed() + 
theme_minimal()



rm(list=ls())
#generate QMP fungal network in gravesoil
otu<-read.csv('QMP_fungi_genera_GS_top100.csv',row.names=1)

cor<- corAndPvalue(t(otu),use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')
r[p > 0.01 | abs(r) < 0.60] = 0

g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) 
g = delete.vertices(g, names(degree(g)[degree(g) == 0]))   

E(g)$corr = E(g)$weight            
E(g)$weight = abs(E(g)$weight)  

tax_all = read.csv('fungi_genera_tax.csv',row.names=1)  
tax<-tax_all[rownames(otu),]
tax<-data.frame(cbind(rownames(otu),tax))
colnames(tax)<-c('Genera','Phylum')
rownames(tax)<-tax[,1]

tax = tax[as.character(V(g)$name), ]                  

V(g)$Phylum = tax$Phylum                              


node_list = data.frame(
  label = names(V(g)),
  phylum = V(g)$Phylum
  )  

head(node_list)
write.csv(node_list, 'network.node_QMP_fungi_GS.csv')   


edge = data.frame(as_edgelist(g))                      
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

edge_list$PN<- cut(edge_list$correlation, breaks = c(-Inf,0, Inf), labels = c("neg","pos"), right=FALSE)

df_QMP_fungi<-data.frame(cbind(paste(edge_list$source,edge_list$target,sep='_'),edge_list$PN))


head(edge_list)
write.csv(edge_list, 'network.edge_list_QMP_fungi_GS.csv')  
write_graph(g, 'network_QMP_fungi_GS.graphml', format = 'graphml') 


#generate RMP fungal network in gravesoil
otu<-read.csv('RMP_fungi_genera_GS_top100.csv',row.names=1)

cor<- corAndPvalue(t(otu),use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')
r[p > 0.01 | abs(r) < 0.60] = 0

g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) 
g = delete.vertices(g, names(degree(g)[degree(g) == 0]))   

E(g)$corr = E(g)$weight            
E(g)$weight = abs(E(g)$weight)  

tax_all = read.csv('fungi_genera_tax.csv',row.names=1)  
tax<-tax_all[rownames(otu),]
tax<-data.frame(cbind(rownames(otu),tax))
colnames(tax)<-c('Genera','Phylum')
rownames(tax)<-tax[,1]

tax = tax[as.character(V(g)$name), ]                  

V(g)$Phylum = tax$Phylum                              


node_list = data.frame(
  label = names(V(g)),
  phylum = V(g)$Phylum
  )  

head(node_list)
write.csv(node_list, 'network.node_RMP_fungi_GS.csv')   


edge = data.frame(as_edgelist(g))                      
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

edge_list$PN<- cut(edge_list$correlation, breaks = c(-Inf,0, Inf), labels = c("neg","pos"), right=FALSE)
df_RMP_fungi<-data.frame(cbind(paste(edge_list$source,edge_list$target,sep='_'),edge_list$PN))


head(edge_list)
write.csv(edge_list, 'network.edge_list_RMP_fungi_GS.csv')  
write_graph(g, 'network_RMP_fungi_GS.graphml', format = 'graphml') 

#generate confusion matrix
df1<-data.frame(df_QMP_fungi[,2])
df2<-data.frame(df_RMP_fungi[,2])

rownames(df1)<-df_QMP_fungi[,]$X1
rownames(df2)<-df_RMP_fungi[,]$X1

colnames(df1)<-'QMP'
colnames(df2)<-'RMP'

merge<-merge(df1,df2,by=0,all=T)
merge<-replace(merge,merge=='2','1Pos')
merge<-replace(merge,merge=='1','2Neg')
merge[is.na(merge)]<-'3NS'

cm<-table(merge$QMP,merge$RMP)

confusion_matrix_df <- as.data.frame.matrix(cm)
confusion_matrix_df $Response <- rownames(confusion_matrix_df )
draw_data <- melt(confusion_matrix_df )

colnames(draw_data)<-c("QMP","RMP","num")
library(ggplot2)
ggplot(draw_data,aes(QMP,RMP,fill = num))+ 
geom_tile() +
geom_text(aes(label = num),col='red',size=6)+
scale_fill_gradientn(colors = rev(hcl.colors(10,"Blues")))+ 
coord_fixed() + 
theme_minimal()





rm(list=ls())
#generate QMP fungal network in tissue
otu<-read.csv('QMP_fungi_genera_MS_top100.csv',row.names=1)

cor<- corAndPvalue(t(otu),use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')
r[p > 0.01 | abs(r) < 0.60] = 0

g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) 
g = delete.vertices(g, names(degree(g)[degree(g) == 0]))   

E(g)$corr = E(g)$weight            
E(g)$weight = abs(E(g)$weight)  

tax_all = read.csv('fungi_genera_tax.csv',row.names=1)  
tax<-tax_all[rownames(otu),]
tax<-data.frame(cbind(rownames(otu),tax))
colnames(tax)<-c('Genera','Phylum')
rownames(tax)<-tax[,1]

tax = tax[as.character(V(g)$name), ]                  

V(g)$Phylum = tax$Phylum                              


node_list = data.frame(
  label = names(V(g)),
  phylum = V(g)$Phylum
  )  

head(node_list)
write.csv(node_list, 'network.node_QMP_fungi_MS.csv')   


edge = data.frame(as_edgelist(g))                      
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

edge_list$PN<- cut(edge_list$correlation, breaks = c(-Inf,0, Inf), labels = c("neg","pos"), right=FALSE)

df_QMP_fungi<-data.frame(cbind(paste(edge_list$source,edge_list$target,sep='_'),edge_list$PN))


head(edge_list)
write.csv(edge_list, 'network.edge_list_QMP_fungi_MS.csv')  
write_graph(g, 'network_QMP_fungi_MS.graphml', format = 'graphml') 


#generate RMP fungal network in tissue
otu<-read.csv('RMP_fungi_genera_MS_top100.csv',row.names=1)

cor<- corAndPvalue(t(otu),use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')
r[p > 0.01 | abs(r) < 0.60] = 0

g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) 
g = delete.vertices(g, names(degree(g)[degree(g) == 0]))   

E(g)$corr = E(g)$weight            
E(g)$weight = abs(E(g)$weight)  

tax_all = read.csv('fungi_genera_tax.csv',row.names=1)  
tax<-tax_all[rownames(otu),]
tax<-data.frame(cbind(rownames(otu),tax))
colnames(tax)<-c('Genera','Phylum')
rownames(tax)<-tax[,1]

tax = tax[as.character(V(g)$name), ]                  

V(g)$Phylum = tax$Phylum                              


node_list = data.frame(
  label = names(V(g)),
  phylum = V(g)$Phylum
  )  

head(node_list)
write.csv(node_list, 'network.node_RMP_fungi_MS.csv')   


edge = data.frame(as_edgelist(g))                      
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

edge_list$PN<- cut(edge_list$correlation, breaks = c(-Inf,0, Inf), labels = c("neg","pos"), right=FALSE)
df_RMP_fungi<-data.frame(cbind(paste(edge_list$source,edge_list$target,sep='_'),edge_list$PN))


head(edge_list)
write.csv(edge_list, 'network.edge_list_RMP_fungi_MS.csv')  
write_graph(g, 'network_RMP_fungi_MS.graphml', format = 'graphml') 

#generate confusion matrix
df1<-data.frame(df_QMP_fungi[,2])
df2<-data.frame(df_RMP_fungi[,2])

rownames(df1)<-df_QMP_fungi[,]$X1
rownames(df2)<-df_RMP_fungi[,]$X1

colnames(df1)<-'QMP'
colnames(df2)<-'RMP'

merge<-merge(df1,df2,by=0,all=T)
merge<-replace(merge,merge=='2','1Pos')
merge<-replace(merge,merge=='1','2Neg')
merge[is.na(merge)]<-'3NS'

cm<-table(merge$QMP,merge$RMP)

confusion_matrix_df <- as.data.frame.matrix(cm)
confusion_matrix_df $Response <- rownames(confusion_matrix_df )
draw_data <- melt(confusion_matrix_df )

colnames(draw_data)<-c("QMP","RMP","num")
library(ggplot2)
ggplot(draw_data,aes(QMP,RMP,fill = num))+ 
geom_tile() +
geom_text(aes(label = num),col='red',size=6)+
scale_fill_gradientn(colors = rev(hcl.colors(10,"Blues")))+ 
coord_fixed() + 
theme_minimal()




