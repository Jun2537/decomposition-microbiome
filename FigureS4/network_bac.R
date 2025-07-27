setwd('E:/abs/network')

library(WGCNA)
library(ggplot2)  
library(vegan)
library(psych)
library(igraph)

data<-read.csv('tax.csv')
mydata<-unique(data)
mydata<-mydata[mydata$Genera!='Unassigned',]
write.csv(mydata,'genera_tax.csv')


rm(list=ls())
otu<-read.csv('AA_genera_GS_top100.csv',row.names=1)



cor<- corAndPvalue(t(otu),use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')

r[p > 0.01 | abs(r) < 0.60] = 0
write.csv(data.frame(r, check.names = FALSE), 'corr.matrix_AA_GS.csv')

g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) 
g = delete.vertices(g, names(degree(g)[degree(g) == 0]))   

E(g)$corr = E(g)$weight            
E(g)$weight = abs(E(g)$weight)  

tax_all = read.csv('genera_tax.csv',row.names=1)  
tax<-tax_all[rownames(otu),]
write.csv(tax,'AA_GS_tax.csv')

tax<-read.csv('AA_GS_tax.csv',row.names=1)

tax = tax[as.character(V(g)$name), ]                  #为节点加上分类信息

V(g)$Phylum = tax$Phylum                              #门


node_list = data.frame(
  label = names(V(g)),
  phylum = V(g)$Phylum
  )  

head(node_list)
write.csv(node_list, 'network.node_AA_GS_bac.csv')   


edge = data.frame(as_edgelist(g))                      #创建边列表
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

edge_list$PN<- cut(edge_list$correlation, breaks = c(-Inf,0, Inf), labels = c("neg","pos"), right=FALSE)

head(edge_list)
write.csv(edge_list, 'network.edge_list_AA_GS_bac.csv')  

write_graph(g, 'network_AA_GS_bac.graphml', format = 'graphml') 




rm(list=ls())
otu<-read.csv('RA_genera_GS_top100.csv',row.names=1)



cor<- corAndPvalue(t(otu),use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')

r[p > 0.01 | abs(r) < 0.60] = 0
write.csv(data.frame(r, check.names = FALSE), 'corr.matrix_RA_GS.csv')

g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) 
g = delete.vertices(g, names(degree(g)[degree(g) == 0]))   

E(g)$corr = E(g)$weight            
E(g)$weight = abs(E(g)$weight)  

tax_all = read.csv('genera_tax.csv',row.names=1)  
tax<-tax_all[rownames(otu),]
write.csv(tax,'RA_GS_tax.csv')

tax<-read.csv('RA_GS_tax.csv',row.names=1)

tax = tax[as.character(V(g)$name), ]                  #为节点加上分类信息

V(g)$Phylum = tax$Phylum                              #门


node_list = data.frame(
  label = names(V(g)),
  phylum = V(g)$Phylum
  )  

head(node_list)
write.csv(node_list, 'network.node_RA_GS_bac.csv')   


edge = data.frame(as_edgelist(g))                      #创建边列表
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

edge_list$PN<- cut(edge_list$correlation, breaks = c(-Inf,0, Inf), labels = c("neg","pos"), right=FALSE)

head(edge_list)
write.csv(edge_list, 'network.edge_list_RA_GS_bac.csv')  

write_graph(g, 'network_RA_GS_bac.graphml', format = 'graphml') 



rm(list=ls())
otu<-read.csv('AA_genera_MS_top100.csv',row.names=1)

cor<- corAndPvalue(t(otu),use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')

r[p > 0.01 | abs(r) < 0.60] = 0
write.csv(data.frame(r, check.names = FALSE), 'corr.matrix_AA_MS.csv')

g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) 
g = delete.vertices(g, names(degree(g)[degree(g) == 0]))   

E(g)$corr = E(g)$weight            
E(g)$weight = abs(E(g)$weight)  

tax_all = read.csv('genera_tax.csv',row.names=1)  
tax<-tax_all[rownames(otu),]
write.csv(tax,'AA_MS_tax.csv')

tax<-read.csv('AA_MS_tax.csv',row.names=1)

tax = tax[as.character(V(g)$name), ]                  #为节点加上分类信息

V(g)$Phylum = tax$Phylum                              #门


node_list = data.frame(
  label = names(V(g)),
  phylum = V(g)$Phylum
  )  

head(node_list)
write.csv(node_list, 'network.node_AA_MS_bac.csv')   


edge = data.frame(as_edgelist(g))                      #创建边列表
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

edge_list$PN<- cut(edge_list$correlation, breaks = c(-Inf,0, Inf), labels = c("neg","pos"), right=FALSE)

head(edge_list)
write.csv(edge_list, 'network.edge_list_AA_MS_bac.csv')  

write_graph(g, 'network_AA_MS_bac.graphml', format = 'graphml') 

rm(list=ls())
otu<-read.csv('RA_genera_MS_top100.csv',row.names=1)

cor<- corAndPvalue(t(otu),use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')

r[p > 0.01 | abs(r) < 0.60] = 0
write.csv(data.frame(r, check.names = FALSE), 'corr.matrix_RA_MS.csv')

g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) 
g = delete.vertices(g, names(degree(g)[degree(g) == 0]))   

E(g)$corr = E(g)$weight            
E(g)$weight = abs(E(g)$weight)  

tax_all = read.csv('genera_tax.csv',row.names=1)  
tax<-tax_all[rownames(otu),]
write.csv(tax,'RA_MS_tax.csv')

tax<-read.csv('RA_MS_tax.csv',row.names=1)

tax = tax[as.character(V(g)$name), ]                  #为节点加上分类信息

V(g)$Phylum = tax$Phylum                              #门


node_list = data.frame(
  label = names(V(g)),
  phylum = V(g)$Phylum
  )  

head(node_list)
write.csv(node_list, 'network.node_RA_MS_bac.csv')   


edge = data.frame(as_edgelist(g))                      #创建边列表
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

edge_list$PN<- cut(edge_list$correlation, breaks = c(-Inf,0, Inf), labels = c("neg","pos"), right=FALSE)

head(edge_list)
write.csv(edge_list, 'network.edge_list_RA_MS_bac.csv')  

write_graph(g, 'network_RA_MS_bac.graphml', format = 'graphml') 






