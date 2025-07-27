#Figure6A
setwd('E:/abs/Git/Figure6')
rm(list=ls())
library(igraph)
data<-read.csv('Genera_bac_fun_QMP.csv',row.names=1)

#Select the bacterial genera that are present in more than three samples
pa<-decostand(data,'pa')
pa<-pa[rowSums(pa)>3,]
otu<-data[rownames(pa),]

cor<- corAndPvalue(t(otu),use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')

r[p > 0.001 | abs(r) < 0.60] = 0

g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) 
g = delete.vertices(g, names(degree(g)[degree(g) == 0]))   

E(g)$corr = E(g)$weight            
E(g)$weight = abs(E(g)$weight)  

#generate network module
M<-cluster_fast_greedy(as.undirected(g))
#Extract the names of the nodes in the module
M[[1]]
M[[2]]
...

#Read taxonomic information
tax = read.csv('tax_anno.csv', row.names=1, header=T)       

tax = tax[as.character(V(g)$name), ]                  
V(g)$Kingdom = tax$Kingdom                            
V(g)$Genera = tax$Genera                              

#node list
node_list = data.frame(
  label = names(V(g)),
  kingdom = V(g)$Kingdom,
  genera=V(g)$Genera)  
  
head(node_list)
write.csv(node_list, 'network.node_list.csv')   
   
#edge list
edge = data.frame(as_edgelist(g))                      
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

head(edge_list)
write.csv(edge_list, 'network.edge_list.csv')  

write_graph(g, 'network.graphml', format = 'graphml') 

#Calculate network parameters
nodes_num = length(V(g))                   
edges_num = length(E(g))                   
positive.cor_num = sum(E(g)$corr>0)        
negative.cor_num = sum(E(g)$corr<0)        
average_degree = mean(degree(g))           
average_path_length = average.path.length(g, directed = FALSE)     
network_diameter = diameter(g, directed = FALSE)                   
network_density = graph.density(g)                                 
clustering_coefficient = transitivity(g)                           




#Figure6B
setwd('E:/abs/Git/Figure6')

mod1<-read.csv('Module1_genera.csv',row.names=1)
mod1<-t(mod1)
mod1<-decostand(mod1,'hellinger')
pca_result<-prcomp(mod1,center=T, scale= F)
pc_mod1<-pca_result$x


mod2<-read.csv('Module2_genera.csv',row.names=1)
mod2<-t(mod2)
mod2<-decostand(mod2,'hellinger')
pca_result<-prcomp(mod2,center=T, scale= F)
pc_mod2<-pca_result$x

mod3<-read.csv('Module3_genera.csv',row.names=1)
mod3<-t(mod3)
mod3<-decostand(mod3,'hellinger')
pca_result<-prcomp(mod3,center=T, scale= F)
pc_mod3<-pca_result$x

mod4<-read.csv('Module4_genera.csv',row.names=1)
mod4<-t(mod4)
mod4<-decostand(mod4,'hellinger')
pca_result<-prcomp(mod4,center=T, scale= F)
pc_mod4<-pca_result$x

mod5<-read.csv('Module5_genera.csv',row.names=1)
mod5<-t(mod5)
mod5<-decostand(mod5,'hellinger')
pca_result<-prcomp(mod5,center=T, scale= F)
pc_mod5<-pca_result$x

mod6<-read.csv('Module6_genera.csv',row.names=1)
mod6<-t(mod6)
mod6<-decostand(mod6,'hellinger')
pca_result<-prcomp(mod6,center=T, scale= F)
pc_mod6<-pca_result$x

pc<-cbind(pc_mod1[,1],pc_mod2[,1],pc_mod3[,1],pc_mod4[,1],pc_mod5[,1],pc_mod6[,1])

metab<-read.csv('class_normal.csv',row.names=1)

data<-cbind(pc,metab)

library(WGCNA)
library(ggplot2)  
library(vegan)
library(psych)
library(igraph)

cor<- corAndPvalue(data,use="pairwise.complete.obs",method = "spearman",alternative='two.sided')
r = cor$cor
p = cor$p                     
p = p.adjust(p, method = 'BH')

r[p > 0.001 | abs(r) < 0.60] = 0

r<-r[-1:-6,1:6]
r<-r[abs(rowSums(r))>0,]
nrow(r)

r[r==0] = NA
write.csv(r,'r.csv',row.names=T)

#Manually correct the name error
r<-read.csv('r.csv',row.names=1)

anno<-read.csv('name_anno.csv',row.names=1)
r_anno<-anno[rownames(r),]

r<-cbind(r,r_anno)
r<-r[order(r[,7]),]



library(ComplexHeatmap)

right_annotation=HeatmapAnnotation(
Superclass=r[,7],which='row',show_annotation_name=FALSE,
col=list(Sperclass=c('Alkaloids and derivatives'='#B9717D',
'Benzenoids'='#F5C8B5',
'Lipids and lipid-like molecules'='#AD93B8',
'Nucleosides, nucleotides, andanalogues'='#7EB4C6',
'Organic acids and derivatives'='#98CA8B',
'Organic oxygen compounds'='#795B34',
'Organoheterocyclic compounds'='#E06E6E',
'Phenylpropanoids and polyketides'='#6E86AC'))
)

anno<-data.frame(r[,7])
colnames(anno)='anno'

mat<-data.frame(r[,-7])

bar_ha=HeatmapAnnotation(df=anno,which='row',show_annotation_name=T,
col=list(anno=c('Alkaloids and derivatives'='#B9717D',
'Benzenoids'='#F5C8B5',
'Lipids and lipid-like molecules'='#AD93B8',
'Nucleosides, nucleotides, and analogues'='#7EB4C6',
'Organic acids and derivatives'='#98CA8B',
'Organic oxygen compounds'='#795B34',
'Organoheterocyclic compounds'='#E06E6E',
'Phenylpropanoids and polyketides'='#6E86AC'))
)


Heatmap(mat,
col=colorRampPalette(c('#240EFF','#B491F9','#ECEAEF','#FFA48C','#FF1A09'))(100),
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        row_names_side =  'right',
        column_title = NULL,
        rect_gp = gpar(col = "white", lwd = 1),
        row_names_gp = gpar(fontsize = 10),
        row_split=anno,
        left_annotation =bar_ha,
        column_names_gp = gpar(fontsize = 10))









