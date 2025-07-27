rm(list=ls())

setwd('E:/abs/Git/Figure7andS6-S12/bac/MS')
library(randomForest)
library(sampling)

#model established and feature selected for Bacteria_QMP in Tissue
comm<-read.csv('MS_bac_genera_QMP.csv',row.names=1)
comm<-data.frame(t(comm))

Day<-read.csv('Group_MS.csv',row.names=1)

set.seed(1)
sub_train<-strata(Day,stratanames=("Day"),size=rep(4,7),method="srswor")
train_Day<-Day[sub_train$ID_unit,]
train_Data<-comm[sub_train$ID_unit,]

test_Day<-Day[-sub_train$ID_unit,]
test_Data<-comm[-sub_train$ID_unit,]


set.seed(1)
rf1 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE)  
print(rf1)


otu_train.cv <- replicate(5, rfcv(train_Data, train_Day, cv.fold=10,step=1.2), simplify = FALSE)
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)

#Fig.S6C
library(ggplot2)
p <- ggplot(otu_train.cv.mean, aes(Group.1, x)) +
    geom_line() +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
    labs(title = '',x = 'Number of ASVs', y = 'Cross-validation error')
p+geom_vline(xintercept = 9,color='red')


imp= as.data.frame(rf1$importance)
imp = imp[order(imp[,1],decreasing = T),]

train_select<-train_Data[,rownames(imp[1:9,])]
test_select<-test_Data[,rownames(imp[1:9,])]

n_train<-nrow(train_select)
n_test<-nrow(test_select)
predict_train<-matrix(NA,nrow=n_train,ncol=100)
predict_test<-matrix(NA,nrow=n_test,ncol=100)
r_squared<-1

for (i in 1:100) {
  set.seed(i)
  rf1 = randomForest(train_Day~.,data=train_select, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train[,i] <- predict(rf1,train_select)
  predict_test[,i] <- predict(rf1,test_select)
  r_squared[i]<-rf1$rsq[length(rf1$rsq)]
  i=i+1
}

mean(r_squared)

(mae_train<-mean(abs(rowMeans(predict_train)-train_Day)))
(sd_train<-sd(abs(rowMeans(predict_train)-train_Day)))

(mae_test<-mean(abs(rowMeans(predict_test)-test_Day)))
(sd_test<-sd(abs(rowMeans(predict_test)-test_Day)))

factor<-c(replicate(n_train,'train set'),replicate(n_test,'test_set'))
predict<-c(rowMeans(predict_train),rowMeans(predict_test))
True_value<-c(train_Day,test_Day)

df<-data.frame(cbind(factor,predict,True_value))

#Fig.7A_bacteria_QMP
library(ggplot2)
mycolor<-c('#E64B35','#3C5488')
p <- ggplot(df, aes(as.numeric(True_value), as.numeric(predict), color=factor)) +
       geom_point(size=2.5,alpha=0.8) +
       scale_color_manual(values = mycolor)+ 
       theme_bw()+
       theme(panel.grid =element_blank())+
       xlab(paste("True value (days)", sep = "")) +
       ylab(paste("Predictive value (days)", sep = "")) +
       geom_abline(intercept = 0, slope = 1)


#resutlt of entire dataset
n_train<-nrow(train_Data)
n_test<-nrow(test_Data)
predict_train_all<-matrix(NA,nrow=n_train,ncol=100)
predict_test_all<-matrix(NA,nrow=n_test,ncol=100)
r_squared_all<-1

for (i in 1:100) {
  set.seed(i)
  rf2 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train_all[,i] <- predict(rf2,train_Data)
  predict_test_all[,i] <- predict(rf2,test_Data)
  r_squared_all[i]<-rf2$rsq[length(rf2$rsq)]
  i=i+1
}

mean(r_squared_all)

(mae_train_all<-mean(abs(rowMeans(predict_train_all)-train_Day)))
(sd_train_all<-sd(abs(rowMeans(predict_train_all)-train_Day)))

(mae_test_all<-mean(abs(rowMeans(predict_test_all)-test_Day)))
(sd_test_all<-sd(abs(rowMeans(predict_test_all)-test_Day)))

wilcox.test(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),alternative='two.sided',paired=T)
wilcox.test(colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day)),alternative='two.sided',paired=T)
wilcox.test(r_squared,r_squared_all,alternative='two.sided',paired=T)


pred<-data.frame(cbind(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),
                 colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day))))


colnames(pred)<-c('train_select','train_all','test_select','test_all')
library(tidyr)

#Fig.S10A_bacteria_QMP
pred<-data.frame(pivot_longer(pred,cols=1:4,names_to='type',values_to='MAE'))
pred<-separate(pred,col='type',into=c('set','data'),sep='_')

ggplot(pred, aes(x=factor(set,levels=c('train','test')), y=MAE, fill=data)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Mean absolute error (days)", sep = ""))+
  theme_bw()


#Fig.S11A_bacteria_QMP
pred2<-data.frame(cbind(r_squared,r_squared_all))
colnames(pred2)<-c('selected','all')
pred2<-data.frame(pivot_longer(pred2,cols=1:2,names_to='type',values_to='R2'))
ggplot(pred2, aes(x=type, y=R2, fill=type)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Coefficient of determination", sep = ""))+
  theme_bw()


#model established and feature selected for Bacteria_RMP in Tissue
comm<-read.csv('MS_bac_genera_RMP.csv',row.names=1)
comm<-data.frame(t(comm))

Day<-read.csv('Group_MS.csv',row.names=1)

set.seed(1)
sub_train<-strata(Day,stratanames=("Day"),size=rep(4,7),method="srswor")
train_Day<-Day[sub_train$ID_unit,]
train_Data<-comm[sub_train$ID_unit,]

test_Day<-Day[-sub_train$ID_unit,]
test_Data<-comm[-sub_train$ID_unit,]


set.seed(1)
rf1 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE)  
print(rf1)

otu_train.cv <- replicate(5, rfcv(train_Data, train_Day, cv.fold=10,step=1.2), simplify = FALSE)
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)

#Fig.S6D
library(ggplot2)
p <- ggplot(otu_train.cv.mean, aes(Group.1, x)) +
    geom_line() +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
    labs(title = '',x = 'Number of ASVs', y = 'Cross-validation error')
p+geom_vline(xintercept = 29,color='red')


imp= as.data.frame(rf1$importance)
imp = imp[order(imp[,1],decreasing = T),]

train_select<-train_Data[,rownames(imp[1:29,])]
test_select<-test_Data[,rownames(imp[1:29,])]

n_train<-nrow(train_select)
n_test<-nrow(test_select)
predict_train<-matrix(NA,nrow=n_train,ncol=100)
predict_test<-matrix(NA,nrow=n_test,ncol=100)
r_squared<-1

for (i in 1:100) {
  set.seed(i)
  rf1 = randomForest(train_Day~.,data=train_select, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train[,i] <- predict(rf1,train_select)
  predict_test[,i] <- predict(rf1,test_select)
  r_squared[i]<-rf1$rsq[length(rf1$rsq)]
  i=i+1
}

mean(r_squared)

(mae_train<-mean(abs(rowMeans(predict_train)-train_Day)))
(sd_train<-sd(abs(rowMeans(predict_train)-train_Day)))

(mae_test<-mean(abs(rowMeans(predict_test)-test_Day)))
(sd_test<-sd(abs(rowMeans(predict_test)-test_Day)))

factor<-c(replicate(n_train,'train set'),replicate(n_test,'test_set'))
predict<-c(rowMeans(predict_train),rowMeans(predict_test))
True_value<-c(train_Day,test_Day)

df<-data.frame(cbind(factor,predict,True_value))

#Fig.7A_bacteria_RMP
library(ggplot2)
mycolor<-c('#E64B35','#3C5488')
p <- ggplot(df, aes(as.numeric(True_value), as.numeric(predict), color=factor)) +
       geom_point(size=2.5,alpha=0.8) +
       scale_color_manual(values = mycolor)+ 
       theme_bw()+
       theme(panel.grid =element_blank())+
       xlab(paste("True value (days)", sep = "")) +
       ylab(paste("Predictive value (days)", sep = "")) +
       geom_abline(intercept = 0, slope = 1)

#resutlt of entire dataset
n_train<-nrow(train_Data)
n_test<-nrow(test_Data)
predict_train_all<-matrix(NA,nrow=n_train,ncol=100)
predict_test_all<-matrix(NA,nrow=n_test,ncol=100)
r_squared_all<-1

for (i in 1:100) {
  set.seed(i)
  rf2 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train_all[,i] <- predict(rf2,train_Data)
  predict_test_all[,i] <- predict(rf2,test_Data)
  r_squared_all[i]<-rf2$rsq[length(rf2$rsq)]
  i=i+1
}

mean(r_squared_all)

(mae_train_all<-mean(abs(rowMeans(predict_train_all)-train_Day)))
(sd_train_all<-sd(abs(rowMeans(predict_train_all)-train_Day)))

(mae_test_all<-mean(abs(rowMeans(predict_test_all)-test_Day)))
(sd_test_all<-sd(abs(rowMeans(predict_test_all)-test_Day)))

wilcox.test(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),alternative='two.sided',paired=T)
wilcox.test(colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day)),alternative='two.sided',paired=T)
wilcox.test(r_squared,r_squared_all,alternative='two.sided',paired=T)


pred<-data.frame(cbind(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),
                 colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day))))


colnames(pred)<-c('train_select','train_all','test_select','test_all')

#Fig.S10A_bacateria_RMP
library(tidyr)
pred<-data.frame(pivot_longer(pred,cols=1:4,names_to='type',values_to='MAE'))
pred<-separate(pred,col='type',into=c('set','data'),sep='_')

ggplot(pred, aes(x=factor(set,levels=c('train','test')), y=MAE, fill=data)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Mean absolute error (days)", sep = ""))+
  theme_bw()

#Fig.S11A_bacateria_RMP
pred2<-data.frame(cbind(r_squared,r_squared_all))
colnames(pred2)<-c('selected','all')
pred2<-data.frame(pivot_longer(pred2,cols=1:2,names_to='type',values_to='R2'))
ggplot(pred2, aes(x=type, y=R2, fill=type)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Coefficient of determination", sep = ""))+
  theme_bw()



rm(list=ls())
setwd('E:/abs/Git/Figure7andS6-S12/fun/MS')

#model established and feature selected for fungi_QMP in Tissue
comm<-read.csv('MS_fun_genera_QMP.csv',row.names=1)
comm<-data.frame(t(comm))

Day<-read.csv('Group_MS.csv',row.names=1)

set.seed(1)
sub_train<-strata(Day,stratanames=("Day"),size=rep(4,7),method="srswor")
train_Day<-Day[sub_train$ID_unit,]
train_Data<-comm[sub_train$ID_unit,]

test_Day<-Day[-sub_train$ID_unit,]
test_Data<-comm[-sub_train$ID_unit,]


set.seed(1)
rf1 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE)  
print(rf1)

otu_train.cv <- replicate(5, rfcv(train_Data, train_Day, cv.fold=10,step=1.2), simplify = FALSE)
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)

#Fig.S7C
library(ggplot2)
p <- ggplot(otu_train.cv.mean, aes(Group.1, x)) +
    geom_line() +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
    labs(title = '',x = 'Number of ASVs', y = 'Cross-validation error')
p+geom_vline(xintercept = 12,color='red')


imp= as.data.frame(rf1$importance)
imp = imp[order(imp[,1],decreasing = T),]

train_select<-train_Data[,rownames(imp[1:12,])]
test_select<-test_Data[,rownames(imp[1:12,])]

n_train<-nrow(train_select)
n_test<-nrow(test_select)
predict_train<-matrix(NA,nrow=n_train,ncol=100)
predict_test<-matrix(NA,nrow=n_test,ncol=100)
r_squared<-1

for (i in 1:100) {
  set.seed(i)
  rf1 = randomForest(train_Day~.,data=train_select, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train[,i] <- predict(rf1,train_select)
  predict_test[,i] <- predict(rf1,test_select)
  r_squared[i]<-rf1$rsq[length(rf1$rsq)]
  i=i+1
}

mean(r_squared)

(mae_train<-mean(abs(rowMeans(predict_train)-train_Day)))
(sd_train<-sd(abs(rowMeans(predict_train)-train_Day)))

(mae_test<-mean(abs(rowMeans(predict_test)-test_Day)))
(sd_test<-sd(abs(rowMeans(predict_test)-test_Day)))

factor<-c(replicate(n_train,'train set'),replicate(n_test,'test_set'))
predict<-c(rowMeans(predict_train),rowMeans(predict_test))
True_value<-c(train_Day,test_Day)

df<-data.frame(cbind(factor,predict,True_value))

#Fig.7A_fungi_QMP
library(ggplot2)
mycolor<-c('#E64B35','#3C5488')
p <- ggplot(df, aes(as.numeric(True_value), as.numeric(predict), color=factor)) +
       geom_point(size=2.5,alpha=0.8) +
       scale_color_manual(values = mycolor)+ 
       theme_bw()+
       theme(panel.grid =element_blank())+
       xlab(paste("True value (days)", sep = "")) +
       ylab(paste("Predictive value (days)", sep = "")) +
       geom_abline(intercept = 0, slope = 1)

#resutlt of entire dataset
n_train<-nrow(train_Data)
n_test<-nrow(test_Data)
predict_train_all<-matrix(NA,nrow=n_train,ncol=100)
predict_test_all<-matrix(NA,nrow=n_test,ncol=100)
r_squared_all<-1

for (i in 1:100) {
  set.seed(i)
  rf2 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train_all[,i] <- predict(rf2,train_Data)
  predict_test_all[,i] <- predict(rf2,test_Data)
  r_squared_all[i]<-rf2$rsq[length(rf2$rsq)]
  i=i+1
}

mean(r_squared_all)

(mae_train_all<-mean(abs(rowMeans(predict_train_all)-train_Day)))
(sd_train_all<-sd(abs(rowMeans(predict_train_all)-train_Day)))

(mae_test_all<-mean(abs(rowMeans(predict_test_all)-test_Day)))
(sd_test_all<-sd(abs(rowMeans(predict_test_all)-test_Day)))

wilcox.test(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),alternative='two.sided',paired=T)
wilcox.test(colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day)),alternative='two.sided',paired=T)
wilcox.test(r_squared,r_squared_all,alternative='two.sided',paired=T)


pred<-data.frame(cbind(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),
                 colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day))))


colnames(pred)<-c('train_select','train_all','test_select','test_all')

#Fig.S10A_fungi_QMP
library(tidyr)
pred<-data.frame(pivot_longer(pred,cols=1:4,names_to='type',values_to='MAE'))
pred<-separate(pred,col='type',into=c('set','data'),sep='_')

ggplot(pred, aes(x=factor(set,levels=c('train','test')), y=MAE, fill=data)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Mean absolute error (days)", sep = ""))+
  theme_bw()

#Fig.S11A_fungi_QMP
pred2<-data.frame(cbind(r_squared,r_squared_all))
colnames(pred2)<-c('selected','all')
pred2<-data.frame(pivot_longer(pred2,cols=1:2,names_to='type',values_to='R2'))
ggplot(pred2, aes(x=type, y=R2, fill=type)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Coefficient of determination", sep = ""))+
  theme_bw()




rm(list=ls())
setwd('E:/abs/Git/Figure7andS6-S12/fun/MS')
#model established and feature selected for fungi_RMP in Tissue
comm<-read.csv('MS_fun_genera_RMP.csv',row.names=1)
comm<-data.frame(t(comm))

Day<-read.csv('Group_MS.csv',row.names=1)

set.seed(1)
sub_train<-strata(Day,stratanames=("Day"),size=rep(4,7),method="srswor")
train_Day<-Day[sub_train$ID_unit,]
train_Data<-comm[sub_train$ID_unit,]

test_Day<-Day[-sub_train$ID_unit,]
test_Data<-comm[-sub_train$ID_unit,]


set.seed(1)
rf1 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE)  
print(rf1)

otu_train.cv <- replicate(5, rfcv(train_Data, train_Day, cv.fold=10,step=1.2), simplify = FALSE)
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)

#Fig.S7D
library(ggplot2)
p <- ggplot(otu_train.cv.mean, aes(Group.1, x)) +
    geom_line() +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
    labs(title = '',x = 'Number of ASVs', y = 'Cross-validation error')
p+geom_vline(xintercept = 15,color='red')


imp= as.data.frame(rf1$importance)
imp = imp[order(imp[,1],decreasing = T),]

train_select<-train_Data[,rownames(imp[1:15,])]
test_select<-test_Data[,rownames(imp[1:15,])]

n_train<-nrow(train_select)
n_test<-nrow(test_select)
predict_train<-matrix(NA,nrow=n_train,ncol=100)
predict_test<-matrix(NA,nrow=n_test,ncol=100)
r_squared<-1

for (i in 1:100) {
  set.seed(i)
  rf1 = randomForest(train_Day~.,data=train_select, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train[,i] <- predict(rf1,train_select)
  predict_test[,i] <- predict(rf1,test_select)
  r_squared[i]<-rf1$rsq[length(rf1$rsq)]
  i=i+1
}

mean(r_squared)

(mae_train<-mean(abs(rowMeans(predict_train)-train_Day)))
(sd_train<-sd(abs(rowMeans(predict_train)-train_Day)))

(mae_test<-mean(abs(rowMeans(predict_test)-test_Day)))
(sd_test<-sd(abs(rowMeans(predict_test)-test_Day)))

factor<-c(replicate(n_train,'train set'),replicate(n_test,'test_set'))
predict<-c(rowMeans(predict_train),rowMeans(predict_test))
True_value<-c(train_Day,test_Day)

df<-data.frame(cbind(factor,predict,True_value))

#Fig.7A_fungi_RMP
library(ggplot2)
mycolor<-c('#E64B35','#3C5488')
p <- ggplot(df, aes(as.numeric(True_value), as.numeric(predict), color=factor)) +
       geom_point(size=2.5,alpha=0.8) +
       scale_color_manual(values = mycolor)+ 
       theme_bw()+
       theme(panel.grid =element_blank())+
       xlab(paste("True value (days)", sep = "")) +
       ylab(paste("Predictive value (days)", sep = "")) +
       geom_abline(intercept = 0, slope = 1)

#resutlt of entire dataset
n_train<-nrow(train_Data)
n_test<-nrow(test_Data)
predict_train_all<-matrix(NA,nrow=n_train,ncol=100)
predict_test_all<-matrix(NA,nrow=n_test,ncol=100)
r_squared_all<-1

for (i in 1:100) {
  set.seed(i)
  rf2 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train_all[,i] <- predict(rf2,train_Data)
  predict_test_all[,i] <- predict(rf2,test_Data)
  r_squared_all[i]<-rf2$rsq[length(rf2$rsq)]
  i=i+1
}

mean(r_squared_all)

(mae_train_all<-mean(abs(rowMeans(predict_train_all)-train_Day)))
(sd_train_all<-sd(abs(rowMeans(predict_train_all)-train_Day)))

(mae_test_all<-mean(abs(rowMeans(predict_test_all)-test_Day)))
(sd_test_all<-sd(abs(rowMeans(predict_test_all)-test_Day)))

wilcox.test(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),alternative='two.sided',paired=T)
wilcox.test(colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day)),alternative='two.sided',paired=T)
wilcox.test(r_squared,r_squared_all,alternative='two.sided',paired=T)


pred<-data.frame(cbind(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),
                 colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day))))


colnames(pred)<-c('train_select','train_all','test_select','test_all')

#Fig.S10A_fungi_RMP
library(tidyr)
pred<-data.frame(pivot_longer(pred,cols=1:4,names_to='type',values_to='MAE'))
pred<-separate(pred,col='type',into=c('set','data'),sep='_')

ggplot(pred, aes(x=factor(set,levels=c('train','test')), y=MAE, fill=data)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Mean absolute error (days)", sep = ""))+
  theme_bw()

#Fig.S11A_fungi_RMP
pred2<-data.frame(cbind(r_squared,r_squared_all))
colnames(pred2)<-c('selected','all')
pred2<-data.frame(pivot_longer(pred2,cols=1:2,names_to='type',values_to='R2'))
ggplot(pred2, aes(x=type, y=R2, fill=type)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Coefficient of determination", sep = ""))+
  theme_bw()


rm(list=ls())
setwd('E:/abs/Git/Figure7andS6-S12/metab')

#model established and feature selected for metabolites in Tissue
comm<-read.csv('metabolites_normal.csv',row.names=1)
Day<-read.csv('Group_MS.csv',row.names=1)

set.seed(1)
sub_train<-strata(Day,stratanames=("Day"),size=rep(4,7),method="srswor")
train_Day<-Day[sub_train$ID_unit,]
train_Data<-comm[sub_train$ID_unit,]

test_Day<-Day[-sub_train$ID_unit,]
test_Data<-comm[-sub_train$ID_unit,]

set.seed(1)
rf1 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE)  
print(rf1)

otu_train.cv <- replicate(5, rfcv(train_Data, train_Day, cv.fold=10,step=1.2), simplify = FALSE)
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)

#Fig.S8
library(ggplot2)
p <- ggplot(otu_train.cv.mean, aes(Group.1, x)) +
    geom_line() +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
    labs(title = '',x = 'Number of ASVs', y = 'Cross-validation error')
p+geom_vline(xintercept = 302,color='red')


imp= as.data.frame(rf1$importance)
imp = imp[order(imp[,1],decreasing = T),]

train_select<-train_Data[,rownames(imp[1:302,])]
test_select<-test_Data[,rownames(imp[1:302,])]

n_train<-nrow(train_select)
n_test<-nrow(test_select)
predict_train<-matrix(NA,nrow=n_train,ncol=100)
predict_test<-matrix(NA,nrow=n_test,ncol=100)
r_squared<-1

for (i in 1:100) {
  set.seed(i)
  rf1 = randomForest(train_Day~.,data=train_select, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train[,i] <- predict(rf1,train_select)
  predict_test[,i] <- predict(rf1,test_select)
  r_squared[i]<-rf1$rsq[length(rf1$rsq)]
  i=i+1
}

mean(r_squared)

(mae_train<-mean(abs(rowMeans(predict_train)-train_Day)))
(sd_train<-sd(abs(rowMeans(predict_train)-train_Day)))

(mae_test<-mean(abs(rowMeans(predict_test)-test_Day)))
(sd_test<-sd(abs(rowMeans(predict_test)-test_Day)))

factor<-c(replicate(n_train,'train set'),replicate(n_test,'test_set'))
predict<-c(rowMeans(predict_train),rowMeans(predict_test))
True_value<-c(train_Day,test_Day)

df<-data.frame(cbind(factor,predict,True_value))

#Fig.7B
library(ggplot2)
mycolor<-c('#E64B35','#3C5488')
p <- ggplot(df, aes(as.numeric(True_value), as.numeric(predict), color=factor)) +
       geom_point(size=2.5,alpha=0.8) +
       scale_color_manual(values = mycolor)+ 
       theme_bw()+
       theme(panel.grid =element_blank())+
       xlab(paste("True value (days)", sep = "")) +
       ylab(paste("Predictive value (days)", sep = "")) +
       geom_abline(intercept = 0, slope = 1)

#resutlt of entire dataset
n_train<-nrow(train_Data)
n_test<-nrow(test_Data)
predict_train_all<-matrix(NA,nrow=n_train,ncol=100)
predict_test_all<-matrix(NA,nrow=n_test,ncol=100)
r_squared_all<-1

for (i in 1:100) {
  set.seed(i)
  rf2 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train_all[,i] <- predict(rf2,train_Data)
  predict_test_all[,i] <- predict(rf2,test_Data)
  r_squared_all[i]<-rf2$rsq[length(rf2$rsq)]
  i=i+1
}

mean(r_squared_all)

(mae_train_all<-mean(abs(rowMeans(predict_train_all)-train_Day)))
(sd_train_all<-sd(abs(rowMeans(predict_train_all)-train_Day)))

(mae_test_all<-mean(abs(rowMeans(predict_test_all)-test_Day)))
(sd_test_all<-sd(abs(rowMeans(predict_test_all)-test_Day)))

wilcox.test(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),alternative='two.sided',paired=T)
wilcox.test(colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day)),alternative='two.sided',paired=T)
wilcox.test(r_squared,r_squared_all,alternative='two.sided',paired=T)


pred<-data.frame(cbind(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),
                 colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day))))


colnames(pred)<-c('train_select','train_all','test_select','test_all')

#Fig.S10C
library(tidyr)
pred<-data.frame(pivot_longer(pred,cols=1:4,names_to='type',values_to='MAE'))
pred<-separate(pred,col='type',into=c('set','data'),sep='_')

ggplot(pred, aes(x=factor(set,levels=c('train','test')), y=MAE, fill=data)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Mean absolute error (days)", sep = ""))+
  theme_bw()

#Fig.S11C
pred2<-data.frame(cbind(r_squared,r_squared_all))
colnames(pred2)<-c('selected','all')
pred2<-data.frame(pivot_longer(pred2,cols=1:2,names_to='type',values_to='R2'))
ggplot(pred2, aes(x=type, y=R2, fill=type)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Coefficient of determination", sep = ""))+
  theme_bw()



rm(list=ls())
setwd('E:/abs/Git/Figure7andS6-S12/multi')
#model established and feature selected for multi-omics in Tissue
comm<-read.csv('multi_MS.csv',row.names=1)
comm<-data.frame(t(comm))
Day<-read.csv('Group_MS.csv',row.names=1)

set.seed(1)
sub_train<-strata(Day,stratanames=("Day"),size=rep(4,7),method="srswor")
train_Day<-Day[sub_train$ID_unit,]
train_Data<-comm[sub_train$ID_unit,]

test_Day<-Day[-sub_train$ID_unit,]
test_Data<-comm[-sub_train$ID_unit,]


set.seed(1)
rf1 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE)  
print(rf1)

otu_train.cv <- replicate(5, rfcv(train_Data, train_Day, cv.fold=10,step=1.2), simplify = FALSE)
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)

#Fig.S9
library(ggplot2)
p <- ggplot(otu_train.cv.mean, aes(Group.1, x)) +
    geom_line() +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
    labs(title = '',x = 'Number of ASVs', y = 'Cross-validation error')
p+geom_vline(xintercept =58,color='red')


imp= as.data.frame(rf1$importance)
imp = imp[order(imp[,1],decreasing = T),]

train_select<-train_Data[,rownames(imp[1:58,])]
test_select<-test_Data[,rownames(imp[1:58,])]

n_train<-nrow(train_select)
n_test<-nrow(test_select)
predict_train<-matrix(NA,nrow=n_train,ncol=100)
predict_test<-matrix(NA,nrow=n_test,ncol=100)
r_squared<-1

for (i in 1:100) {
  set.seed(i)
  rf1 = randomForest(train_Day~.,data=train_select, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train[,i] <- predict(rf1,train_select)
  predict_test[,i] <- predict(rf1,test_select)
  r_squared[i]<-rf1$rsq[length(rf1$rsq)]
  i=i+1
}

mean(r_squared)

(mae_train<-mean(abs(rowMeans(predict_train)-train_Day)))
(sd_train<-sd(abs(rowMeans(predict_train)-train_Day)))

(mae_test<-mean(abs(rowMeans(predict_test)-test_Day)))
(sd_test<-sd(abs(rowMeans(predict_test)-test_Day)))

factor<-c(replicate(n_train,'train set'),replicate(n_test,'test_set'))
predict<-c(rowMeans(predict_train),rowMeans(predict_test))
True_value<-c(train_Day,test_Day)

df<-data.frame(cbind(factor,predict,True_value))

#Fig.7C
library(ggplot2)
mycolor<-c('#E64B35','#3C5488')
p <- ggplot(df, aes(as.numeric(True_value), as.numeric(predict), color=factor)) +
       geom_point(size=2.5,alpha=0.8) +
       scale_color_manual(values = mycolor)+ 
       theme_bw()+
       theme(panel.grid =element_blank())+
       xlab(paste("True value (days)", sep = "")) +
       ylab(paste("Predictive value (days)", sep = "")) +
       geom_abline(intercept = 0, slope = 1)


#resutlt of entire dataset
n_train<-nrow(train_Data)
n_test<-nrow(test_Data)
predict_train_all<-matrix(NA,nrow=n_train,ncol=100)
predict_test_all<-matrix(NA,nrow=n_test,ncol=100)
r_squared_all<-1

for (i in 1:100) {
  set.seed(i)
  rf2 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train_all[,i] <- predict(rf2,train_Data)
  predict_test_all[,i] <- predict(rf2,test_Data)
  r_squared_all[i]<-rf2$rsq[length(rf2$rsq)]
  i=i+1
}

mean(r_squared_all)

(mae_train_all<-mean(abs(rowMeans(predict_train_all)-train_Day)))
(sd_train_all<-sd(abs(rowMeans(predict_train_all)-train_Day)))

(mae_test_all<-mean(abs(rowMeans(predict_test_all)-test_Day)))
(sd_test_all<-sd(abs(rowMeans(predict_test_all)-test_Day)))

wilcox.test(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),alternative='two.sided',paired=T)
wilcox.test(colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day)),alternative='two.sided',paired=T)
wilcox.test(r_squared,r_squared_all,alternative='two.sided',paired=T)


pred<-data.frame(cbind(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),
                 colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day))))


colnames(pred)<-c('train_select','train_all','test_select','test_all')

#Fig.S10D
library(tidyr)
pred<-data.frame(pivot_longer(pred,cols=1:4,names_to='type',values_to='MAE'))
pred<-separate(pred,col='type',into=c('set','data'),sep='_')

ggplot(pred, aes(x=factor(set,levels=c('train','test')), y=MAE, fill=data)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Mean absolute error (days)", sep = ""))+
  theme_bw()

#Fig.S11D
pred2<-data.frame(cbind(r_squared,r_squared_all))
colnames(pred2)<-c('selected','all')
pred2<-data.frame(pivot_longer(pred2,cols=1:2,names_to='type',values_to='R2'))
ggplot(pred2, aes(x=type, y=R2, fill=type)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Coefficient of determination", sep = ""))+
  theme_bw()




rm(list=ls())
setwd('E:/abs/Git/Figure7andS6-S12/bac/GS')

#model established and feature selected for Bacteria_QMP in Grave soil
comm<-read.csv('AA_genera_GS.csv',row.names=1)
comm<-data.frame(t(comm))
Day<-read.csv('Group_GS.csv',row.names=1)

set.seed(1)
sub_train<-strata(Day,stratanames=("Day"),size=rep(4,7),method="srswor")
train_Day<-Day[sub_train$ID_unit,]
train_Data<-comm[sub_train$ID_unit,]

test_Day<-Day[-sub_train$ID_unit,]
test_Data<-comm[-sub_train$ID_unit,]


set.seed(1)
rf1 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE)  
print(rf1)

otu_train.cv <- replicate(5, rfcv(train_Data, train_Day, cv.fold=10,step=1.2), simplify = FALSE)
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)

#Fig.S6A
library(ggplot2)
p <- ggplot(otu_train.cv.mean, aes(Group.1, x)) +
    geom_line() +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
    labs(title = '',x = 'Number of ASVs', y = 'Cross-validation error')
p+geom_vline(xintercept =233,color='red')


imp= as.data.frame(rf1$importance)
imp = imp[order(imp[,1],decreasing = T),]

train_select<-train_Data[,rownames(imp[1:233,])]
test_select<-test_Data[,rownames(imp[1:233,])]

n_train<-nrow(train_select)
n_test<-nrow(test_select)
predict_train<-matrix(NA,nrow=n_train,ncol=100)
predict_test<-matrix(NA,nrow=n_test,ncol=100)
r_squared<-1

for (i in 1:100) {
  set.seed(i)
  rf1 = randomForest(train_Day~.,data=train_select, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train[,i] <- predict(rf1,train_select)
  predict_test[,i] <- predict(rf1,test_select)
  r_squared[i]<-rf1$rsq[length(rf1$rsq)]
  i=i+1
}

mean(r_squared)

(mae_train<-mean(abs(rowMeans(predict_train)-train_Day)))
(sd_train<-sd(abs(rowMeans(predict_train)-train_Day)))

(mae_test<-mean(abs(rowMeans(predict_test)-test_Day)))
(sd_test<-sd(abs(rowMeans(predict_test)-test_Day)))

factor<-c(replicate(n_train,'train set'),replicate(n_test,'test_set'))
predict<-c(rowMeans(predict_train),rowMeans(predict_test))
True_value<-c(train_Day,test_Day)

df<-data.frame(cbind(factor,predict,True_value))

#Fig.S12_bacteria_QMP
library(ggplot2)
mycolor<-c('#E64B35','#3C5488')
p <- ggplot(df, aes(as.numeric(True_value), as.numeric(predict), color=factor)) +
       geom_point(size=2.5,alpha=0.8) +
       scale_color_manual(values = mycolor)+ 
       theme_bw()+
       theme(panel.grid =element_blank())+
       xlab(paste("True value (days)", sep = "")) +
       ylab(paste("Predictive value (days)", sep = "")) +
       geom_abline(intercept = 0, slope = 1)

#resutlt of entire dataset
n_train<-nrow(train_Data)
n_test<-nrow(test_Data)
predict_train_all<-matrix(NA,nrow=n_train,ncol=100)
predict_test_all<-matrix(NA,nrow=n_test,ncol=100)
r_squared_all<-1

for (i in 1:100) {
  set.seed(i)
  rf2 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train_all[,i] <- predict(rf2,train_Data)
  predict_test_all[,i] <- predict(rf2,test_Data)
  r_squared_all[i]<-rf2$rsq[length(rf2$rsq)]
  i=i+1
}

mean(r_squared_all)

(mae_train_all<-mean(abs(rowMeans(predict_train_all)-train_Day)))
(sd_train_all<-sd(abs(rowMeans(predict_train_all)-train_Day)))

(mae_test_all<-mean(abs(rowMeans(predict_test_all)-test_Day)))
(sd_test_all<-sd(abs(rowMeans(predict_test_all)-test_Day)))

wilcox.test(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),alternative='two.sided',paired=T)
wilcox.test(colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day)),alternative='two.sided',paired=T)
wilcox.test(r_squared,r_squared_all,alternative='two.sided',paired=T)


pred<-data.frame(cbind(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),
                 colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day))))


colnames(pred)<-c('train_select','train_all','test_select','test_all')
library(tidyr)

#Fig.S10B_bacteria_QMP
pred<-data.frame(pivot_longer(pred,cols=1:4,names_to='type',values_to='MAE'))
pred<-separate(pred,col='type',into=c('set','data'),sep='_')

ggplot(pred, aes(x=factor(set,levels=c('train','test')), y=MAE, fill=data)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Mean absolute error (days)", sep = ""))+
  theme_bw()


#Fig.S11B_bacteria_QMP
pred2<-data.frame(cbind(r_squared,r_squared_all))
colnames(pred2)<-c('selected','all')
pred2<-data.frame(pivot_longer(pred2,cols=1:2,names_to='type',values_to='R2'))
ggplot(pred2, aes(x=type, y=R2, fill=type)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Coefficient of determination", sep = ""))+
  theme_bw()



rm(list=ls())
setwd('E:/abs/Git/Figure7andS6-S12/bac/GS')
#model established and feature selected for Bacteria_RMP in Grave soil
comm<-read.csv('RA_genera_GS.csv',row.names=1)
comm<-data.frame(t(comm))
Day<-read.csv('Group_GS.csv',row.names=1)

set.seed(1)
sub_train<-strata(Day,stratanames=("Day"),size=rep(4,7),method="srswor")
train_Day<-Day[sub_train$ID_unit,]
train_Data<-comm[sub_train$ID_unit,]

test_Day<-Day[-sub_train$ID_unit,]
test_Data<-comm[-sub_train$ID_unit,]


set.seed(1)
rf1 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE)  
print(rf1)

otu_train.cv <- replicate(5, rfcv(train_Data, train_Day, cv.fold=10,step=1.2), simplify = FALSE)
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)

#Fig.S6B
library(ggplot2)
p <- ggplot(otu_train.cv.mean, aes(Group.1, x)) +
    geom_line() +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
    labs(title = '',x = 'Number of ASVs', y = 'Cross-validation error')
p+geom_vline(xintercept =66,color='red')


imp= as.data.frame(rf1$importance)
imp = imp[order(imp[,1],decreasing = T),]

train_select<-train_Data[,rownames(imp[1:66,])]
test_select<-test_Data[,rownames(imp[1:66,])]

n_train<-nrow(train_select)
n_test<-nrow(test_select)
predict_train<-matrix(NA,nrow=n_train,ncol=100)
predict_test<-matrix(NA,nrow=n_test,ncol=100)
r_squared<-1

for (i in 1:100) {
  set.seed(i)
  rf1 = randomForest(train_Day~.,data=train_select, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train[,i] <- predict(rf1,train_select)
  predict_test[,i] <- predict(rf1,test_select)
  r_squared[i]<-rf1$rsq[length(rf1$rsq)]
  i=i+1
}

mean(r_squared)

(mae_train<-mean(abs(rowMeans(predict_train)-train_Day)))
(sd_train<-sd(abs(rowMeans(predict_train)-train_Day)))

(mae_test<-mean(abs(rowMeans(predict_test)-test_Day)))
(sd_test<-sd(abs(rowMeans(predict_test)-test_Day)))

factor<-c(replicate(n_train,'train set'),replicate(n_test,'test_set'))
predict<-c(rowMeans(predict_train),rowMeans(predict_test))
True_value<-c(train_Day,test_Day)

df<-data.frame(cbind(factor,predict,True_value))

#Fig.S12_bacteria_RMP
library(ggplot2)
mycolor<-c('#E64B35','#3C5488')
p <- ggplot(df, aes(as.numeric(True_value), as.numeric(predict), color=factor)) +
       geom_point(size=2.5,alpha=0.8) +
       scale_color_manual(values = mycolor)+ 
       theme_bw()+
       theme(panel.grid =element_blank())+
       xlab(paste("True value (days)", sep = "")) +
       ylab(paste("Predictive value (days)", sep = "")) +
       geom_abline(intercept = 0, slope = 1)

#resutlt of entire dataset
n_train<-nrow(train_Data)
n_test<-nrow(test_Data)
predict_train_all<-matrix(NA,nrow=n_train,ncol=100)
predict_test_all<-matrix(NA,nrow=n_test,ncol=100)
r_squared_all<-1

for (i in 1:100) {
  set.seed(i)
  rf2 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train_all[,i] <- predict(rf2,train_Data)
  predict_test_all[,i] <- predict(rf2,test_Data)
  r_squared_all[i]<-rf2$rsq[length(rf2$rsq)]
  i=i+1
}

mean(r_squared_all)

(mae_train_all<-mean(abs(rowMeans(predict_train_all)-train_Day)))
(sd_train_all<-sd(abs(rowMeans(predict_train_all)-train_Day)))

(mae_test_all<-mean(abs(rowMeans(predict_test_all)-test_Day)))
(sd_test_all<-sd(abs(rowMeans(predict_test_all)-test_Day)))

wilcox.test(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),alternative='two.sided',paired=T)
wilcox.test(colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day)),alternative='two.sided',paired=T)
wilcox.test(r_squared,r_squared_all,alternative='two.sided',paired=T)


pred<-data.frame(cbind(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),
                 colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day))))


colnames(pred)<-c('train_select','train_all','test_select','test_all')
library(tidyr)

#Fig.S10B_bacteria_RMP
pred<-data.frame(pivot_longer(pred,cols=1:4,names_to='type',values_to='MAE'))
pred<-separate(pred,col='type',into=c('set','data'),sep='_')

ggplot(pred, aes(x=factor(set,levels=c('train','test')), y=MAE, fill=data)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Mean absolute error (days)", sep = ""))+
  theme_bw()


#Fig.S11B_bacteria_RMP
pred2<-data.frame(cbind(r_squared,r_squared_all))
colnames(pred2)<-c('selected','all')
pred2<-data.frame(pivot_longer(pred2,cols=1:2,names_to='type',values_to='R2'))
ggplot(pred2, aes(x=type, y=R2, fill=type)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Coefficient of determination", sep = ""))+
  theme_bw()


rm(list=ls())
setwd('E:/abs/Git/Figure7andS6-S12/fun/GS')

#model established and feature selected for fungi_QMP in Grave soil
comm<-read.csv('GS_fun_genera_QMP.csv',row.names=1)
comm<-data.frame(t(comm))
Day<-read.csv('Group_GS.csv',row.names=1)

set.seed(1)
sub_train<-strata(Day,stratanames=("Day"),size=rep(4,7),method="srswor")
train_Day<-Day[sub_train$ID_unit,]
train_Data<-comm[sub_train$ID_unit,]

test_Day<-Day[-sub_train$ID_unit,]
test_Data<-comm[-sub_train$ID_unit,]


set.seed(1)
rf1 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE)  
print(rf1)

otu_train.cv <- replicate(5, rfcv(train_Data, train_Day, cv.fold=10,step=1.2), simplify = FALSE)
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)

#Fig.S7A
library(ggplot2)
p <- ggplot(otu_train.cv.mean, aes(Group.1, x)) +
    geom_line() +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
    labs(title = '',x = 'Number of ASVs', y = 'Cross-validation error')
p+geom_vline(xintercept =28,color='red')


imp= as.data.frame(rf1$importance)
imp = imp[order(imp[,1],decreasing = T),]

train_select<-train_Data[,rownames(imp[1:28,])]
test_select<-test_Data[,rownames(imp[1:28,])]

n_train<-nrow(train_select)
n_test<-nrow(test_select)
predict_train<-matrix(NA,nrow=n_train,ncol=100)
predict_test<-matrix(NA,nrow=n_test,ncol=100)
r_squared<-1

for (i in 1:100) {
  set.seed(i)
  rf1 = randomForest(train_Day~.,data=train_select, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train[,i] <- predict(rf1,train_select)
  predict_test[,i] <- predict(rf1,test_select)
  r_squared[i]<-rf1$rsq[length(rf1$rsq)]
  i=i+1
}

mean(r_squared)

(mae_train<-mean(abs(rowMeans(predict_train)-train_Day)))
(sd_train<-sd(abs(rowMeans(predict_train)-train_Day)))

(mae_test<-mean(abs(rowMeans(predict_test)-test_Day)))
(sd_test<-sd(abs(rowMeans(predict_test)-test_Day)))

factor<-c(replicate(n_train,'train set'),replicate(n_test,'test_set'))
predict<-c(rowMeans(predict_train),rowMeans(predict_test))
True_value<-c(train_Day,test_Day)

df<-data.frame(cbind(factor,predict,True_value))

#Fig.S12_fungi_QMP
library(ggplot2)
mycolor<-c('#E64B35','#3C5488')
p <- ggplot(df, aes(as.numeric(True_value), as.numeric(predict), color=factor)) +
       geom_point(size=2.5,alpha=0.8) +
       scale_color_manual(values = mycolor)+ 
       theme_bw()+
       theme(panel.grid =element_blank())+
       xlab(paste("True value (days)", sep = "")) +
       ylab(paste("Predictive value (days)", sep = "")) +
       geom_abline(intercept = 0, slope = 1)

#resutlt of entire dataset
n_train<-nrow(train_Data)
n_test<-nrow(test_Data)
predict_train_all<-matrix(NA,nrow=n_train,ncol=100)
predict_test_all<-matrix(NA,nrow=n_test,ncol=100)
r_squared_all<-1

for (i in 1:100) {
  set.seed(i)
  rf2 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train_all[,i] <- predict(rf2,train_Data)
  predict_test_all[,i] <- predict(rf2,test_Data)
  r_squared_all[i]<-rf2$rsq[length(rf2$rsq)]
  i=i+1
}

mean(r_squared_all)

(mae_train_all<-mean(abs(rowMeans(predict_train_all)-train_Day)))
(sd_train_all<-sd(abs(rowMeans(predict_train_all)-train_Day)))

(mae_test_all<-mean(abs(rowMeans(predict_test_all)-test_Day)))
(sd_test_all<-sd(abs(rowMeans(predict_test_all)-test_Day)))

wilcox.test(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),alternative='two.sided',paired=T)
wilcox.test(colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day)),alternative='two.sided',paired=T)
wilcox.test(r_squared,r_squared_all,alternative='two.sided',paired=T)


pred<-data.frame(cbind(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),
                 colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day))))


colnames(pred)<-c('train_select','train_all','test_select','test_all')
library(tidyr)

#Fig.S10B_fungi_QMP
pred<-data.frame(pivot_longer(pred,cols=1:4,names_to='type',values_to='MAE'))
pred<-separate(pred,col='type',into=c('set','data'),sep='_')

ggplot(pred, aes(x=factor(set,levels=c('train','test')), y=MAE, fill=data)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Mean absolute error (days)", sep = ""))+
  theme_bw()


#Fig.S11B_fungi_QMP
pred2<-data.frame(cbind(r_squared,r_squared_all))
colnames(pred2)<-c('selected','all')
pred2<-data.frame(pivot_longer(pred2,cols=1:2,names_to='type',values_to='R2'))
ggplot(pred2, aes(x=type, y=R2, fill=type)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Coefficient of determination", sep = ""))+
  theme_bw()



rm(list=ls())
setwd('E:/abs/Git/Figure7andS6-S12/fun/GS')

#model established and feature selected for fungi_RMP in Grave soil
comm<-read.csv('GS_fun_genera_RMP.csv',row.names=1)
comm<-data.frame(t(comm))
Day<-read.csv('Group_GS.csv',row.names=1)

set.seed(1)
sub_train<-strata(Day,stratanames=("Day"),size=rep(4,7),method="srswor")
train_Day<-Day[sub_train$ID_unit,]
train_Data<-comm[sub_train$ID_unit,]

test_Day<-Day[-sub_train$ID_unit,]
test_Data<-comm[-sub_train$ID_unit,]


set.seed(1)
rf1 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE)  
print(rf1)

otu_train.cv <- replicate(5, rfcv(train_Data, train_Day, cv.fold=10,step=1.2), simplify = FALSE)
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)

#Fig.S7B
library(ggplot2)
p <- ggplot(otu_train.cv.mean, aes(Group.1, x)) +
    geom_line() +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
    labs(title = '',x = 'Number of ASVs', y = 'Cross-validation error')
p+geom_vline(xintercept =8,color='red')


imp= as.data.frame(rf1$importance)
imp = imp[order(imp[,1],decreasing = T),]

train_select<-train_Data[,rownames(imp[1:8,])]
test_select<-test_Data[,rownames(imp[1:8,])]

n_train<-nrow(train_select)
n_test<-nrow(test_select)
predict_train<-matrix(NA,nrow=n_train,ncol=100)
predict_test<-matrix(NA,nrow=n_test,ncol=100)
r_squared<-1

for (i in 1:100) {
  set.seed(i)
  rf1 = randomForest(train_Day~.,data=train_select, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train[,i] <- predict(rf1,train_select)
  predict_test[,i] <- predict(rf1,test_select)
  r_squared[i]<-rf1$rsq[length(rf1$rsq)]
  i=i+1
}

mean(r_squared)

(mae_train<-mean(abs(rowMeans(predict_train)-train_Day)))
(sd_train<-sd(abs(rowMeans(predict_train)-train_Day)))

(mae_test<-mean(abs(rowMeans(predict_test)-test_Day)))
(sd_test<-sd(abs(rowMeans(predict_test)-test_Day)))

factor<-c(replicate(n_train,'train set'),replicate(n_test,'test_set'))
predict<-c(rowMeans(predict_train),rowMeans(predict_test))
True_value<-c(train_Day,test_Day)

df<-data.frame(cbind(factor,predict,True_value))

#Fig.S12_fungi_RMP
library(ggplot2)
mycolor<-c('#E64B35','#3C5488')
p <- ggplot(df, aes(as.numeric(True_value), as.numeric(predict), color=factor)) +
       geom_point(size=2.5,alpha=0.8) +
       scale_color_manual(values = mycolor)+ 
       theme_bw()+
       theme(panel.grid =element_blank())+
       xlab(paste("True value (days)", sep = "")) +
       ylab(paste("Predictive value (days)", sep = "")) +
       geom_abline(intercept = 0, slope = 1)

#resutlt of entire dataset
n_train<-nrow(train_Data)
n_test<-nrow(test_Data)
predict_train_all<-matrix(NA,nrow=n_train,ncol=100)
predict_test_all<-matrix(NA,nrow=n_test,ncol=100)
r_squared_all<-1

for (i in 1:100) {
  set.seed(i)
  rf2 = randomForest(train_Day~.,data=train_Data, ntree=1000, proximity=TRUE, importance=TRUE) 
  predict_train_all[,i] <- predict(rf2,train_Data)
  predict_test_all[,i] <- predict(rf2,test_Data)
  r_squared_all[i]<-rf2$rsq[length(rf2$rsq)]
  i=i+1
}

mean(r_squared_all)

(mae_train_all<-mean(abs(rowMeans(predict_train_all)-train_Day)))
(sd_train_all<-sd(abs(rowMeans(predict_train_all)-train_Day)))

(mae_test_all<-mean(abs(rowMeans(predict_test_all)-test_Day)))
(sd_test_all<-sd(abs(rowMeans(predict_test_all)-test_Day)))

wilcox.test(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),alternative='two.sided',paired=T)
wilcox.test(colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day)),alternative='two.sided',paired=T)
wilcox.test(r_squared,r_squared_all,alternative='two.sided',paired=T)


pred<-data.frame(cbind(colMeans(abs(predict_train-train_Day)),colMeans(abs(predict_train_all-train_Day)),
                 colMeans(abs(predict_test-test_Day)),colMeans(abs(predict_test_all-test_Day))))


colnames(pred)<-c('train_select','train_all','test_select','test_all')
library(tidyr)

#Fig.S10B_fungi_RMP
pred<-data.frame(pivot_longer(pred,cols=1:4,names_to='type',values_to='MAE'))
pred<-separate(pred,col='type',into=c('set','data'),sep='_')

ggplot(pred, aes(x=factor(set,levels=c('train','test')), y=MAE, fill=data)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Mean absolute error (days)", sep = ""))+
  theme_bw()


#Fig.S11B_fungi_RMP
pred2<-data.frame(cbind(r_squared,r_squared_all))
colnames(pred2)<-c('selected','all')
pred2<-data.frame(pivot_longer(pred2,cols=1:2,names_to='type',values_to='R2'))
ggplot(pred2, aes(x=type, y=R2, fill=type)) +
  geom_boxplot()+
  xlab(paste(''))+
  ylab(paste("Coefficient of determination", sep = ""))+
  theme_bw()





