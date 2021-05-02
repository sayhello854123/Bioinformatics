library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(patchwork)
pFilter=0.05 
gene="Risk score" 
immuneFile="./cibersoft/CIBERSORTx_Job1_Results.txt"
expFile="./cox/tcgaRisk.txt" 
immune=read.table(immuneFile,sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])
immune <- immune[,-c(1:3)]
immune <- immune[,-c(2:10)]
immune <- immune[,-c(5:10)]
rt=read.table(expFile,sep="\t",header=T,check.names=F,row.names = 1)
rownames(rt) = substr(rownames(rt),1,12)
exp <- rt[,-c(1:19)]
sameSample=intersect(row.names(immune),row.names(exp))
immune1=immune[sameSample,]
exp1=exp[sameSample,]
colnames(exp1)[1]=c('Risk score')

#相关性检验
#outTab=data.frame()
x=as.numeric(exp1[,1])
y=as.numeric(immune1[,1])
if(sd(y)>0.001){
  df1=as.data.frame(cbind(x,y))
  corT=cor.test(x,y,method="spearman")
  cor=corT$estimate
  pValue=corT$p.value}
  p1=ggplot(df1, aes(x, y)) + ggtitle("TCGA") +
    theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
    ylab(colnames(immune1)[1])+xlab(gene)+
    geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
    stat_cor(method = 'spearman', aes(x =x, y =y))
###
  y=as.numeric(immune1[,2])
  if(sd(y)>0.001){
    df1=as.data.frame(cbind(x,y))
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pValue=corT$p.value}
  p2=ggplot(df1, aes(x, y)) + ggtitle("TCGA") +
    theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
    ylab(colnames(immune1)[2])+xlab(gene)+
    geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
    stat_cor(method = 'spearman', aes(x =x, y =y))
###
  y=as.numeric(immune1[,3])
  if(sd(y)>0.001){
    df1=as.data.frame(cbind(x,y))
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pValue=corT$p.value}
  p3=ggplot(df1, aes(x, y)) + ggtitle("TCGA") +
    theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
    ylab(colnames(immune1)[3])+xlab(gene)+
    geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
    stat_cor(method = 'spearman', aes(x =x, y =y))
  ###
  y=as.numeric(immune1[,4])
  if(sd(y)>0.001){
    df1=as.data.frame(cbind(x,y))
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pValue=corT$p.value}
  p4=ggplot(df1, aes(x, y)) + ggtitle("TCGA") +
    theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
    ylab(colnames(immune1)[4])+xlab(gene)+
    geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
    stat_cor(method = 'spearman', aes(x =x, y =y))

###METABRIC
  pFilter=0.05 
  gene="Risk score" 
  immuneFile="./cibersoft/CIBERSORTx_Job1_Results_metabric.csv"
  expFile="./cox/metabricRisk.txt" 
  immune=data.table::fread(immuneFile,data.table = F)
  rownames(immune) <- immune[,1]
  immune <- immune[,-1]
  immune=immune[immune[,"P-value"]<pFilter,]
  immune=as.matrix(immune[,1:(ncol(immune)-3)])
  immune <- immune[,-c(1:3)]
  immune <- immune[,-c(2:10)]
  immune <- immune[,-c(5:10)] 
  rt=read.table(expFile,sep="\t",header=T,check.names=F,row.names = 1)
  exp <- rt[,-c(1:19)]
  sameSample=intersect(row.names(immune),row.names(exp))
  immune1=immune[sameSample,]
  exp1=exp[sameSample,]
  colnames(exp1)[1]=c('Risk score')
##### 
   x=as.numeric(exp1[,1])
  y=as.numeric(immune1[,1])
  if(sd(y)>0.001){
    df1=as.data.frame(cbind(x,y))
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pValue=corT$p.value}
  p5=ggplot(df1, aes(x, y)) + ggtitle("METABRIC") +
    theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
    ylab(colnames(immune1)[1])+xlab(gene)+
    geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
    stat_cor(method = 'spearman', aes(x =x, y =y))
  ####
  y=as.numeric(immune1[,2])
  if(sd(y)>0.001){
    df1=as.data.frame(cbind(x,y))
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pValue=corT$p.value}
  p6=ggplot(df1, aes(x, y)) + ggtitle("METABRIC") +
    theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
    ylab(colnames(immune1)[2])+xlab(gene)+
    geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
    stat_cor(method = 'spearman', aes(x =x, y =y))
  ####
  y=as.numeric(immune1[,3])
  if(sd(y)>0.001){
    df1=as.data.frame(cbind(x,y))
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pValue=corT$p.value}
  p7=ggplot(df1, aes(x, y)) + ggtitle("METABRIC") +
    theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
    ylab(colnames(immune1)[3])+xlab(gene)+
    geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
    stat_cor(method = 'spearman', aes(x =x, y =y))
  ####
  y=as.numeric(immune1[,4])
  if(sd(y)>0.001){
    df1=as.data.frame(cbind(x,y))
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pValue=corT$p.value}
  p8=ggplot(df1, aes(x, y)) + ggtitle("METABRIC") +
    theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
    ylab(colnames(immune1)[4])+xlab(gene)+
    geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
    stat_cor(method = 'spearman', aes(x =x, y =y))
  P <- (p2|p3)/(p4|p1)/(p6|p7)/(p8|p5)
  ggplot2::ggsave(filename = './cor/COR.pdf',plot = P, width =10, height =16)  
  