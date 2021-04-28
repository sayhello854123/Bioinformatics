

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library("limma")                                                      #引用limma包
setwd("D:\\biowolf\\80geneFilter\\07.uniq")               #设置工作目录
rt=read.table("symbol.txt",sep="\t",header=T,check.names=F)           #读取文件

#对数据进行处理
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp)，colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]

#输出uniq基因文件
uniq=rbind(ID=colnames(data),data)
write.table(uniq,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)        #输出文件


