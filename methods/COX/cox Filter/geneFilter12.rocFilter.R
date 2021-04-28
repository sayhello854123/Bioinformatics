

#install.packages("survivalROC")

library(survivalROC)

setwd("D:\\biowolf\\80geneFilter\\12.RocFilter")                      #工作目录（需修改）
rocFilter=0.6                                                                     #ROC过滤值
rt=read.table("indepSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)      #读取输入文件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
	   roc=survivalROC(Stime=rt$futime, 
	                   status=rt$fustat, 
	                   marker = rt[,i], 
	                   predict.time =1, 
	                   method="KM")
	   if(roc$AUC>rocFilter){
	       sigGenes=c(sigGenes,i)
	       outTab=rbind(outTab,cbind(gene=i,AUC=roc$AUC))
	   }
}
write.table(outTab,file="ROC.xls",sep="\t",row.names=F,quote=F)    #输出基因和p值表格文件
rocSigExp=rt[,sigGenes]
rocSigExp=cbind(id=row.names(rocSigExp),rocSigExp)
write.table(rocSigExp,file="rocSigExp.txt",sep="\t",row.names=F,quote=F)

