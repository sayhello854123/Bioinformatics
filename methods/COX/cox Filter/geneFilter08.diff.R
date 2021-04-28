

setwd("D:\\biowolf\\80geneFilter\\08.diffFilter")            #设置工作目录
inputFile="uniq.symbol.txt"                                              #输入文件
pFilter=0.05                                                      #p值临界值
logFCfilter=1                                                     #logFC临界值
conNum=32                                                         #normal组样品数目
treatNum=375                                                      #tumor组样品数目

#读取输入文件
outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))
data=read.table(inputFile,sep="\t",header=T,check.names=F,row.names=1)
data=as.matrix(data)

#差异分析
for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  pvalue=wilcoxTest$p.value
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
	if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		  outTab=rbind(outTab，cbind(gene=i, conMean=conGeneMeans, treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
	 }
}

#输出所有基因的差异情况
write.table(outTab,file="all.xls",sep="\t",row.names=F,quote=F)

#输出差异表格
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$pValue))<pFilter),]
write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)

#输出差异基因表达
heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="diffGeneExp.txt",sep="\t",col.names=F,quote=F)

