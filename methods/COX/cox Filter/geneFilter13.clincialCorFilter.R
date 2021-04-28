

setwd("D:\\biowolf\\80geneFilter\\13.clinicalCorFilter")                #修改工作目录
expFile="rocSigExp.txt"                                                             #表达数据文件
clinicalFile="clinical.txt"                                                         #临床数据文件

exp=read.table(expFile,sep="\t",header=T,check.names=F,row.names=1)                 #读取表达数据文件
cli=read.table(clinicalFile,sep="\t",header=T,check.names=F,row.names=1)            #读取临床数据文件
exp=exp[,3:ncol(exp)]
samSample=intersect(row.names(exp),row.names(cli))
exp=exp[samSample,]
cli=cli[samSample,]
pFilter=0.05   

#临床相关性分析，输出表格
outTab=c()
outTab=rbind(outTab,c("id",colnames(cli),"SigNum"))
colnames(outTab)=c("id",colnames(cli),"SigNum")
for(i in colnames(exp)){
		clinicalPvalVector=c()
		sigSum=0
		for(clinical in colnames(cli)){
		    tab1=table(cli[,clinical])
		    labelNum=length(tab1)
			  rt1=cbind(expression=exp[,i],clinical=cli[,clinical])
			  if(labelNum==2){
			      cliTest<-wilcox.test(expression ~ clinical, data=rt1)
			  }else{
			      cliTest<-kruskal.test(expression ~ clinical, data = rt1)}
			  pValue=cliTest$p.value
			  clinicalPvalVector=c(clinicalPvalVector,pValue)
			  if(pValue<pFilter){
			     sigSum=sigSum+1
			  }
		}
		geneClinical=c(i,clinicalPvalVector,sigSum)
		outTab=rbind(outTab,geneClinical)
}
write.table(outTab,file="clinicalCor.xls",sep="\t",col.names=F,row.names=F,quote=F)

