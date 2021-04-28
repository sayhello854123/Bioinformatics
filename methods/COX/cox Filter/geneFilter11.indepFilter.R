

#install.packages('survival')

library(survival)
setwd("D:\\biowolf\\80geneFilter\\11.indepFilter")                     #修改工作目录
expFile="surSigExp.txt"                                                            #表达数据文件
clinicalFile="clinical.txt"                                                        #临床数据文件

exp=read.table(expFile,sep="\t",header=T,check.names=F,row.names=1)                #读取表达数据文件
cli=read.table(clinicalFile,sep="\t",header=T,check.names=F,row.names=1)           #读取临床数据文件
samSample=intersect(row.names(exp),row.names(cli))
exp=exp[samSample,]
cli=cli[samSample,]
pFilter=0.05   

#独立预后分析，输出表格
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(exp[,3:ncol(exp)])){
    rt=cbind(exp[,1:2],cli,gene=exp[,i])
    multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
    coxSummary=summary(multiCox)
	  if(coxSummary$coefficients["gene","Pr(>|z|)"]<pFilter){
	      sigGenes=c(sigGenes,i)
	      outTab=rbind(outTab,
	               cbind(id=i,
	               HR=coxSummary$conf.int["gene","exp(coef)"],
	               HR.95L=coxSummary$conf.int["gene","lower .95"],
	               HR.95H=coxSummary$conf.int["gene","upper .95"],
	               pvalue=coxSummary$coefficients["gene","Pr(>|z|)"])
	      )
	  }
}
write.table(outTab,file="indep.xls",sep="\t",row.names=F,quote=F)
indepSigExp=exp[,sigGenes]
indepSigExp=cbind(id=row.names(indepSigExp),indepSigExp)
write.table(indepSigExp,file="indepSigExp.txt",sep="\t",row.names=F,quote=F)

