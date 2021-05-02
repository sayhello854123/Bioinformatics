rm(list=ls())
options(stringsAsFactors = F)

library(survival)
risk=read.table('./cox/tcgaRisk.txt',sep = '\t',header = T,check.names = F,row.names = 1)
clinical=read.table('./bulk data/01.data/cli_tcga.txt',sep = '\t',header = T,check.names = F,row.names = 1)
sameSample=intersect(row.names(clinical),row.names(risk))
risk=risk[sameSample,]
clinical=clinical[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],clinical,riskScore=risk[,(ncol(risk)-1)])

#单因素独立预后分析
uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="./cox/unicox/tcga.uniCox.txt",sep="\t",row.names=F,quote=F)

#多因素独立预后分析
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="./cox/unicox/tcga.multiCox.txt",sep="\t",row.names=F,quote=F)

risk=read.table('./cox/metabricRisk.txt',sep = '\t',header = T,check.names = F,row.names = 1)
cli1 <- read.table('./bulk data/brca_metabric/123.txt',sep = '\t',header = T,check.names = F,row.names = 1)
cli2 <- read.table('./bulk data/brca_metabric/234.txt',sep = '\t',header = T,check.names = F,row.names = 1)
same <- intersect(row.names(cli1),row.names(cli2))
clinical <- cbind(cli1[same,],cli2[same,])[,-c(1:2)]
sameSample=intersect(row.names(clinical),row.names(risk))
risk=risk[sameSample,]
clinical=clinical[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],clinical,riskScore=risk[,(ncol(risk)-1)])

uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="./cox/unicox/metabric.uniCox.txt",sep="\t",row.names=F,quote=F)

#多因素独立预后分析
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="./cox/unicox/metabric.multiCox.txt",sep="\t",row.names=F,quote=F)


bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  pdf(file=forestFile, width = 6,height = 4.3)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
}
bioForest(coxFile="./cox/unicox/tcga.uniCox.txt",forestFile="./cox/unicox/tcga.uniForest.pdf",forestCol="green")
bioForest(coxFile="./cox/unicox/tcga.multiCox.txt",forestFile="./cox/unicox/tcga.multiForest.pdf",forestCol="red")
bioForest(coxFile="./cox/unicox/metabric.uniCox.txt",forestFile="./cox/unicox/metabric.uniForest.pdf",forestCol="green")
bioForest(coxFile="./cox/unicox/metabric.multiCox.txt",forestFile="./cox/unicox/metabric.multiForest.pdf",forestCol="red")
