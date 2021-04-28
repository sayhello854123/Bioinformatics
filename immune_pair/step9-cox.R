

###选取生存相关的gene pair
rm(list = ls())
load("./cox/cox_input.RData")
library(survival)
library(glmnet)
library(survivalROC)
pFilter=0.01
rt$futime=rt$futime/365
outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    diff=survdiff(Surv(futime, fustat) ~rt[,gene],data = rt)
    pValue=1-pchisq(diff$chisq,df=1)
    if(pValue<pFilter){
      sigGenes=c(sigGenes,gene)
      outTab=rbind(outTab,
                   cbind(gene=gene,
                         #KM=pValue,
                         HR=coxSummary$conf.int[,"exp(coef)"],
                         HR.95L=coxSummary$conf.int[,"lower .95"],
                         HR.95H=coxSummary$conf.int[,"upper .95"],
                         coxPvalue=coxP) )
    }
  }
}
rt=rt[,sigGenes]

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 10000)
cvfit=cv.glmnet(x, y, family="cox", maxit = 1000)
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]

trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore))

predictTime=3 #1年的ROC曲线，需要做3年或5年改成相应的数值
roc=survivalROC(Stime=outTab$futime, status=outTab$fustat, marker = outTab$riskScore, predict.time =predictTime, method="KM")
sum=roc$TP-roc$FP
cutOp=roc$cut.values[which.max(sum)]
cutTP=roc$TP[which.max(sum)]
cutFP=roc$FP[which.max(sum)]
pdf(file="./cox/ROC.pdf",width=6,height=6)
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1.1),col="black", 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.2, cex.lab=1.2, cex.axis=1.2,font=1.2)
points(cutFP,cutTP, pch=20, col="red",cex=1.5)
text(cutFP+0.1,cutTP-0.05,paste0("Cutoff:",sprintf("%0.3f",cutOp)))
dev.off()

risk=as.vector(ifelse(trainScore>cutOp,"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="./cox/tcgaRisk.txt",sep="\t",quote=F,row.names=F)

load("./cox/metabric_input.RData")
rt$futime=rt$futime/365
rt[,3:ncol(rt)][rt[,3:ncol(rt)]<0]=0
testFinalGeneExp=rt[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(testScore>cutOp,"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="./cox/metabricRisk.txt",sep="\t",quote=F,row.names=F)
