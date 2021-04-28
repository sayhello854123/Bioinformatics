rt <- read.table('minyi ID.txt',header=T,sep="\t",row.names=1,check.names=F)
rt1 <- read.table('tcgaImmuneExp.share.txt',header=T,sep="\t",row.names=1,check.names=F)
same <- intersect(row.names(rt),row.names(rt1))
data <- rt1[same,]
time <- read.table('time.txt',header=T,sep="\t",row.names=1,check.names=F)
data1 <- t(data)
rownames(data1)<- substr(rownames(data1),1,12)
sameID <- intersect(row.names(data1),row.names(time))
input_data <- cbind(time[sameID,],data1[sameID,])
out=rbind(ID=colnames(input_data),input_data)
write.table(out,file="surtest_input.txt",sep="\t",quote=F,col.names=F)


###获取21个免疫基因的风险值
library(survival)                                                
rt=read.table("surtest_input.txt",header=T,sep="\t",check.names=F,row.names=1)   
rt$futime=rt$futime/365


multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)


outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)


riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
            file="risk.txt",
            sep="\t",
            quote=F,
            row.names=F)
###多重ROC比较##
rt=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)  
rt$IRGrisk=rt$riskScore
rt1 <- read.table("tcgaRisk.txt",header=T,sep="\t",check.names=F,row.names=1)
rt1$IRGPIrisk=rt1$riskScore
same <- intersect(row.names(rt),row.names(rt1))
rt <- rt[,-c(3:4)]
rt1 <- rt1[,-c(3:18)]
roc_input <- cbind(rt[same,],rt1[same,])
roc_input <- roc_input[,-c(4:5)]
out=rbind(ID=colnames(roc_input),roc_input)
write.table(out,file="roc_input.txt",sep="\t",quote=F,col.names=F)


library(survivalROC)
rt=read.table("roc_input.txt",header=T,sep="\t",check.names=F,row.names=1)    
rocCol=rainbow(ncol(rt)-2)
aucText=c()

pdf(file="multiROC.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$IRGPIrisk, predict.time =5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("IRGPIrisk"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)


j=1
for(i in colnames(rt[,3:(ncol(rt)-1)])){
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =5, method="KM")
  j=j+1
  aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
}
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()
