##对测试集进行筛选###
pFilter=0.05                           #定义单因素显著性###直接定义0.005可获取8个miRNA的结果

library(survival)                                                  #引用包
rt=read.table("train.txt",header=T,sep="\t",check.names=F,row.names=1)       #读取输入文件
###筛选32个miRNA与患者的总生存期相关
sigGenes=c("futime","fustat")
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
  sigGenes=c(sigGenes,i)
  outTab=rbind(outTab,
                     cbind(id=i,
                           HR=coxSummary$conf.int[,"exp(coef)"],
                           HR.95L=coxSummary$conf.int[,"lower .95"],
                           HR.95H=coxSummary$conf.int[,"upper .95"],
                           pvalue=coxSummary$coefficients[,"Pr(>|z|)"])) 
  }                                     
}
###这里为了能清楚的演示我保留了32个单因素的结果，其实可以直接定义0.005的结果
write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)

####选择了8个miRNA###的
rt=read.table("uniSigExp.txt",header=T,sep="\t",check.names=F,row.names=1) 
pFilter=0.005 
sigGenes=c("futime","fustat")
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])) 
  }                                     
}
write.table(outTab,file="uniCox1.xls",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp1.txt",sep="\t",row.names=F,quote=F)


####Kaplan-Meier方法8个挑选3个与生存相关，一个被证明与事实不相符合##
pFilter=0.05

rt=read.table("uniSigExp1.txt",header=T,sep="\t",check.names=F,row.names=1) 
sigGenes=c(sigGenes,gene)
for(gene in colnames(rt[,3:ncol(rt)])){
  if(sd(rt[,gene])<0.1){
    next}
  a=rt[,gene]<=median(rt[,gene])
  #km方法
  diff=survdiff(Surv(futime, fustat) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  fit=survfit(Surv(futime, fustat) ~ a, data = rt)
  #cox方法
  #cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  #coxSummary = summary(cox)
  #coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(pValue<pFilter){
    sigGenes=c(sigGenes,gene)
    if(pValue<0.001){
      pValue=signif(pValue,4)
      pValue=format(pValue, scientific = TRUE)
    }else{
      pValue=round(pValue,3)
    }
    pdf(file=paste(gene,".survival.pdf",sep=""),
        width = 5.5,            #图片的宽度
        height =5,              #图片的高度	      
    )
    plot(fit, 
         lwd=2,
         col=c("red","blue"),
         xlab="Time (year)",
         mark.time=T,
         ylab="Survival rate",
         ylim=c(0,1.09),
         main=paste(gene,"(p=", pValue ,")",sep="") )
    legend("topright", 
           c(paste(gene," high expression",sep=""), 
             paste(gene," low expression",sep="") ), 
           lwd=2, 
           col=c("red","blue"))
    dev.off()
  }
}
###删除hsa-miR-485-5p##
rt <- rt[,-7]##hsa-miR-485-5p在第7列
uniSigExp=cbind(id=row.names(rt),rt)
write.table(uniSigExp,file="mucox_input.txt",sep="\t",row.names=F,quote=F)



####多元cox分析
library(survival)
library(survminer)

rt=read.table("mucox_input.txt",header=T,sep="\t",check.names=F,row.names=1)  #读取train输入文件


#使用train组构建COX模型
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

#输出模型相关信息
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



#输出train组风险文件
outTab <- as.matrix(outTab)
rownames(outTab) <- outTab[,1]
sigGenes=c("futime","fustat")
sigGenes1=rt[,sigGenes]
trainFinalGeneExp=rt[,rownames(outTab)]
same <- intersect(row.names(sigGenes1),row.names(trainFinalGeneExp))
trainFinal=cbind(sigGenes1[same,],trainFinalGeneExp[same,])
myFun=function(x){crossprod(as.numeric(x),multiCoxSum$coefficients[,"coef"])}
trainScore=apply(trainFinalGeneExp,1,myFun)
risktrian=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(trainFinal,riskScore=as.vector(trainScore),risktrian)
out=rbind(ID=colnames(outTab),outTab)
write.table(out,file="trainRisk.txt",sep="\t",quote=F,col.names=F)###生成的数据最后一列名因子错误，我不想改了自行打开添加risk的标题


#输出test组风险文件
Test=read.table("test.txt",header=T,sep="\t",check.names=F,row.names=1)          #读取test输入文件
TestFinalGeneExp=Test[,names(trainFinalGeneExp)]
testScore=apply(TestFinalGeneExp,1,myFun)     #利用train得到模型预测test样品风险
riskTest=as.vector(ifelse(testScore>median(testScore),"high","low"))
sigGenes2=Test[,sigGenes]
same1<- intersect(row.names(sigGenes2),row.names(TestFinalGeneExp))
testFinal=cbind(sigGenes2[same1,],TestFinalGeneExp[same1,])
outTab1=cbind(testFinal,riskScore=as.vector(testScore),riskTest)
out1=rbind(ID=colnames(outTab1),outTab1)
write.table(out1,file="testRisk.txt",sep="\t",quote=F,col.names=F)###生成的数据最后一列名因子错误，我不想改了自行打开添加risk的标题

###生成的trainRisk和testRisk 最后一行标题名是NA自己手动改一下

###在train组和test组以及整个组中预测生存的5个miRNA
library(survival)

#绘制train组生存曲线
rt1=read.table("trainRisk.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt1)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt1)
summary(fit)    #查看五年生存率
pdf(file="survivalTrain.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival curve (p=", pValue ,")",sep=""),
     mark.time=T)
legend("topright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()

#绘制test组生存曲线
rt2=read.table("testRisk.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt2)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt2)
summary(fit)    #查看五年生存率
pdf(file="survivalTest.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival curve (p=", pValue ,")",sep=""),
     mark.time=T)
legend("topright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()

rt <- rbind(rt1,rt2)
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
summary(fit)    #查看五年生存率
pdf(file="survivalall.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival curve (p=", pValue ,")",sep=""),
     mark.time=T)
legend("topright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()

#####ROC####
library(survivalROC)
time=3    ####预测时间可以修改成1年 或3年或5年 ###此处文章为3年
rt1=read.table("trainRisk.txt",header=T,sep="\t",check.names=F,row.names=1)
pdf(file="rocTrain.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc1=survivalROC(Stime=rt1$futime, status=rt1$fustat, marker = rt1$riskScore, 
                predict.time =time, method="KM")
plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(roc1$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

rt2=read.table("testRisk.txt",header=T,sep="\t",check.names=F,row.names=1)
pdf(file="rocTest.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc2=survivalROC(Stime=rt2$futime, status=rt2$fustat, marker = rt2$riskScore, 
                predict.time =time, method="KM")
plot(roc2$FP, roc2$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(roc2$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

rt <- rbind(rt1,rt2)
pdf(file="rocall.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                predict.time =time, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

###survStat###
#绘制train组生存状态图
rt1=read.table("trainRisk.txt",header=T,sep="\t",check.names=F,row.names=1)       #读取train输入文件
rt1=rt1[order(rt1$riskScore),]
riskClass=rt1[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
color=as.vector(rt1$fustat)
color[color==1]="red"
color[color==0]="green"
pdf(file="survStatTrain.pdf",width = 12,height = 5)
plot(rt1$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
dev.off()

#绘制test组生存状态图
rt2=read.table("testRisk.txt",header=T,sep="\t",check.names=F,row.names=1)
rt2=rt2[order(rt2$riskScore),]
riskClass=rt2[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
color=as.vector(rt2$fustat)
color[color==1]="red"
color[color==0]="green"
pdf(file="survStatTest.pdf",width = 12,height = 5)
plot(rt2$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
dev.off()


#绘制全部组生存状态图
rt <- rbind(rt1,rt2)
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
color=as.vector(rt$fustat)
color[color==1]="red"
color[color==0]="green"
pdf(file="survStatALL.pdf",width = 12,height = 5)
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
dev.off()
