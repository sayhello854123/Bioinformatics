###选取你们做的差异基因的表达文件#
##这个矩阵的列名为TCGA的样本##行名为你们的差异基因
##获取49样本中的肿瘤样本####
rt <- read.table('49_diffmRNAExp.txt',sep = '\t',header = T,check.names = F,row.names = 1)
group_list=ifelse(as.numeric(substr(colnames(rt ),14,15)) < 02,'tumor','normal')####substr 挑选字符串按所在特定位置结果
lusc_tumor<- rt[,group_list == "tumor" ]###选出所有的的肿瘤样本
out=rbind(ID=colnames(lusc_tumor),lusc_tumor)
write.table(out,file="lusc_tumor_49.txt",sep="\t",quote=F,col.names=F)

###制作临床数据###
##https://blog.csdn.net/sayhello1025/article/details/103474816###是获取下面输入文件的教程######
pheo <- read.table('clinical.txt',sep = '\t',header = T,check.names = F)
rownames(pheo) <- pheo$bcr_patient_barcode

##生存时间##
pheo[, 9][is.na(pheo[,9])] = 0
pheo[, 10][is.na(pheo[, 10])] = 0
pheo$futime <- as.numeric(pheo[, 9]) + as.numeric(pheo[, 10])

###生存状态##
pheo$status <- ifelse(pheo$vital_status == "Alive", 0, 1)###0是活着1是死了

##age##
###下面是为了验证一下是否相等
a <- pheo$age_at_initial_pathologic_diagnosis
b <- floor(abs(pheo$days_to_birth/365) ) ##绝对值向下取整
c <- a==b
table(c)
##感觉差不多可以直接用
pheo$age_at_initial_pathologic_diagnosis[is.na(pheo$age_at_initial_pathologic_diagnosis)] <- '0'
pheo$age <- as.numeric(pheo$age_at_initial_pathologic_diagnosis) 
## stage
library(stringr) 
pheo$stage <- str_split_fixed(pheo$stage_event, 'T', 2)[,1]
##后面的T,M,N分期由于无法定位具体的位置我手动分的

###临床数据##
clinical <- pheo[,c("futime","status","age","gender","stage")]
out=rbind(ID=colnames(clinical),clinical)
write.table(out,file="clinical.txt",sep="\t",quote=F,col.names=F)
time_input <- clinical[,-c(3:5)]
out=rbind(ID=colnames(time_input),time_input)
write.table(out,file="time.txt",sep="\t",quote=F,col.names=F)



##合并临床数据##
rt <- read.table('lusc_tumor_49.txt',sep = '\t',header = T,check.names = F,row.names = 1)
colnames(rt) <- substr(colnames(rt),1,12)
rt <- t(rt)
rt1 <- read.table('time.txt',sep = '\t',header = T,check.names = F,row.names = 1)
same=intersect(row.names(rt),row.names(rt1))
sur_input <- cbind(rt1[same,],rt[same,])
out=rbind(ID=colnames(sur_input),sur_input)
write.table(out,file="sur_input.txt",sep="\t",quote=F,col.names=F)

library(survival)
library(survminer)
pFilter=0.05                                                #由于图形太多，只对p小于0.05的基因绘图
rt=read.table("sur_input.txt",header=T,sep="\t",check.names=F,row.names = 1)    #读取输入文件
rt$futime=rt$futime/365                                     #如果以月为单位，除以30；以年为单位，除以365
outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  if(sd(rt[,gene])<0.1){
    next}
  a=rt[,gene]<=median(rt[,gene])
  #km方法
  diff=survdiff(Surv(futime, fustat) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  fit=survfit(Surv(futime, fustat) ~ a, data = rt)
  #cox方法
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if((pValue<pFilter) & (coxP<pFilter)){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab, 
                 cbind(gene=gene,
                       KM=pValue,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       coxPvalue=coxP) )
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
write.table(outTab,file="survival.xls",sep="\t",row.names=F,quote=F)    #输出基因和p值表格文件
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="surSigExp1.txt",sep="\t",row.names=F,quote=F)
