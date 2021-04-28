rt <- read.table('LUAD miRNAmatrix.txt',sep = '\t',header = T,check.names = F,row.names = 1)
group=sapply(strsplit(colnames(rt),"-"),'[',4)### strsplit 函数是为了将字符串按照固定的格式分开
A <- ifelse(group=='01A','tumor','NORMAL')
table(A)
luad_tumor <- rt[,A=='tumor']
luad_tumor <- t(luad_tumor)
str(luad_tumor)
rt1 <- read.table('LUADGene.txt',sep = '\t',header = T,check.names = F,row.names = 1)
same <- intersect(row.names(luad_tumor),row.names(rt1))
survival_input <- rt1[same,]
rownames(survival_input)<- substr(rownames(survival_input),1,12)
time <- read.table('time.txt',sep = '\t',header = T,check.names = F,row.names = 1)
same_id <- intersect(row.names(survival_input),row.names(time))
data=cbind(time[same_id,],survival_input[same_id,])
data <- data[,-4]
input_data<-na.omit(data)
out=rbind(ID=colnames(input_data),input_data)
write.table(out,file="sur_input.txt",sep="\t",quote=F,col.names=F)

library(survival)
pFilter=0.05 #显著性过滤条件

rt=read.table("sur_input.txt",header=T,sep="\t",check.names=F,row.names=1)     #读取输入文件
rt$futime=rt$futime/365     
gene="hsa-mir-625"
a=rt[,gene]<=median(rt[,gene])
diff=survdiff(Surv(futime, fustat) ~a,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
}else{
  pValue=round(pValue,3)
}

fit <- survfit(Surv(futime, fustat) ~ a, data = rt)

pdf(file="survival.pdf",
    width=6,
    height=6)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     mark.time=T,
     ylab="Survival rate",
     main=paste(gene,"(p=", pValue ,")",sep="") )
legend("topright", 
       c("High expression","Low expression"), 
       lwd=2, 
       col=c("red","blue"))
dev.off()
summary(fit)
library(survminer)
surPlot=ggsurvplot(fit,  size = 1,  # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palettes
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           xlim = c(0,20),
           xlab = "Time in years",
           break.time.by = 2, 
           #axes.offset=F,
           title = "hsa-mir-625",
           legend.title ="express",
           risk.table.y.text.col = T, # Risk table color by groups
           legend.labs = c("high","low"), # Change legend labels
           risk.table.height = 0.25,# Useful to change when you have multiple groups
           conf.int.style = "step",  # customize style of confidence intervals
           surv.median.line = "hv",  # add the median survival pointer.
           ggtheme = theme_light()  # Change ggplot2 theme
)
pdf(file="ggsurvival.pdf",width = 6.5,height =5.5)
print(surPlot)
dev.off()
