

#install.packages("survival")

library(survival)

setwd("D:\\biowolf\\80geneFilter\\10.survivalFilter")            #工作目录（需修改）
pFilter=0.05                                                                 #显著性过滤条件
rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1)     #读取输入文件
rt$futime=rt$futime/365                                                      #如果以月为单位，除以30；以年为单位，除以365

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
	   if(sd(rt[,gene])<0.1){
	     next}
	   a=rt[,gene]<=median(rt[,gene])
	  
		 rt1=rt[a,]
	   b=setdiff(rownames(rt),rownames(rt1))
	   rt2=rt[b,]
	   surTab1=summary(survfit(Surv(futime, fustat) ~ 1, data = rt1))
	   surTab2=summary(survfit(Surv(futime, fustat) ~ 1, data = rt2))
	   survivalTab1=cbind(time=surTab1$time, surv=surTab1$surv,lower=surTab1$lower,upper=surTab1$upper)
	   survivalTab1=survivalTab1[survivalTab1[,"time"]<5,]
	   if(class(survivalTab1)=="matrix"){
	     survivalTab1=survivalTab1[nrow(survivalTab1),]
	   }
	   survivalTab2=cbind(time=surTab2$time, surv=surTab2$surv,lower=surTab2$lower,upper=surTab2$upper)
	   survivalTab2=survivalTab2[survivalTab2[,"time"]<5,]
	   if(class(survivalTab2)=="matrix"){
	     survivalTab2=survivalTab2[nrow(survivalTab2),]
	   }
	   fiveYearsDiff=abs(survivalTab1["surv"]-survivalTab2["surv"])

     #km方法
	   diff=survdiff(Surv(futime, fustat) ~a,data = rt)
	   pValue=1-pchisq(diff$chisq,df=1)
	   fit=survfit(Surv(futime, fustat) ~ a, data = rt)
	   #cox方法
	   cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
	   coxSummary = summary(cox)
	   coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	
	   if((pValue<pFilter) & (coxP<pFilter) & (fiveYearsDiff>0.15)){
	       sigGenes=c(sigGenes，gene)
	       outTab=rbind(outTab, 
	                    cbind(gene=gene,
	                          KM=pValue,
	                          HR=coxSummary$conf.int[,"exp(coef)"],
	                          HR.95L=coxSummary$conf.int[,"lower .95"],
	                          HR.95H=coxSummary$conf.int[,"upper .95"],
			                      coxPvalue=coxP) )
		 }
}
write.table(outTab,file="survival.xls",sep="\t",row.names=F,quote=F)    #输出基因和p值表格文件
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="surSigExp.txt",sep="\t",row.names=F,quote=F)

