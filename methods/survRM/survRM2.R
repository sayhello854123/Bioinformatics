
library(survRM2)
require(survival)
library(survRM2)
library(Hmisc)
library(survcomp)
library(prodlim)

setwd("C:\\Users\\Administrator\\Desktop\\survRM")
data <- read.table("test1.txt",header = T,row.names = 1)


compareRiskModal=function(os,ev,risks,max_time){

  tau <- max_time 
  n.grid <- 100 

  OS <-Surv(os,ev)
  marker.pp <- seq(from=0,to=1,length=n.grid)
  mdat=cbind()
  c.all=list()
  c.indexs=c()
  
  for(i in 1:ncol(risks)){
    marker <-as.numeric(risks[,i])
    marker.qq<-quantile(marker,marker.pp)
    fitdat.df<-data.frame(marker=marker)
    newdat.df<-data.frame(marker=marker.qq)
    
    cox.model1<-coxph(OS~marker,data=fitdat.df)
    rms.calc <-summary(survfit(cox.model1,newdata=newdat.df),rmean=tau)
    rms.mean <-rms.calc$table[,"*rmean"]
    mdat=cbind(mdat,rms.mean)
    c1=concordance.index(x=as.numeric(risks[,i]), surv.time=os, surv.event=ev,
                      method="noether")
    c.all=c(c.all,list(c1))
    c.indexs=c(c.indexs,paste0(colnames(risks)[i],':',round(c1$c.index,2),',95%CI(',round(c1$lower,2),'-',round(c1$upper,2),'),p=',signif(c1$p.value,2)))
  }
  
  mt=ceiling(max(mdat))
  
  cls=rainbow(ncol(risks))
  if(ncol(risks)>1){
    plot(marker.pp,as.numeric(mdat[,1]),pch=20,col = cls[1], cex=1.3,ylim=c(0,mt),xlab="",ylab="")
    if(ncol(risks)>2){
      for(i in 2:(ncol(risks)-1)){
        par(new=TRUE)
        plot(marker.pp,as.numeric(mdat[,i]),pch=20,col = cls[i], cex=1.3,ylim=c(0,mt),xlab="",ylab="")  
      }
    }
    par(new=TRUE)
    plot(marker.pp,as.numeric(mdat[,ncol(risks)]),pch=20,col = cls[ncol(risks)],
         cex=1.3,ylim=c(0,mt),xlab="Percentile of Score",ylab="RMS(months)")
  }else{
    plot(marker.pp,as.numeric(mdat[,1]),pch=20,col = "red", cex=1.3,ylim=c(0,mt),xlab="Percentile of Score",ylab="RMS(months)")
  }
  trimws(c.indexs, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  legend("topright",c.indexs, pch=20,col=cls, title= "C-index", inset = .05,cex = 1.2)
  
  if(ncol(risks)>1){
    all.cp=c()
  for(i in 1:(ncol(risks)-1)){
    c1=concordance.index(x=risks[,i], surv.time=os, surv.event=ev,
                         method="noether")
    for(j in (i+1):ncol(risks)){
      c2=concordance.index(x=risks[,j], surv.time=os, surv.event=ev,
                           method="noether")
      p=min(cindex.comp(c1, c2)$p.value,cindex.comp(c2, c1)$p.value)
      all.cp=c(all.cp,paste0(colnames(risks)[i],'-vs-',colnames(risks)[j],':p=',signif(p,2)))
    }
    
    trimws(all.cp, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
    legend("bottomleft",all.cp, pch=20,title= "Compare C-index", inset = .05,cex = 1.2)
  }
  }
}

pdf('surv.pdf',width = 10,height = 7)
compareRiskModal(data$futime,
                 data$fustat,
                 data.frame(ICPI=data$ICPIriskScore,PAIR=data$PAIRriskScore),
                 60)
dev.off()

