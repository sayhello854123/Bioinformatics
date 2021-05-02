library(survival)
library(survminer)
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(face = "bold",hjust=0.5)
    )
}
inputFile='./cox/tcgaRisk.txt'
inputFile='./cox/metabricRisk.txt'
#读取输入文件
rt=read.table(inputFile,header=T,sep="\t")
#比较高低风险组生存差异，得到显著性p值
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%0.3f",pValue))
  }
fit1 <- survfit(Surv(futime, fustat) ~ risk, data = rt)
fit2 <- survfit(Surv(futime, fustat) ~ risk, data = rt)
a <- coxph(Surv(futime, fustat) ~risk,data = rt)
tcga = summary(a)
tcga$conf.int[,"exp(coef)"]
tcga$conf.int[,"lower .95"]
tcga$conf.int[,"upper .95"]
surPlot <- list()
#绘制生存曲线
surPlot[[1]]=ggsurvplot(fit1, 
                     data=rt,
                     title = "Training Survival Curves \nTCGA",ggtheme=custom_theme(),
                     pval=pValue,
                     pval.size=6,
                     risk.table=TRUE,
                     legend = c(0.83,0.82),
                     legend.labs=c("High risk", "Low risk"), risk.table.col = "strata",
                     legend.title="Risk",
                     xlab="Time(years)",
                     break.time.by = 3,
                     risk.table.title="",
                     palette=c("red", "blue"),
                     risk.table.height=.25)

surPlot[[2]]=ggsurvplot(fit2, 
                        data=rt,
                        title = "Validation Survival Curves\n METABRIC",
                        ggtheme=custom_theme(),
                        pval=pValue,
                        pval.size=6,
                        risk.table=TRUE,
                        legend = c(0.83,0.82),
                        legend.labs=c("High risk", "Low risk"), risk.table.col = "strata",
                        legend.title="Risk",
                        xlab="Time(years)",
                        break.time.by = 3,
                        risk.table.title="",
                        palette=c("red", "blue"),
                        risk.table.height=.25)
pdf(file='./cox/survival.pdf',onefile = FALSE,width = 11,height =5)
arrange_ggsurvplots(surPlot, print = TRUE,
                    ncol = 2, nrow = 1, risk.table.height = 0.3)
dev.off()
