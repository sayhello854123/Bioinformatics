library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
set.seed(12345)
durg_plot <- function(x,y,z,a){
  rt = read.table(x, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  data=data[rowMeans(data)>0.5,]
  
  #删掉正常样品
  group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
  group=sapply(strsplit(group,""), "[", 1)
  group=gsub("2","1",group)
  data=data[,group==0]
  dim(data)
  
  train <- read.table(y, header=T, sep="\t", check.names=F,row.names = 1)
  test <-  read.table(z, header=T, sep="\t", check.names=F,row.names = 1)
  drug <-  read.table(a)
  drug <- drug$V1
  picDir='./drug_plot'
  dir.create(picDir)
  
  for(i in drug){
    senstivity=pRRopheticPredict(data, i, selection=1)
    senstivity=senstivity[senstivity!="NaN"]
    #读取风险输入文件
    risk=rbind(train,test)
    
    #风险文件和药物敏感性合并
    sameSample=intersect(row.names(risk), names(senstivity))
    risk=risk[sameSample, "group",drop=F]
    senstivity=senstivity[sameSample]
    rt=cbind(risk, senstivity)
    rt$risk=factor(rt$group, levels=c("Low", "High"))
    type=levels(factor(rt[,"risk"]))
    comp=combn(type, 2)
    my_comparisons=list()
    for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
    
    #绘制箱线图
    boxplot=ggboxplot(rt, x="risk", y="senstivity", fill="risk",
                      xlab="Risk",
                      ylab=paste0(i, " senstivity (IC50)"),
                      legend.title="Risk",
                      palette=c("green", "red")
    )+ 
      stat_compare_means(comparisons=my_comparisons)
    filename <- paste0(i, ".pdf",sep='')
    outfile <- paste(picDir,filename,sep="/")
    pdf(file=outfile, width=5, height=4.5)
    print(boxplot)
    dev.off()
  }
}
x <- './durg/03protein_coding_combat_exp.txt'
y <- './durg/training.txt'
z <- './durg/validation.txt'
a <- './durg/durg.txt'
durg_plot(x,y,z,a)