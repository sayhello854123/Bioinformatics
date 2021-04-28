pFilter=0.05  
rt=data.table::fread("./cibersoft/CIBERSORTx_Job1_Results_metabric.csv ",data.table = F)
rownames(rt) <- rt[,1]
rt <- rt[,-1]
data=rt[rt[,"P-value"]<pFilter,]
data=data[,1:(ncol(rt)-3)]
risk=read.table("./cox/metabricRisk.txt",header=T,sep="\t",row.names=1,check.names=F)
sameSample=intersect(row.names(data),row.names(risk))
data=data[sameSample,]
risk=risk[sameSample,]
riskHigh=risk[risk$risk=="high",]
riskLow=risk[risk$risk=="low",]
dataHigh=data[row.names(riskHigh),]
dataLow=data[row.names(riskLow),]
newData=rbind(cbind(dataHigh,risk="High-risk"),cbind(dataLow,risk="Low-risk"))

risk="risk"
dotCol=c("red","blue")
outTab=data.frame()
for(cell in colnames(newData[,1:(ncol(newData)-1)])){
  rt1=rbind(cell=newData[,cell],risk=newData[,risk])
  rt1=as.matrix(t(rt1))
  wilcoxTest=wilcox.test(cell ~ risk, data=rt1)
  pValue=wilcoxTest$p.value
  sig=ifelse(pValue<0.001,"***",ifelse(pValue<0.01,"**",ifelse(pValue<0.05,"*"," ")))
  outTab=rbind(outTab,cbind(cell=cell,pValue=pValue,sig))
  pval=0
  if(pValue<0.001){
    pval=signif(pValue,4)
    pval=format(pval, scientific = TRUE)
  }else{
    pval=sprintf("%0.3f",pValue)
  }
  
  if(pValue<0.05){	  
    b = boxplot(cell ~ risk, data = rt1,outline = FALSE, plot=F) 
    yMin=min(b$stats)
    yMax = max(b$stats/5+b$stats)
    ySeg = max(b$stats/10+b$stats)
    ySeg2 = max(b$stats/12+b$stats)
    n = ncol(b$stats)
    
    pdf(file=paste0(cell,".pdf"),width=6,height=5)
    par(mar = c(4,7,3,3))
    boxplot(cell ~ risk, data = rt1, ylab = cell,col=dotCol,xlab="",names=c("High risk","Low risk"),
            cex.main=1.2, cex.lab=1, cex.axis=1, ylim=c(yMin,yMax), outline = FALSE)
    segments(1,ySeg, n,ySeg);segments(1,ySeg, 1,ySeg2);segments(n,ySeg, n,ySeg2)
    text((1+n)/2,ySeg,labels=paste0("p=",pval),cex=1,pos=3)
    dev.off()
  }
}
write.table(outTab,file="./cibersoft/METABRIC/corStat.txt",sep="\t",row.names=F,quote=F)

#输出雷达图需要的输入文件
fmsb=rbind(colMeans(dataHigh),colMeans(dataLow))
row.names(fmsb)=c("High","Low")
write.table(cbind(id=rownames(fmsb),fmsb),file="./cibersoft/METABRIC/fmsbInput.txt",sep="\t",quote=F,row.names=F)

library(fmsb) 

data=read.table("./cibersoft/METABRIC/fmsbInput.txt",header=T,sep="\t",row.names=1,check.names=F)   #读取输入文件
data=rbind(rep(max(data),22),rep(0,22),data)
#定义图形顺序
sortCellNames= c("B cells memory",
                 "B cells naive",              
                 "Dendritic cells activated",
                 "Dendritic cells resting",     
                 "Eosinophils",
                 "Macrophages M0",              
                 "Macrophages M1",
                 "Macrophages M2",              
                 "Mast cells activated",
                 "NK cells resting",            
                 "Plasma cells",
                 "Mast cells resting", 
                 "T cells CD4 memory activated",
                 "T cells CD4 naive", 
                 "T cells CD4 memory resting",        
                 "Monocytes",
                 "Neutrophils",                 
                 "NK cells activated",           
                 "T cells CD8","T cells follicular helper",   
                 "T cells gamma delta",
                 "T cells regulatory (Tregs)")
data=data[,sortCellNames]
#定义颜色
colors=c("red","blue")
#定义显著性
corStat=read.table("./cibersoft/METABRIC/corStat.txt",header=T,sep="\t",row.names=1,check.names=F)
corStat=corStat[sortCellNames,]
colnames(data)=paste0(colnames(data),corStat$sig)

#输出结果
pdf(file="./cibersoft/METABRIC/radar.pdf",height=8,width=8)
radarchart( data  , axistype=1 , 
            pcol=colors,                 #设置颜色
            plwd=2 ,                     #线条粗线
            plty=1,                      #虚线，实线
            cglcol="grey",               #背景线条颜色
            cglty=1,                     #背景线条虚线，实线 
            caxislabels=seq(0,1,0.05),   #坐标刻度
            cglwd=0.8,                   #背景线条粗细
            axislabcol="grey",           #刻度颜色
            vlcex=0.75                   #字体大小
)
legend("topright",legend=rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors ,cex=1.2, pt.cex=2)
dev.off()