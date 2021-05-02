library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
pFilter=0.05 
gene="Risk score" 
immuneFile="./cibersoft/CIBERSORTx_Job1_Results_metabric.csv"
expFile="./cox/metabricRisk.txt" 
immune=data.table::fread(immuneFile,data.table = F)
rownames(immune) <- immune[,1]
immune <- immune[,-1]
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

rt=read.table(expFile,sep="\t",header=T,check.names=F,row.names = 1)
exp <- rt[,-c(1:19)]
sameSample=intersect(row.names(immune),row.names(exp))
immune1=immune[sameSample,]
exp1=exp[sameSample,]
colnames(exp1)[1]=c('Risk score')

#相关性检验
outTab=data.frame()
x=as.numeric(exp1[,1])
#按免疫细胞循环
for(j in colnames(immune1)[1:22]){
  y=as.numeric(immune1[,j])
  if(sd(y)>0.001){
    df1=as.data.frame(cbind(x,y))
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pValue=corT$p.value
    p1=ggplot(df1, aes(x, y)) + ggtitle("METABRIC") +
      theme(plot.title = element_text(size=18,  family="serif",face = "bold",hjust = 0.5)) +
      ylab(j)+xlab(gene)+
      geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
      stat_cor(method = 'spearman', aes(x =x, y =y))
    if(pValue<pFilter){
      pdf(file=paste0(j,".pdf"),width=5,height=5)
      print(p1)
      dev.off()
      outTab=rbind(outTab,cbind(Cell=j,pValue))
    }
  }
}
write.table(outTab,file="./cor/METABRIC/immuneCor.result.txt",sep="\t",row.names=F,quote=F)
