library(ggplot2)
library(survminer)
rt <- read.table('CIBERSORT-Results.txt',header = T,sep = '\t',check.names = F,row.names = 1)
group <- ifelse(as.numeric(substr(rownames(rt),14,15)) < 10,'tumor','normal')
table(group)
normal <- rt[group=='normal',]
tumor <- rt[group=='tumor',]
Data=rbind(cbind(normal,type="normal"),cbind(tumor,type="tumor"))
rt1 <- rbind(cell=Data[,1],type=Data[,length(Data)])
rt1=as.data.frame(t(rt1))
rt1$subtype <- ifelse(rt1$type==1,'normal','tumor')
rt1 <- rt1[,-2]
p <- ggplot(rt1,aes(subtype,cell))+geom_jitter(position=position_jitter(0.1) ,cex=3,aes(colour = factor(subtype)))
p+scale_color_manual(values = c("#6A5ACD", "#CD3333"))+theme_classic()+ theme(legend.position="none")+geom_boxplot()

 
a <- ggplot(rt1, aes(x=subtype, y=cell,color=subtype)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.1),cex=2.5,aes(colour = factor(subtype)))+
 theme_classic()+theme(legend.position="none")
  a+scale_color_manual(values=c('#436EEE', '#CD0000'))+#stat_compare_means(method = "anova", label.y = 0.4)+
  labs(x="subtype", y = "CIBERSORT fractions",title="B cells naive")

