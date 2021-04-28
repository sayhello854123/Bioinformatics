rt=read.table("symbol.txt",sep="\t",header=T,check.names=F,row.names = 1)
rownames(rt) <- rt[,1]
qw=read.table("cli.txt",sep="\t",header=T,check.names=F)  
library(TCGAbiolinks)
subtypes <- PanCancerAtlas_subtypes()
a <- ifelse(subtypes$cancer.type == "BRCA",'brca','other')
b <- subtypes[a=='brca',]
c <- ifelse(b$Subtype_mRNA=="Basal",'TNBC','xx')
tnbc <- b[c=='TNBC',]
rt1 <- rt[,colnames(rt) %in% tnbc$pan.samplesID]
save(tnbc,file = 'tcga-tnbc.RData')
group_list=ifelse(as.numeric(substr(colnames(data1),14,15)) < 10,'tumor','normal')
tcga_normal <- data1[, group_list == 'normal']
rt2=subtypes[match(colnames(tcga_normal),subtypes$pan.samplesID),] 
rt3=rt2[rt2$Subtype_mRNA%in%c("Normal"),]
tcga_normal_need=tcga_normal[,colnames(tcga_normal)%in%rt3$pan.samplesID]
sameGene=intersect(row.names(rt1),row.names(tcga_normal_need))
input_data=cbind(tcga_normal_need[sameGene,],rt1[sameGene,])
save(input_data,file = 'TNBC-input.RData')
