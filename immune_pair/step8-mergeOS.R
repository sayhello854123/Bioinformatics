rm(list = ls())
data1=data.table::fread('./bulk data/01.data/tcgapair.csv',data.table = F)
rownames(data1)=data1[,1]
data1=data1[,-1]
colnames(data1)[1047]=c("TCGA-BH-A0E0-01A")
table(apply(data1,1,mean))
data1=data1[apply(data1,1,mean)>0.3,]
data1=data1[apply(data1,1,mean)<0.7,]

data3=data.table::fread('./bulk data/01.data/metabricpair.csv',data.table = F)
rownames(data3)=data3[,1]
data3=data3[,-1]
colnames(data3)[1903]=c("MB-7299")
table(apply(data3,1,mean))
data3=data3[apply(data3,1,mean)>0.3,]
data3=data3[apply(data3,1,mean)<0.7,]

same=intersect(row.names(data1),row.names(data3))
data1=data1[same,]
data3=data3[same,]

#####合并生存时间##
dir.create('cox')
cli <- data.table::fread('./bulk data/02.tcga/TCGA-BRCA.survival.tsv/TCGA-BRCA.survival.tsv',data.table = F)
rownames(cli) <- cli[,1]
group <- ifelse(cli$OS.time>1,'T','F')
table(group)
cli <- cli[group=='T',]
cli$futime <- cli$OS.time
cli$fustat <- cli$OS
same=intersect(row.names(cli),row.names(t(data1)))
rt=cbind(cli[same,],t(data1)[same,])[,-c(1:4)]
save(rt,file = './cox/cox_input.RData')

cli <- read.table('./bulk data/brca_metabric/123.txt',sep = '\t',header = T,check.names = F,row.names = 1)
metabric_cli=cli [,c(1:2)]
metabric_cli=na.omit(metabric_cli)
metabric_cli$futime <- metabric_cli$futime*30
group <- ifelse(metabric_cli$futime>1,'T','F')
table(group)
metabric_cli<- metabric_cli[group=='T',]
same=intersect(row.names(metabric_cli),row.names(t(data3)))
rt=cbind(metabric_cli[same,],t(data3)[same,])
save(rt,file = './cox/metabric_input.RData')


