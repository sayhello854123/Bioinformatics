###1 整理metabric的数据
rm(list = ls())
library('limma')
exp <- data.table::fread('./bulk data/brca_metabric/data_expression_median.txt',data.table =  F)
rt=exp[,-2]
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=na.omit(data)
gene <- read.table("./immune single cell/cell_identify/mark.txt", header=F, check.names=F, sep="\t")
sameGene=intersect(as.vector(gene[,1]),rownames(data))
metabric_expimuune <- data[sameGene,]
save(metabric_expimuune,file = './bulk data/03.metabric/metabric_exp.RData')

##2.筛选符合的样本##
rm(list = ls())
load('./bulk data/03.metabric/metabric_exp.RData')
load('./bulk data/brca_metabric/metabric_cli.RData')
metabric_cli=metabric_cli[,c(1:2)]
metabric_cli=na.omit(metabric_cli)
metabric_cli$futime <- metabric_cli$futime*12
group <- ifelse(metabric_cli$futime>1,'T','F')
table(group)
metabric_cli<- metabric_cli[group=='T',]
same <- intersect(row.names(metabric_cli),row.names(t(metabric_expimuune)))
metabricexp<- metabric_expimuune[,same]
dim(metabricexp)
save(metabricexp,file = './bulk data/03.metabric/metabric_exp2.RData')
