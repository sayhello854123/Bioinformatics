 ###1.处理TCGA_brca数据###
library('limma')
symbol <- data.table::fread('./bulk data/02.tcga/gencode.v22.annotation.gene.probeMap',data.table = F)
rownames(symbol) <- symbol[,1]
tcga_fpkm <- data.table::fread('./bulk data/02.tcga/TCGA-BRCA.htseq_fpkm.tsv/TCGA-BRCA.htseq_fpkm.tsv',data.table = F)
rownames(tcga_fpkm) <- tcga_fpkm[,1]
same <- intersect(row.names(symbol),row.names(tcga_fpkm))
tcga_B <- cbind(symbol[same,],tcga_fpkm[same,])
tcga_B <- tcga_B[,-1]
rt<- tcga_B[,-c(2:6)]
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data[1:4][1:4]
save(data,file='./bulk data/02.tcga/tcga_brca.RData')

###2.获取TCGA_brca 肿瘤的免疫基因####
rm(list = ls())
load('./bulk data/02.tcga/tcga_brca.RData')
group=sapply(strsplit(colnames(data),"-"),'[',4)
table(group)
data1 <- data[,group=='01A']
group_list1=ifelse(as.numeric(substr(colnames(data1),14,15)) < 02,'tumor','normal')
tcga_brca <- data1[,group_list1 == "tumor" ]
dim(tcga_brca)
gene <- read.table("./immune single cell/cell_identify/mark.txt", header=F, check.names=F, sep="\t")
sameGene=intersect(as.vector(gene[,1]),rownames(tcga_brca))
tcga_brca_geneExp=tcga_brca[sameGene,]
save(tcga_brca_geneExp,file='./bulk data/02.tcga/tcgaexp.RData')

####3.筛选符合条件的样本
rm(list = ls())
load('./bulk data/02.tcga/tcgaexp.RData')
cli <- data.table::fread('./bulk data/02.tcga/TCGA-BRCA.survival.tsv/TCGA-BRCA.survival.tsv',data.table = F)
rownames(cli) <- cli[,1]
group <- ifelse(cli$OS.time>1,'T','F')
table(group)
cli <- cli[group=='T',]
same <- intersect(as.vector(cli[,1]),row.names(t(tcga_brca_geneExp)))
tcga_brca_immune <- tcga_brca_geneExp[,same]
cli <- data.table::fread('./bulk data/02.tcga/TCGA-BRCA.GDC_phenotype.tsv/TCGA-BRCA.GDC_phenotype.tsv',data.table = F)
group <- cli$sample_type.samples
table(cli$sample_type.samples)
cli <- cli[group=='Primary Tumor',]
dim(cli)
same <- intersect(as.vector(cli[,1]),row.names(t(tcga_brca_immune)))
tcga_brca_immune <- tcga_brca_immune[,same]
save(tcga_brca_immune,file ='./bulk data/02.tcga/tcgaexp1.RData')
