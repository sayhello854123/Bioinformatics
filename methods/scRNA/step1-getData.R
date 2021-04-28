rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(GEOquery)
library(dplyr)
library(Seurat)
library(patchwork)
library(clusterProfiler)
library('org.Hs.eg.db')
library(stringr) 
'''
获取单细胞表达矩阵
'''
rt <- data.table::fread('./GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt',data.table = F)
gset <- getGEO('GSE75688', destdir=".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F) 
phenoDat <- gset[[1]]
phenoDat <- pData(phenoDat)
'''
https://static-content.springer.com/esm/art%3A10.1038%2Fncomms15081/MediaObjects/41467_2017_BFncomms15081_MOESM467_ESM.xlsx
经过原文QC data的过滤我们保留了515个样本
'''
load('metadata.RData')

###处理样本基因
rownames(rt) <- rt[,1]
rt1 <- rt[,-c(1:3)]
same <- intersect(metadata[,1],colnames(rt1))
rt2 <- rt1[,same]
"""
First, TPM values <1 were considered unreliable and substituted with zero
"""
rt2[rt2 < 1] <- 0
"""
Second, TPM values were log2-transformed after adding a value of one
"""
rt3 <- log2(rt2+1)
"""
Third, genes expressed in <10% of all tumour groups were removed
"""
c <- rowSums(rt3 != 0)/length(rt3)
group <- ifelse(c>=0.1,'T','F')
table(group)
rt4 <- rt3[group=='T',]
"""
最终获取14647个基因，515个合格样本
我们的目标是研究原发乳腺癌
"""
table(metadata$`cell type:ch1`)
"""
metastatic breast cancer    primary breast cancer 
                     105                      410 
继续提取410个原发样本
"""
group <- ifelse(metadata$`cell type:ch1`=='primary breast cancer','t','f')
table(group)
metadata <- metadata[group=='t',]
same <- intersect(metadata[,1],colnames(rt4))
rt5 <- rt4[,same]
ids=data.frame(ensembl_id=str_split(rownames(rt5),
                                    '[.]',simplify = T)[,1],
               median=apply(rt5,1,median)
)
g2s=unique(toTable(org.Hs.egSYMBOL))
g2e=unique(toTable(org.Hs.egENSEMBL)) 
s2e=merge(g2e,g2s,by='gene_id')
table(ids$ensembl_id %in% s2e$ensembl)
ids$symbol=s2e[match(ids$ensembl_id,s2e$ensembl),3]
length(unique(ids$symbol))
head(ids) 
ids=ids[order(ids$symbol,ids$median,decreasing = T),]#把ids$symbol按照ids$median排序
ids=ids[!duplicated(ids$symbol),]#取出不重复的ids$symbol
dim(ids) 
rt5= rt5[rownames(ids),]#取出表达矩阵中ids有的行
same <- intersect(row.names(ids),row.names(rt5))
rt5 <- cbind(ids[same,],rt5[same,])
rt5 <- rt5[-10499,]
rownames(rt5) <- rt5[,3]
rt5 <- rt5[,-c(1:3)]
save(metadata,file = 'metadata.RData')
save(rt5,file = 'input.RData')
