rm(list=ls())
options(stringsAsFactors = F)
library(limma)
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
rt <- data.table::fread('./single cell data/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt',data.table = F)
rt <- rt[,-1]
rt <- rt[,-2]
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data <- as.data.frame(data)


gset <- getGEO('GSE75688', destdir=".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F) 
phenoDat <- gset[[1]]
phenoDat <- pData(phenoDat)
table(phenoDat$`cell type:ch1`)
primary <- phenoDat[phenoDat$`cell type:ch1`=='primary breast cancer',]
sample <- data.table::fread('./single cell data/GSE75688_final_sample_information.txt/GSE75688_final_sample_information.txt',data.table = F)
table(sample$type)
sample <- sample[sample$type=='SC',]
table(sample$index)
sample1 <- sample[sample$index=='Tumor',]
sample2 <- sample[sample$index2=='Immune',]
same1 <- intersect(sample1[,1],colnames(data))
tumor <- data[,same1]
same2 <- intersect(colnames(tumor),primary[,1])
tumor <- tumor[,same2]
same3 <- intersect(sample2[,1],colnames(data))
immune <- data[,same3]
rt <- read.table('./single cell data/phenoDat.txt',header = T,check.names = F,sep = '\t')
rownames(rt) <- rt[,2]
same1 <- intersect(row.names(rt),colnames(immune))
meta <- rt[same1,]
rownames(meta) <- meta[,1]
save(meta,file = './single cell data/metadata_immune.RData')
save(immune,file = './single cell data/immune.RData')
same2 <- intersect(row.names(rt),colnames(tumor))
metadata <- rt[same2,]
rownames(metadata) <- metadata[,1]
save(metadata,file = './single cell data/metadata_tumor.RData')
save(tumor,file ='./single cell data/tumor.RData')
