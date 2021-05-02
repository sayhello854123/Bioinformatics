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
library(ggplot2)

load("./single cell data/metadata_tumor.RData")
load("./single cell data/tumor.RData")
raw.data=tumor
# Find ERCC's, compute the percent ERCC, and drop them from the data.
erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)
percent.ercc <- Matrix::colSums(raw.data[erccs, ])/Matrix::colSums(raw.data)
fivenum(percent.ercc)
ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = FALSE)
raw.data <- raw.data[-ercc.index,]
dim(raw.data)
"""
一共获取了92个ERCC基因，基因数目从55860   281 到55768   281
"""
folders=c('BC01','BC02','BC03','BC04','BC05','BC06','BC07','BC08','BC10','BC11')
BRCA <- CreateSeuratObject(counts = raw.data, project = folders, min.cells = 3, min.features = 50)#这一步基因数目从55768删除到24071
row.names(metadata) <- metadata[,1]
BRCA <- AddMetaData(object = BRCA, metadata = metadata)
BRCA <- AddMetaData(object = BRCA, percent.ercc, col.name = "percent.ercc")
head(BRCA@meta.data)
# Calculate percent ribosomal genes and add to metadata
ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = BRCA@assays$RNA@data), value = TRUE)
percent.ribo <- Matrix::colSums(BRCA@assays$RNA@counts[ribo.genes, ])/Matrix::colSums(BRCA@assays$RNA@data)
fivenum(percent.ribo)
BRCA <- AddMetaData(object = BRCA, metadata = percent.ribo, col.name = "percent.ribo")
rownames(BRCA)[grepl('^mt-',rownames(BRCA),ignore.case = T)]
rownames(BRCA)[grepl('^Rp[sl]',rownames(BRCA),ignore.case = T)]
#seurat包的函数：PercentageFeatureSet 可以用来计算线粒体基因含量
BRCA[["percent.mt"]] <- PercentageFeatureSet(BRCA, pattern = "^MT-")
fivenum(BRCA[["percent.mt"]][,1])
rb.genes <- rownames(BRCA)[grep("^RP[SL]",rownames(BRCA),ignore.case = T)]
C<-GetAssayData(object = BRCA, slot = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
fivenum(percent.ribo)
A<- VlnPlot(object = BRCA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
BRCA <- subset(BRCA,subset = nFeature_RNA > 500 & nFeature_RNA < 8200 & percent.mt < 25)
B<-VlnPlot(object = BRCA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save(A,B,file = './single cell data/picture/A.RData')
dim(BRCA)
ggplot2::ggsave(filename = './single cell data/QC/pearplot_before_qc.png',plot = A, width = 12, height = 5)
ggplot2::ggsave(filename = './single cell data/QC/vlnplot_after_qc.png',plot = B, width = 12, height = 5)
saveRDS(BRCA, file="./single cell data/scRNA.rds")
