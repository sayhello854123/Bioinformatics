options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if(!require("data.table")) install.packages("data.table",update = F,ask = F)
if(!require("GSVA")) BiocManager::install("GSVA",update = F,ask = F)

### 数据下载：https://pan.baidu.com/s/1fah6gNKqQnwL57A_EtzHWw 
### 提取码： 9jvu
### 需要表达量数据和marker数据
### marker数据的获取
### Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and Predictors of Response to Checkpoint Blockade
### https://www.cell.com/cms/10.1016/j.celrep.2016.12.019/attachment/f353dac9-4bf5-4a52-bb9a-775e74d5e968/mmc3.xlsx

rm(list = ls())
### 1.准备细胞marker
cellMarker <- data.table::fread("cellMarker.csv",data.table = F)
colnames(cellMarker)[2] <- "celltype"

cellMarker <- lapply(split(cellMarker,cellMarker$celltype), function(x){
  dd = x$Metagene
  unique(dd)
})

save(cellMarker,file = "cellMarker_ssGSEA.Rdata")

### 2.准备表达量矩阵
### 行是基因，列是样本
expr <- data.table::fread("exprMat.txt",data.table = F)
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.matrix(expr)

### 3.使用ssGSEA量化免疫浸润
library(GSVA)
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")

### 简单作图看一下
library(pheatmap)
pheatmap(gsva_data, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = TRUE,#列聚类，可以看出样本之间的区分度
         ##annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,# 显示行名
         ##scale = "row", #以行来标准化
         color =colorRampPalette(c("blue", "white","red"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 50, cellheight = 10,# 格子比例
         fontsize = 10)
