rm(list=ls())
options(stringsAsFactors = F)
library(ggplot2)
library(Seurat)
library(GEOquery)
library(dplyr)
library(Seurat)
library(patchwork)
dir.create("./immune single cell/PCA")
scRNA <- readRDS("./immune single cell/scRNA.rds")
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) 
top10 <- head(VariableFeatures(scRNA), 10) 
plot1 <- VariableFeaturePlot(scRNA) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom") 
ggplot2::ggsave(filename = './immune single cell/PCA/VariableFeatures.pdf',plot = plot, width = 12, height = 5)
#数据中心化
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
'''
上一步找到的高变基因，常常会包含一些细胞周期相关基因。它们会导致细胞聚类发生一定的偏移，
即相同类型的细胞在聚类时会因为细胞周期的不同而分开。
如在Spatially and functionally distinct subclasses of breast cancer-associated fibroblasts revealed by single cell RNA sequencing 
中 vCAFs和cCAFs中只有cell cycle genes的表达差异。
查看我们选择的高变基因中有哪些细胞周期相关基因
'''
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNA))
#细胞周期评分
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNA))
scRNA <- CellCycleScoring(object=scRNA,  g2m.features=g2m_genes,  s.features=s_genes)
#查看细胞周期基因对细胞聚类的影响
scRNAa <- RunPCA(scRNA, features = c(s_genes, g2m_genes))
p1 <- DimPlot(scRNAa, reduction = "pca", group.by = "Phase")
scRNAb <- ScaleData(scRNAa, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(scRNA))
scRNAc <- RunPCA(scRNAb, features = c(s_genes, g2m_genes))
p2 <- DimPlot(scRNAc, reduction = "pca", group.by = "Phase")
plotA <- CombinePlots(plots = list(p1, p2),legend="bottom")
ggplot2::ggsave(filename = './immune single cell/PCA/cellcycle_pca.pdf',plot = plotA, width = 8, height = 5)

dim(scRNAb)
#PCA降维并提取主成分
scRNA <- RunPCA(scRNAb, features = VariableFeatures(scRNA)) 
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident") 
plot2 <- ElbowPlot(scRNA, ndims=20, reduction="pca") 
plotc <- plot1+plot2
ggplot2::ggsave(filename = './immune single cell/PCA/pca.pdf',plot = plotc, width = 8, height = 4)
DimHeatmap(scRNA, dims = 1:6, cells = 100, balanced = TRUE)
scRNA <- JackStraw(scRNA, num.replicate = 100)
scRNA <- ScoreJackStraw(scRNA, dims = 1:18)
E1 = JackStrawPlot(scRNA, dims = 1:18)
ggplot2::ggsave(filename = './immune single cell/PCA/pca_empirical.pdf',plot = E1, width = 8, height = 7)
saveRDS(scRNA, file="./immune single cell/scRNA1.rds")
