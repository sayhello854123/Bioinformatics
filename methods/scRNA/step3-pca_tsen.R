rm(list=ls())
options(stringsAsFactors = F)
library(ggplot2)
library(Seurat)
library(GEOquery)
library(dplyr)
library(Seurat)
library(patchwork)

BRCA <- readRDS("scRNA1.rds")
##PCA分析
all.genes <- VariableFeatures(BRCA)
BRCA <- ScaleData(BRCA, features = all.genes) #PCA降维之前的标准预处理步骤
'''
上一步找到的高变基因，常常会包含一些细胞周期相关基因。它们会导致细胞聚类发生一定的偏移，
即相同类型的细胞在聚类时会因为细胞周期的不同而分开。
如在Spatially and functionally distinct subclasses of breast cancer-associated fibroblasts revealed by single cell RNA sequencing 
中 vCAFs和cCAFs中只有cell cycle genes的表达差异。
查看我们选择的高变基因中有哪些细胞周期相关基因
'''
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(BRCA))
#细胞周期评分
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(BRCA))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(BRCA))
BRCA <- CellCycleScoring(object=BRCA,  g2m.features=g2m_genes,  s.features=s_genes)
#查看细胞周期基因对细胞聚类的影响
scRNAa <- RunPCA(BRCA, features = c(s_genes, g2m_genes))
p <- DimPlot(scRNAa, reduction = "pca", group.by = "Phase")
ggsave("cellcycle_pca.png", p, width = 8, height = 6)
#需要消除细胞周期的影响
scRNA <- ScaleData(BRCA, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(BRCA))

#PCA降维并提取主成分
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident") 
plot2 <- ElbowPlot(scRNA, ndims=20, reduction="pca") 
plotc <- plot1+plot2
ggsave("pca.pdf", plot = plotc, width = 8, height = 4) 
ggsave("pca.png", plot = plotc, width = 8, height = 4)
scRNAb <- JackStraw(scRNA, num.replicate = 100)
scRNAb <- ScoreJackStraw(scRNAb, dims = 1:20)
D <- JackStrawPlot(scRNAb, dims = 1:20)
ggsave("pca1.png",D, width = 8, height = 6)

##发现前16个PC很关键纳入后续分析

##TSNE聚类分析
pc.num =1:16
scRNA <- FindNeighbors(scRNA, dims = pc.num)                #计算邻接距离
scRNA <- FindClusters(scRNA, resolution = 0.5)                 #对细胞分组,优化标准模块化
table(scRNA@meta.data$seurat_clusters)
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)

#非线性降维
scRNA = RunTSNE(scRNA, dims = pc.num)
embed_tsne <- Embeddings(scRNA, 'tsne')
write.csv(embed_tsne,'embed_tsne.csv')
plot1 = DimPlot(scRNA, reduction = "tsne") 
ggsave("tSNE.pdf", plot = plot1, width = 8, height = 7)
ggsave("tSNE.png", plot = plot1, width = 8, height = 7)

#UMAP
scRNA <- RunUMAP(scRNA, dims = pc.num)
embed_umap <- Embeddings(scRNA, 'umap')
write.csv(embed_umap,'embed_umap.csv') 
plot2 = DimPlot(scRNA, reduction = "umap") 
ggsave("UMAP.pdf", plot = plot2, width = 8, height = 7)
ggsave("UMAP.png", plot = plot2, width = 8, height = 7)

#合并tSNE与UMAP
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
ggsave("tSNE_UMAP.png", plot = plotc, width = 10, height = 5)
saveRDS(scRNA, file="scRNA2.rds")
