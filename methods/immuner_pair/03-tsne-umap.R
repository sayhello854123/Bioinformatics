
scRNA <- readRDS("./immune single cell/scRNA1.rds")
dir.create("./immune single cell/cluster")
##TSNE聚类分析
pc.num =1:11
scRNA <- FindNeighbors(scRNA, dims = pc.num)                #计算邻接距离
scRNA <- FindClusters(scRNA, resolution = 0.8)                 #对细胞分组,优化标准模块化
table(scRNA@meta.data$seurat_clusters)
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)

#非线性降维
scRNA = RunTSNE(scRNA, dims = pc.num)
embed_tsne <- Embeddings(scRNA, 'tsne')
write.csv(embed_tsne,'./immune single cell/embed_tsne.csv')
plot1 = DimPlot(scRNA, reduction = "tsne") 
ggsave("./immune single cell/cluster/tSNE.pdf", plot = plot1, width = 8, height = 7)
ggsave("./immune single cell/cluster/tSNE.png", plot = plot1, width = 8, height = 7)
#UMAP
scRNA <- RunUMAP(scRNA, dims = pc.num)
embed_umap <- Embeddings(scRNA, 'umap')
write.csv(embed_umap,'./immune single cell/cluster/embed_umap.csv') 
plot2 = DimPlot(scRNA, reduction = "umap") 
ggsave("./immune single cell/cluster/UMAP.pdf", plot = plot2, width = 8, height = 7)
ggsave("./immune single cell/cluster/UMAP.png", plot = plot2, width = 8, height = 7)

plotc <- plot1+plot2+ plot_layout(guides = 'collect')
ggsave("./immune single cell/cluster/tSNE_UMAP.pdf", plot = plotc, width = 10, height = 5)
ggsave("./immune single cell/cluster/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)
##保存数据
saveRDS(scRNA, file="./immune single cell/scRNA2.rds")

