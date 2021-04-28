##发现前16个PC很关键纳入后续分析
scRNA <- readRDS("./single cell data/scRNA1.rds")
dir.create("./single cell data/cluster")
##TSNE聚类分析
pc.num =1:14
scRNA <- FindNeighbors(scRNA, dims = pc.num)                #计算邻接距离
scRNA <- FindClusters(scRNA, resolution = 0.8)                 #对细胞分组,优化标准模块化
table(scRNA@meta.data$seurat_clusters)
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)

#非线性降维
scRNA = RunTSNE(scRNA, dims = pc.num)
embed_tsne <- Embeddings(scRNA, 'tsne')
write.csv(embed_tsne,'./single cell data/embed_tsne.csv')
plot1 = DimPlot(scRNA, reduction = "tsne") 
ggsave("./single cell data/cluster/tSNE.pdf", plot = plot1, width = 8, height = 7)
ggsave("./single cell data/cluster/tSNE.png", plot = plot1, width = 8, height = 7)
#UMAP
scRNA <- RunUMAP(scRNA, dims = pc.num)
embed_umap <- Embeddings(scRNA, 'umap')
write.csv(embed_umap,'./single cell data/cluster/embed_umap.csv') 
plot2 = DimPlot(scRNA, reduction = "umap") 
ggsave("./single cell data/cluster/UMAP.pdf", plot = plot2, width = 8, height = 7)
ggsave("./single cell data/cluster/UMAP.png", plot = plot2, width = 8, height = 7)
E5 <- plot1+plot2
E5 <- plot1+plot2+ plot_layout(guides = 'collect')
save(E5,file = './single cell data/picture/c.RData')
ggsave("./single cell data/cluster/tSNE_UMAP.pdf", plot = plotc, width = 10, height = 5)
ggsave("./single cell data/cluster/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)
##保存数据
saveRDS(scRNA, file="./single cell data/scRNA2.rds")

