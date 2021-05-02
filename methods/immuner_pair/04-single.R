rm(list=ls())
options(stringsAsFactors = F)
library(ggplot2)
library(Seurat)
library(GEOquery)
library(dplyr)
library(Seurat)
library(patchwork)
#cluster差异基因
dir.create("./immune single cell/cell_identify")
scRNA <- readRDS("./immune single cell/scRNA2.rds")
diff.wilcox = FindAllMarkers(scRNA)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "./immune single cell/cell_identify/diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "./immune single cell/cell_identify/top10_diff_genes_wilcox.csv", row.names = F)
top10_genes <- read.csv("./immune single cell/cell_identify/top10_diff_genes_wilcox.csv")
top10_genes = CaseMatch(search = as.vector(top10_genes$gene), match = rownames(scRNA)) 
plot1 = DoHeatmap(scRNA, features = top10_genes, group.by = "seurat_clusters", group.bar = T, size = 4)
ggsave("./immune single cell/cell_identify/top10_markers.pdf", plot=plot1, width=8, height=6) 
ggsave("./immune single cell/cell_identify/top10_markers.png", plot=plot1, width=8, height=6)

#挑选部分基因
select_genes <- c('C1QB','SH2D1A','CD79A','TCF4')
#vlnplot展示
p1 <- VlnPlot(scRNA, features = select_genes, pt.size=0, group.by="celltype", ncol=2)
ggsave("./immune single cell/cell_identify/selectgenes_VlnPlot.png", p1, width=6 ,height=8)
#featureplot展示
p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label=T, ncol=2)
ggsave("./immune single cell/cell_identify/selectgenes_FeaturePlot.png", p2, width=8 ,height=12)
p3=p1|p2
ggsave("./immune single cell/cell_identify/selectgenes.png", p3, width=10 ,height=8)


library(SingleR)
refdata <- HumanPrimaryCellAtlasData()
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),
                  'celltype'] <- celltype$celltype[i]
}
p1 = DimPlot(scRNA, group.by="celltype", label=T, label.size=5, reduction='tsne')
p2 = DimPlot(scRNA, group.by="celltype", label=T, label.size=5, reduction='umap')
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
ggsave("./immune single cell/cell_identify/tSNE_celltype.pdf", p1, width=7 ,height=6)
ggsave("./immune single cell/cell_identify/UMAP_celltype.pdf", p2, width=7 ,height=6)
ggsave("./immune single cell/cell_identify/celltype.pdf", p3, width=10 ,height=5)
ggsave("./immune single cell/cell_identify/celltype.png", p3, width=10 ,height=5)

metadata<- scRNA@meta.data
alltyple <- cbind(scRNA@meta.data$orig.ident,scRNA@meta.data$seurat_clusters)
alltyple <- as.data.frame(alltyple)
alltyple <- cbind(alltyple,scRNA@meta.data$orig.ident)
alltyple$sample=as.character(alltyple$`scRNA@meta.data$orig.ident`)
alltyple$clusters=as.character(scRNA@meta.data$seurat_clusters)
ggplot(alltyple,aes(x=clusters)) + geom_bar(aes(fill=factor(sample)))
saveRDS(scRNA, file="./immune single cell/scRNA3.rds")
