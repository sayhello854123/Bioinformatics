rm(list=ls())
options(stringsAsFactors = F)
library(ggplot2)
library(Seurat)
library(GEOquery)
library(dplyr)
library(Seurat)
library(patchwork)
#cluster差异基因
dir.create("cell_identify")
scRNA <- readRDS("scRNA2.rds")
diff.wilcox = FindAllMarkers(scRNA)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(all.markers, "cell_identify/diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "cell_identify/top10_diff_genes_wilcox.csv", row.names = F)

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
ggsave("cell_identify/tSNE_celltype.pdf", p1, width=7 ,height=6)
ggsave("cell_identify/UMAP_celltype.pdf", p2, width=7 ,height=6)
ggsave("cell_identify/celltype.pdf", p3, width=10 ,height=5)
ggsave("cell_identify/celltype.png", p3, width=10 ,height=5)
saveRDS(scRNA, file="scRNA3.rds")
