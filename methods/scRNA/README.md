# 单细胞处理<br> 
## Single-cell RNA-seq enables comprehensive tumour and immune cell profiling in primary breast cancer<br>     

接下来我们看一下数据的预处理<br> 
  To remove cells with low-quality sequencing values, four filtering criteria were applied: (1) number of total reads; (2) mapping rate; (3) number of detected genes; and (4) portion of intergenic region. 

[数据的前处理](https://github.com/sayhello854123/R/blob/master/scRNA/Seurat/step1-getData.R)<br> 
##Seurat 处理流程<br>
导入数据发现，进行ERCC基因的筛选<br>
由于每个单细胞都是独特的，不可能开展重复实验并评估噪音。因此，必须采取一些质量控制手段，以确保数据的可靠性。专家建议，向每个细胞裂解液中加入已知序列和数量的合成mRNA，如外源RNA对照联盟（ERCC）开发的加标RNA。这些RNA的读数将提供样本间差异的信息。

```
load("input.RData")
load("metadata.RData")
raw.data=rt5
# Find ERCC's, compute the percent ERCC, and drop them from the data.
erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)
```
![](https://github.com/sayhello854123/R/blob/master/scRNA/Seurat/pictuer/2.png)  
得到的结果显示没有外援基因  
创建Seurat对象
```
BRCA <- CreateSeuratObject(counts = raw.data,project = "scRNA", min.cells = 3, min.features = 50)
row.names(metadata) <- metadata[,1]
BRCA <- AddMetaData(object = BRCA, metadata = metadata)
```
计算线粒体基因的百分比
```
BRCA[["percent.mt"]] <- PercentageFeatureSet(BRCA, pattern = "^MT-")
A<- VlnPlot(object = BRCA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Calculate percent ribosomal genes 
ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = BRCA@assays$RNA@data), value = TRUE)
percent.ribo <- Matrix::colSums(BRCA@assays$RNA@counts[ribo.genes, ])/Matrix::colSums(BRCA@assays$RNA@data)
fivenum(percent.ribo)
BRCA <- AddMetaData(object = BRCA, metadata = percent.ribo, col.name = "percent.ribo")
```
![](https://github.com/sayhello854123/R/blob/master/scRNA/Seurat/pictuer/03.png)
通过上图我们发现nFeature_RNA中大部细胞的基因数目2000-6500，第二个图纵坐标是测到的所有基因序列是数目，第三个图是线粒体的百分比，这里面我们发现没有线粒体
```
## QC
B <- FeatureScatter(BRCA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
BRCA <- subset(BRCA, subset = nFeature_RNA > 2000 & nFeature_RNA < 6500 )
```
![](https://github.com/sayhello854123/R/blob/master/scRNA/Seurat/pictuer/05.jpg)
我们发现我们测序深度和基因数目相关性很高  

对数据进行标准化
```
BRCA <- NormalizeData(object = BRCA, normalization.method = "LogNormalize", scale.factor = 10000)  
```
Identification of highly variable features
```
BRCA<- FindVariableFeatures(BRCA, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(BRCA), 20)
plot1 <- VariableFeaturePlot(BRCA)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
```
![](https://github.com/sayhello854123/R/blob/master/scRNA/Seurat/pictuer/4.jpg)
 
PCA分析
```
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
```
上一步找到的高变基因，常常会包含一些细胞周期相关基因。它们会导致细胞聚类发生一定的偏移，
即相同类型的细胞在聚类时会因为细胞周期的不同而分开。
如在Spatially and functionally distinct subclasses of breast cancer-associated fibroblasts revealed by single cell RNA sequencing 
中 vCAFs和cCAFs中只有cell cycle genes的表达差异。
查看我们选择的高变基因中有哪些细胞周期相关基因  
```
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(BRCA))
```
细胞周期评分  
```
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(BRCA))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(BRCA))
BRCA <- CellCycleScoring(object=BRCA,  g2m.features=g2m_genes,  s.features=s_genes)
```  
查看细胞周期基因对细胞聚类的影响  
```
scRNAa <- RunPCA(BRCA, features = c(s_genes, g2m_genes))
p <- DimPlot(scRNAa, reduction = "pca", group.by = "Phase")
ggsave("cellcycle_pca.png", p, width = 8, height = 6)
```
![](https://github.com/sayhello854123/R/blob/master/scRNA/Seurat/pictuer/cellcycle_pca.png)
 需要消除细胞周期的影响
```
scRNA <- ScaleData(BRCA, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(BRCA))
```  
PCA降维并提取主成分  
```
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
```
![](https://github.com/sayhello854123/R/blob/master/scRNA/Seurat/pictuer/pca.png)<br>
![](https://github.com/sayhello854123/R/blob/master/scRNA/Seurat/pictuer/pca1.png)<br>

发现前16个PC很关键纳入后续分析


##TSNE聚类分析<br>
```
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
```
![](https://github.com/sayhello854123/R/blob/master/scRNA/Seurat/pictuer/tSNE.png)<br>

```
#UMAP
scRNA <- RunUMAP(scRNA, dims = pc.num)
embed_umap <- Embeddings(scRNA, 'umap')
write.csv(embed_umap,'embed_umap.csv') 
plot2 = DimPlot(scRNA, reduction = "umap") 
ggsave("UMAP.pdf", plot = plot2, width = 8, height = 7)
ggsave("UMAP.png", plot = plot2, width = 8, height = 7)
```
![](https://github.com/sayhello854123/R/blob/master/scRNA/Seurat/pictuer/UMAP.png)<br>
```
合并tSNE与UMAP
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
ggsave("tSNE_UMAP.png", plot = plotc, width = 10, height = 5)
saveRDS(scRNA, file="scRNA2.rds")
```
![](https://github.com/sayhello854123/R/blob/master/scRNA/Seurat/pictuer/tSNE_UMAP.png)<br>
#细胞类型鉴定<br>
通过marker基因鉴定细胞类型目前仍是最常用的方法，应用此方法的基础是找到各个cluster的显著高表达的基因，这就需要差异分析提供候选基因列表。Seura提供多种差异分析的方法，默认wilcox方法，MASK方法和DESeq2方法也可以留意一下。
```
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
```
SingleR鉴定细胞类型
```
library(SingleR)
refdata <- HumanPrimaryCellAtlasData()
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
```
与CellMarker：http://biocc.hrbmu.edu.cn/CellMarker/index.jsp<br>
以及PanglaoDB：https://panglaodb.se/index.html<br>
结果一致<br>
```
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
```
![](https://github.com/sayhello854123/R/blob/master/scRNA/Seurat/cell_identify/tSNE_celltype.pdf)<br>
![](https://github.com/sayhello854123/R/blob/master/scRNA/Seurat/cell_identify/UMAP_celltype.pdf)<br>
![](https://github.com/sayhello854123/R/blob/master/scRNA/Seurat/cell_identify/celltype.png)<br>
结果发现除了0簇为Monocyte 其余全为Epithelial cell<br>
#进一步判断不同组分来源样本
```
metadata<- scRNA@meta.data
alltyple <- cbind(scRNA@meta.data$orig.ident,scRNA@meta.data$seurat_clusters)
alltyple <- as.data.frame(alltyple)
alltyple <- cbind(alltyple,scRNA@meta.data$orig.ident)
alltyple$sample=as.character(alltyple$`scRNA@meta.data$orig.ident`)
alltyple$clusters=as.character(scRNA@meta.data$seurat_clusters)
ggplot(alltyple,aes(x=clusters)) + geom_bar(aes(fill=factor(sample)))
```
![](https://github.com/sayhello854123/R/blob/master/scRNA/Seurat/cell_identify/xxl.png)<br>
