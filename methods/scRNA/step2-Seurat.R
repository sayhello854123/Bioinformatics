
load("input.RData")
load("metadata.RData")
raw.data=rt5
# Find ERCC's, compute the percent ERCC, and drop them from the data.
erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)

BRCA <- CreateSeuratObject(counts = raw.data,project = "scRNA", min.cells = 3, min.features = 50)
row.names(metadata) <- metadata[,1]
BRCA <- AddMetaData(object = BRCA, metadata = metadata)

#计算线粒体基因的百分比
BRCA[["percent.mt"]] <- PercentageFeatureSet(BRCA, pattern = "^MT-")
A<- VlnPlot(object = BRCA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Calculate percent ribosomal genes 
ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = BRCA@assays$RNA@data), value = TRUE)
percent.ribo <- Matrix::colSums(BRCA@assays$RNA@counts[ribo.genes, ])/Matrix::colSums(BRCA@assays$RNA@data)
fivenum(percent.ribo)
BRCA <- AddMetaData(object = BRCA, metadata = percent.ribo, col.name = "percent.ribo")



## QC
B <- FeatureScatter(BRCA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
BRCA <- subset(BRCA, subset = nFeature_RNA > 2000 & nFeature_RNA < 6500 )
#对数据进行标准化
BRCA <- NormalizeData(object = BRCA, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features
BRCA<- FindVariableFeatures(BRCA, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(BRCA), 20)
plot1 <- VariableFeaturePlot(BRCA)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
saveRDS(BRCA, file="scRNA1.rds")

