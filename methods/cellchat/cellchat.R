library(CellChat)
library(tidyverse)
library(ggalluvial)
library(Seurat)
library(patchwork)
library(ggpubr)
library(NMF)
rm(list=ls())
options(stringsAsFactors = FALSE)
picDir <- './01.single cell/02.cellchat/'
scRNA <- readRDS('00.data/01.single cell/scRNA2.rds')
#scRNA <- UpdateSeuratObject(scRNA)
#数据整理
cell <- c("TAM","Malignant cell","B cell","CAF" ,"TEC","T cell","HPC-like")
scRNA@meta.data$need_type <-ifelse(scRNA@meta.data$Type  %in% cell ,'need','non-need')
phe=scRNA@meta.data
table(phe$need_type)
cells.use <- row.names(scRNA@meta.data)[which(phe$need_type=='need')]
scRNA<-subset(scRNA,cells=cells.use) 
#创建cellchat对象
scRNA <- readRDS('00.data/01.single cell/scRNA6.rds')
data.input = scRNA@assays$RNA@data
meta = scRNA@meta.data


cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Type")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use ="Type" ) 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#可视化由多个配体受体或信号通路调节的细胞通信
cell <- c("B cell","CAF","HPC-like","Malignant cell","T cell","TAM","TEC")
for(i in cell){
  p <- netVisual_bubble(cellchat, sources.use = which(cell==i), targets.use = c(1:7), remove.isolate = FALSE)
 filename <- paste0(i, ".pdf",sep='')
  outfile <- paste(picDir,filename,sep="\\")
  pdf(file=outfile,
      width = 8,         #图片的宽度
      height =6)         #图片的高度
  print(p)
  dev.off()
}
##识别传入和传出信号功能贡献最大的细胞
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
pdf('./01.single cell/02.cellchat/heatmap1.pdf',height=7,width=6)
print(ht1)
dev.off()
pdf('./01.single cell/02.cellchat/heatmap2.pdf',height=7,width=6)
print(ht2)
dev.off()
# 细胞通讯模式和信号网络
selectK(cellchat, pattern = "outgoing")
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 3)
p1 = netAnalysis_river(cellchat, pattern = "outgoing")
ggsave('./01.single cell/02.cellchat/out_river.pdf', p1, width = 12, height = 6)
p2 =netAnalysis_dot(cellchat, pattern = "outgoing")
ggsave('./01.single cell/02.cellchat/out_dimplot.pdf', p2, width = 12, height = 6)
selectK(cellchat, pattern = "incoming")
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = 5)
p3 = netAnalysis_river(cellchat, pattern = "incoming")
ggsave('./01.single cell/02.cellchat/incom_river.pdf', p3, width = 12, height = 6)
p4 <- netAnalysis_dot(cellchat, pattern = "incoming")
ggsave('./01.single cell/02.cellchat/incom_dimplot.pdf', p4, width = 12, height = 6)
# 按功能相似性聚类
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
# Visualization in 2D-space
p5 = netVisual_embedding(cellchat, type = "functional")
p5 <- p5+ggtitle("Identify signaling groups based on their functional similarity") +theme(plot.title = element_text(hjust = 0.5))
ggsave("./01.single cell/02.cellchat/custer_pathway_function.pdf", p5, width = 9, height = 6)
#p = netVisual_embeddingZoomIn(cellchat, type = "functional")
#ggsave("custer_pathway_function2.png", p, width = 8, height = 6)

# 按结构相似性聚类
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
p8 = netVisual_embedding(cellchat, type = "structural")
p8 <- p8+ggtitle("Identify signaling groups based on structure similarity") +theme(plot.title = element_text(hjust = 0.5)) 
ggsave("./01.single cell/02.cellchat/custer_pathway_structure.pdf", p8, width = 9, height = 6)

##plot
P <- (p0|p1)/p2/p+plot_annotation(tag_levels = 'A')
ggsave("./01.single cell/cellchat.png", plot =P , width =14, height =12) 
##
p6 <-  netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:7), remove.isolate = FALSE)
p7 <-  netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1:7), remove.isolate = FALSE)
library(jpeg)
library(ggpubr)
A1 <- readJPEG('./01.single cell/02.cellchat/heatmap.jpg')
p0<-ggplot()+
  background_image(A1)+
  theme_void()
P1 <- p6+p7
ggsave('./01.single cell/02.cellchat/picture1.pdf', P1, width = 12, height = 6)
P2 <- p1+p3
ggsave('./01.single cell/02.cellchat/picture2.pdf', P2, width = 16, height = 6)
P3 <- p2+p4
ggsave('./01.single cell/02.cellchat/picture3.pdf', P3, width = 12, height = 6)
P4 <- p5+p8
ggsave('./01.single cell/02.cellchat/picture4.pdf', P4, width = 12, height = 6)
A2 <- readJPEG('./01.single cell/02.cellchat/picture1.jpg')
p1<-ggplot()+
  background_image(A1)+
  theme_void()
A3 <- readJPEG('./01.single cell/02.cellchat/picture2.jpg')
p2<-ggplot()+
  background_image(A1)+
  theme_void()
P <- p0/(P1)/(P2)/(P3)+plot_annotation(tag_levels = 'A')
ggsave("./01.single cell/02.cellchat/picture.png", plot =P , width =15, height = 16) 
