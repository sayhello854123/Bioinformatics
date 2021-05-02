library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle)
library(clusterProfiler)
library(org.Hs.eg.db)
rm(list=ls())
dir.create("enrich")
scRNA <- readRDS("./immune single cell/scRNA3.rds")

#比较cluster0和cluster1的差异表达基因
dge.cluster <- FindMarkers(scRNA,ident.1 = 2,ident.2 = 3)
sig_dge.cluster <- subset(dge.cluster, p_val_adj<0.01&abs(avg_log2FC)>1)
#比较B_cell和T_cells的差异表达基因
dge.celltype <- FindMarkers(scRNA, ident.1 = 'T_cells', ident.2 = 'Monocyte', group.by = 'celltype')
sig_dge.celltype <- subset(dge.celltype, p_val_adj<0.01&abs(avg_log2FC)>1)
ego_CC <- enrichGO(gene          = row.names(sig_dge.celltype),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.01)
ego_MF <- enrichGO(gene          = row.names(sig_dge.celltype),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.01)
ego_BP <- enrichGO(gene          = row.names(sig_dge.celltype),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.01)   

p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("barplot for Biological process")
p_CC <- barplot(ego_CC,showCategory = 10) + ggtitle("barplot for Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("barplot for Molecular function")
plotc <- p_BP/p_CC/p_MF
ggsave('./enrich/enrichGO1.png', plotc, width = 9,height = 12)


genelist <- bitr(row.names(sig_dge.celltype), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
p1 <- barplot(ekegg, showCategory=10)+ ggtitle("barplot for KEGG")
p2 <- dotplot(ekegg, showCategory=10)+ggtitle("dotplot for KEGG")
plotB = p1/p2
plotA <- plotc/plotB 
ggsave("./enrich/enrichKEGG.png", plotB, width =9, height =12 )
ggsave("./enrich/enrich.png", plot =plotA, width =12, height =25)

