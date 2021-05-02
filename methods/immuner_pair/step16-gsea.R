library(limma)
library(org.Hs.eg.db)
library(fgsea)
library(ggplot2)
library(data.table)
library(clusterProfiler)
library(patchwork)
dir.create('GSEA')
riskFile="./cox/tcgaRisk.txt" 
load('./bulk data/02.tcga/tcga_brca.RData')
group=sapply(strsplit(colnames(data),"-"),'[',4)
table(group)
data1 <- data[,group=='01A']
group_list1=ifelse(as.numeric(substr(colnames(data1),14,15)) < 02,'tumor','normal')
tcga_brca <- data1[,group_list1 == "tumor" ]
dim(tcga_brca)

#读取risk文件
risk=read.table(riskFile,header=T,sep="\t",row.names=1,check.names=F)
sameSample=intersect(colnames(tcga_brca),row.names(risk))
data=tcga_brca[,sameSample]
risk=risk[sameSample,]

#low risk和high risk组比较
riskLow=risk[risk$risk=="low",]
riskHigh=risk[risk$risk=="high",]
dataLow=data[,row.names(riskLow)]
dataHigh=data[,row.names(riskHigh)]
meanLow=rowMeans(dataLow)
meanHigh=rowMeans(dataHigh)
meanLow[meanLow<0.00001]=0.00001
meanHigh[meanHigh<0.00001]=0.00001
logFC=log2(meanHigh/meanLow)
logFC=sort(logFC)
genes=names(logFC)

#输出GSEA需要的输入文件
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    #找出基因对应的id
id=cbind(genes,entrezIDs,logFC)
write.table(file="./GSEA/TCGA/id.txt",id,sep="\t",quote=F,row.names=F)

gmtFile="./GSEA/c2.cp.kegg.v7.1.entrez.gmt"                #GSEA数据集文件

rt=read.table("./GSEA/TCGA/id.txt",sep="\t",header=T,check.names=F)#读取输入文件

rt=rt[is.na(rt[,"entrezIDs"])==F,]
#write.table(rt,file="xid.txt",sep="\t",quote=F,row.names = F) 
geneList=rt$logFC
names(geneList) = rt$entrezIDs
geneList = sort(geneList, decreasing = TRUE)
#读取gmt文件
hallmarks <- read.gmt(gmtFile)
# 需要网络，大家会很拥挤，但是速度很快
y <- GSEA(geneList,TERM2GENE =hallmarks)
                nperm=10000)
cnetplot(y,foldChange = geneList)
y2<- setReadable(y,OrgDb=org.Hs.eg.db, keyType = "ENTREZID")
cnetplot(y2,showCategory = 4,
         foldChange = geneList,
         colorEdge = T)
c <- dotplot(y,showCategory=12,split=".sign")+facet_grid(~.sign)
ggplot2::ggsave(filename = './GSEA/TCGA/dotplot.pdf',plot = c, width =8, height =5.5)

yd <- data.frame(y)
library(enrichplot)
gseaplot2(y,"HALLMARK_E2F_TARGETS",color = "red",pvalue_table = T)
gseaplot2(y,11,color = "red",pvalue_table = T)
ridgeplot(y)
gseaplot2(y, geneSetID = 1:3)


### 自己添加文字加文字
index <- "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY"
gseaplot2(y,index ,color = "green")

anno <- yd[index , c("enrichmentScore","NES", "p.adjust")]
colnames(anno)[1] <- "ES"
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

P1 <- gseaplot2(y,index,color = "green")+
  annotate("text",0.6, 0.8, label = lab, hjust=0, vjust=0,size = 5)+ggtitle("Natural killer cells mediated cytotoxicity") +
  theme(plot.title = element_text(size=18,  family="serif",face = "bold",hjust = 0.5))
#ggplot2::ggsave(filename = './GSEA/TCGA/NATURAL_KILLER.pdf',plot = P, width =8, height = 7)

###
index <- "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY"
gseaplot2(y,index ,color = "green")

anno <- yd[index , c("enrichmentScore","NES", "p.adjust")]
colnames(anno)[1] <- "ES"
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

P2 <- gseaplot2(y,index,color = "green")+
  annotate("text",0.6, 0.8, label = lab, hjust=0, vjust=0,size = 5)+ggtitle("T cell receptor signaling pathway") +
  theme(plot.title = element_text(size=18,  family="serif",face = "bold",hjust = 0.5))
#ggplot2::ggsave(filename = './GSEA/TCGA/T cell signaling.pdf',plot = P, width =8, height = 7)

###
index <- "KEGG_JAK_STAT_SIGNALING_PATHWAY"
gseaplot2(y,index ,color = "green")

anno <- yd[index , c("enrichmentScore","NES", "p.adjust")]
colnames(anno)[1] <- "ES"
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

P3 <- gseaplot2(y,index,color = "green")+
  annotate("text",0.6, 0.8, label = lab, hjust=0, vjust=0,size = 5)+ggtitle("JAK-STAT signaling pathway") +
  theme(plot.title = element_text(size=18,  family="serif",face = "bold",hjust = 0.5))
#ggplot2::ggsave(filename = './GSEA/TCGA/JAK-STAT.pdf',plot = P, width =8, height = 7)

###
index <- "KEGG_CHEMOKINE_SIGNALING_PATHWAY"
gseaplot2(y,index ,color = "green")

anno <- yd[index , c("enrichmentScore","NES", "p.adjust")]
colnames(anno)[1] <- "ES"
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

P4 <- gseaplot2(y,index,color = "green")+
  annotate("text",0.6, 0.8, label = lab, hjust=0, vjust=0,size = 5)+ggtitle("	Chemokine signaling pathway") +
  theme(plot.title = element_text(size=18,  family="serif",face = "bold",hjust = 0.5))
#ggplot2::ggsave(filename = './GSEA/TCGA/Chemokine signaling.pdf',plot = P, width =8, height = 7)

###
index <- "KEGG_HEMATOPOIETIC_CELL_LINEAGE"
gseaplot2(y,index ,color = "green")

anno <- yd[index , c("enrichmentScore","NES", "p.adjust")]
colnames(anno)[1] <- "ES"
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

P5 <- gseaplot2(y,index,color = "green")+
  annotate("text",0.6, 0.8, label = lab, hjust=0, vjust=0,size = 5)+ggtitle("Hematopoietic cell lineage") +
  theme(plot.title = element_text(size=18,  family="serif",face = "bold",hjust = 0.5))
#ggplot2::ggsave(filename = './GSEA/TCGA/Hematopoietic cell.pdf',plot = P, width =8, height = 7)

###
index <- "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"
gseaplot2(y,index ,color = "green")

anno <- yd[index , c("enrichmentScore","NES", "p.adjust")]
colnames(anno)[1] <- "ES"
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

P6 <- gseaplot2(y,index,color = "green")+
  annotate("text",0.6, 0.8, label = lab, hjust=0, vjust=0,size = 5)+ggtitle("	Cytokine-cytokine receptor interaction") +
  theme(plot.title = element_text(size=18,  family="serif",face = "bold",hjust = 0.5))
#ggplot2::ggsave(filename = './GSEA/TCGA/Cytokine-cytokine.pdf',plot = P, width =8, height = 7)

P <-(P1|P2)/(P3|P4)/(P5|P6)
ggplot2::ggsave(filename = './GSEA/GSEA.pdf',plot = P, width =12, height =14)
