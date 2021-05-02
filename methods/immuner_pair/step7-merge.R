load('./bulk data/02.tcga/tcgaexp1.RData')
load('./bulk data/03.metabric/metabric_exp2.RData')
dim(tcga_brca_immune)
dim(metabricexp)
data1=tcga_brca_immune[apply(tcga_brca_immune,1,mad)>0.5,]
dim(data1)
data2=metabricexp[apply(metabricexp,1,mad)>0.5,]
dim(data2)
same <- intersect(row.names(data1),row.names(data2))
tcga <- data1[same,]
metabric <- data2[same,]

gene2=read.table("./bulk data/01.data/Geneappend3.txt", header=F, check.names=F, sep="\t")
sameGene=intersect(intersect(as.vector(gene[,1]),rownames(tcga)),row.names(metabric)) 

tcga_immune <- tcga[sameGene,]
metabric_immune <- metabric[sameGene,]
save(tcga,file = './bulk data/01.data/tcga_single_exp.RData')
save(metabric,file = './bulk data/01.data/metabriic_single_exp.RData')
save(tcga_immune,file = './bulk data/01.data/tcga_immune_single_exp.RData')
save(metabric_immune,file = './bulk data/01.data/metabric_immune_single_exp.RData')

out=rbind(ID=colnames(tcga_immune),tcga_immune)
write.table(out,file="./bulk data/01.data/tcga.txt",sep="\t",quote=F,col.names=F)

out=rbind(ID=colnames(metabric_immune),metabric_immune)
write.table(out,file="./bulk data/01.data/metabric.txt",sep="\t",quote=F,col.names=F)


A=rownames(data1)##将每个数据集中的deg赋值给A/B/C，主要是为了方面
B=rownames(data2)
C=gene2$V1
D <- gene$V1
venn.plot <- venn.diagram(
  list(TCGA=A,METABRIC=B,IMMPORT= C,single_cell = D),
  filename = "./bulk data/1.tiff",##韦恩图的名字
  lty = 1,
  lwd = 1,
  col = "black",  ##圈的颜色
  fill = c("#DC143C", "#00BFFF", "#7CFC00",'#FFD700'),##对应每个圈的颜色，有几个数据集，就需要有相应数量的颜色
  alpha = 0.60,
  cat.col = "black",##此处设置每个数据集的名称颜色，也可以使用c（）函数输入三种颜色
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.07,
  cex = 0.8
)
