library('limma')
dir.create('01.data')
symbol <- data.table::fread('./01.data/gencode.v22.annotation.gene.probeMap',data.table = F)
rownames(symbol) <- symbol[,1]
tcga_fpkm <- data.table::fread('./TCGA-LUAD.htseq_fpkm.tsv/TCGA-LUAD.htseq_fpkm.tsv',data.table = F)
rownames(tcga_fpkm) <- tcga_fpkm[,1]
same <- intersect(row.names(symbol),row.names(tcga_fpkm))
tcga_B <- cbind(symbol[same,],tcga_fpkm[same,])
tcga_B <- tcga_B[,-1]
rt<- tcga_B[,-c(2:6)]
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data[1:4][1:4]
save(data,file = './01.data/TCGA_luad.RData')
