#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")

library("org.Hs.eg.db")          #引用包
rt=read.table("symbol.txt",sep="\t",check.names=F,header=T)    #读取文件
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    #找出基因对应的id
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)

rt1=read.table("id.txt",sep="\t",header=T,check.names=F)
##is.na函数会将NA判断为T,数值为F
rt1=rt1[is.na(rt[,"entrezID"])==F,]