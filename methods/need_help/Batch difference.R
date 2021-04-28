library(limma)
library(edgeR)
foldChange=1 
padj=0.01

setwd("G:/")
result='result'  
dir.create(result)###创造文件夹
outPath <- "G:/result"


a = list.files("text") ##将文件夹内的文件生成一个list
dir = paste("./text/",a,sep="") ##设置读入文件的路径
n = length(dir) 
dir1 <- paste(a,sep="")###设置输出文件的名字 
out_fileName <- sapply(dir1,function(x){
  paste(x, "", sep='')}) ##格式这里我的文件名中有txt 
out_filePath  <- sapply(out_fileName, function(x){
  paste(outPath ,x,sep='/')}) ##输出路径名
                                     #用paste命令构建路径变量dir
                             



for (i in 1:n){
  rt = read.table(file = dir[i], header=T, sep="\t")
  rownames(rt)=rt[,1]
  rt = rt[,-1]
  exp=rt[,1:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  data=data[rowMeans(data)>0,]
  
  #group=c("normal","tumor","tumor","normal","tumor") 
  group=c("normal","tumor")                     
  #design <- model.matrix(~group)
  y <- DGEList(counts=data,group=group)
  
  bcv = 0.4                                    # For human bcv=0.4,other bcv=0.1 
  et <- exactTest(y, dispersion=bcv^2) 
  results = et$table
  results$q = p.adjust(results$PValue, method = 'fdr')
  
  topTags(et)
  ordered_tags <- topTags(et, n=100000)
  
  allDiff=ordered_tags$table
  allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
  diff=allDiff
  diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
  write.table(diffSig,file=out_filePath[i],sep = '\t',quote = F)
}
