rm(list = ls())
options(stringsAsFactors = F)
# 首先读取两个文件
sig_matrix <-"LM22-ref.txt"  # cibersoft 内置数据库挖掘
mixture_file <- "mRNA2.txt"  # 约80M，TCGA数据库

# 两个表达矩阵需要取交集
if(!file.exists( 'input.Rdata')){
  
  #read in data
  X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
  Y <- read.table(mixture_file, header=T, sep="\t", check.names=F)
  Y <- Y[!duplicated(Y[,1]),]
  rownames(Y)<-Y[,1]
  Y<-Y[,-1]
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  Y[1:4,1:4]
  X[1:4,1:4]
  dim(X)
  dim(Y)
  
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  
  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}
  
  QN = F
  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  dim(X)
  dim(Y)
  
  
  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  Y[1:4,1:4]
  X[1:4,1:4]
  boxplot(X[,1:4])
  save(X,Y,file = 'input.Rdata')
}

