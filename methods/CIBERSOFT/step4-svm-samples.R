rm(list = ls())
options(stringsAsFactors = F)
load(file = 'input.Rdata')
Y[1:4,1:4]
X[1:4,1:4]
dim(X)
dim(Y)
library(preprocessCore)
library(parallel)
library(e1071)
source("CIBERSORT.R") 

header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
print(header)
load(file = 'nulldist_perm_1000.Rdata')

output <- matrix()
itor <- 1
mix <- dim(Y)[2]
pval <- 9999

P=1000
# 表达矩阵的每个样本，都需要计算一下LM22的比例
#iterate through mix
while(itor <= mix){
  
  ##################################
  ## Analyze the first mixed sample
  ##################################
  
  y <- Y[,itor]
  
  #标准化样本数据集
  y <- (y - mean(y)) / sd(y)
  
  #执行SVR核心算法
  result <- CoreAlg(X, y)
  
  #获得结果
  w <- result$w
  mix_r <- result$mix_r
  mix_rmse <- result$mix_rmse
  
  #计算p-value
  if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
  
  #输出output
  out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
  if(itor == 1) {output <- out}
  else {output <- rbind(output, out)}
  itor <- itor + 1
  
}
head(output)

#save results
write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)

#return matrix object containing all results
obj <- rbind(header,output)
obj <- obj[,-1]
obj <- obj[-1,]
obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
rownames(obj) <- colnames(Y)
colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
obj
save(obj,file = 'output_obj.Rdata')

