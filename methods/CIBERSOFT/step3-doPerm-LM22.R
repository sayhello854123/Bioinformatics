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

itor <- 1
Ylist <- as.list(data.matrix(Y))
dist <- matrix()
# 就是把 CoreAlg 函数运行1000次
perm=1000
while(itor <= perm){
  print(itor) # 打印进度
  
  #random mixture
  yr <- as.numeric(Ylist[ sample(length(Ylist),dim(X)[1]) ])
  
  #standardize mixture
  yr <- (yr - mean(yr)) / sd(yr)
  
  #run CIBERSORT core algorithm
  result <- CoreAlg(X, yr)
  
  mix_r <- result$mix_r
  
  #store correlation
  if(itor == 1) {dist <- mix_r}
  else {dist <- rbind(dist, mix_r)}
  
  itor <- itor + 1
}
newList <- list("dist" = dist)
nulldist=sort(newList$dist) 
# 这个nulldist 主要是用来计算P值
if(F){
  
  P=perm
  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}
  print(nulldist)
}

header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
print(header)
save(nulldist,file = 'nulldist_perm_1000.Rdata')
load(file = 'nulldist_perm_1000.Rdata')
print(nulldist)
fivenum(print(nulldist))


