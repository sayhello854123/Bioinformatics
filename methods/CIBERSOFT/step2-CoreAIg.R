rm(list = ls())
options(stringsAsFactors = F)
library(preprocessCore)
library(parallel)
library(e1071)
source("CIBERSORT.R") 

load(file = 'input.Rdata')
Y[1:4,1:4]
X[1:4,1:4]
dim(X)
dim(Y)
# 下面的演示是为了搞清楚 CoreAlg 函数
# 并不需要保存任何信息

# 从表达矩阵Y里面，随机挑选LM22矩阵基因数量的表达量值
Ylist <- as.list(data.matrix(Y))
yr <- as.numeric(Ylist[ sample(length(Ylist),dim(X)[1]) ])
# yr 这个时候是一个假设的样本
#standardize mixture,就是scale 函数
yr <- (yr - mean(yr)) / sd(yr)
boxplot(yr)
# 每次随机挑选的yr，都是需要走后面的流程

# 一切都是默认值的支持向量机
# 这里的X是LM22矩阵，不同的免疫细胞比例组合成为不同的yr
# 这里的yr是随机的，反推免疫细胞比例
out=svm(X,yr)
out
out$SV
# SVM 需要自行学习哦

# 需要修改的参数包括：type="nu-regression",kernel="linear",nu=nus,scale=F
svn_itor <- 3

y=yr
res <- function(i){
  if(i==1){nus <- 0.25}
  if(i==2){nus <- 0.5}
  if(i==3){nus <- 0.75}
  model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
  model
}

#Execute In a parallel way the SVM
if(Sys.info()['sysname'] == 'Windows') {
  
  out <- mclapply(1:svn_itor, res, mc.cores=1) 
}else {
  out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
}
# 运行了Support Vector Machines，函数是 svm {e1071}
out
#Initiate two variables with 0
nusvm <- rep(0,svn_itor)
corrv <- rep(0,svn_itor)

 
t <- 1
while(t <= svn_itor) {
  
  # 得到两个向量之间矩阵乘法的权重，此时应该只得到一个数字。
  # 这样做是乘以系数
  
  # 支持向量是数据集的点，它们靠近分隔类别的平面
  # 现在的问题是，我没有任何类别（离散变量，例如“运动”、“电影”），但我有一个连续变量
  mySupportVectors <- out[[t]]$SV
  
  # 系数定义
  myCoefficients <- out[[t]]$coefs
  weights = t(myCoefficients) %*% mySupportVectors
  
  # 设置权重和相关性
  weights[which(weights<0)]<-0
  w<-weights/sum(weights)
   
  # 根据对应的权重与参考集相乘
  u <- sweep(X,MARGIN=2,w,'*')
  
  # 统计每行总和
  k <- apply(u, 1, sum)
  nusvm[t] <- sqrt((mean((k - y)^2))) 
  corrv[t] <- cor(k, y)
  t <- t + 1
}
#pick best model
rmses <- nusvm
corrv
mn <- which.min(rmses)
mn  
model <- out[[mn]]
# 从nus为0.25,0.5,0.75的3个模型里面挑选一个即可

#get and normalize coefficients
 
q <- t(model$coefs) %*% model$SV 

q[which(q<0)]<-0
# w 就是计算后的22种免疫细胞的比例
w <- (q/sum(q))

mix_rmse <- rmses[mn]
mix_r <- corrv[mn]

# 会返回这个随机的y的免疫细胞组成情况，就是权重w
newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
newList

# 根据对应的权重与参考集相乘
u <- sweep(X,MARGIN=2,w,'*') 
k <- apply(u, 1, sum)
plot(y,k)
sqrt((mean((k - y)^2))) 
cor(k, y)
# 通常这个预测结果都惨不忍睹





# 每次把表达矩阵通过去卷积拆分成为LM22的免疫细胞比例结果

# 并且包装成为函数，如下：

#' CIBERSORT R script v1.03 (last updated 07-10-2015)
#' Note: Signature matrix construction is not currently available; use java version for full functionality.
#' Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#Core algorithm
CoreAlg <- function(X, y){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-e1071::svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- parallel::mclapply(1:svn_itor, res, mc.cores=1) else
    out <- parallel::mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}

