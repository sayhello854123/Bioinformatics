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
# 运行CIBERSORT官方函数，获得三个官方算法定义的Function：
# CIBERSORT，CoreAlg和doPerm。

sig_matrix = 'LM22-ref.txt'
mixture_file = 'test.txt'
# 主要是 read.table 的参数可以随意修改，根据你自己的表达矩阵txt文件适应性调整即可
X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
Y <- read.table(mixture_file, header=T, sep="\t", check.names=F)

# 运行 CIBERSORT 会写出 CIBERSORT-Results.txt 文件覆盖掉之前的结果
# 务必看清楚 test.txt 这个表达矩阵格式哦！！！
# CIBERSORT(sig_matrix, mixture_file, perm=10, QN=TRUE)
# 后面可以对 CIBERSORT-Results.txt 里面的结果进行各式各样的可视化