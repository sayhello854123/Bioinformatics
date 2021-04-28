###  http://xena.ucsc.edu/
### 处理TCGA数据
rm(list = ls())
dd <- data.table::fread("tcga_RSEM_gene_tpm",data.table = F)
test <- dd[1:100,1:100]

### 清洁数据三部曲
### 1.基因注释 2.行列转置 3.增加分组

### 获取临床信息,之后只保留有临床信息的样本
clin <- data.table::fread("Survival_SupplementalTable_S1_20171025_xena_sp",data.table = F)

colnames(clin)[2:3] <- c("TCGA_id","type")
### 同时有临床信息的样本
index <- clin$sample %in% colnames(dd)
### 获取有临床信息的样本
clin <- clin[index,]
dd <- dd[,c("sample",clin$sample)]
colnames(dd)[1] <- "gene_id"

test <- dd[1:10,1:10]
### 探针转换
### 下载后先解压gencode.v23.annotation.gtf.gz
### 读入需要时间，耐心等待
### gtf1 <- rtracklayer::import('gencode.v23.annotation.gtf')
### gtf_df <- as.data.frame(gtf1)
### save(gtf_df,file = "gtf_df.Rdata")
load("gtf_df.Rdata")

### 提取编码mRNA
library(dplyr)
library(tidyr)
library(tibble)

tcga_panmRNA_expr <- gtf_df %>% 
  #筛选gene,和编码指标
  filter(type=="gene",gene_type=="protein_coding") %>%
  #选取两列
  dplyr::select(c(gene_name,gene_id)) %>% 
  #和表达量数据合并
  inner_join(dd,by ="gene_id") %>% 
  ## 去掉多余的列
  dplyr::select(-gene_id) %>% 
  ## 去掉基因名称中可能的NA，可有可无
  filter(gene_name!="NA") %>% 
  ## 去掉重复
  distinct(gene_name,.keep_all = T) %>% 
  ## 列名转为行名
  column_to_rownames("gene_name")

test <- tcga_panmRNA_expr[1:10,1:10]

### 防止你的机子崩溃掉，可以保存一下数据
save(clin,tcga_panmRNA_expr,file = "tcga_panmRNA_expr.Rdata")
#################################################################################
## 看到教程中的rm应该感到高兴
rm(list=ls())
### 用ssGSEA来量化浸润水平
### 1.加载marker
load(file = "cellMarker_ssGSEA.Rdata")
### 2.加载表达量
load(file = "tcga_panmRNA_expr.Rdata")
expr <- tcga_panmRNA_expr
test <- expr[1:10,1:10]
expr <- as.matrix(expr)
library(GSVA)
### 挺耗时间的，调用了12个线程，17:12开始, 17:37结束
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")
## save(gsva_data,file = "gsva_data_TCGA.Rdata")
## load(file = "gsva_data_TCGA.Rdata")
test <- gsva_data[1:10,1:10]
tcga_gsva <- as.data.frame(t(gsva_data))
test <- tcga_gsva[1:10,1:10]
## 添加分组信息
tcga_gsva <- cbind(clin,subtype=substring(rownames(tcga_gsva),14,15),tcga_gsva)
save(tcga_gsva,file = "pancancer_tcga_gsva.Rdata")

#################################################################################
#### 此时可以把表达量的数据变成清洁数据
## 加载表达量数据
load(file = "tcga_panmRNA_expr.Rdata")
test <- tcga_panmRNA_expr[1:10,1:10]
tcga_panmRNA_expr <- as.data.frame(t(tcga_panmRNA_expr))
test <- tcga_panmRNA_expr[1:10,1:10]

tcga_panmRNA_expr <- cbind(clin,subtype=substring(rownames(tcga_panmRNA_expr),14,15),tcga_panmRNA_expr)
test <- tcga_panmRNA_expr[1:10,1:10]
save(tcga_panmRNA_expr,file = "tcga_panmRNA_expr_with_clin.Rdata")

### 量化结束，可以做点分析了
