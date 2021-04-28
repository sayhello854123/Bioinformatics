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
load(file = 'output_obj.Rdata')

# Step3:将CIBERSORT_Result挑选整理，去除没有差异表达的细胞。
library(dplyr)
library(tidyr)
library(tidyverse)
cibersort_raw <- read.table("CIBERSORT-Results.txt",header = T,sep = '\t') %>%
  rename("Patients" = "Mixture") %>%
  select(-c("P.value","Correlation","RMSE"))
# 通过管道符一步步先将CIBERSORT_Results读入R语言中，并将其第一列列名“Mixture”修改为“Patiens”。
#并赋值给cibersort_raw。

cibersort_tidy <- cibersort_raw %>%
  remove_rownames() %>%
  column_to_rownames("Patients")
# 将cibersort_raw第一列变为列名后赋值给cibersort_tidy。

flag <- apply(cibersort_tidy,2,function(x) sum(x == 0) < 
                dim(cibersort_tidy)[1]/2)
# 筛选出0值太多的一些细胞。

cibersort_tidy <- cibersort_tidy[,which(flag)] %>%
  as.matrix() %>%
  t()
# 留下在大部分样本中有所表达的细胞。

bk <- c(seq(0,0.2,by = 0.01),seq(0.21,0.85,by=0.01))
# breaks用来定义数值和颜色的对应关系。

# Step4:将CIBERSORT_Result进行可视化
#1）热图
library(pheatmap)
library(RColorBrewer)
pheatmap(
  cibersort_tidy,
  breaks = bk,
  cluster_cols = T,
  scale = "row",
  cluster_row = T,
  border_color = NA,
  show_colnames = F,
  show_rownames = T,
  color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
            colorRampPalette(colors = c("white","red"))(length(bk)/2)
  ))
#调整参数让热图更加美观。

#柱状图可视化细胞占比预测
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
cibersort_barplot <- cibersort_raw %>%
  gather(key = Cell_type,value = Proportion,2:23)
#使用RColorBrewer包配置需要的色彩方案，使用gather函数中的key-value对应关系重建细胞名称和比例的对应关系并赋值给cibersort_barplot

#cibersort_barplot$Patient1 <- factor(cibersort_barplot$Patient,
#                                   levels = str_sort(unique(cibersort_barplot$Patient),
#                                                      numeric = T))

ggplot(cibersort_barplot,aes(Patients,Proportion,fill = Cell_type)) + 
  geom_bar(position = "stack",stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + theme_bw() +
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(23))
#调整参数让柱状图更加美观。

#直观箱线图
ggplot(cibersort_barplot,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,coulour = "black") + theme_bw() + 
  labs(x = "", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
  scale_fill_manual(values = mypalette(23))
#调整参数让柱状图更加美观。
