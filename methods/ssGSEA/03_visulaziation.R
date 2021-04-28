### 特定癌种中，癌和癌旁的差异
rm(list = ls())
load(file = "pancancer_tcga_gsva.Rdata")
test2 <- tcga_gsva[1:10,1:40]

table(tcga_gsva$type)
library(dplyr)
dd <- tcga_gsva %>% 
  filter(type=="BRCA") %>% 
  select(-c(1:2,4:34)) %>% 
  mutate(subtype = as.character(subtype)) %>% 
  filter(subtype %in% c("01","11")) %>% 
  mutate(sample = ifelse(subtype =="01","tumor","normal")) %>% 
  select(sample,everything())

### 调整数据
library(dplyr)
library(tidyr)
dd1 <- dd %>% 
  pivot_longer(cols=4:31,
               names_to= "celltype",
               values_to = "NES")

library(ggplot2)
library(ggpubr)
### 箱线图
ggplot(data =dd1, aes(x = celltype, y = NES))+
  geom_boxplot(aes(fill = sample),outlier.shape = NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=sample), label = "p.signif")

### 小提琴
ggplot(data =dd1, aes(x = celltype, y = NES))+
  geom_violin(aes(fill = sample),position = position_dodge(1),scale = "width")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=sample), label = "p.signif")

### 混合叠加
ggplot(data =dd1, aes(x = celltype, y = NES))+
  geom_boxplot(aes(fill = sample),position = position_dodge(1),width=.3,outlier.shape = NA)+
  geom_violin(aes(colour = sample),position = position_dodge(1),scale = "width",fill=NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=sample), label = "p.signif")
