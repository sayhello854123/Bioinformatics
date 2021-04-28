rm(list = ls())
load(file = "tcga_panmRNA_expr_with_clin.Rdata")
test1 <- tcga_panmRNA_expr[1:10,1:40]

load(file = "pancancer_tcga_gsva.Rdata")
test2 <- tcga_gsva[1:10,1:40]
### 提取乳腺癌的数据

library(dplyr)
## 获取表达量数据
expr_data <- tcga_panmRNA_expr %>% 
  filter(type=="BRCA") %>% 
  dplyr::select(-c(1:35))
test <- expr_data[1:10,1:10]

## 获取免疫数据
immu_data <- tcga_gsva %>% 
  filter(type=="BRCA") %>% 
  dplyr::select(-c(1:35))

gene <- "KLF5"
y <- as.numeric(expr_data[,gene])

cor_data <- do.call(rbind,lapply(colnames(immu_data),function(x){
  dd <- cor.test(as.numeric(immu_data[,x]),y,method ="spearman",exact=FALSE)
  data.frame(cell=x,cor=dd$estimate,p.value=dd$p.value)
}))


### 画图展示全貌
library(dplyr)
library(ggplot2)
cor_data %>% 
  filter(p.value <0.05) %>% 
  ggplot(aes(cor,forcats::fct_reorder(cell,cor)))+
  geom_segment(aes(xend=0,yend=cell))+
  geom_point(aes(col=p.value,size=abs(cor)))+
  scale_colour_gradientn(colours=c("#7fc97f","#984ea3"))+
  #scale_color_viridis_c(begin = 0.5, end = 1)+
  scale_size_continuous(range =c(2,8))+
  theme_bw()+
  ylab(NULL)

### 筛选p值有意义的细胞
imucells <- cor_data %>% 
  filter(p.value <0.05) %>% 
  arrange(desc(cor)) %>% 
  pull(cell) %>% 
  as.vector()

imucells

library(ggplot2)
corr_eqn <- function(x,y,digits=2) {
  test <- cor.test(x,y,method="spearman")
  paste(paste0("n = ",length(x)),
        paste0("r = ",round(test$estimate,digits),"(pearson)"),
        paste0("p.value= ",round(test$p.value,digits)),
        sep = ", ")
}
imucell <- "CD56dim natural killer cell"
## 提取数据
plot_df <- data.frame(
  gene = expr_data[,gene],
  imucell = immu_data[,imucell]
)
## 作图
plot_df %>% 
  ggplot(aes(gene,imucell))+
  geom_point(col="#984ea3")+
  geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#fdc086")+
  geom_rug(col="#7fc97f")+
  theme_minimal()+
  xlab(paste0(gene," (log2(TPM))"))+
  ylab(paste0(imucell," (NES)"))+
  labs(title = paste0(corr_eqn(plot_df$gene,plot_df$imucell)))+
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm"))
