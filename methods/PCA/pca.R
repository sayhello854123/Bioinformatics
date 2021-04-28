library("FactoMineR")
library("factoextra")
library("tidyverse")
data <- read.table("SUM.txt", header=T, row.names=NULL,sep="\t")
rownames(data) <- data[,1]
data <- data[,-1]
data <- data[rowSums(data)!=0,]### 删除所有列值为0的样本
pca <- prcomp(t(data),scale = TRUE,center = TRUE)
var_explained  <-  pca$sdev^2/sum(pca$sdev^2)
write.table(pca$rotation,file="PC.xls",quote=F,sep="\t")   #输出特征向量
write.table(predict(pca),file="newTab.xls",quote=F,sep="\t")   #输出新表
pca.sum=summary(pca)
write.table(pca.sum$importance,file="importance.xls",quote=F,sep="\t")#输出PC比重
png("PCA1.png",width = 800,height =650)
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 50))
dev.off()
png("PCA2.1.png",width = 800,height = 600)
fviz_pca_ind(pca, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE)
dev.off()
png("PCA2.2.png",width = 800,height = 600)
fviz_pca_ind(pca, col.ind = "cos2", pointsize = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) 
dev.off()
png("PCA2.3.png",width = 800,height = 600)
fviz_cos2(pca, choice = "ind") + coord_flip()
dev.off()
png("PCA2.png",width = 800,height =650)
fviz_pca_ind(pca,col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE )
dev.off()
png("PCA3.png",width = 1050,height =1450)
pca$x %>% 
  as.data.frame %>%
  rownames_to_column("continent_letter") %>%
  separate(continent_letter,c("continent")) %>%
  ggplot(aes(x=PC1,y=PC2)) + geom_point(aes(color=continent),size=4) +
  theme_bw(base_size=32) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top")
dev.off()
