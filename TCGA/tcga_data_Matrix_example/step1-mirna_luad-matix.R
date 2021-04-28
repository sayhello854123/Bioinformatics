#install.packages("rjson") 
library("rjson")
result <- fromJSON(file = "metadata.cart.2020-09-09.json")
name <- unlist(lapply(result, function(x){x[[6]]}))
cid <- unlist(lapply(result, function(x){x[[11]][[1]][[4]]}))
tcga_id <- unlist(lapply(result, function(x){x[[12]][[3]][[1]][[7]]}))
id <- data.frame(name=name,cid=cid,tcga_id=tcga_id)
rownames(id) <- id[,1]
dir <- './gdc_download_20200909_022932.149639'
fi <- list.files(path = dir,pattern = '*.quantification.txt',recursive = T)
miRNA <- lapply(fi,function(x){
                             result <- read.table(file = file.path(dir,x),sep = '\t')[,1:2]
                              return(result)
                           })
mi_df <- t(do.call(cbind,miRNA))
colnames(mi_df)=mi_df[1,]
mi_df=mi_df[seq(2,nrow(mi_df),by=2),]
group=sapply(strsplit(fi,"/"),'[',2)
fi[1]
group[1]
rownames(mi_df)=group
same <- intersect(row.names(id),row.names(mi_df))
miA <- cbind(id[same,],mi_df[same,])
rownames(miA) <- miA[,3]
miA <- miA[,-c(1:4)]
group1=sapply(strsplit(rownames(miA),"_"),'[',1)
rownames(miA)=group1
miRNA <- t(miA)  

group2=sapply(strsplit(colnames(miRNA),"-"),'[',4) 
table(group2)
group_list1=ifelse(as.numeric(substr(colnames(miRNA),14,15)) < 10,'tumor','normal')
tumor <- miRNA[,group_list1 == "tumor" ]
normal <- miRNA[,group_list1 == "normal" ]
miRNA <- cbind(normal,tumor)
