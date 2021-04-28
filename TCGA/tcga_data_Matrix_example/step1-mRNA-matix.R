library("rjson")
result <- fromJSON(file = "metadata.cart.2020-09-27.json")
file_name <- unlist(lapply(result, function(x){x$file_name}))##result中file_name的位置来确定如何写
case_id <- unlist(lapply(result, function(x){x$associated_entities[[1]]$case_id}))##原理同上
TCGA_id <- unlist(lapply(result, function(x){x$associated_entities[[1]]$entity_submitter_id}))
id <- data.frame(file_name=file_name,case_id=case_id,TCGA_id=TCGA_id)
rownames(id) <- id[,1]
dir <- './gdc'
file <- list.files(path = dir,pattern = '*.htseq.counts',recursive = T)##tcga 文件后缀
Count <- lapply(file,function(x){
  result <- read.table(file = file.path(dir,x),sep = '\t')[,1:2]
  return(result)
})
need_matrix <- t(do.call(cbind,Count))
colnames(need_matrix)=need_matrix[1,]
need_matrix=need_matrix[seq(2,nrow(need_matrix),by=2),]
rownames(need_matrix)=file
rownames(id) <- substr(rownames(id),1,49)
same <- intersect(row.names(id),row.names(need_matrix))
miA<- cbind(id[same,],need_matrix[same,])
miA <- miA[,-c(1:2)]
rownames(miA) <- miA[,1]
miA <- miA[,-1]
SKCM_mRNA<- t(miA)
group_list1=ifelse(as.numeric(substr(colnames(SKCM_mRNA),14,15)) < 10,'tumor','normal')
tumor <- SKCM_mRNA[,group_list1 == "tumor" ]
normal <- SKCM_mRNA[,group_list1 == "normal" ]
SKCM_mRNA_matrix<- cbind(normal,tumor)
