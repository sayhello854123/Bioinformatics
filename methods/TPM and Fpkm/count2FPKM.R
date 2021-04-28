"""
GTF 文件下载地址：ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.chr.gtf.gz
"""
options(stringsAsFactors = F)
expMatrix <- read.csv("easy_input.csv",
                      row.names = 1, header = TRUE, as.is = T)
#查看前三个基因的read count
expMatrix[1:3,]

for(i in 1:4){
  colnames(expMatrix)[i]=gsub("\\.","-",colnames(expMatrix)[i])
  
}

# 1 载入GTF文件 ----------------------------------------------------------------------

library(GenomicFeatures)

txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.90.chr.gtf",format="gtf")




# 通过exonsBy获取每个gene上的所有外显子的起始位点和终止位点，然后用reduce去除掉重叠冗余的部分，最后计算长度 -----------


exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})

exons_gene_2_lens=data.frame(t(data.frame(exons_gene_lens)))

names(exons_gene_2_lens)="length"

head(exons_gene_2_lens)[1]

write.csv(exons_gene_2_lens,"eff2_length.csv",row.names = T)

# 获得每个基因下所有外显子的总长度后，就可以利用下述公式计算FPKM了 --------------------------------------

#FPKM= total exon reads/ (mapped reads (Millions) * exon length(KB))



# 2 -----------------------------------------------------------------------
eff_length2 <-read.csv("eff2_length.csv", row.names = 1, header = T)

eff_length2$gene_id <- rownames(eff_length2)
# 从输入数据里提取基因名
feature_ids <- rownames(expMatrix)

# 检查gtf文件和表达量输入文件里基因名的一致性
if (! all(feature_ids %in% rownames(eff_length2))){
  tbl <- table(feature_ids %in% rownames(eff_length2))
  msg1 <- sprintf("%i gene is shared, %i gene is specified", tbl[[2]],tbl[[1]])
  warning(msg1)
} 

if (! identical(feature_ids, rownames(eff_length2))){
  msg2 <- sprintf("Given GTF file only contain %i gene, but experssion matrix has %i gene", nrow(eff_length2), nrow(expMatrix))
  warning(msg2)
}

# trim the expression matrix and effetive gene length
expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2),]
mm <- match(rownames(expMatrix), rownames(eff_length2))
eff_length2 <- eff_length2[mm, ]


#identical : The safe and reliable way to test two objects for being exactly  equal. 
#It returns TRUE in this case, FALSE in every other case.--------

if (identical(rownames(eff_length2), rownames(expMatrix))){
  print("GTF and expression matix now have the same gene and gene in same order")
}


countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}


fpkms <- apply(expMatrix, 2, countToFpkm, effLen = eff_length2$length)

fpkms.m<-data.frame(fpkms)
colnames(fpkms.m)<-colnames(expMatrix)
dim(fpkms.m)

write.csv(fpkms.m,"fpkm.csv",row.names = T)
