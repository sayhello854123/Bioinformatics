"""
我上传的两个转换都要看，不要只看这一个，不然无法完成
"""
options(stringsAsFactors = F)
expMatrix <- read.csv("easy_input.csv",
                      row.names = 1, header = TRUE, as.is = T)
# 1 首先要保证表达矩阵的行名和存放基因长度向量的名字一致, 这一步非常重要-----------------------------------------------------------------------

eff_length2 <-read.csv("eff2_length.csv", row.names = 1, header = T)
eff_length2$gene_id <- rownames(eff_length2)

for(i in 1:4){
  colnames(expMatrix)[i]=gsub("\\.","-",colnames(expMatrix)[i])
  
}

# 1 载入GTF文件 ----------------------------------------------------------------------

library(GenomicFeatures)

txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.90.chr.gtf",format="gtf") #使用最新的gtf文件




# 通过exonsBy获取每个gene上的所有外显子的起始位点和终止位点，然后用reduce去除掉重叠冗余的部分，最后计算长度 -----------


exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})

exons_gene_2_lens=data.frame(t(data.frame(exons_gene_lens)))

names(exons_gene_2_lens)="length"

head(exons_gene_2_lens)[1]

#write.csv(exons_gene_2_lens,"eff2_length.csv",row.names = T)

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



# 2 -----------------------------------------------------------------------

# trim the expression matrix and effetive gene length
expMatrix <- expMatrix[feature_ids %in%rownames(eff_length2) ,]
mm <- match(rownames(expMatrix), rownames(eff_length2))
############################
a=c(1,2,3)

b=c(2,3,4,5)
match(a,b)
b[match(a,b)]
############################
max(mm)
eff_length2 <- eff_length2[mm, ]

if (identical(rownames(eff_length2), rownames(expMatrix))){
  print("GTF and expression matix now have the same gene and gene in same order")
}


# 3 -----------------------------------------------------------------------
#计算公式
#TPMi=(Ni/Li)*1000000/sum(Ni/Li+……..+ Nm/Lm)


x <- expMatrix/eff_length2$length
y=t(x)
expMatrix_tpm <- t( t(x) / colSums(x) ) * 1e6 
#查看前三个基因的TPM值
expMatrix_tpm[1:3,]

write.csv(expMatrix_tpm,"expMatrix_tpm.csv",row.names = T)

# 4 -----------------------------------------------------------------------

#计算基因的长度
#如果只是想粗略了解一下表达情况，可以简单把基因在染色体上的起始位置和结束位置之差用作标准化的长度。
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)

#查看基因组参数
mart = useMart('ensembl')
listDatasets(mart)

#你需要哪个基因组，就复制它在dataset列里的词，放在下面这行的`dataset = `参数里
#此处以人类为例
bmart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# 从输入数据里提取基因名
feature_ids <- rownames(expMatrix)
str(feature_ids)
attributes = c(
  "ensembl_gene_id",
  #"hgnc_symbol",
  "chromosome_name",
  "start_position",
  "end_position"
)
filters = "ensembl_gene_id"

feature_info <- biomaRt::getBM(attributes = attributes, 
                               filters = filters, 
                               values = feature_ids[1], 
                               mart = bmart)

mm <- match(feature_ids, feature_info[[filters]])
feature_info_full <- feature_info[mm, ]
rownames(feature_info_full) <- feature_ids

# 计算基因的有效长度
eff_length <- abs(feature_info_full$end_position - feature_info_full$start_position)
names(eff_length) <- feature_info_full$ensembl_gene_id
