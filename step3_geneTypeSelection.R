#设置路径
setwd("E:/肿瘤信息学")

#清楚缓存，释放内存
rm(list=ls())
gc()

#加载包
# BiocManager::install("rtracklayer")#,force = TRUE)
# remove.packages("IRanges") ##or remove mannually at C:\Users\lenovo\Documents\R\win-library\3.5
# packageVersion("IRanges")
# update.packages('IRanges')
# BiocManager::install("BiocParallel",force = TRUE)
ps <- c('rtracklayer', 'dplyr','forcats','magrittr','BiocParallel','IRanges') #假如你需要批量导入的包是这四个
for(i in ps){library(i, character.only = T)}; rm(i) # 移除了名为i的对象，删除垃圾对象

#读取数据
expSet1=read.table("gene_expression_matrix_final.txt",row.names=1,header =T,sep='\t',stringsAsFactors=F, comment.char = "!",check.names = F)
expSet1$gene_id=rownames(expSet1)
test <- expSet1[1:50,];#View(test)
tail(expSet1$gene_id,10) # usually last 5 rows are not gene names
expSet1 <- expSet1[1:(length(expSet1$gene_id)-5),]

expSet2 <- expSet1 %>% 
  tidyr::separate(gene_id,into = c("gene_id","drop"),sep="\\.") %>% 
  dplyr::select(-drop)
write.csv(expSet2,file = "gene_expression_matrix_nopoint_final.csv",quote=F)


gtf1 <- rtracklayer::import('Homo_sapiens.GRCh38.89.chr.gtf')
gtf_df <- as.data.frame(gtf1)
save(gtf_df,file = "gtf_df.Rda")
test <- gtf_df[1:50,];#View(test)
mRNA_expSet <- gtf_df %>% # get mRNA 
  dplyr::filter(type=="gene",gene_biotype=="protein_coding") %>% 
  dplyr::select(c(gene_name,gene_id)) %>% 
  dplyr::inner_join(expSet2,by ="gene_id") 
#save(mRNA_expSet,file = "mRNA_expSet.Rda")
write.csv(mRNA_expSet,file = "mRNA_expSet.csv",row.names=F,quote=F)

