#设置路径
setwd("E:/肿瘤信息学")

#清楚缓存，释放内存
rm(list=ls())
gc()

#读入基因表达数据和样本信息
raw=read.table("gene_expression_matrix.txt",row.names=1,header =T,sep='\t',
               stringsAsFactors=F, comment.char = "!",check.names = F)
sample=read.table("gdc_sample_sheet.2023-03-08.tsv",header =T,
                  sep='\t',stringsAsFactors=F, comment.char = "!")

#
id1=colnames(raw)
id2=sample$File.Name
sample.sorted=sample[match(id1,id2),c(2,7,8)]
barcode=sample.sorted[,2]
typeID=substr(barcode,14,1000)
colnames(raw)=barcode

#group=ifelse(substr(typeID,1,1)=='0',"tumor","normal")
group=rep(1,length(typeID))*NA
group[substr(typeID,1,1)=='0']="tumor"
group[substr(typeID,1,1)=='1']="normal"
sample.sorted$condition=group

expSet=raw[,order(group)]
group_inf=sample.sorted[order(group),-1]
write.table(expSet,"gene_expression_matrix_final.txt",quote=FALSE,col.names=NA,row.names=T,sep='\t')
write.table(group_inf,"group_information.txt",quote=FALSE,row.names=F,col.names=T,sep='\t')
