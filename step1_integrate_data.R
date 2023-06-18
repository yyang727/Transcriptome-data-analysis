#设置路径
setwd("E:/肿瘤信息学")

#清除缓存
rm(list=ls())
gc()

#读取文件夹
first_category_name = list.files("E:/肿瘤信息学/gdc_download_20230308_234819.474455")            
#list.files命令得到"路径”文件夹下所有文件夹的名称
dir = paste("E:/肿瘤信息学/gdc_download_20230308_234819.474455/",first_category_name,sep="")   
#用paste命令构建路径变量dir,第一级目录(不包括文件名）的详细路径
n = length(dir)   

data3=NULL
for(i in 1:n){
  temp=list.files(path=dir[i],pattern = "*.tsv")#获取存放数据文件的文件名
  for(j in 1:length(temp)){
    data1=read.table(paste0(dir[i],"/",temp[j]),header =F,sep='\t',quote="",stringsAsFactors=F, comment.char = "#")
    #用paste0形成一个完整的能打开文件的路径，且不需要将第一行设置为表头以及避免注释内容的干扰
    fileName=temp[j]
    data1=data1[-(1:5),]  #去掉前五行的统计信息
    data3=cbind(data3,data1[,4])  #提取每个文件夹中的第四列数据
    colnames(data3)[ncol(data3)]=fileName  #将列名替换为每个样本独特的名（文件夹名）
  }
}
rownames(data3)=data1[,1]  #把每一行的基因ENSEMBL_ID赋值给行名

#输出得到的表达矩阵（此时行名为文件夹名）
write.table(data3,"gene_expression_matrix.txt",quote=FALSE,col.names=NA,row.names=T,sep='\t')