# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("ggplot2")
# loading R packages 
#library(tidyverse)
library(magrittr)
library(glue)
library(DESeq2)
library(ggplot2)
# for more information, please refer to:
# http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

rm(list=ls())

raw=read.csv("mRNA_expSet.csv",header =T,stringsAsFactors=F, comment.char = "!",check.names = F)
rownames(raw)=raw$gene_id
geneIDset=raw[,c(1,2)]
raw=raw[-c(1,2)]
group_inf=read.table("group_information.txt",header =T,sep='\t',stringsAsFactors=F, comment.char = "!")
group_inf$condition=factor(group_inf$condition,levels=c("normal","tumor"))

##dir.create("DESeq2")
dds <- DESeqDataSetFromMatrix(countData =as.matrix(raw), colData=group_inf,design = ~ condition)
mcols(dds)=geneIDset

#########上面的语句也可以用下面的代码替换###########
# se <- SummarizedExperiment(assays = list(counts = as.matrix(raw)))
# se$condition <- factor(group)
# dds <- DESeqDataSet(se, design = ~ condition)
####################################################

# filter genes with low expression level
######这个方法也可以#############
# rs <- rowMeans(raw)
# geneKeep <- rs > quantile(rs, 0.2)
# sum(geneKeep)

######这个方法也可以#############
geneKeep <- rowSums(counts(dds)) >= 10
dds <- dds[geneKeep, ]
head(dds)
View(dds)
# save vsd profile for visualization later
##获得归一化的基因表达矩阵：a matrix of transformed, normalized counts，可用于做热图
vsd <- vst(dds, blind = T)
saveRDS(vsd, "DESeq2/vsd.rds")

#expSet.normalized=as.data.frame(counts(dds,normalized=T))
#normalized:logical indicating whether or not to divide the counts by the size factors or normalization factors

##################################

# perform DEG test
dds <- DESeq(dds)
resultsNames(dds)


res=results(dds,contrast=c("condition","tumor","normal"),tidy=T)
res$gene_name=mcols(dds)@listData[["gene_name"]]

resOrdered=res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
summary(res)
write.table(resOrdered,"DESeq2/allResult.txt", row.names = F, quote = F, sep = "\t")
##write.table里col.names = NA表示行名输出时行名那一列的列名为空

# ####对log FC进行压缩,但p-value q-value不变，因此影响DEG数量，一般变少######
# DEGres <- list()
# DEGres$tumor_vs_normal <- lfcShrink(dds, coef = "condition_tumor_vs_normal", type = "apeglm")
# #lfcShrink:Adds shrunken log2 fold changes (LFC) and SE to a results table from DESeq run without LFC shrinkage.
# iwalk(DEGres, ~ write.csv(.x, glue("DESeq2/{.y}.DEG.csv")))
# #Apply a function to each element of a vector, and its index
saveRDS(dds, "DESeq2/dds.rds")#增加
# ##########################################

#plotMA(res)
res.dat <- data.frame(res)
res.up <- res.dat[which(res.dat$log2FoldChange> 1 & res.dat$padj < 0.01 ),]
res.up <- res.up[order(res.up$log2FoldChange, -log10(res.up$pvalue), decreasing = T),]
res.down <- res.dat[which(res.dat$log2FoldChange <(-1) & res.dat$padj < 0.01),]
res.down <- res.down[order(abs(res.down$log2FoldChange), -log10(res.down$pvalue), decreasing = T),]
write.table(res.up, "DESeq2/res_DEG_up.txt", row.names = F, quote = F, sep = "\t")
write.table(res.down, "DESeq2/res_DEG_down.txt", row.names = F, quote = F, sep = "\t")
