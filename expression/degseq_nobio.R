# Title     : Identify DEGs for no bio repeat samples
# Created by: stl23
# Created on: 2021/8/5

library(edgeR)
library(DEGseq)
args <- commandArgs(trailingOnly = TRUE)
samples <- args[1]
groups <- args[2]
path <- args[3]
prefix <- args[4]
if (!dir.exists(paths = path)) {
      stop("Directory provided does not exist")
    }
input_file <- file.path(path,'../All.read.count.for_r.txt')
sample_name <- unlist(strsplit(samples,','))
group_name <- unlist(strsplit(groups,','))
group_compare <- group_name[!duplicated(group_name)]
sample1 <- sample_name[which(group_name %in% group_compare[1])]
sample2 <- sample_name[which(group_name %in% group_compare[2])]
out_gene <- paste0(path,'/',prefix)
out_DE_gene <- paste0(path,'/',prefix,'.RNAseq_different_expression_genes_results.tsv')
out_up_gene <-  paste0(path,'/',prefix,'.RNAseq_UP_genes_results.tsv')
out_down_gene <- paste0(path,'/',prefix, '.RNAseq_DOWN_genes_results.tsv')

## Load data (Geneid Chr Start End Strand Length file_name1 file_name2 ...)
genecountData <- read.table(file = input_file, sep="\t", header=TRUE,row.names=1,check.names=FALSE) ## https://zhuanlan.zhihu.com/p/338940362
group <- group_name
genelist <- DGEList(counts=genecountData, group=group)
#genelist <- DGEList(counts=genecountData[2:length(colnames(genecountData))], group=group)
## filter
keep <- rowSums(cpm(genelist)>1) >= 1   # at least 1 sample cpm > 1
genelist <- genelist[keep, , keep.lib.sizes=FALSE]
## TMM normalisation
TMM <- calcNormFactors(genelist, method="TMM")

## predict DE genes by edgeR
#bcv <- 0.4  # human:0.4; genetically identical model organisms:0.1 ; technical replicates: 0.01
#et <- exactTest(TMM, dispersion = bcv ^ 2)
#DEgenes <- decideTestsDGE(et, p.value = 0.05, lfc = 0)
## Identify differently expressed genes by DEGseq
TMM <- estimateCommonDisp(TMM)
pseudo_data <- TMM$counts
pseudo_data2 <- cbind(rownames(pseudo_data),apply(pseudo_data,2,as.character))
geneExpMatrix1 <- cbind(rownames(pseudo_data),pseudo_data2[,sample1])
geneExpMatrix2 <- cbind(rownames(pseudo_data),pseudo_data2[,sample2])
colnames(geneExpMatrix1) <- c('GeneID',sample1)
colnames(geneExpMatrix2) <- c('GeneID',sample2)
exp1 <- length(sample1)+1
exp2 <- length(sample2)+1
DE_data <- DEGexp(geneExpMatrix1=geneExpMatrix1, geneCol1=1, expCol1=c(2:exp1), groupLabel1=group_compare[1],
geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2:exp2), groupLabel2=group_compare[2],method="LRT",
pValue=0.005,qValue=0.005,outputDir=path)

## Filter DE genes by |log2FC| >1 && qvalue < 0.005
DE_data <- read.table(paste0(path,'/output_score.txt'), header=TRUE)
names(DE_data) <- c("GeneID",group_compare[1],group_compare[2],"log2FoldChange","log2FoldChange_normalized","p_value","q_value_Benjamini_1995","q_value_Storey_2003","Signature_p_value_0.005")

#if (DE_data$Signature_p_value_0.005 == "TRUE"){
#    if (DE_data[,2] > DE_data[,3]){
#    DE_data$log2_FC [is.na(DE_data$log2_FC)] <- 'Inf'
#    }
#    else {
#    DE_data$log2_FC [is.na(DE_data$log2_FC)] <- '-Inf'
#    }
#}
## log2Fold_Change(A/B) if A is 0(B not 0) then replace 'NA' to '-Inf'; else if B is 0 then replae 'NA' to 'Inf'
DE_data$log2FoldChange <- ifelse(DE_data[,2] == 0,
                          DE_data$log2FoldChange[is.na(DE_data$log2FoldChange)] <- '-Inf',
                          ifelse(DE_data[,3] == 0 ,DE_data$log2_FC[is.na(DE_data$log2FoldChange)] <- 'Inf',DE_data$log2FoldChange)
                          )
DE_data_sig <- subset(DE_data,Signature_p_value_0.005 == "TRUE" & (log2FoldChange >1 | log2FoldChange < (-1)))
DE_data_sig_up <- subset(DE_data,Signature_p_value_0.005 == "TRUE" & log2FoldChange >1)
DE_data_sig_down <- subset(DE_data,Signature_p_value_0.005 == "TRUE" & log2FoldChange < (-1))
write.table(DE_data_sig, file = out_DE_gene,sep="\t", quote=FALSE, row.names=FALSE)
write.table(DE_data_sig_up, file = out_up_gene,sep="\t", quote=FALSE, row.names=FALSE)
write.table(DE_data_sig_down, file = out_down_gene,sep="\t", quote=FALSE, row.names=FALSE)
