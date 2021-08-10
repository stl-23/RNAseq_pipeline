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
input_file <- file.path(path,'All.read.count.txt')
sample_name <- unlist(strsplit(samples,','))
group_name <- unlist(strsplit(groups,','))
out_gene <- paste0(path,'/',prefix)
out_DE_gene <- paste0(path,'/',prefix,'.RNAseq_different_expression_genes_results.tsv')
out_up_gene <-  paste0(path,'/',prefix,'.RNAseq_UP_genes_results.tsv')
out_down_gene <- paste0(path,'/',prefix, '.RNAseq_DOWN_genes_results.tsv')

## Load data (Geneid Chr Start End Strand Length file_name1 file_name2 ...)
genecountData <- read.table(file = input_file, sep="\t", header=TRUE,skp=1)
## Data format (Geneid sample1 sample2 ...)
data <- genecountData[,-2:-6]
names(data[,-1]) <- sample_name ## change file names to sample names
data <- as.matrix(data)
rownames(data) <- data$Geneid
## Filter
#new_data <- new_data[rowMeans(new_data)>1,] ## Filter low-abundance genes

group <- group_name
genelist <- DGEList(counts=data, group=group)
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
pseudo_data <- TMM$pseudo.counts
pseudo_data <- cbind(rownames(data),apply(pseudo_data,2,as.character))
geneExpMatrix1 <- pseudo_data[,1:2]
geneExpMatrix2 <- pseudo_data[,-2]
DE_data <- DEGexp(geneExpMatrix1=geneExpMatrix1, geneCol1=1, expCol1=c(2), groupLabel1=group_name[1],
geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(3), groupLabel2=group_name[2],method="LRT",
pValue=0.005,qValue=0.005,outputDir={out_gene})

## Filter DE genes by |log2FC| >1 && qvalue < 0.005
DE_data <- read.table(paste0(out_gene,output_score.txt), header=TRUE)
names(DE_data) <- c("GeneID",sample_name,"log2FoldChange","log2FoldChange_normalized","p_value","q_value_Benjamini_1995",
"q_value_Storey_2003","Signature_p_value_0.005")
if (DE_data$Signature_p_value_0.005 == "TRUE"){
    if (DE_data$value1 > DE_data$value2){
    DE_data$log2FoldChange [is.na(DE_data$log2FoldChange)] <- 'Inf'
    }
    else {
    DE_data$log2FoldChange [is.na(DE_data$log2FoldChange)] <- '-Inf'
    }
DE_data_sig <- DE_data[DE_data$q_value_Storey_2003 < 0.005 & (DE_data$log2FoldChange >1 | DE_data$log2FoldChange < -1)]
DE_data_sig_up <- DE_data[DE_data$q_value_Storey_2003 < 0.005 & DE_data$log2FoldChange >1]
DE_data_sig_down <- DE_data[DE_data$q_value_Storey_2003 < 0.005 & DE_data$log2FoldChange <-1]
write.table(DE_data_sig, file = out_DE_gene,sep="\t", quote=FALSE, row.names=FALSE)
write.table(DE_data_sig_up, file = out_up_gene,sep="\t", quote=FALSE, row.names=FALSE)
write.table(DE_data_sig_down, file = out_down_gene,sep="\t", quote=FALSE, row.names=FALSE)
