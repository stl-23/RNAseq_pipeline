# Title     : Identify DEGs for bio repeat samples (n<4)
# Created by: stl23
# Created on: 2021/8/5

library(DESeq2)
args <- commandArgs(trailingOnly = TRUE)
groups <- args[1]
path <- args[2]
prefix <- args[3]
if (!dir.exists(paths = path)) {
      stop("Directory provided does not exist")
    }
input_file <- file.path(path,'../All.read.count.for_r.txt')
group_name <- unlist(strsplit(groups,','))
out_gene_csv <- paste0(path,'/',prefix,'.RNAseq_gene_results.csv')
out_DE_gene <- paste0(path,'/',prefix,'.RNAseq_different_expression_genes_results.tsv')
out_up_gene <- paste0(path,'/',prefix,'.RNAseq_UP_genes_results.tsv')
out_down_gene <- paste0(path,'/',prefix,'.RNAseq_DOWN_genes_results.tsv')

# Load gene data(Geneid file_name1 file_name2...)
genecountData <- read.table(file = input_file, sep="\t", row.names = 1, header=T)
## Data format (Geneid sample1 sample2 ...)
#genecountData <- genecountData[,-2:-6]
#names(genecountData[,-1]) <- sample_name ## change file names to sample names
#genecountData <- as.matrix(genecountData)
#rownames(genecountData) <- genecountData$Geneid
#colnames(genecountData)[1] <- 'GeneID'
condition <- factor(group_name,levels = unique(group_name))  ## order limited: c("Control","Treatment") --> Treatment vs. Control

# Get target matrix
genecolData <- data.frame(row.names=colnames(genecountData)[1:length(colnames(genecountData))], condition)
genecountData <- genecountData[, rownames(genecolData)]
genedds <- DESeqDataSetFromMatrix(countData = as.matrix(round(genecountData)), colData = genecolData, design = ~ condition)
# Identify signficant differently expressed Transcripts/genes
genedataset <- DESeq(genedds)
generes <- as.data.frame(results(genedataset))
generes <- cbind(rownames(generes),generes)
colnames(generes)[1] <- 'GeneID'

# Sort by q value
generesordered <- generes[order(generes$padj),]

write.csv(as.data.frame(generesordered),file = out_gene_csv)

diff_gene_deseq2 <-subset(generesordered, padj < 0.05)

up_regulated_genes <- subset(generesordered, padj < 0.05 & log2FoldChange > 0)
down_regulated_genes <- subset(generesordered, padj < 0.05 & log2FoldChange < 0)

write.table(diff_gene_deseq2,file = out_DE_gene,sep="\t", quote=FALSE, row.names=FALSE)
write.table(up_regulated_genes,file = out_up_gene,sep="\t", quote=FALSE, row.names=FALSE)
write.table(down_regulated_genes, file = out_down_gene,sep="\t", quote=FALSE, row.names=FALSE)
