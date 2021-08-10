# Title     : identify DEGs for bio repeat samples (n>=4)
# Created by: stl23
# Created on: 2021/8/5

library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

args <- commandArgs(trailingOnly = TRUE)
samples <- args[1]
groups <- args[2]
inpath <- args[3]
outpath <- args[4]
prefix <- args[5]
if (!dir.exists(paths = inpath) & !dir.exists(paths = outpath)) {
      stop("Directory provided does not exist")
    }
out_gene_csv <- paste0(outpath,'/',prefix,'.RNAseq_gene_results.csv')
out_transcript_csv <- paste0(outpath,'/',prefix,'.RNAseq_transcript_results.csv')
out_gene_tsv <- paste0(outpath,'/',prefix,'.RNAseq_gene_results.tsv')
out_transcript_tsv <- paste0(outpath,'/',prefix,'.RNAseq_transcript_results.tsv')
out_DE_transcript <- paste0(outpath,'/',prefix,'.RNAseq_different_expression_transcripts_results.tsv')
out_DE_gene <- paste0(outpath,'/',prefix,'.RNAseq_different_expression_genes_results.tsv')

sample_name <- unlist(strsplit(samples,','))
group_name <- unlist(strsplit(groups,','))
pheno_data <- data.frame(ids=sample_name,type=group_name)  ## load from a file if prefer to
bg <- ballgown(dataDir = inpath, samplePattern = "RNAseq", pData=pheno_data)

# Filter low-abundance genes
bg_filt <- subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Identify signficant differently expressed Transcripts/genes
results_transcripts <- stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes <- stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM") ## covariate is the group type

## Add names
results_transcripts <- data.frame(transcriptNames=ballgown::transcriptNames(bg_filt),transcriptIDs=ballgown::transcriptIDs(bg_filt),geneNames=ballgown::geneNames(bg_filt),geneIDs=ballgown::geneIDs(bg_filt),results_transcripts)
tmp <- unique(results_transcripts[,c("geneNames", "geneIDs")])
results_genes <- merge(results_genes, tmp, by.x=c("id"), by.y=c("geneIDs"), all.x=TRUE)

# Sort from the smallest P value to largest
results_transcripts <- arrange(results_transcripts,pval)
results_genes <- arrange(results_genes,pval)

# Output as CSV
write.csv(results_transcripts,file = out_transcript_csv,row.names=FALSE)
write.csv(results_genes,file = out_gene_csv,row.names=FALSE)

# Output as TSV
write.table(results_transcripts,file = out_transcript_tsv,sep="\t", quote=FALSE, row.names=FALSE)
write.table(results_genes, file = out_gene_tsv,sep="\t", quote=FALSE, row.names=FALSE)
# Identify genes with q value < 0.05
DE_transcripts <- subset(results_transcripts,results_transcripts$qval<0.05)
DE_genes <- subset(results_genes,results_genesqval<0.05)
write.table(DE_transcripts,file = out_DE_transcript,sep="\t", quote=FALSE, row.names=FALSE)
write.table(DE_genes,file = out_DE_gene,sep="\t", quote=FALSE, row.names=FALSE)
