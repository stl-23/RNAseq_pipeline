
import os
import rnaseq.getmyconfig as getmyconfig
preDE = getmyconfig('Transcript','preDE').strip("'")


def makeballgown(input_dir,samples,groups,out_dir,prefix):
    ### make ballgown R script (recommonded large sample sizes: n>=4 per group)
    ## Note that Ballgown’s statistical test is a standard linear model-based comparison.
    ## For small sample sizes (n < 4 per group), it is often better to perform regularization.
    ## This can be done using the limma,Deseq2,edgeR packages
    # https://github.com/alyssafrazee/ballgown
    # https://pmbio.org/module-06-rnaseq/0006/03/01/RNAseq_Differential_Expression/

    sample_name = ','.join(['"'+sample+'"' for sample in samples])
    group_name = ','.join(['"'+group+'"' for group in groups])
    out_transcript_csv = os.path.join(out_dir, prefix+".RNAseq_transcript_results.csv")
    out_gene_csv = os.path.join(out_dir, prefix+".RNAseq_gene_results.csv")
    out_transcript_tsv = os.path.join(out_dir,prefix+".RNAseq_transcript_results.tsv")
    out_gene_tsv = os.path.join(out_dir,prefix+".RNAseq_gene_results.tsv")
    out_DE_transcript = os.path.join(out_dir,prefix+".RNAseq_different_expression_transcripts_results.tsv")
    out_DE_gene = os.path.join(out_dir,prefix+".RNAseq_different_expression_genes_results.tsv")

    cmd = f"""library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

pheno_data <- data.frame(ids=c({sample_name}),type=c({group_name}))  ## load from a file if prefer to 
bg <- ballgown(dataDir = "{input_dir}", samplePattern = "RNAseq", pData=pheno_data)

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
write.csv(results_transcripts,"{out_transcript_csv}",row.names=FALSE)
write.csv(results_genes,"{out_gene_csv}",row.names=FALSE)

# Output as TSV
write.table(results_transcripts,"{out_transcript_tsv}",sep="\t", quote=FALSE, row.names=FALSE)
write.table(results_genes,"{out_gene_tsv}",sep="\t", quote=FALSE, row.names=FALSE)
# Identify genes with q value < 0.05
DE_transcripts <- subset(results_transcripts,results_transcripts$qval<0.05)
DE_genes <- subset(results_genes,results_genesqval<0.05)
write.table(DE_transcripts,"{out_DE_transcript}",sep="\t", quote=FALSE, row.names=FALSE)
write.table(DE_genes,"{out_DE_gene}",sep="\t", quote=FALSE, row.names=FALSE)
""".format(**locals())
#    elif biorepeat == 'N':
#        cmd += f"""
# Identify genes with q value < 0.005 and |log2(FC)| >1 
#DE_transcripts < - subset(results_transcripts, results_transcripts$qval<0.005 & (results_transcripts$fc>1 | results_transcripts$fc< -1 ))
#DE_genes <- subset(results_genes,results_genes$qval<0.005 & (results_genes$fc>1 | results_genes$fc< -1))
#write.table(DE_transcripts,"{out_DE_transcript}",sep="\t", quote=FALSE, row.names=FALSE)
#write.table(DE_genes,"{out_DE_genes}",sep="\t", quote=FALSE, row.names=FALSE)
#""".format(**locals())

    return cmd

def makedeseq2(input_dir,samples,groups,out_dir,prefix,read_length):
    ### make deseq2 R script (allowed small sample sizes: n < 4 per group)
    sample_name = ','.join(['"'+sample+'"' for sample in samples])
    group_name = ','.join(['"'+group+'"' for group in groups])
    #out_transcript_csv = os.path.join(out_dir, prefix+".RNAseq_transcript_results.csv")
    out_gene_csv = os.path.join(out_dir, prefix+".RNAseq_gene_results.csv")
    #out_DE_transcript = os.path.join(out_dir,prefix+".RNAseq_different_expression_transcripts_results.tsv")
    out_DE_gene = os.path.join(out_dir,prefix+".RNAseq_different_expression_genes_results.tsv")
    out_up_gene = os.path.join(out_dir,prefix+".RNAseq_UP_genes_results.tsv")
    out_down_gene = os.path.join(out_dir,prefix+".RNAseq_DOWN_genes_results.tsv")

    gtf = [os.path.join(input_dir,sample,sample+'out.gtf') for sample in samples]
    with open(os.path.join(out_dir,'sample.list'),'w') as fh:
        for index,sample in enumerate(samples):
            fh.write(sample+' '+gtf[index]+'\n')
    cmd_make_matrix = "python {preDE} -i {out_dir}/sample.list -l {read_len}" \
                      " -g {out_dir}/{prefix}.gene_count_matrix.csv " \
                      "-t {out_dir}/{prefix}.transcript_count_matrix.csv".format(preDE=preDE,out_dir=out_dir,
                                                                                 prefix=prefix,read_len=read_length)
    #### prepDE.py first column is geneid|gene name, try remove gene names



    cmd_make_deseq2 = f"""library(DESeq2)
# Load gene count matrix
genecountData <- as.matrix(read.csv("{out_dir}/{prefix}.gene_count_matrix.csv"))
colnames(genecountData)[1] <- 'GeneID'
condition <- factor(c({group_name}))  ## order limited: c("Control","Treatment") --> Treatment vs. Control
# Get target matrix
genecolData <- data.frame(row.names=colnames(genecountData), condition)

genecountData <- genecountData[, rownames(genecolData)]

genedds <- DESeqDataSetFromMatrix(countData = genecountData, colData = genecolData, design = ~ condition)
# Identify signficant differently expressed Transcripts/genes
genedataset <- DESeq(genedds)
generes <- results(genedataset)

# Sort by q value
generesordered <- generes[order(generes$padj),]

write.csv(as.data.frame(generesordered),file="{out_gene_csv}")

diff_gene_deseq2 <-subset(generesordere, padj < 0.05)

up_regulated_genes <- subset(generesordere, padj < 0.05 & log2FoldChange > 0)
down_regulated_genes <- subset(generesordere, padj < 0.05 & log2FoldChange < 0)

write.table(diff_gene_deseq2,"{out_DE_gene}",sep="\t", quote=FALSE, row.names=FALSE)
write.table(up_regulated_genes,"{out_up_gene}",sep="\t", quote=FALSE, row.names=FALSE)
write.table(down_regulated_genes,"{out_down_gene}",sep="\t", quote=FALSE, row.names=FALSE)
""".format(**locals())

    return cmd_make_matrix,cmd_make_deseq2

