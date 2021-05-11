import os
import rnaseq.utils as utils
import rnaseq.getmyconfig as getmyconfig

featurecount = getmyconfig('Transcript','featurecount').strip("'")

def makefeaturecount(bams,gtf,out_dir,level='exon',attrType='gene_id',threads=8,ref_only=False,paired=True,
                     other_parameters=None):
    bam = ' '.join(bams)
    if not ref_only:
        cmd_count = f"""{featurecount} -T {threads} -a {out_dir}/merged.gtf -o {out_dir}/All.read.count.txt
        -t {level} -g {attrType} {bam}""".format(**locals())
    else:
        cmd_count = f"""{featurecount} -T {threads} -a {gtf} -o {out_dir}/All.read.count.txt
        -t {level} -g {attrType} {bam}""".format(**locals())

    if paired:
        cmd_count += " -p"

    if other_parameters:
        cmd_count += other_parameters

    return cmd_count


def makeDEscript(samples,groups,out_dir,prefix):
    sample_name = ','.join(['"'+sample+'"' for sample in samples])
    group_name = ','.join(['"'+group+'"' for group in groups])
    #out_transcript_csv = os.path.join(out_dir, prefix+".RNAseq_transcript_results.csv")
    out_gene = os.path.join(out_dir, prefix)
    #out_DE_transcript = os.path.join(out_dir,prefix+".RNAseq_different_expression_transcripts_results.tsv")
    out_DE_gene = os.path.join(out_dir,prefix+".RNAseq_different_expression_genes_results.tsv")
    utils.makedir(out_gene)

    cmd_identify = f"""library(edgeR)
library(DEGseq)
## Load data (Geneid Chr Start End Strand Length file_name1 file_name2 ...)
genecountData <- read.table("{out_dir}/All.read.count.txt", sep="\t", header=TRUE,skp=1))
## Data format (Geneid sample1 sample2 ...)
data <- genecountData[,-2:-6]
names(data[,-1]) <- c({sample_name}) ## change file names to sample names
data <- as.matrix(data)
rownames(data) <- data$Geneid
## Filter
#new_data <- new_data[rowMeans(new_data)>1,] ## Filter low-abundance genes

group <- c({group_name})
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
DEGexp(geneExpMatrix1=geneExpMatrix1, geneCol1=1, expCol1=c(2), groupLabel1={group_name}[1],
geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(3), groupLabel2={group_name}[2],method="LRT",
pValue=0.005,qValue=0.005,outputDir=out_gene)
""".format(**locals())
    cmd_filter = """
## Filter DE genes by |log2FC| >1 && qvalue < 0.005
DE_data <- read.table("%s/output_score.txt", header=TRUE)
names(DE_data) <- c("GeneNames","%s","log2FoldChange","log2FoldChange normalized","p-value","q-value_Benjamini_1995","q-value_Storey_2003","Signature_p-value_0.005")
if (DE_data$ignature_p-value_0.005 == "TRUE"){
    if (DE_data$value1 > DE_data$value2)){
    DE_data$log2FoldChange [is.na(DE_data$log2FoldChange)] <- Inf
    }
    else {
    DE_data$log2FoldChange [is.na(DE_data$log2FoldChange)] <- -Inf
    }
DE_data_sig <- DE_data[DE_data$q-value_Storey_2003. < 0.005 & (DE_data$log2FoldChange >1 | DE_data$log2FoldChange < -1)]
write.table(DE_data_sig, "%s",sep="\t", quote=FALSE, row.names=FALSE)
"""%(out_gene,sample_name, out_DE_gene)
    return cmd_identify, cmd_filter