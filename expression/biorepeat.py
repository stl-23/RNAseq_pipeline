import os
from utils import getConfig
preDE = getConfig('Transcript','preDE')
Rscript = getConfig('Transcript','Rscript')

degseq_bio = os.path.join(os.path.dirname(os.path.abspath(__file__)),'degseq_bio.R')
ballgown =  os.path.join(os.path.dirname(os.path.abspath(__file__)),'ballgown.R')

def makeballgown(input_dir,samples,groups,out_dir,prefix):
    ### make ballgown R script (recommonded large sample sizes: n>=4 per group)
    ## Note that Ballgown’s statistical test is a standard linear model-based comparison.
    ## For small sample sizes (n < 4 per group), it is often better to perform regularization.
    ## This can be done using the limma,Deseq2,edgeR packages
    # https://github.com/alyssafrazee/ballgown
    # https://pmbio.org/module-06-rnaseq/0006/03/01/RNAseq_Differential_Expression/

    sample_names = ','.join(samples)
    group_names = ','.join(groups)
    cmd = "{Rscript} {ballgown} {samples} {groups} {inpath} {outpath} {prefix}".format(Rscript=Rscript,
                                                                              ballgown=ballgown,
                                                                              samples=sample_names,
                                                                              groups=group_names,
                                                                              inpath=input_dir,
                                                                              outpath=out_dir,
                                                                              prefix=prefix)
    return cmd

def makedeseq2(input_dir,samples,groups,out_dir,prefix):
    ### make deseq2 R script (allowed small sample sizes: n < 4 per group)
    group_names = ','.join(groups)
    #out_transcript_csv = os.path.join(out_dir, prefix+".RNAseq_transcript_results.csv")
    #out_gene_csv = os.path.join(out_dir, prefix+".RNAseq_gene_results.csv")
    #out_DE_transcript = os.path.join(out_dir,prefix+".RNAseq_different_expression_transcripts_results.tsv")
    #out_DE_gene = os.path.join(out_dir,prefix+".RNAseq_different_expression_genes_results.tsv")
    #out_up_gene = os.path.join(out_dir,prefix+".RNAseq_UP_genes_results.tsv")
    #out_down_gene = os.path.join(out_dir,prefix+".RNAseq_DOWN_genes_results.tsv")

    #gtf = [os.path.join(input_dir,sample,sample+'.out.gtf') for sample in samples]
    #with open(os.path.join(out_dir,'sample.list'),'w') as fh:
    #    for index,sample in enumerate(samples):
    #        fh.write(sample+' '+gtf[index]+'\n')
    #cmd_make_matrix = "python {preDE} -i {out_dir}/sample.list -l {read_len}" \
    #                  " -g {out_dir}/{prefix}.gene_count_matrix.csv " \
    #                  "-t {out_dir}/{prefix}.transcript_count_matrix.csv".format(preDE=preDE,out_dir=out_dir,
    #                                                                             prefix=prefix,read_len=read_length)
    #### prepDE.py first column is geneid|gene name, try remove gene names

    cmd_make_deseq2 = "{Rscript} {degseq_bio} {groups} {path} {prefix}".format(Rscript=Rscript,
                                                                              degseq_bio=degseq_bio,
                                                                              groups=group_names,
                                                                              path=out_dir,
                                                                              prefix=prefix)
   
    return cmd_make_deseq2

