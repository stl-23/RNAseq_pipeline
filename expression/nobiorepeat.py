import os
import utils
Rscript = utils.getConfig('Transcript','Rscript')

degseq_nobio = os.path.join(os.path.dirname(os.path.abspath(__file__)),'degseq_nobio.R')
def makeDEscript(samples,groups,out_dir,prefix):
    sample_name = ','.join(samples)
    group_name = ','.join(groups)
    #out_transcript_csv = os.path.join(out_dir, prefix+".RNAseq_transcript_results.csv")
    out_gene = os.path.join(out_dir, prefix)
    #out_DE_transcript = os.path.join(out_dir,prefix+".RNAseq_different_expression_transcripts_results.tsv")
    #out_DE_gene = os.path.join(out_dir,prefix+".RNAseq_different_expression_genes_results.tsv")
    #out_up_gene = os.path.join(out_dir,prefix+".RNAseq_UP_genes_results.tsv")
    #out_down_gene = os.path.join(out_dir,prefix+".RNAseq_DOWN_genes_results.tsv")

    utils.makedir(out_gene)
    cmd = "{Rscript}  {degseq_nobio} {samples} {groups} {path} {prefix}".format(Rscript=Rscript,
                                                                     degseq_nobio=degseq_nobio,
                                                                     samples=sample_name,
                                                                     groups=group_name,
                                                                     path=out_dir,
                                                                     prefix=prefix)
    return cmd
