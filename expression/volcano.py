import os
import utils
Rscript = utils.getConfig('Transcript','Rscript')

volcano_bio = os.path.join(os.path.dirname(os.path.abspath(__file__)),'volcano_bio.R')
volcano_nobio = os.path.join(os.path.dirname(os.path.abspath(__file__)),'volcano_nobio.R')

def volcano(group1,group2,path,biorepeat=True):
    ## bio replicates: padj < 0.05 && log2FC >0 (up regulated); padj < 0.05 && log2FC < 0 (down regulated)
    ## no bio replicates: padj < 0.005 && log2FC >1 (up regulated); padj < 0.005 && log2FC < -1 (down regulated)
    # code from https://www.yunbios.net/R-language-volcano-mapping.html
    if biorepeat:
        cmd = "{Rscript} {volcano_bio} {path} {group1} {group2}".format(Rscript=Rscript,
                                                                         volcano_bio=volcano_bio,
                                                                         path=path,
                                                                         group1=group1,
                                                                         group2=group2)
    else:
        cmd = "{Rscript} {volcano_nobio} {path} {group1} {group2}".format(Rscript=Rscript,
                                                                         volcano_nobio=volcano_nobio,
                                                                         path=path,
                                                                         group1=group1,
                                                                         group2=group2)

    return cmd
