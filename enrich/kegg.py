import os
from utils import getConfig
Rscript = getConfig('Transcript','Rscript')
kegg = os.path.join(os.path.dirname(os.path.abspath(__file__)),'kegg.R')

def makeKEGG(file,out_dir,prefix,org='hsa'):
    cmd = "{Rscript} {kegg} {file} {out} {prefix} {org}".format(Rscript=Rscript,
                                                                 kegg=kegg,
                                                                 file=file,
                                                                 out=out_dir,
                                                                 prefix=prefix,
                                                                 org=org)
    return cmd
