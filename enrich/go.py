import os
from utils import getConfig
Rscript = getConfig('Transcript','Rscript')
go = os.path.join(os.path.dirname(os.path.abspath(__file__)),'go.R')
def makeGO(file,out_dir,prefix,org='hsa'):
    ## code from https://blog.csdn.net/sinat_30623997/article/details/79250940?utm_source=blogxgwz1
    if org == 'hsa':
        OrgDb='org.Hs.eg.db'
        cmd = "{Rscript} {go} {file} {out} {prefix} {org}".format(Rscript=Rscript,
                                                              go=go,
                                                              file=file,
                                                              out=out_dir,
                                                              prefix=prefix,
                                                              org=OrgDb)
    return cmd
