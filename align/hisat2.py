#!/usr/bin/env python

from utils import getConfig

#hisat2 = getConfig('Transcript','hisat2')
#samtools = getConfig('Transcript','samtools')

def align(fq1,fq2,ref,out,threads=2):
    hisat2 = getConfig('Transcript','hisat2')
    samtools = getConfig('Transcript','samtools')
    if fq2:  ## paired-end
        cmd = ("{hisat2} --dta -p {threads} -x {ref} -1 {fq1} -2 {fq2} --summary-file {out}.sum.txt" 
               " | {samtools} sort -@ {threads} -o {out}.bam").format(**locals())
    else:    ## single-end
        cmd = ("{hisat2} --dta -p {threads} -x {ref} -U {fq1} --summary-file {out}.sum.txt" 
               " | {samtools} sort -@ {threads} -o {out}.bam").format(**locals())
    return cmd
