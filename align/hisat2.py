#!/usr/bin/env python

import rnaseq.getmyconfig as getmyconfig

hisat2 = getmyconfig('Transcript','hisat2').strip("'")
samtools = getmyconfig('Transcript','samtools').strip("'")

def align(fq1,fq2,ref,out,threads=2):
    if fq2:
        cmd = ("{hisat2} --dta -p {threads} -x {ref} -1 {fq1} -2 {fq2} --summary-file {out}.sum.txt" 
               " | {samtools} sort -@ {threads} -o {out}.bam").format(**locals())
    else:
        cmd = ("{hisat2} --dta -p {threads} -x {ref} -U {fq1} --summary-file {out}.sum.txt" 
               " | {samtools} sort -@ {threads} -o {out}.bam").format(**locals())
    return cmd
