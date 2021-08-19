import os
from utils import getConfig

featurecount = getConfig('Transcript', 'featurecount')
format_file =  os.path.join(os.path.dirname(os.path.abspath(__file__)),'format_featurecount.py')
def makefeaturecount(bams,gtf,out_dir,paired=True,threads=8,levels='exon',attrType='gene_id',other_parameters=None):
    bam = ' '.join(bams)
    cmd_count = f"""{featurecount} -T {threads} -a {gtf} -o {out_dir}/All.read.count.txt -t {levels} -g {attrType} {bam}""".format(featurecount=featurecount,
    threads=threads,
    gtf=gtf,
    out_dir=out_dir,
    levels=levels,
    attrType=attrType,
    bam=bam)

    if paired:
        cmd_count += " -p"

    if other_parameters:
        cmd_count += other_parameters
    cmd_count = cmd_count+ "\n"+ """python {format_file} {out_dir}/All.read.count.txt {out_dir}/All.read.count.for_r.txt""".format(format_file=format_file,out_dir=out_dir)

    return cmd_count
