from utils import getConfig

featurecount = getConfig('Transcript', 'featurecount')
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

    return cmd_count
