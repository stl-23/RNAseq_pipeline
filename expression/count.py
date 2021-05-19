import getmyconfig

featurecount = getmyconfig.getConfig('Transcript', 'featurecount').strip("'")
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