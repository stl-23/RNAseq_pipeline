import os
from utils import getConfig,makedir

stringtie = getConfig('Transcript','stringtie')

def assemblyandquantitation(samples,bams,gtf,out_dir,threads=8,merge=True):
    dic = {}
    gtfs = []
    for index,sample in enumerate(samples):
        bam = bams[index]
        out = os.path.join(out_dir, sample)
        out_name = os.path.join(out, sample)
        makedir(out)
        if merge: ## reference and assembled transcripts
            ### transcripts prediction for each sample
            cmd_run = f"{stringtie} {bam} -p {threads} -G {gtf} -l {sample} -o {out_name}.gtf".format(**locals())
            gtfs.append("{out}.gtf".format(out=out_name))
            ### Estimate merged transcript abundances and create table counts for Ballgown or Deseq2
            cmd_count = f"{stringtie} {bam} -e -B -p {threads} -G {out_dir}/merged.gtf -o {out_name}.out.gtf " \
                    f"-A {out_name}.abundance.tab ".format(**locals())
            dic[sample] = [cmd_run,cmd_count]
        else:     ## reference only
            cmd_count = f"{stringtie} {bam} -e -B -p {threads} -G {gtf} -l {sample} -o {out_name}.out.gtf " \
                        f"-A {out_name}.abundance.tab".format(**locals())
            dic[sample] = [cmd_count]
    if merge:
        ### merge reference and assembly gtf files
        gtf_files = ' '.join(gtfs)
        cmd_merge = f"{stringtie} --merge -p {threads} -G {gtf} -o {out_dir}/merged.gtf {gtf_files}".format(**locals())
        dic['merge'] = cmd_merge
    return dic
