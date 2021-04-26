import os
import rnaseq.getmyconfig as getmyconfig

stringtie = getmyconfig('Transcript','stringtie').strip("'")

def assembly(samples,bams,gtf,out_dir,threads=8):
    dic = {}
    gtfs = []
    for index,sample in enumerate(samples):
        bam = bams[index]
        out = os.path.join(out_dir, sample)
        ### transcripts prediction for each sample
        cmd_run = f"{stringtie} {bam} -p {threads} -G {gtf} -l {sample} -o {out}.gtf".format(**locals())
        gtfs.append("{out}.gtf".format(out=out))
        ### Estimate transcript abundances and create table counts for Ballgown
        cmd_count = f"{stringtie} {bam} -e -b {out_dir} -p {threads} -G {out_dir}/merged.gtf -o {out}.count.gtf -A " \
                    f"{out}.abundance.txt ".format(**locals())
        dic[sample] = [cmd_run,cmd_count]
    gtf_list = ' '.join(gtfs)
    cmd_merge = f"{stringtie} --merge -p {threads} -G {gtf} -o {out_dir}/merged.gtf {gtf_list}".format(**locals())

    dic['merge'] = cmd_merge

    return dic