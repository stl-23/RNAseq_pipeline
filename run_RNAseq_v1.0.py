#!/usr/bin/env python
import os
import argparse
import glob
import utils
import subprocess

from align import hisat2
from splicing import splicing
from assembly import stringtie
from expression import (count,biorepeat,nobiorepeat,volcano)
from enrich import (go,kegg)

def checksamples(inputs_dir,samples_dic):
    new_dic = {}
    paired = True
    for group,samples in samples_dic.items():
        for sample in samples:
            sample_f_list = glob.glob(os.path.join(os.path.abspath(inputs_dir).rstrip('/'),sample+'_1*'))
            sample_r_list = glob.glob(os.path.join(os.path.abspath(inputs_dir).rstrip('/'),sample+'_2*'))
            if sample_f_list and sample_r_list:
                sample_f = sample_f_list[0]
                sample_r = sample_r_list[0]
                if group in new_dic:
                    new_dic[group].append(sample+':'+sample_f+';'+sample_r)
                else:
                    new_dic[group] = [sample+':'+sample_f+';'+sample_r]
            elif sample_f_list and not sample_r_list:
                paired = False
                sample_f = sample_f_list[0]
                if group in new_dic:
                    new_dic[group].append(sample + ':' + sample_f)
                else:
                    new_dic[group] = [sample + ':' + sample_f]
            else:
                raise Exception('Can not find files,please check and name files as samplexx_1/2.fq/.fastq(.gz)')

    return new_dic,paired

def mapping(dic,ref_prefix,out_path):
    map_dic = {}
    for group, samples in dic.items():
        for sample in samples:
            sample_name = sample.split(':')[0]
            sample_list = sample.split(':')[1].split(';')
            out_name = os.path.join(out_path, sample_name)
            if len(sample_list) == 2: # paired end
                sample_f,sample_r = sample_list[0],sample_list[1]
                map_dic[sample_name] = hisat2.align(sample_f,sample_r,ref_prefix,out_name,mthreads)
            elif len(sample_list) == 1: # single end
                sample_f = sample_list[0]
                map_dic[sample_name] = hisat2.align(sample_f,'',ref_prefix,out_name,mthreads)
            else:
                raise Exception('Load no input files')

    return map_dic

def obtainDEG(groups,samples_dic,DEG_path):
    DEG_dic = {}
    for group in groups:
        A,B = group.split(':')[0],group.split(':')[1]
        two_group_samples = samples_dic[A] + samples_dic[B]
        two_group_names = [A] * len(samples_dic[A]) + [B] * len(samples_dic[B])
        two_group_out_path = os.path.join(DEG_path,A + '_vs_' + B)
        utils.makedir(two_group_out_path)
        #DE_file = os.path.join(two_group_out_path, A + '_vs_' + B + ".RNAseq_different_expression_genes_results.tsv")
        if len(samples_dic[A]) < 2 or len(samples_dic[B]) < 2: ## no bio repeats
            nobio_cmds = nobiorepeat.makeDEscript(two_group_samples,two_group_names,two_group_out_path,A+'_vs_'+B)
            volcano_cmd = volcano.volcano(A, B, two_group_out_path, False)
            DEG_dic[A + '_vs_' + B] = [nobio_cmds, volcano_cmd]
        elif len(samples_dic[A]) >= 2 and len(samples_dic[B]) >= 2: ## bio repeats
            bio_cmd = biorepeat.makedeseq2(assembly_path,two_group_samples,two_group_names,two_group_out_path,A+'_vs_'+B)
            volcano_cmd = volcano.volcano(A, B, two_group_out_path, True)
            DEG_dic[A + '_vs_' + B] = [bio_cmd, volcano_cmd]

    return DEG_dic

def obtainAS(groups,samples_dic,map_result_path,AS_path,gtf,paired_or_single,readlength=50,theads=8):
    AS_dic = {}
    for group in groups:
        A,B = group.split(':')[0],group.split(':')[1]
        two_group_out_path = os.path.join(AS_path,A + '_vs_' + B)
        temp_path = os.path.join(two_group_out_path,'temp/')
        utils.makedir(two_group_out_path)
        utils.makedir(temp_path)
        A_paths = [os.path.join(map_result_path,i+'.bam') for i in samples_dic[A]]
        B_paths = [os.path.join(map_result_path,i+'.bam') for i in samples_dic[B]]
        with open(os.path.join(temp_path,'b1.txt'),'w') as fa, open(os.path.join(temp_path,'b2.txt'),'w') as fb:
            fa.write(','.join(A_paths))
            fb.write(','.join(B_paths))
        as_cmd = splicing.alternative_splicing(two_group_out_path,gtf,paired_or_single,readlength,theads)
        AS_dic[A + '_vs_' + B] = as_cmd

    return AS_dic

def obtainEnrich(groups,DEG_path,org):
    Enrich_dic = {}
    for group in groups:
        A,B = group.split(':')[0],group.split(':')[1]
        two_group_out_path = os.path.join(DEG_path,A + '_vs_' + B)
        DE_file = os.path.join(two_group_out_path, A + '_vs_' + B + ".RNAseq_different_expression_genes_results.tsv")
        go_cmd = go.makeGO(DE_file,two_group_out_path,A+'_vs_'+B,org)
        kegg_cmd = kegg.makeKEGG(DE_file, two_group_out_path, A + '_vs_' + B, org)
        Enrich_dic[A+'_vs_'+B] = [go_cmd,kegg_cmd]
     
    return Enrich_dic

if __name__ == '__main__':
    examplelog = """EXAMPLES:
    python3 
    """
    parser = argparse.ArgumentParser(description='RNAseq pipeline v1.0',
                                     epilog=examplelog,
                                     add_help=False)

    general = parser.add_argument_group(title='General options')
    general.add_argument('-h', '--help', action="help",
                         help="show the help and exit")
    general.add_argument('-i', '--input', type=str,
                         help="The input directory of clean reads")
    general.add_argument('-o', '--output', type=str,
                         help="The output directory")
    general.add_argument('-r', '--reference', type=str,
                         help="The fasta file of reference")
    general.add_argument('-g', '--gtf', type=str,default='',
                         help="The gene annotation file of reference,gtf format")
    general.add_argument('-v','--build_version',type=str,default='',choices=['hg19','hg38'],
                          help='Human genome build version,if used,do not set -r and -g')
    general.add_argument('--samples',type=str,
                         help="Sample names. The sample names should be accordant with that in input directory"
                              "e.g.  sample1,sample2,sample3:sample4,sample5,sample6:sample7..."
                              "a group of samples should be seperated by comma")
    general.add_argument('--groups',type=str,
                         help="Group names. Should be consistent with the order of samples names."
                              "e.g. group1:group2:group3...")
    general.add_argument('--compare',type=str,
                         help="DEG(Differentially Expressed Genes) group pairs. group1:group2,group1:group3...")
    general.add_argument('--merge', action='store_true',
                         help="Merge novel transcripts with reference gtf file if set True,otherwise use reference gtf only")
    general.add_argument('--readlength',type=int,default=50,
                         help="average read length for alternative splicing by rMATs.")
    general.add_argument('--script',action='store_true',
                         help='Only generate shell scripts')
    general.add_argument('--jobs',type=int,default=1,
                         help='The maximum jobs when run in local at the same time')
    general.add_argument('--mthreads',type=int,default=5,
                         help="Maximum Threads")
    enrich = parser.add_argument_group(title='Enrichment options')
    enrich.add_argument('--org',type=str,default='hsa',
                        help="Species for GO/KEGG enrichment")

    args = parser.parse_args()
    inputs_dir = args.input
    outputs_dir = args.output
    ref = args.reference
    gtf = args.gtf
    buildver = args.build_version
    samples = args.samples
    groups = args.groups
    compare = args.compare
    script = args.script
    mthreads = args.mthreads
    merge = args.merge
    org = args.org
    jobs = args.jobs
    readlength = args.readlength
### parse parameters ###
    index_shell = os.path.join(os.path.dirname(os.path.abspath(__file__)),'align','index.sh')
    statistics_shell = os.path.join(os.path.dirname(os.path.abspath(__file__)),'align','statistics.sh')
    softwares_list = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'softwares.config')
    ref_prefix = ''
    if buildver and not (ref or gtf):
        if buildver == 'hg19':
            genomicsdb = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'genomicsdb','hg19/')
            ref = os.path.join(genomicsdb,'hg19.fa')
            gtf = os.path.join(genomicsdb,'hg19.gtf')
            ref_prefix = os.path.splitext(ref)[0]
        elif buildver == 'hg38':
            genomicsdb = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'genomicsdb','hg38/')
            ref = os.path.join(genomicsdb,'hg38.fa')
            ref_prefix = os.path.splitext(ref)[0]
            gtf = os.path.join(genomicsdb,'hg38.gtf')
    elif not buildver and (ref and gtf):
        #ref_prefix = os.path.splitext(os.path.basename(os.path.abspath(ref)))[0]
        ref_prefix = os.path.splitext(os.path.abspath(ref))[0]
        #if maptool == 'hisat2':
        hisat2_suffix = ['.1.ht2','.2.ht2','.3.ht2','.4.ht2','.5.ht2','.6.ht2','.7.ht2','.8.ht2']
        if not all([os.path.exists(ref_prefix+i) for i in hisat2_suffix]):
            subprocess.check_call(['sh', index_shell, ref, 'hisat2', ref_prefix,softwares_list])
            #subprocess.check_call(['sh', statistics_shell, ref, ref_prefix,softwares_list])
        ref = ref_prefix+'.aline.fa'
        #gtf = ref_prefix+'.gtf'

  ## {group1:[sample1,sample2,sample3],group2:[sample4,sample5,sample6]}
    samples_dic = { i[0]:i[1].split(",") for i in zip(groups.strip().split(":"),samples.strip().split(":")) }
  ## main work
    map_result_path = os.path.join(os.path.abspath(outputs_dir), '01.mapping_results')
    assembly_path = os.path.join(os.path.abspath(outputs_dir), '02.quantitation_results')
    AS_path = os.path.join(os.path.abspath(outputs_dir), '03.alternative_splicing_results')
    DEG_path = os.path.join(os.path.abspath(outputs_dir), '04.DEG_results')
    
    utils.makedir(map_result_path)
    utils.makedir(assembly_path)
    utils.makedir(AS_path)
    utils.makedir(DEG_path)
    # mapping
    map_dic = {}
    if samples_dic:
        full_path_sample_dic,paired = checksamples(inputs_dir, samples_dic)
        map_dic = mapping(full_path_sample_dic,ref_prefix,map_result_path)
    if script:
        for sample in map_dic:
            utils.out_cmd('s1_'+sample+'.sh', map_dic[sample])
    else:
        print("Step1:Mapping...")
        utils.multi_run(utils.run_shell_cmd,list(map_dic.values()),jobs)

    # assemble and quantitation by stringtie
    sample_names = map_dic.keys()
    bams = [os.path.join(map_result_path, sample_name + '.bam') for sample_name in sample_names]
    assemble_quantity_dic = stringtie.assemblyandquantitation(sample_names,bams,gtf,assembly_path,mthreads,merge)

    if script:
        if 'merge' in assemble_quantity_dic:
            assemble_quantity_dic_tmp = assemble_quantity_dic
            merge_gtf_cmd = assemble_quantity_dic_tmp.pop('merge')
            for sample in assemble_quantity_dic_tmp:
                utils.out_cmd('s2.1_'+sample+'.assemble.sh', assemble_quantity_dic_tmp[sample][0])
                utils.out_cmd('s2.2_'+sample+'.count.sh', assemble_quantity_dic_tmp[sample][1])
            #merge_gtf_cmd = assemble_quantity_dic['merge']
            utils.out_cmd('s2.1.5_merge_gtf.sh', merge_gtf_cmd)
        else:
            for sample in assemble_quantity_dic:
                utils.out_cmd('s2_' + sample + '.count.sh', assemble_quantity_dic[sample][0])
    else:
        if 'merge' in assemble_quantity_dic:
            print("Transcripts quantitation...")
            ## predict novo transcripts
            new_trans_cmd = [i[0] for i in assemble_quantity_dic.values()]
            utils.multi_run(utils.run_shell_cmd, new_trans_cmd, jobs)
            ## merge novo and ref transcripts
            utils.multi_run(utils.run_shell_cmd, assemble_quantity_dic['merge'], jobs)
            ## generate tab
            tab_cmd = [i[1] for i in assemble_quantity_dic.values()]
            utils.multi_run(utils.run_shell_cmd, tab_cmd, jobs)
        else:
            ## generate tab
            print("Step2:Transcripts quantitation...")
            utils.multi_run(utils.run_shell_cmd, list(assemble_quantity_dic.values()), jobs)

    # read count
    if merge:
        merged_gtf_file = os.path.join(assembly_path, 'merged.gtf')
        count_matrix_cmd = count.makefeaturecount(bams, merged_gtf_file, DEG_path, paired, mthreads)
    else:
        count_matrix_cmd = count.makefeaturecount(bams, gtf, DEG_path, paired, mthreads)
    if script:
        utils.out_cmd('s3.1_all_count_matrix.sh', count_matrix_cmd)

    else:
        print("Step3.1:Read count table generating...")
        utils.multi_run(utils.run_shell_cmd, count_matrix_cmd, jobs)
    ## FPKM and TPM
    fpkm_tpm_cmd = count.fpkm_tpm_sum(samples,groups,assembly_path,DEG_path)
    if script:
        utils.out_cmd('s3.2_all_fpkm_tpm.sh', fpkm_tpm_cmd)
    else:
        print("Step3.2:FPKM/TPM")
        utils.multi_run(utils.run_shell_cmd,fpkm_tpm_cmd, 1)
    

    ## get differentially alternative splicing
    compare_groups = compare.split(',')
    if paired:
        paired_or_single = "paired"
    else:
        paired_or_single = "single"
    AS_dic = obtainAS(compare_groups,samples_dic,map_result_path,AS_path,gtf,paired_or_single,readlength,mthreads)
    if script:
        for group in AS_dic:
            utils.out_cmd('s4.1_'+group+'.AS.sh', AS_dic[group])
    else:
        print("Step4: %s differentially alternative splicing identifying..." %(','.join(AS_dic.keys())))
        as_cmds = []
        for group in AS_dic:
            as_cmds.append(AS_dic[group])
        utils.multi_run(utils.run_shell_cmd,as_cmds,jobs)
    # get DEGs
    DEG_dic = obtainDEG(compare_groups,samples_dic,DEG_path)

    if script:
        for group in DEG_dic:
            A,B = group.split('_vs_')[0],group.split('_vs_')[1]
            if len(samples_dic[A]) < 2 or len(samples_dic[B]) < 2: ## no bio repeats
                utils.out_cmd('s5.1_'+group+'.DEG.sh', DEG_dic[group][0])
                utils.out_cmd('s5.2_'+ group+'.volcano.sh', DEG_dic[group][1])
            elif len(samples_dic[A]) >= 2 and len(samples_dic[B]) >= 2: ## bio repeats
                utils.out_cmd('s5.1_'+group+'.DEG.sh', DEG_dic[group][0])
                utils.out_cmd('s5.2_'+group+'.volcano.sh', DEG_dic[group][1])
    else:
        print("Step5: %s differently expressed genes(DEGs) identifying..." %(','.join(DEG_dic.keys())))
        cal_matrix_cmds = []
        DEG_cmds = []
        volcano_cmds = []
        for group in DEG_dic:
            A,B = group.split('_vs_')[0],group.split('_vs_')[1]
            if len(samples_dic[A]) < 2 or len(samples_dic[B]) < 2: ## no bio repeats
                DEG_cmds.append(DEG_dic[group][0])
                volcano_cmds.append(DEG_dic[group][1])
            elif len(samples_dic[A]) >= 2 and len(samples_dic[B]) >= 2: ## bio repeats
                #cal_matrix_cmds.append(DEG_dic[group][0])
                DEG_cmds.append(DEG_dic[group][0])
                volcano_cmds.append(DEG_dic[group][1])
        #if cal_matrix_cmds:
        #    utils.multi_run(utils.run_shell_cmd, cal_matrix_cmds, jobs)
        utils.multi_run(utils.run_shell_cmd, DEG_cmds, jobs)
        utils.multi_run(utils.run_shell_cmd, volcano_cmds, jobs)
    
    # GO and KEGG enrichment      
    Enrich_dic = obtainEnrich(compare_groups,DEG_path,org)
    if script:
        for group in Enrich_dic:
            utils.out_cmd('s6_'+group+'.All.GO.sh',Enrich_dic[group][0])
            utils.out_cmd('s6_'+group+'.All.KEGG.sh',Enrich_dic[group][1])
    else:
        print("Step6:Enrichment analysis...")
        go_cmds = [i[0] for i in Enrich_dic.values()]
        kegg_cmds = [i[1] for i in Enrich_dic.values()]
        utils.multi_run(utils.run_shell_cmd, go_cmds, jobs)
        utils.multi_run(utils.run_shell_cmd, kegg_cmds, jobs)

