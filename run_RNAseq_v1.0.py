#!/usr/bin/env python
import os
import argparse
import glob
import utils
import subprocess

from align import hisat2
from assembly import stringtie
from expression import (count,biorepeat,nobiorepeat,volcano)
from enrich import (go,kegg)

def checksamples(inputs_dir,samples_dic):
    new_dic = {}
    paired = True
    for group,samples in samples_dic.items():
        for sample in samples:
            sample_f = glob.glob(os.path.join(os.path.abspath(inputs_dir).rstrip('/'),sample+'_1*'))
            sample_r = glob.glob(os.path.join(os.path.abspath(inputs_dir).rstrip('/'),sample+'_2*'))
            if sample_f and sample_r:
                if group in new_dic:
                    new_dic[group].append(sample+':'+sample_f+';'+sample_r)
                else:
                    new_dic[group] = [sample+':'+sample_f+';'+sample_r]
            elif sample_f and not sample_r:
                paired = False
                if group in new_dic:
                    new_dic[group].append(sample + ':' + sample_f)
                else:
                    new_dic[group] = [sample + ':' + sample_f]
            else:
                print('Can not find files,please check and use names as sample_1/2.fq/.fastq(.gz)')

    return new_dic,paired

def mapping(dic,ref_prefix,out_path):
    map_dic = {}
    for group, samples in dic.items():
        for sample in samples:
            sample_name = sample.split(':')[0]
            sample_f,sample_r = sample.split(':')[1].split(';')
            out_name = os.path.join(out_path, sample_name)
            if sample_f and sample_r: # paired end
                map_dic[sample_name] = hisat2.align(sample_f,sample_r,ref_prefix,out_name,mthreads)
            elif not sample_r: # single end
                map_dic[sample_name] = hisat2.align(sample_f,ref_prefix,out_name,mthreads)

    return map_dic



if __name__ == 'main':
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
    general.add_argument('-s','--samples',type=str,
                         help="Samples and groups names. The sample names should be accordant with that in input directory"
                              "e.g.  group1:sample1,sample2,sample3;group2:sample4,sample5,sample6;..."
                              "a group of samples should be seperated by comma,and groups should be seperated "
                              "by semicolon")
    general.add_argument('-c','--compare',type=str,
                         help="DEG(Differentially Expressed Genes) group pairs. group1:group2;group1:group3...")
    general.add_argument('--readlength',type=int,
                         help="Average read length")
    general.add_argument('--merge',type=str,default=True,choices=[True,False],
                         help="Merge novel transcripts with reference gtf file if choose True,otherwise use reference gtf only")
    general.add_argument('--script',type=str,default=True,choices=[True,False],
                         help='Only generate shell scripts')
    mapping = parser.add_argument_group(title='Mapping options')
    mapping.add_argument('--mthreads',type=int,
                         help="Threads for mapping (hisat2) ")
    assembly = parser.add_argument_group(title='Assembly options')
    assembly.add_argument('--dthreads', type=int,
                         help="Threads for assembly and quantitation (stringtie) ")
    enrich = parser.add_argument_group(title='Enrichment options')
    enrich.add_argument('--go_orgdb',type=str,default='org.Hs.eg.db',
                        help="OrgDb database for GO enrichment")
    enrich.add_argument('--kegg_species',type=str,default='hsa',
                        help="Species database for KEGG enrichment")

    args = parser.parse_args()
    inputs_dir = args.input
    outputs_dir = args.output
    ref = args.reference
    gtf = args.gtf
    buildver = args.build_version
    samples = args.samples
    compare = args.compare
    script = args.script
    mthreads = args.mthreads
    dthreads = args.dthreads
    merge = args.merge
    readlength = args.readlength
    go_orgdb = args.go_orgdb
    kegg_species = args.kegg_species

### parse parameters ###
    index_shell = os.path.join(os.path.abspath(os.path.dirname(__file__)),'align','index.sh')
    statistics_shell = os.path.join(os.path.abspath(os.path.dirname(__file__)),'align','statistics.sh')
    if buildver and not ref and not gtf:
        if buildver == 'hg19':
            genomicsdb = os.path.join(os.path.dirname(os.path.abspath(inputs_dir).rstrip('/')), 'genomicsdb','hg19/')
            ref = os.path.join(genomicsdb,'hg19.fa')
            gtf = os.path.join(genomicsdb,'hg19.gtf')
            ref_prefix = os.path.splitext(ref)[0]
        elif buildver == 'hg38':
            genomicsdb = os.path.join(os.path.dirname(os.path.abspath(inputs_dir).rstrip('/')), 'genomicsdb','hg38/')
            ref = os.path.join(genomicsdb,'hg38.fa')
            ref_prefix = os.path.splitext(ref)[0]
            gtf = os.path.join(genomicsdb,'hg38.gtf')
    elif not buildver and (ref and gtf):
        ref_prefix = os.path.splitext(os.path.basename(os.path.abspath(ref)))[0]
        #if maptool == 'hisat2':
        subprocess.check_call(['sh', index_shell, ref, 'hisat2', ref_prefix])
        subprocess.check_call(['sh', statistics_shell, ref, ref_prefix])
        ref = ref_prefix+'.fa'
        gtf = ref_prefix+'.gtf'

  ## {group1:[sample1,sample2,sample3],group2:[sample4,sample5,sample6]}
    try:
        samples_dic = {i.split(':')[0]:i.split(':')[1].split(',') for i in samples.strip().split(";")}
    except IOError:
        print('Wrong group:sample format')
        exit()
    finally:
        samples_dic = {}

  ## main work
    map_result_path = os.path.join(os.path.abspath(outputs_dir), 'mapping_results')
    assembly_path = os.path.join(os.path.abspath(outputs_dir), 'AQ_results')
    DEG_path = os.path.join(os.path.abspath(outputs_dir), 'DEG_results')

    paired = True
    # mapping
    if samples_dic:
        full_path_sample_dic,paired = checksamples(inputs_dir, samples_dic)
        map_dic = mapping(full_path_sample_dic,ref_prefix,map_result_path)
    if script:
        for sample in map_dic:
            with open('s1_'+sample+'.sh','w') as fh:
                fh.write(map_dic[sample])
    else:
        for sample in map_dic:
            utils.run_shell_cmd(map_dic[sample])

    # assemble and quantitation by stringtie
    sample_names = map_dic.keys()
    bams = [os.path.join(map_result_path, sample_name + '.bam') for sample_name in sample_names]
    assemble_quantity_dic = stringtie.assemblyandquantitation(sample_names,bams,gtf,assembly_path,dthreads,merge)
    if 'merge' in assemble_quantity_dic:
        merge_gtf_cmd = assemble_quantity_dic.pop('merge')
        with open('s2.1.5_merge_gtf.sh','w') as fh:
            fh.write(merge_gtf_cmd)
    if script:
        for sample in assemble_quantity_dic:
            with open('s2.1_'+sample+'.assemble.sh','w') as fh:
                fh.write(assemble_quantity_dic[sample][0])
            with open('s2.2_'+sample+'.count.sh','w') as fh:
                fh.write(assemble_quantity_dic[sample][1])
    else:
        for sample in assemble_quantity_dic:
            utils.run_shell_cmd(assemble_quantity_dic[sample][0])
        utils.run_shell_cmd(merge_gtf_cmd)
        for sample in assemble_quantity_dic:
            utils.run_shell_cmd(assemble_quantity_dic[sample][1])

    # get DEGs
    if merge:
        merged_gtf_file = os.path.join(assembly_path, 'merged.gtf')
        count_matrix_cmd = count.featurecount(bams, merged_gtf_file, DEG_path, paired)
    else:
        count_matrix_cmd = count.featurecount(bams, gtf, DEG_path, paired)
    with open('s3.1_all_count_matrix.sh', 'w') as fh:
        fh.write(count_matrix_cmd)

    compare_groups = compare.split(';')
    for group in compare_groups:
        A,B = group.split(':')[0],group.split(':')[1]
        two_group_samples = samples_dic[A] + samples_dic[B]
        two_group_names = [A] * len(samples_dic[A]) + [B] * len(samples_dic[B])
        DE_file = os.path.join(DEG_path, A + 'vs' + B, ".RNAseq_different_expression_genes_results.tsv")
        two_group_out_path = os.path.join(DEG_path,A + 'vs' + B)
        if len(samples_dic[A]) < 2 or len(samples_dic[B]) < 2: ## no bio repeats
            TMM_cmd, identify_cmd = nobiorepeat.makeDEscript(two_group_samples,two_group_names,DEG_path,A+'vs'+B)
            with open('s3.2.1_'+A+'vs'+B+'.TMM.sh','w') as fh:
                fh.write(TMM_cmd)
            with open('s3.2.2_'+A+'vs'+B+'.DEG.sh','w') as fh:
                fh.write(identify_cmd)

            volcano_cmd = volcano.volcano(DE_file,A,B,False)
            with open('s3_3_'++A+'vs'+B+'.volcano.sh','w') as fh:
                fh.write(volcano_cmd)

        elif len(samples_dic[A]) >= 2 and len(samples_dic[B]) >= 2: ## bio repeats
            cal_matrix_cmd, identify2_cmd = biorepeat.makedeseq2(assembly_path,two_group_samples,two_group_names,DEG_path,A+'vs'+B,readlength)
            with open('s3.2.1_' + A + 'vs' + B + '.cal.matrix.sh', 'w') as fh:
                fh.write(cal_matrix_cmd)
            with open('s3.2.2_'+A+'vs'+B+'.DEG.sh','w') as fh:
                fh.write(identify2_cmd)

            volcano_cmd = volcano.volcano(DE_file,A,B,True)
            with open('s3_3_'++A+'vs'+B+'.volcano.sh','w') as fh:
                fh.write(volcano_cmd)

    # GO and KEGG enrichment
        go_cmd = go.makeGO(DE_file,two_group_out_path,A+'vs'+B,go_orgdb)
        with open('s4_'+A+'vs'+B+'All.GO.sh','w') as fh:
            fh.write(go_cmd)
        kegg_cmd = kegg.makeKEGG(DE_file,two_group_out_path,A+'vs'+B,kegg_species)
        with open('s4_'+A+'vs'+B+'All.KEGG.sh','w') as fh:
            fh.write(kegg_cmd)