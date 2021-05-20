#!/usr/bin/env python
import os
import argparse
import sys
import logging
import glob
import utils
import getmyconfig
from align import hisat2
from assembly import stringtie
from expression import (count,biorepeat,nobiorepeat)

def checksamples(samples_dic):
    new_dic = {}
    if samples_dic:
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
                    if group in new_dic:
                        new_dic[group].append(sample + ':' + sample_f)
                    else:
                        new_dic[group] = [sample + ':' + sample_f]
                else:
                    print('Can not find files,please check and use names as sample_1/2.fq/.fastq(.gz)')

    return new_dic

def mapping(dic,ref_prefix):
    out_path = os.path.join(os.path.abspath(outputs_dir), 'mapping_results')
    map_dic = {}
    for group, samples in dic.items():
        for sample in samples:
            sample_name = sample.split(':')[0]
            sample_f,sample_r = sample.split(':')[1].split(';')
            out_name = os.path.join(out_path, sample_name)
            if sample_f and sample_r: # paired end
                map_dic[sample_name] = hisat2.align(sample_f,sample_r,ref_prefix,out_name,mthreads)
            elif not sample_r: # single end
                map_dic[sample_name] = hisat2.align(sample_f,ref_prefix,out_name,threads)

    return map_dic

def quantitation(sample_names,gtf,out_path,):
    input_path = os.path.join(os.path.abspath(outputs_dir), 'mapping_results')
    out_path = os.path.join(os.path.abspath(outputs_dir), 'DEG_results')
    bams = [ os.path.join(input_path, sample_name+'.bam') for sample_name in sample_names]
    assembly = []
    quantitation = []
    assemble_quantity_dic = stringtie.assemblyandquantitation(sample_names,bams,gtf,out_path,dthreads)
    for i in assemble_quantity_dic:
        if i == 'merge':
            merge_gtf = assemble_quantity_dic[i]
        else:
            assembly.append(assemble_quantity_dic[i][0])
            quantitation.append(assemble_quantity_dic[i][1])

    return assembly,merge_gtf,quantitation




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
                         help="The gtf file of reference for variation sites annotation")
    general.add_argument('-v','--build_version',type=str,default='',choices=['hg19','hg38'],
                          help='Human genome build version,if used,do not set -r and -g')
    general.add_argument('-s','--samples',type=str,
                         help="Samples and groups names: group1:sample1,sample2,sample3;group2:sample4,sample5,sample6;..."
                              "a group of samples should be seperated by comma,and groups should be seperated "
                              "by semicolon")
    general.add_argument('-c','--compare',type=str,
                         help="DEG(Differentially Expressed Genes) group pairs. group1:group2;group1:group3...")
    general.add_argument('--script',type=str,default=True,choices=[True,False],
                         help='Only generate shell scripts')
    mapping = parser.add_argument_group(title='Mapping options')
    mapping.add_argument('--mthreads',type=int,
                         help="Threads for mapping (hisat2) ")


    args = parser.parse_args()
    inputs_dir = args.input
    outputs_dir = args.output
    ref = args.reference
    gtf = args.gtf
    buildver = args.build_version
    samples = args.samples
    compare = args.compare
    script = args.script

### parse parameters ###
    index_shell = os.path.abspath(os.path.dirname(__file__))+'/align/index.sh'
    statistics_shell = os.path.abspath(os.path.dirname(__file__))+'/align/statistics.sh'
    if buildver and not ref and not gtf:
        if buildver == 'hg19':
            genomicsdb = os.path.join(os.path.dirname(os.path.abspath(inputs_dir).rstrip('/')), 'genomicsdb/hg19/')
            ref = os.path.join(genomicsdb,'hg19.fa')
            gtf = os.path.join(genomicsdb,'hg19.gtf')
            ref_prefix = os.path.splitext(ref)[0]
        elif buildver == 'hg38':
            genomicsdb = os.path.join(os.path.abspath(inputs_dir), 'genomicsdb/hg38/')
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

    ## make shell script or run
    # mapping
    map_dic = mapping(checksamples(samples_dic),ref_prefix)
    if script:
        for sample in map_dic:
            with open('s1_'+sample+'.sh','w') as fh:
                fh.write(map_dic[sample])
    else:
        for sample in map_dic:
            utils.run_shell_cmd(map_dic[sample])

    # get DEGs

