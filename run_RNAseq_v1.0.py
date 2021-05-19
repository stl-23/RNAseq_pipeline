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
                sample1 = glob.glob(os.path.join(os.path.abspath(inputs_dir).rstrip('/'),sample+'_1*'))
                sample2 = glob.glob(os.path.join(os.path.abspath(inputs_dir).rstrip('/'),sample+'_2*'))
                if sample1 and sample2:
                    new_dic[group] = [sample1,sample2]
                elif sample1 and not sample2:
                    new_dic[group] = [sample1]
                else:
                    print('Can not find files,please check and use names as sample_1/2.fq/.fastq(.gz)')

    return new_dic

def mapping(dic,ref_prefix):
    out_path = os.path.join(os.path.abspath(outputs_dir), 'mapping_results')
    new_dic = checksamples(dic)
    for group, samples in new_dic.items():
        if len(samples) == 2:
            cmd = hisat2.align(samples[0],samples[1],ref_prefix,out_path,threads)
        elif len(samples) == 1:
            cmd = hisat2.align(samples[0],ref_prefix,out_path,threads)

    return cmd

def DEG():
    





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
    general.add_argument('-g', '--gff3', type=str,default='',
                         help="The gff3 file of reference for variation sites annotation")
    general.add_argument('-v','--build_version',type=str,default='',choices=['hg19','hg38'],
                          help='Human genome build version,if used,do not set -r and -g')
    general.add_argument('-s','--samples',type=str,
                         help="Samples and groups names: group1:sample1,sample2,sample3;group2:sample4,sample5,sample6;..."
                              "a group of samples should be seperated by comma,and groups should be seperated "
                              "by semicolon")
    general.add_argument('-c','--compare',type=str,
                         help="DEG(Differentially Expressed Genes) group pairs. group1:group2;group1:group3...")
    mapping = parser.add_argument_group(title='Mapping options')

    args = parser.parse_args()
    inputs_dir = args.input
    outputs_dir = args.output
    ref = args.reference
    gff3 = args.gff3
    buildver = args.build_version
    samples = args.samples
    compare = args.compare

### parse parameters ###
    index_shell = os.path.abspath(os.path.dirname(__file__))+'/align/index.sh'
    statistics_shell = os.path.abspath(os.path.dirname(__file__))+'/align/statistics.sh'
    if buildver and not ref and not gff3:
        if buildver == 'hg19':
            genomicsdb = os.path.join(os.path.dirname(os.path.abspath(inputs_dir).rstrip('/')), 'genomicsdb/hg19/')
            ref = os.path.join(genomicsdb,'hg19.fa')
            gff3 = os.path.join(genomicsdb,'hg19.gff3')
            ref_prefix = os.path.splitext(ref)[0]
        elif buildver == 'hg38':
            genomicsdb = os.path.join(os.path.abspath(inputs_dir), 'genomicsdb/hg38/')
            ref = os.path.join(genomicsdb,'hg38.fa')
            ref_prefix = os.path.splitext(ref)[0]
            gff3 = os.path.join(genomicsdb,'hg38.gff3')
    elif not buildver and (ref and gff3):
        ref_prefix = os.path.splitext(os.path.basename(os.path.abspath(ref)))[0]
        #if maptool == 'hisat2':
        subprocess.check_call(['sh', index_shell, ref, 'hisat2', ref_prefix])
        subprocess.check_call(['sh', statistics_shell, ref, ref_prefix])
        ref = ref_prefix+'.fa'
        gff3 = ref_prefix+'.gff3'

  ## {group1:[sample1,sample2,sample3],group2:[sample4,sample5,sample6]}
    try:
        samples_dic = {i.split(':')[0]:i.split(':')[1].split(',') for i in samples.strip().split(";")}
    except IOError:
        print('Wrong group:sample format')
        exit()
    finally:
        samples_dic = {}

    ##
