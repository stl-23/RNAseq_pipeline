#!/usr/bin/env python
import os
import sys
import glob
#import collections
if len(sys.argv) != 5:
    print("python fpkm_tpm.py sample1,sample2,sample3:sample4,sample5... group1:group2... directory_of_strintie_result output_path")
    exit(0)
samples = sys.argv[1]
groups = sys.argv[2]
count_path = sys.argv[3]
outpath = sys.argv[4]
fw1 = open(outpath+'/All_fpkm.list','w')
fw2 = open(outpath+'/All_tpm.list','w')
fw3 = open(outpath+'/All_fpkm_rowmeans.list','w')
fw4 = open(outpath+'/All_tpm_rowmeans.list','w')
dic = { i[0]:i[1].split(",") for i in zip(groups.strip().split(":"),samples.strip().split(":")) }
sample_lst = [ y for x in samples.strip().split(":") for y in x.split(",")]
#dic_sample_all = collections.OrderedDict()
dic_sample_all = {}
fw1.write('GeneID'+'\t'+'\t'.join(sample_lst)+'\n')
fw2.write('GeneID'+'\t'+'\t'.join(sample_lst)+'\n')
fw3.write('GeneID'+'\t'+'\t'.join(dic.keys())+'\n')
fw4.write('GeneID'+'\t'+'\t'.join(dic.keys())+'\n')
for i in sample_lst:
    dic_sample = {}
    lst = glob.glob('/'.join([count_path,i,i+".abundance.tab"]))
    if not lst:
        raise Exception('No such file:%s' %('/'.join([count_path,i,i+".abundance.tab"])))
    else:
        with open(lst[0]) as fh:
            for lines in fh:
                if lines.startswith("Gene"):continue
                line = lines.strip().split('\t')
                geneid,start,stop,fpkm,tpm = line[0],line[4],line[5],line[-2],line[-1]
                genelength = int(stop)-int(start)+1
                if geneid in dic_sample:  ## alternative transcripts, select the longest transcript's fpkm and tpm
                    pre_length = int(dic_sample[geneid].split('_')[1])
                    pre_fpkm_tpm = '_'.join(dic_sample[geneid].split('_')[2:])
                    max_length = max(genelength,pre_length)
                    if max_length == pre_length:
                        dic_sample[geneid] = '_'.join([i,str(pre_length),pre_fpkm_tpm])
                    elif max_length == genelength:
                        dic_sample[geneid] = '_'.join([i,str(genelength),fpkm,tpm])
                else:
                    dic_sample[geneid] = '_'.join([i,str(genelength),fpkm,tpm])
        for gene in dic_sample:
            if gene in dic_sample_all:
                dic_sample_all[gene].append(dic_sample[gene])
            else:
                dic_sample_all[gene] = [dic_sample[gene]]

## {gene1:[sample1_genelength_fpkm_tpm,sample2_genelength_fpkm_tpm...],gene2:[sample1_genelength_fpkm_tpm,sample2_genelength_fpkm_tpm...]}  --> dic_sample_all
for k,v in dic_sample_all.items():
    fpkm = []
    tpm = []
    fpkm_means = []
    tpm_means = []
    for x,s in enumerate(sample_lst):
        #if len(sample_lst) > len(v):
        #    if s != v[x].split('_',1)[0]:
        #        v[x] = '_'.join([s,'0','0','0'])
        fpkm.append(v[x].split('_')[-2])
        tpm.append(v[x].split('_')[-1])
    fw1.write(k+'\t'+'\t'.join(fpkm)+'\n')
    fw2.write(k+'\t'+'\t'.join(tpm)+'\n')
    for group,samples in dic.items():
        fpkm_group = []
        tpm_group = []
        for sample in samples:
            for info in v:
                lst = info.split('_')
                if sample == lst[0]:
                    fpkm_group.append(float(lst[-2]))
                    tpm_group.append(float(lst[-1]))
        ave_fpkm = float(sum(fpkm_group))/len(samples)
        ave_tpm = float(sum(tpm_group))/len(samples)
        fpkm_means.append(ave_fpkm)
        tpm_means.append(ave_tpm)
    fw3.write(k+'\t'+'\t'.join([str(i) for i in fpkm_means])+'\n')
    fw4.write(k+'\t'+'\t'.join([str(i) for i in tpm_means])+'\n')

fw1.close()
fw2.close()
fw3.close()
fw4.close()
