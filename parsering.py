
import os

def parse_short_read_dir(inputs, outs, seq_type='PE'):
    input_path = os.path.abspath(inputs)
    out_path = os.path.abspath(outs)
    lst = os.listdir(input_path)
    dic = {}
    seq_suffix = [seq + zz for seq in ['.fq', '.fastq'] for zz in ['', '.gz', '.bz2', '.zip']]
    lst = [file for file in lst for s in seq_suffix if file.endswith(s)]
    try:
        if lst:
            samples = [i.split('_1')[0] for i in lst if '_1' in i]
            suffix = [i.split('_1')[1] for i in lst if '_1' in i]
    except IOError:
        print('No such file or directory or wrong file suffix(e.q. sample_1.fq.gz,sample.fq.gz)')

    if samples:
        if seq_type == 'PE':
            for index, sample in enumerate(samples):
                dic[sample] = [os.path.join(input_path, sample + '_1' + suffix[index]),
                               os.path.join(input_path, sample + '_2' + suffix[index]),
                               os.path.join(out_path, sample + '_1' + suffix[index]),
                               os.path.join(out_path, sample + '_2' + suffix[index]),
                               os.path.join(out_path,sample)]
        elif seq_type == 'SE':
            for index, sample in enumerate(samples):
                dic[sample] = [os.path.join(input_path,sample+'_1'+suffix[index]),
                               os.path.join(out_path,sample+'_1'+suffix[index]),
                               os.path.join(out_path,sample)]

    return dic

def parse_long_read_dir(inputs,outs):
    input_path = os.path.abspath(inputs)
    out_path = os.path.abspath(outs)
    lst = os.listdir(input_path)
    dic = {}
    tgs_seq_suffix =  [seq + zz for seq in ['.fa', '.fasta', '.fastq', '.fq'] for zz in ['','.gz', '.bz2', '.zip']]
    tgs_lst = [file for file in lst for suffix in tgs_seq_suffix if file.endswith(suffix)]
    try:
        if tgs_lst:
            samples = [os.path.splitext(os.path.splitext(i)[0])[0] for i in tgs_lst]
            suffix = [os.path.splitext(os.path.splitext(i)[0])[-1]+os.path.splitext(i)[-1] for i in tgs_lst]
    except IOError:
        print('No such file or directory or wrong file suffix(e.q. sample.fa)')

    for index, sample in enumerate(samples):
        dic[sample] = [os.path.join(input_path, sample+suffix[index]),
                       os.path.join(out_path, sample+suffix[index]),
                       os.path.join(out_path,sample)]
    return dic

