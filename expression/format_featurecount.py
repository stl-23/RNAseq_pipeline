import os
import sys
if len(sys.argv) != 3:
    print("python format_featurecount.py All.read.count.txt output.file")
    exit(0)
with open(sys.argv[1]) as fh,open(sys.argv[2],'w') as fw:
    data = []
    for lines in fh:
        if lines.startswith('#'):
            continue
        else:
            line = lines.strip().split('\t')
            if lines.startswith('Geneid'):
                samples = line[6:]
            else:
                data.append(line[0]+'\t'+'\t'.join(line[6:]))
    new_samples = [i.strip().split('/')[-1].rstrip('.bam') for i in samples]
    fw.write('Geneid'+'\t'+'\t'.join(new_samples)+'\n'+'\n'.join(data))
