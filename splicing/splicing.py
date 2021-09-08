import os
from utils import getConfig,makedir

python3 = getConfig('Transcript','python3')
rmats = getConfig('Transcript','rmats')

def alternative_splicing(path,gtf,paired_or_single="paired",readlength=50,theads=8):
    tmp = os.path.join(path,'temp')
    b1 = os.path.join(tmp,'b1.txt')
    b2 = os.path.join(tmp,'b2.txt')
    cmd = f"source ~/.bashrc && {python3} {rmats} --b1 {b1} --b2 {b2} -t {paired_or_single} --gtf {gtf} --readLength {readlength} --variable-read-length \
    --od {path} --tmp {tmp} --nthread {theads}".format(**locals())

    return cmd
