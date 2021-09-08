# RNAseq_pipelinev1.0

# Dependencies download and usage
1. Pull our docker image from Docker Hub
```
docker pull stl23/rnaseq:v1.5
```
2 Run the pipeline
```
docker run --rm -v  /mnt/data/stl/rna_seq/rnaseq/:/scripts -v /mnt/data/stl/rna_seq/test/cleandata/:/inputs -v /mnt/data/stl/rna_seq/test/RNA_results/:/outputs -v /mnt/data/stl/single_cell/refdata-gex-GRCh38-2020-A/:/genomedb stl23/rnaseq:v1.5 bash -c 'cd /outputs/ && python3 /scripts/run_RNAseq_v1.0.py -i /inputs -o /outputs -r /genomedb/fasta/genome.fa -g /genomedb/genes/genes.gtf --samples SRR8075135,SRR8075136,SRR8075137:SRR8075138,SRR8075139,SRR8075140 --groups ctl:case --compare ctl:case --org hsa --readlength 50 --jobs 3 --mthreads 5'
```
# Parameters
```
usage: run_RNAseq_v1.0.py [-h] [-i INPUT] [-o OUTPUT] [-r REFERENCE] [-g GTF]
                          [-v {hg19,hg38}] [--samples SAMPLES]
                          [--groups GROUPS] [--compare COMPARE] [--merge]
                          [--readlength READLENGTH] [--script] [--jobs JOBS]
                          [--mthreads MTHREADS] [--org ORG]

RNAseq pipeline v1.0

General options:
  -h, --help            show the help and exit
  -i INPUT, --input INPUT
                        The input directory of clean reads
  -o OUTPUT, --output OUTPUT
                        The output directory
  -r REFERENCE, --reference REFERENCE
                        The fasta file of reference
  -g GTF, --gtf GTF     The gene annotation file of reference,gtf format
  -v {hg19,hg38}, --build_version {hg19,hg38}
                        Human genome build version,if used,do not set -r and
                        -g
  --samples SAMPLES     Sample names. The sample names should be accordant
                        with that in input directorye.g. sample1,sample2,sampl
                        e3:sample4,sample5,sample6:sample7...a group of
                        samples should be seperated by comma
  --groups GROUPS       Group names. Should be consistent with the order of
                        samples names.e.g. group1:group2:group3...
  --compare COMPARE     DEG(Differentially Expressed Genes) group pairs.
                        group1:group2,group1:group3...
  --merge               Merge novel transcripts with reference gtf file if set
                        True,otherwise use reference gtf only
  --readlength READLENGTH
                        average read length for alternative splicing by rMATs.
  --script              Only generate shell scripts
  --jobs JOBS           The maximum jobs when run in local at the same time
  --mthreads MTHREADS   Maximum Threads

Enrichment options:
  --org ORG             Species for GO/KEGG enrichment
```
