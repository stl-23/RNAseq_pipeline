# RNAseq_pipelinev1.0

# Dependencies download and usage
1. Pull our docker image from Docker Hub
···
docker pull stl23/rnaseq:v1.4
···
2 Run the pipeline
···
docker run --rm -v  /mnt/data/stl/rna_seq/rnaseq/:/scripts -v /mnt/data/stl/rna_seq/test/cleandata/:/inputs -v /mnt/data/stl/rna_seq/test/RNA_results/:/outputs -v /mnt/data/stl/single_cell/refdata-gex-GRCh38-2020-A/:/genomedb stl23/rnaseq:v1.4 bash -c 'cd /outputs/ && python3 /scripts/run_RNAseq_v1.0.py -i /inputs -o /outputs -r /genomedb/fasta/genome.fa -g /genomedb/genes/genes.gtf --samples SRR8075135,SRR8075136,SRR8075137:SRR8075138,SRR8075139,SRR8075140 --groups ctl:case --compare ctl:case --org hsa --jobs 3 --mthreads 5'
···
