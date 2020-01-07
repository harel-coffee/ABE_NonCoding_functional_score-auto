

# Get started

This folder provides inputs and outputs for gRNA deep sequencing, which is the first step for all other analyses.

`Undetermined.fastq.gz` was demultiplexed according to: https://hemtools.readthedocs.io/en/latest/content/Bioinformatics_tools/demultiplexing2.html

gRNA counts and significants qualitification were done using MaGeCK, following : https://hemtools.readthedocs.io/en/latest/content/NGS_pipelines/crispr_seq.html


## MaGeck RRA Commands


`mageck count --output-prefix ABE -l gRNA_lib_mageck_format.csv --fastq *`

`mageck test --normcounts-to-file --count-table ABE.count.txt --norm-method control --output-prefix ABE_RRA_results --control-sgrna control.list -t ABE_HBF_HIGH_R1.fastq.gz,ABE_HBF_HIGH_R2.fastq.gz -c ABE_HBF_LOW_R1.fastq.gz,ABE_HBF_LOW_R2.fastq.gz`

