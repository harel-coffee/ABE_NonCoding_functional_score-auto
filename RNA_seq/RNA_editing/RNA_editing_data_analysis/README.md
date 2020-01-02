

# Get started

Our pipeline follows the methods from:

Grünewald, J., Zhou, R., Iyer, S. et al. CRISPR DNA base editors with reduced RNA off-target and self-editing activities. Nat Biotechnol 37, 1041–1048 (2019) doi:10.1038/s41587-019-0236-6

## Public Data

GSE129894

#### 1. Sampling RNA-seq data such that the number of uniquely mapped reads is 30M.

RNA-seq data is first mapped to the human genome using STAR. Then, based on the mapping rate, we down-sample RNA-seq fastq files using `seqtk` (see `down_sample.sh`):

`seqtk sample -s100 $1 sampling_ratio > $2; gzip $2`

sampling ratio is provided in `sampling_ratio.tsv`. The formula to calculate sampling ratio is: let `X` denote the desired number of uniquely mapped reads (i.e., X=30,000), `m_i`, `T_i` denote the uniquely mapped rate and total reads for sample `s_i`, then we need to find the sampling ratio `r_i` `r_j` for sample `s_i` `s_j`, such that `T_i * m_i * r_i` == `T_j * m_j * r_j`.

#### 2. RNA-seq variant identification using GATK.

GATK RNA-seq best practice pipeline is using the code from HemTools:https://hemtools.readthedocs.io/en/latest/content/NGS_pipelines/rna_seq_variant_call.html


#### 3. Run bam_count for control samples

https://github.com/genome/bam-readcount

The following control samples are used: 
- Hudep2 WT rep 1
- Hudep2 WT rep 2
- GSM3724237 HEK293T-HEKsite2-control

Code to run bam count is provided in `bam_count.lsf`. This is a shell script, syntax is provided in https://hemtools.readthedocs.io/en/latest/content/Gallery/run_lsf.html

This step uses bam count to get a coverage for each genomic position, which are then filtered (`parse_count.py`) based on: (1) read depth >= 10; (2) reference allele frequency >= 0.99.

#### 4. summarize RNA-editing vcf

Once we have a list of non-edited positions (from control samples), we can then find out that out of those positions, how many positions are edited in our treatment samples. The treatment samples are:
- Hudep2 ABEmax rep 1
- Hudep2 ABEmax rep 2
- GSM3724238	HEK293T-HEKsite2-ABEmax-rep2

The actual table was generated for:
- using Hudep2 WT rep 1 as control, how many RNA-edits are happened in WT rep 2?
- using Hudep2 WT rep 2 as control, how many RNA-edits are happened in WT rep 1?
- using Hudep2 WT rep 1 as control, how many RNA-edits are happened in Hudep2 ABEmax rep 1?
- using Hudep2 WT rep 2 as control, how many RNA-edits are happened in Hudep2 ABEmax rep 2?
- using GSM3724237 (HEK293T-HEKsite2-control) control, how many RNA-edits are happened in HEK293T GFP?
- using GSM3724237 (HEK293T-HEKsite2-control) control, how many RNA-edits are happened in HEK293T HEKsite2 ABEmax?

These results were the csv files in `*_AG.csv`, which were used to generate the jitter plot.

The code is `subset_vcf_given_positions.py`


