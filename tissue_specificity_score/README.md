


# Get started


This folder provides details about tissue specificity z-score calculation.

## Input

All the 13147 Adenines were used to calculate z-score. Each A was extended +- 10bp.

`blood_ATAC_bw_mean.tsv` was calculated from https://hemtools.readthedocs.io/en/latest/content/Machine_learning/bw_over_bed.html


Reference ATAC-seq data is obtained from:
- https://www.ncbi.nlm.nih.gov/pubmed/27526324
- https://www.ncbi.nlm.nih.gov/pubmed/31189107

ATAC-seq data was processed using `HemTools atac-seq` : https://hemtools.readthedocs.io/en/latest/content/NGS_pipelines/atac_seq.html

## Algorithm

Average ATAC-seq signal was calculated for each A across 16 different cell types. See group definition in `calculate_z_score.ipynb`.

Z-score calculation is similar to the entropy scores presented in: https://www.nature.com/articles/s41597-019-0071-0 . Basically, average signal (i.e., average read counts) was multiplied by 1,000,000 to get the reads per million value (i.e., RPM). Then, relative accessibility in a tissue type i is calculated as the fraction of reads in total read counts across all cell types. Lastly, z-score was calculated across each A.

## Output

`Ery_z_score.csv`





