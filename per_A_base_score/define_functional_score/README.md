
# Get started

Empirical Brown’s method (EBM) is a non-parametric method to combine dependent p-values. 

## Usage

`python get_functional_score.py`

## Input

`ABE_RRA_results.normalized.txt` is the gRNA count table, EBM uses to calculate covariance matrix.

`ABE_high_vs_low_mageck_RRA_results.sgrna_summary.txt` contains the FDRs from MaGeCK RRA test.

`gRNA_all_A.bed` contains every A's coverage on each gRNA. This file is used to get all gRNAs that is overlaped with a particular A.

## Algorithm

The steps are straightforward. For each unique A, get its overlaped gRNAs (i.e., positions, FDRs, counts), which are the input to EBM. EBM outputs combined p-value and we take -np.log10.
