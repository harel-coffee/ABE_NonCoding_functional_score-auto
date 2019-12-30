
# Get started

Empirical Brown’s method (EBM) is a non-parametric method to combine dependent p-values. 

## Usage

See `locus-level score.ipynb`

## Input

`ABE_RRA_results.normalized.txt` is the gRNA count table, EBM uses to calculate covariance matrix.

`ABE_high_vs_low_mageck_RRA_results.sgrna_summary.txt` contains the FDRs from MaGeCK RRA test.

`gRNA_in_locus.bed` contains all 307 loci and their covered gRNAs.

## Algorithm

The steps are straightforward. For each locus, get its covered gRNAs (i.e., positions, FDRs, counts), which are the input to EBM. EBM outputs combined p-value and we take -np.log10.
