

# Get started

See code in `average_model/get_average_frequency.py`.

# Usage

`cd average_model`

`python get_average_frequency.py`

The code is hard-coded to get A-G conversion rate.

## Input

The input is a list of output dirs from CrisprEsso2. The actual file used is `Quantification_window_nucleotide_percentage_table`, although in the `input.list` it used `Alleles_frequency_table_around_sgRNA_NNNNNNNN`. gRNA informatin is needed because a revcomp gRNA could be in the nucleotide frequency table. In that case, T-C is used.

The data folder only contains the minimal working example, which we only provided the nucleotide frequency table from CrisprEsso2, no other data.

# Algorithm

The nucleotide frequency table is a `5 * len(gRNA)` matrix, which provides the observed frequency from amplicon sequencing. To get the A-G conversion rate, we just need to find every A in the gRNA and then get the value at ['A','G']. 

A cutoff of 1% is used. If the observed frequency is below 1%, then that frequency is considered to be a noisy value and thus not used in the following analysis. Illumina sequencing error is 0.1%.















