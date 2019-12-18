

# Get started

See code in `average_model/get_average_frequency.py`.

Code is written in python2. Can be used in Python3 if you change `string.maketrans` to `bytes.maketrans`.

# Algorithm

The input list is a list of files that contain gRNA information. I only uploaded the actual data that is used in our code in the `data` folder.

Code is hard-coded to get A-G conversion rate.

The actual data used is the nucleotide frequency table from CrisprEsso2.

gRNA informatin is needed because a revcomp gRNA could be in the nucleotide frequency table. In that case, T-C is used.

The nucleotide frequency table gives you the nucleotide conversion rate at each gRNA position. So our algorithm is very simple, just combine the table from all exp and take the average.

# Usage

`cd average_model`

`python get_average_frequency.py`
















