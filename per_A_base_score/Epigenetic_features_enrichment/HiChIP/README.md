# How to run get_number_interaction_pos_neg.py

Make sure you login to a compute node.

### Input bed

Prepare your input bed as 4-col tsv, chr, start, end, name (should be unique for each line)

### Usage

```
module load conda3/202011

source activate /home/yli11/.conda/envs/captureC

python get_number_interaction_pos_neg.py input.bed

```

### Output

output file is named as {input_name}_degree.csv

the second last column is the number of loops overlapped with your query.
