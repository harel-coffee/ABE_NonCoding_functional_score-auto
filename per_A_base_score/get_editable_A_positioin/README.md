
# Get started

gRNA.all.bed contains 6194 gRNA locations.

# Usage

1. Get number of editable A.

Based on our amplicon sequencing result, positions from 1 to 11 (0-index) can be edited with at least 1% editing frequency

```
python main.py 1 11
Editing Window starts at:1, length is: 11
The following output should be an empty dataframe
Empty DataFrame
Columns: [0, 1, flag]
Index: []
('editable_A:', 9112, 'editable_C:', 10103, 'Total:', 19215)
```

The parameters are start position and length. Namely, which region in the gRNA is used to define editable A.

2. Get per A bed file.

```
python get_editable_bp_per_gRNA.py
```


This code iterates each gRNA and extract its A's genomic position. Same position could occur multiple times if it covers by multiple gRNAs.

If you remove duplicated A, you will have 13147 unique positions for A.

However, the gRNA - A information is important for subsequent analysis, so we will not remove duplicated A here.

Use `df = df.drop_duplicates(6)` to remove duplicates. And this is how we get the vcf file in order to extract GERP, CADD, and DeepSEA scores.



