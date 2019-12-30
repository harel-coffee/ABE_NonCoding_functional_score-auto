# Get started

## Usage

`python feature_extraction.py`

This will output `ML_data.csv`, which contains motif-footprint features, DeepSEA, CADD, and DNA sequences (for gkm-SVM).

`python main_classification.py`

## Feature table

`ML_data.csv` is the feature table.

## Machine learning setup

Motifs are downloaded from JASPAR, Homer, CIS-BP, and ENCODE-motif databases.

Key erythroid motifs are used, including GATA1,ZBTB7A,KLF1,BCL11A,NFIX,E2F,NFYA,NFE2. CTCF is also included.

Motifs were mapped to DNA sequences and the top 5 scores were used. Adjusted motif mapping score (i.e., motif-footprint score) was calculated by multiplying the mappint score with footprint score.

Footprint score was calculated as the difference between mean(+- 2bp flanking nucleotides) - mean(top motif match covered nucleotides) 

These motif footprint scores were used to train a random forest model.

For DeepSEA and CADD scores, no training process was used. auROC and auPRC were calculated for each cross-validation fold.

For gkm-SVM, the extracted DNA sequences were the features that were directly inputs to `gkm_svm_train`.





















### Motifs used in this study

```
HSAPIENS-ENCODE-MOTIFS-ZBTB7A_KNOWN2
JASPA_MA0469.2_E2F3
HSAPIENS-ENCODE-MOTIFS-ZBTB7A_KNOWN3
HSAPIENS-ENCODE-MOTIFS-ZBTB7A_KNOWN4
JASPA_MA0864.1_E2F2
JASPA_MA0024.3_E2F1
JASPA_MA0750.2_ZBTB7A
CIS-B_M4453_1.02_BCL11A
CIS-B_M4640_1.02_ZBTB7A
Homer_known_GSE104676_BCL11A
JASPA_MA0060.3_NFYA
JASPA_MA0470.1_E2F4
Homer-Erythroblasts-ZBTB7A-GSE74977
JASPA_MA0671.1_NFIX
JASPA_MA0139.1_CTCF
JASPA_MA0035.3_GATA1
JASPA_MA0493.1_KLF1
JASPA_MA0841.1_NFE2
JASPA_MA0150.2_NFE2L2
HSAPIENS-ENCODE-MOTIFS-ZBTB7A_KNOWN1
```



