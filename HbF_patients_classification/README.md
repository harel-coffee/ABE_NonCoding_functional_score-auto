
# Get started

Input: `snp_data.tsv`
Output: `*.pdf`
Code: `main_classification.py`

## Usage

`python main_classification.py`

## Machine learning setup

SCD patients classification using genetic variants.

Feature set 1: common genetic variants from GWAS study.

Feature set 2: feature set 1 + common genetic variants around (+- 20bp) top HbF scores.

Positive class: HbF patients with HbF > mean + 1*std

Negative class: HbF patients with mean + 1*std >= HbF >= mean - 1*std





