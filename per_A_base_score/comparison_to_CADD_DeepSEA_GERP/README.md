
# Get started

This folder provides analysis comparing our score vs. DeepSEA, CADD, GERP.

## Usage

`python main.py` to combine every score by A's position and make the violin plots

`python peak_vs_neighbors.py` to get the signal plots (similar to DeepTools plotHeatmap)

## Algorithm

Violin plots were made by comparing HbF high (HbF score >=50) vs. HbF low (HbF score == 0).

For signal plots, top 500 scores and their +- 100 neighbors 
(scores are ranked, a neighbor is defined as its closest 
ranked score) were plotted for HbFBase, DeepSEA, and CADD.

## Directories

Raw scores were stored in `DeepSEA_CADD_GERP`.

Some intermediate codes were in `jupyter_notebooks`
