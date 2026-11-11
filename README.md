# R scripts for imaging and proteomics data analysis

Scripts from my PhD work on lysosomal proteomics (Universität Osnabrück).
Code only: no raw or unpublished data are included.

## Contents

| Script | Purpose |
|---|---|
| `colocalization_violin_boxplot.R` | Violin + box plot of Pearson's r from colocalization imaging, with pairwise t-tests against a reference condition |

## Usage

1. Install the packages: `install.packages(c("ggplot2", "ggpubr", "readxl"))`
2. Prepare an Excel sheet with the columns `Condition` and `Value`.
3. Edit the settings block at the top of the script (file path, sheet name, condition order).
4. Run the script in R or RStudio. The plot is saved to `output/` as PNG and PDF.
