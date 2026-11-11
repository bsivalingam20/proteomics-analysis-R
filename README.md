# R scripts for proteomics and imaging data analysis

R scripts for quantitative proteomics and imaging data analysis, written during
my PhD on lysosomal proteomics (Universität Osnabrück). The repository contains
code only: no raw data, sample information or unpublished results.

## Contents

### `proteomics/`: analysis and visualization of mass spectrometry data

| Script | Purpose |
|---|---|
| `volcano_plot_highlight_list.R` | Volcano plot (log2 ratio vs. -log10 p-value) with a highlighted reference protein list and labelled proteins |
| `go_enrichment_diverging_barplot.R` | GO over-representation analysis (clusterProfiler) of up- and downregulated proteins; diverging bar plot of the top terms |
| `go_enrichment_by_ratio_bins.R` | GO enrichment of significantly regulated proteins grouped into log2 ratio bins; one bar plot and one results sheet per bin |
| `heatmap_clustering.R` | Clustered heatmap (pheatmap, Ward clustering) of log2 ratios from several comparisons; exports the cluster of each protein |
| `mfuzz_soft_clustering.R` | Soft (fuzzy c-means) clustering of protein profiles with Mfuzz; cluster plots and a table of membership scores |
| `zscore_correlation_plot.R` | Scatter plot of Z-scores from two analyses with Pearson correlation, linear fit and highlighted protein list |
| `log2ratio_correlation_plot.R` | Scatter plot of log2 ratios from two datasets with Pearson correlation and linear fit |
| `venn_diagram_gene_sets.R` | Venn diagram of several gene lists; exports genes per set, genes found in exactly k sets and genes unique to each set |
| `upset_plot_gene_sets.R` | UpSet plot of several gene lists, plus the membership table (gene x set) |

### `imaging/`: quantification of imaging and Western blot data

| Script | Purpose |
|---|---|
| `colocalization_violin_boxplot.R` | Violin + box plot of Pearson's r from colocalization imaging, with pairwise t-tests against a reference condition |
| `western_blot_quantification_boxplot.R` | Box plots of quantified Western blot intensities with replicates; control vs. treated per sample group with t-tests |

## Usage

1. Install R packages:
   ```r
   install.packages(c("ggplot2", "ggpubr", "readxl", "openxlsx", "dplyr",
                      "pheatmap", "ggVennDiagram", "UpSetR"))
   BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "Mfuzz"))
   ```
2. Prepare the input file described in the header of each script.
3. Edit the settings block at the top of the script (file paths, sheet and group names, cut-offs).
4. Run the script in R or RStudio. Plots and tables are saved to an `output/` folder.

## Notes

- Input tables are expected in the format exported by MaxQuant / Perseus (for example the column names `Gene names` and `Difference`).
- Statistical tests are t-tests without multiple-testing correction unless stated otherwise in the script.
- The scripts use generic placeholder names for files, sheets and conditions; replace them with your own.
