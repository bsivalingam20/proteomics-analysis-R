# ------------------------------------------------------------------------------
# Correlation plot of Z-scores from two analyses
#
# Merges two result tables by gene name, plots the Z-scores of analysis 1
# against analysis 2, adds the Pearson correlation coefficient and a linear
# fit, and highlights proteins from a reference list.
#
# Input : - Two Excel files, each with a column 'Gene names' and one Z-score
#           column (names set below)
#         - Optional Excel file with a list of proteins to highlight
#           (column 'Gene names')
# Output: PNG and PDF of the plot and the merged table (Excel) in 'output/'
#
# Requires: readxl, openxlsx
# ------------------------------------------------------------------------------

library(readxl)
library(openxlsx)

# ---- Settings (edit here) ----------------------------------------------------
file_1      <- "path/to/analysis_1.xlsx"
sheet_1     <- "Sheet1"
zscore_col1 <- "Zscore_analysis1"      # Z-score column in file 1 (x-axis)

file_2      <- "path/to/analysis_2.xlsx"
sheet_2     <- "Sheet1"
zscore_col2 <- "Zscore_analysis2"      # Z-score column in file 2 (y-axis)

highlight_file  <- "path/to/highlight_list.xlsx"   # column 'Gene names'; NULL = none
highlight_sheet <- "Sheet1"
extra_genes     <- c("GENE1")          # additional genes to label (exact names)

x_label     <- "Z-score, analysis 1"
y_label     <- "Z-score, analysis 2"
output_dir  <- "output"
output_name <- "zscore_correlation_plot"
# ------------------------------------------------------------------------------

# Read and merge --------------------------------------------------------------
d1 <- read_excel(file_1, sheet = sheet_1)
d2 <- read_excel(file_2, sheet = sheet_2)

for (check in list(list(d1, zscore_col1, "file 1"), list(d2, zscore_col2, "file 2"))) {
  missing_cols <- setdiff(c("Gene names", check[[2]]), names(check[[1]]))
  if (length(missing_cols) > 0) {
    stop("Missing column(s) in ", check[[3]], ": ", paste(missing_cols, collapse = ", "))
  }
}

merged <- merge(d1[, c("Gene names", zscore_col1)],
                d2[, c("Gene names", zscore_col2)],
                by = "Gene names")
names(merged) <- c("Gene names", "x", "y")

# Keep only proteins with finite Z-scores in both analyses
valid <- merged[is.finite(merged$x) & is.finite(merged$y), ]
cat("Proteins in both analyses with valid Z-scores:", nrow(valid), "\n")

correlation <- cor(valid$x, valid$y)   # Pearson
fit         <- lm(y ~ x, data = valid)

# Proteins to highlight
highlight_genes <- extra_genes
if (!is.null(highlight_file)) {
  ref <- read_excel(highlight_file, sheet = highlight_sheet)
  highlight_genes <- unique(c(ref$`Gene names`, extra_genes))
}
marked <- valid[valid$`Gene names` %in% highlight_genes, ]

# Plot function --------------------------------------------------------------
draw_plot <- function() {
  lim_x <- max(abs(valid$x)) * 1.1
  lim_y <- max(abs(valid$y)) * 1.1
  plot(valid$x, valid$y,
       xlab = x_label, ylab = y_label,
       main = paste0("Z-score correlation (r = ", round(correlation, 2), ")"),
       pch = 16, col = "black",
       xlim = c(-lim_x, lim_x), ylim = c(-lim_y, lim_y))
  abline(fit, col = "blue")

  # Highlighted proteins in magenta
  points(marked$x, marked$y, col = "magenta", pch = 19)

  # Label the additional genes of interest
  labelled <- valid[valid$`Gene names` %in% extra_genes, ]
  text(labelled$x, labelled$y, labels = labelled$`Gene names`,
       pos = 3, cex = 1)
}

# Show and save ---------------------------------------------------------------
draw_plot()

dir.create(output_dir, showWarnings = FALSE)
png(file.path(output_dir, paste0(output_name, ".png")),
    width = 2000, height = 2000, res = 300)
draw_plot()
invisible(dev.off())

pdf(file.path(output_dir, paste0(output_name, ".pdf")), width = 7, height = 7)
draw_plot()
invisible(dev.off())

write.xlsx(merged, file.path(output_dir, paste0(output_name, "_merged_data.xlsx")),
           rowNames = FALSE)

# Optional: label points by clicking on them in the plot window (interactive).
# After running draw_plot(), uncomment and run:
# identify(valid$x, valid$y, labels = valid$`Gene names`)
