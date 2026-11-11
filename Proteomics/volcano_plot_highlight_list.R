# ------------------------------------------------------------------------------
# Volcano plot with a highlighted protein list
#
# Plots log2 ratio against -log10(p-value) from a differential-expression table
# (e.g. exported from Perseus), colours the proteins found in a reference list
# (for example a list of proteins of interest from the literature) and labels
# the listed proteins that are not significantly enriched. Selected genes can
# be labelled in addition.
#
# Input : - Excel file with one sheet of results and the columns
#             Gene names  - gene symbols
#             pvalue      - -log10 of the t-test p-value
#             Difference  - log2 ratio (condition vs. control)
#         - Excel file with the reference list: column 'Gene names'
# Output: PNG and PDF of the plot and an Excel table of the labelled proteins
#         in 'output/'
#
# Requires: readxl, openxlsx
# ------------------------------------------------------------------------------

library(readxl)
library(openxlsx)

# ---- Settings (edit here) ----------------------------------------------------
input_file      <- "path/to/results.xlsx"
input_sheet     <- "sheetname"
reference_file  <- "path/to/reference_list.xlsx"   # column 'Gene names'
reference_sheet <- "Sheet1"
highlight_genes <- c("GENE1")      # extra genes to label (exact gene names)

lfc_threshold   <- 2               # log2 ratio cut-off (vertical line)
p_threshold     <- 1.3             # -log10(p-value) cut-off (horizontal line); 1.3 = p 0.05
x_limits        <- c(-8, 8)
y_limits        <- c(0, 7)
x_label         <- "log2 ratio (condition vs. control)"
y_label         <- "-log10 (P-value t-test)"
legend_label    <- "Reference list"

output_dir      <- "output"
output_name     <- "volcano_plot"
# ------------------------------------------------------------------------------

# Read and prepare data -------------------------------------------------------
res <- read_excel(input_file, sheet = input_sheet)
ref <- read_excel(reference_file, sheet = reference_sheet)

required_cols <- c("Gene names", "pvalue", "Difference")
missing_cols  <- setdiff(required_cols, names(res))
if (length(missing_cols) > 0) {
  stop("Missing column(s) in results sheet: ", paste(missing_cols, collapse = ", "))
}
if (!"Gene names" %in% names(ref)) stop("Reference list needs a 'Gene names' column.")

colnames(res)[colnames(res) == "pvalue"] <- "-log10(P-Value)"
dat <- res[, c("Gene names", "-log10(P-Value)", "Difference")]

in_reference <- dat$`Gene names` %in% ref$`Gene names`

# Reference proteins to label: not significant (below p cut-off), or significant
# but below the log2 ratio cut-off
subset_below <- dat[in_reference & dat$`-log10(P-Value)` < p_threshold, ]
subset_above <- dat[in_reference & dat$Difference < lfc_threshold &
                      dat$`-log10(P-Value)` > p_threshold, ]
labelled_table <- merge(subset_above, subset_below, by = "Gene names", all = TRUE)

# Plot function --------------------------------------------------------------
draw_volcano <- function() {
  plot(dat$Difference, dat$`-log10(P-Value)`,
       xlab = x_label, ylab = y_label,
       pch = 16, col = "grey", cex = 1.5,
       xlim = x_limits, ylim = y_limits)

  # Reference proteins in green
  points(dat$Difference[in_reference], dat$`-log10(P-Value)`[in_reference],
         pch = 16, cex = 1.5, col = "green")

  abline(v = lfc_threshold, col = "black", lty = 2, lwd = 2)
  abline(h = p_threshold,   col = "black", lty = 2, lwd = 2)
  legend("topleft", legend = legend_label, col = "green", pch = 16, cex = 0.8)

  # Labels for reference proteins below the cut-offs
  text(subset_below$Difference, subset_below$`-log10(P-Value)`,
       labels = subset_below$`Gene names`, pos = 4, cex = 0.8, col = "black")
  text(subset_above$Difference, subset_above$`-log10(P-Value)`,
       labels = subset_above$`Gene names`, pos = 4, cex = 0.8, col = "black")

  # Additional genes of interest (black, labelled above the point)
  sel <- dat$`Gene names` %in% highlight_genes
  if (any(sel)) {
    points(dat$Difference[sel], dat$`-log10(P-Value)`[sel],
           pch = 16, cex = 1.5, col = "black")
    text(dat$Difference[sel], dat$`-log10(P-Value)`[sel],
         labels = dat$`Gene names`[sel], pos = 3, cex = 1)
  }
}

# Show and save ---------------------------------------------------------------
draw_volcano()

dir.create(output_dir, showWarnings = FALSE)
png(file.path(output_dir, paste0(output_name, ".png")),
    width = 2000, height = 2000, res = 300)
draw_volcano()
invisible(dev.off())

pdf(file.path(output_dir, paste0(output_name, ".pdf")), width = 7, height = 7)
draw_volcano()
invisible(dev.off())

write.xlsx(labelled_table,
           file.path(output_dir, paste0(output_name, "_labelled_proteins.xlsx")),
           rowNames = FALSE)

# Optional: label points by clicking on them in the plot window (interactive).
# After running draw_volcano(), uncomment and run:
# identify(dat$Difference, dat$`-log10(P-Value)`, labels = dat$`Gene names`)
