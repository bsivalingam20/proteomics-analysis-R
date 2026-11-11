# ------------------------------------------------------------------------------
# Colocalization quantification: violin + box plot with pairwise t-tests
#
# Plots per-cell Pearson's correlation coefficients (r) from an imaging
# colocalization analysis for several conditions and compares each condition
# against a reference condition (the first listed condition) with a t-test.
#
# Input : Excel file (.xlsx) with one sheet per experiment and two columns:
#           Condition  - condition name (e.g. A, B)
#           Value      - Pearson's r for one cell / image
# Output: PNG and PDF of the plot in the folder 'output/'
#
# Requires: ggplot2, ggpubr, readxl
# ------------------------------------------------------------------------------

library(ggplot2)
library(ggpubr)
library(readxl)

# ---- Settings (edit here) ----------------------------------------------------
input_file  <- NULL            # path to the .xlsx file; NULL opens a file dialog
sheet_name  <- "sheetname"     # sheet that holds the data
cond_levels <- c("A", "B", "C", "D")  # plot order; first = reference
plot_title  <- "X-Y colocalization"
y_label     <- "Pearson's correlation coefficient (r)"
output_dir  <- "output"
output_name <- "coloc_violin_boxplot"
# ------------------------------------------------------------------------------

# Read data -------------------------------------------------------------------
if (is.null(input_file)) input_file <- file.choose()

dat <- read_excel(input_file, sheet = sheet_name)

required_cols <- c("Condition", "Value")
missing_cols  <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing column(s) in sheet '", sheet_name, "': ",
       paste(missing_cols, collapse = ", "))
}

dat$Condition <- factor(dat$Condition, levels = cond_levels)

# Pairwise comparisons of every condition vs. the reference (first level)
reference   <- cond_levels[1]
comparisons <- lapply(cond_levels[-1], function(x) c(reference, x))

# Plot ------------------------------------------------------------------------
# Note: p-values are unadjusted t-tests. Add p.adjust.method (e.g. "BH") to
# stat_compare_means() if you want multiple-testing correction.
p <- ggplot(dat, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = "gray30") +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5,
               fill = "white") +
  stat_compare_means(comparisons = comparisons,
                     method = "t.test",
                     label = "p.signif",   # stars: * ** *** ****, ns
                     hide.ns = FALSE,
                     bracket.size = 0.6,
                     tip.length = 0.02) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(-0.25, 1.25)) +
  theme_minimal(base_size = 14) +
  labs(x = NULL, y = y_label, title = plot_title) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"))

print(p)

# Save ------------------------------------------------------------------------
dir.create(output_dir, showWarnings = FALSE)
ggsave(file.path(output_dir, paste0(output_name, ".png")), p,
       width = 5, height = 5, dpi = 300)
ggsave(file.path(output_dir, paste0(output_name, ".pdf")), p,
       width = 5, height = 5)
