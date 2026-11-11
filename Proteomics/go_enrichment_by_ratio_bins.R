# ------------------------------------------------------------------------------
# GO enrichment of proteins binned by log2 ratio
#
# Purpose: find out from which log2 ratio onwards the most significant GO terms
# are of interest (e.g. lysosome-related), to justify a stringent ratio cut-off.
#
# Significantly regulated proteins (up or down) are grouped into bins of the
# absolute log2 ratio. A GO over-representation analysis (clusterProfiler) is
# run for every bin and, optionally, for all regulated proteins together as a
# reference. Saves one bar plot per bin, one Excel file with a results sheet per
# bin (including gene symbols) and a summary sheet comparing the bins.
#
# Input : tab-delimited text file, e.g. exported from Perseus, with columns
#           Gene.names  - gene symbols (one symbol per row)
#           Difference  - log2 ratio
#           pvalue      - -log10 of the t-test p-value
#         Decimal commas are converted to decimal points automatically.
# Output: PNG plots and an Excel file with the GO results in 'output/'
#
# Requires: clusterProfiler, org.Hs.eg.db, ggplot2, openxlsx
# ------------------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)   # human annotation; use another OrgDb for other species
library(ggplot2)
library(openxlsx)

# ---- Settings (edit here) ----------------------------------------------------
input_file    <- "path/to/differential_expression.txt"
direction     <- "down"                # "down" (Difference < 0) or "up" (> 0)
p_threshold   <- 1.3                   # -log10(p-value) cut-off; 1.3 = p 0.05
bin_edges     <- c(0, 2, 3, 4, 5, 6)   # bins of |log2 ratio|: [0,2), [2,3), ...
include_all   <- TRUE                  # also run all regulated proteins together
ontology      <- "MF"                  # "BP", "MF" or "CC"
show_terms    <- 20                    # number of terms in each bar plot
term_pattern  <- "lysosom"             # GO terms matching this text are flagged
output_dir    <- "output"
output_name   <- "GO_terms_binned"
# ------------------------------------------------------------------------------

# Read and clean data ---------------------------------------------------------
data <- read.table(input_file, sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE)

required_cols <- c("Gene.names", "Difference", "pvalue")
missing_cols  <- setdiff(required_cols, names(data))
if (length(missing_cols) > 0) {
  stop("Missing column(s) in input file: ", paste(missing_cols, collapse = ", "))
}

# Convert to numeric first (handles decimal commas), then filter
to_numeric <- function(x) as.numeric(gsub(",", ".", as.character(x)))
data$Difference <- to_numeric(data$Difference)
data$pvalue     <- to_numeric(data$pvalue)

data_filter <- data[data$pvalue >= p_threshold, ]
data_filter <- if (direction == "down") {
  data_filter[data_filter$Difference < 0, ]
} else {
  data_filter[data_filter$Difference > 0, ]
}
data_filter$abs_ratio <- abs(data_filter$Difference)
cat("Proteins passing the filters (", direction, "):", nrow(data_filter), "\n")

# Define the bins -------------------------------------------------------------
bins <- list()
if (include_all) bins[["all"]] <- c(0, Inf)
for (i in seq_len(length(bin_edges) - 1)) {
  bins[[paste0(bin_edges[i], "to", bin_edges[i + 1])]] <-
    c(bin_edges[i], bin_edges[i + 1])
}

# GO enrichment per bin -------------------------------------------------------
dir.create(output_dir, showWarnings = FALSE)
results <- list()
summary_rows <- list()

for (bin_label in names(bins)) {
  lo <- bins[[bin_label]][1]
  hi <- bins[[bin_label]][2]

  genes <- data_filter$Gene.names[data_filter$abs_ratio >= lo &
                                    data_filter$abs_ratio < hi]
  cat("Bin", bin_label, ":", length(genes), "genes\n")
  if (length(genes) == 0) next

  # Gene symbols -> Entrez IDs
  entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

  # Note: no custom background (universe) is set, so all annotated genes of the
  # OrgDb are used. The detected proteins can be passed to enrichGO(universe = ...)
  # as a more suitable background for proteomics data.
  ego <- enrichGO(gene          = entrez$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = ontology,
                  pAdjustMethod = "BH")

  df <- if (is.null(ego)) data.frame() else as.data.frame(ego)
  if (nrow(df) == 0) {
    cat("  no enriched GO terms in this bin\n")
    summary_rows[[bin_label]] <- data.frame(
      bin = bin_label, n_genes = length(genes), n_terms = 0,
      top_term = NA, top_p_adjust = NA,
      best_matching_term = NA, best_matching_rank = NA,
      best_matching_p_adjust = NA)
    next
  }

  # Results table with gene symbols instead of Entrez IDs
  df <- df[order(df$p.adjust), ]
  df$geneSymbol <- as.data.frame(setReadable(ego, org.Hs.eg.db, "ENTREZID"))[
    rownames(df), "geneID"]
  results[[bin_label]] <- df

  # Summary: top term and best term matching 'term_pattern'
  hit <- grep(term_pattern, df$Description, ignore.case = TRUE)
  summary_rows[[bin_label]] <- data.frame(
    bin = bin_label, n_genes = length(genes), n_terms = nrow(df),
    top_term = df$Description[1], top_p_adjust = df$p.adjust[1],
    best_matching_term = if (length(hit) > 0) df$Description[hit[1]] else NA,
    best_matching_rank = if (length(hit) > 0) hit[1] else NA,
    best_matching_p_adjust = if (length(hit) > 0) df$p.adjust[hit[1]] else NA)

  # Bar plot
  p <- barplot(ego, showCategory = show_terms,
               title = paste0("GO term enrichment (", ontology, "), |log2 ratio| ",
                              bin_label, ", ", direction)) +
    scale_fill_gradient(low = "steelblue", high = "skyblue")
  ggsave(file.path(output_dir, paste0(output_name, "_", bin_label, ".png")),
         p, width = 8, height = 6, dpi = 300)
}

# Save results ----------------------------------------------------------------
if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  print(summary_df)
  write.xlsx(c(list(summary = summary_df), results),
             file.path(output_dir, paste0(output_name, ".xlsx")))
} else {
  warning("No proteins found in any bin.")
}
