# ------------------------------------------------------------------------------
# Heatmap of log2 ratios from several comparisons, with hierarchical clustering
#
# Reads one sheet per comparison from an Excel file (e.g. exported from Perseus),
# merges them by gene name, draws a clustered heatmap (pheatmap) and splits the
# proteins into clusters.
#
# Input : Excel file with one sheet per comparison. Each sheet needs
#           Gene names  - gene symbols
#           <value_col> - the value to plot (e.g. log2 ratio; name set below)
#         Alternatively, an already merged Excel sheet (first column = gene
#         names, one numeric column per comparison), see 'merged_file'.
# Output: PNG and PDF of the heatmap, the merged table (Excel) and the cluster
#         assignment (CSV) in 'output/'
#
# Requires: readxl, openxlsx, pheatmap
# ------------------------------------------------------------------------------

library(readxl)
library(openxlsx)
library(pheatmap)

# ---- Settings (edit here) ----------------------------------------------------
input_file    <- "path/to/results.xlsx"
sheet_names   <- c("sheet1", "sheet2", "sheet3", "sheet4", "sheet5")
column_labels <- sheet_names          # labels shown below the heatmap
value_col     <- "Difference"         # column to plot in every sheet

merged_file   <- NULL                 # optional: already merged Excel file;
merged_sheet  <- "Sheet1"             # if set, the sheets above are not read

n_clusters      <- 10                 # number of row clusters
center_at_zero  <- FALSE              # TRUE: white = 0, symmetric colour scale
cell_width      <- 100
cell_height     <- 0.15
output_dir      <- "output"
output_name     <- "heatmap"
# ------------------------------------------------------------------------------

# Read and merge --------------------------------------------------------------
if (!is.null(merged_file)) {
  merged_data <- as.data.frame(read_excel(merged_file, sheet = merged_sheet))
} else {
  read_one <- function(sheet, label) {
    d <- read_excel(input_file, sheet = sheet)
    missing_cols <- setdiff(c("Gene names", value_col), names(d))
    if (length(missing_cols) > 0) {
      stop("Missing column(s) in sheet '", sheet, "': ",
           paste(missing_cols, collapse = ", "))
    }
    d <- as.data.frame(d[, c("Gene names", value_col)])
    names(d)[2] <- label
    d
  }
  tables <- mapply(read_one, sheet_names, column_labels, SIMPLIFY = FALSE)
  # Inner join: only proteins present in all sheets are kept
  merged_data <- Reduce(function(a, b) merge(a, b, by = "Gene names"), tables)
}
cat("Proteins in the heatmap:", nrow(merged_data), "\n")

# Convert to matrix -----------------------------------------------------------
data_matrix <- as.matrix(merged_data[, -1])
rownames(data_matrix) <- make.unique(as.character(merged_data[[1]]))

# Clustering cannot handle missing values
complete <- complete.cases(data_matrix)
if (any(!complete)) {
  cat("Removing", sum(!complete), "protein(s) with missing values\n")
  data_matrix <- data_matrix[complete, , drop = FALSE]
}

# Colour scale: blue = down, red = up ----------------------------------------
colors <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- NA
if (center_at_zero) {
  lim    <- max(abs(data_matrix))
  breaks <- seq(-lim, lim, length.out = 101)
}

# Heatmap ---------------------------------------------------------------------
make_heatmap <- function(filename = NA) {
  pheatmap(data_matrix,
           color = colors,
           breaks = breaks,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "ward.D2",
           scale = "none",            # keep the original values
           fontsize_row = 8,
           fontsize_col = 10,
           main = "Log2 ratio heatmap",
           cellwidth = cell_width,
           cellheight = cell_height,
           legend = TRUE,
           filename = filename)
}

heatmap_result <- make_heatmap()      # show on screen

dir.create(output_dir, showWarnings = FALSE)
invisible(make_heatmap(file.path(output_dir, paste0(output_name, ".png"))))
invisible(make_heatmap(file.path(output_dir, paste0(output_name, ".pdf"))))

# Cluster assignment ----------------------------------------------------------
row_clusters <- cutree(heatmap_result$tree_row, k = n_clusters)
protein_clusters <- data.frame(Protein = rownames(data_matrix),
                               Cluster = row_clusters)
write.csv(protein_clusters,
          file.path(output_dir, paste0(output_name, "_clusters.csv")),
          row.names = FALSE)

write.xlsx(merged_data,
           file.path(output_dir, paste0(output_name, "_merged_data.xlsx")),
           rowNames = FALSE)
