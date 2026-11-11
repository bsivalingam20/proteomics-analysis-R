# ------------------------------------------------------------------------------
# UpSet plot of several gene sets
#
# Reads one gene list per sheet from an Excel file, builds a membership table
# (gene x set, 0/1) and draws an UpSet plot (UpSetR) of the overlaps. An UpSet
# plot stays readable with more sets than a Venn diagram.
#
# Input : Excel file with one sheet per gene list and a column 'Gene.names'
#         (as read by openxlsx; the header 'Gene names' becomes 'Gene.names')
# Output: PNG and PDF of the UpSet plot and the membership table (Excel)
#         in 'output/'
#
# Requires: openxlsx, UpSetR
# ------------------------------------------------------------------------------

library(openxlsx)
library(UpSetR)

# ---- Settings (edit here) ----------------------------------------------------
input_file  <- "path/to/gene_lists.xlsx"
sheet_names <- c("sheet1", "sheet2", "sheet3", "sheet4", "sheet5", "sheet6")
set_labels  <- c("Set1", "Set2", "Set3", "Set4", "Set5", "Set6")  # names in the plot
output_dir  <- "output"
output_name <- "upset_plot"
# ------------------------------------------------------------------------------

# Read gene lists (duplicates, NA and empty entries removed) -------------------
read_genes <- function(sheet) {
  x <- read.xlsx(input_file, sheet = sheet)
  if (!"Gene.names" %in% names(x)) {
    stop("Column 'Gene.names' not found in sheet '", sheet, "'")
  }
  genes <- as.character(x$Gene.names)
  unique(genes[!is.na(genes) & genes != ""])
}

sets <- lapply(sheet_names, read_genes)
names(sets) <- set_labels
print(lengths(sets))

# Membership table: one row per gene, one 0/1 column per set -------------------
all_genes  <- unique(unlist(sets))
membership <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)
for (s in set_labels) membership[[s]] <- as.numeric(all_genes %in% sets[[s]])

# UpSet plot ------------------------------------------------------------------
draw_upset <- function() {
  print(upset(membership, sets = set_labels, nsets = length(set_labels),
              order.by = "freq"))
}

draw_upset()

dir.create(output_dir, showWarnings = FALSE)
png(file.path(output_dir, paste0(output_name, ".png")),
    width = 2400, height = 1600, res = 300)
draw_upset()
invisible(dev.off())

pdf(file.path(output_dir, paste0(output_name, ".pdf")), width = 10, height = 7)
draw_upset()
invisible(dev.off())

# Export membership table -----------------------------------------------------
write.xlsx(membership,
           file.path(output_dir, paste0(output_name, "_membership.xlsx")),
           overwrite = TRUE)
