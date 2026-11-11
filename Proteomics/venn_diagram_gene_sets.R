# ------------------------------------------------------------------------------
# Venn diagram of several gene sets
#
# Reads one gene list per sheet from an Excel file, draws a Venn diagram
# (ggVennDiagram) and exports the genes in each set, the genes found in exactly
# k sets and the genes unique to each set.
#
# Input : Excel file with one sheet per gene list and a column 'Gene.names'
#         (as read by openxlsx; the header 'Gene names' becomes 'Gene.names')
# Output: PNG and PDF of the Venn diagram and an Excel file with the gene lists
#         in 'output/'
#
# Requires: openxlsx, ggVennDiagram, ggplot2
# Note: ggVennDiagram supports up to 7 sets.
# ------------------------------------------------------------------------------

library(openxlsx)
library(ggVennDiagram)
library(ggplot2)

# ---- Settings (edit here) ----------------------------------------------------
input_file  <- "path/to/gene_lists.xlsx"
sheet_names <- c("sheet1", "sheet2", "sheet3", "sheet4", "sheet5", "sheet6")
set_labels  <- c("Set1", "Set2", "Set3", "Set4", "Set5", "Set6")  # names in the plot
output_dir  <- "output"
output_name <- "venn_diagram"
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
n_sets <- length(sets)
print(lengths(sets))

# Venn diagram ----------------------------------------------------------------
venn_sets <- sets
names(venn_sets) <- paste0(set_labels, " (n=", lengths(sets), ")")  # sizes counted automatically
p <- ggVennDiagram(venn_sets) +
  scale_fill_gradient(low = "white", high = "blue")
print(p)

dir.create(output_dir, showWarnings = FALSE)
ggsave(file.path(output_dir, paste0(output_name, ".png")),
       p, width = 9, height = 7, dpi = 300)
ggsave(file.path(output_dir, paste0(output_name, ".pdf")),
       p, width = 9, height = 7)

# Export gene lists -----------------------------------------------------------
all_genes     <- unique(unlist(sets))
in_set        <- sapply(sets, function(g) all_genes %in% g)   # genes x sets (TRUE/FALSE)
sets_per_gene <- rowSums(in_set)

sheet_label <- function(x) substr(x, 1, 31)   # Excel sheet names: max. 31 characters

lists <- list()
for (s in set_labels) lists[[sheet_label(paste0("All_", s))]] <- sets[[s]]
for (k in 2:n_sets) {
  lists[[sheet_label(paste0("Exactly_", k, "_sets"))]] <- all_genes[sets_per_gene == k]
}
for (s in set_labels) {
  lists[[sheet_label(paste0("Unique_", s))]] <- all_genes[in_set[, s] & sets_per_gene == 1]
}

write.xlsx(lists, file.path(output_dir, paste0(output_name, "_lists.xlsx")),
           overwrite = TRUE)
