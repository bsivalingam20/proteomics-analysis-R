# ------------------------------------------------------------------------------
# GO enrichment of up- and downregulated proteins: diverging bar plot
#
# Splits a differential-expression table into up- and downregulated proteins,
# runs GO over-representation analysis for each group (clusterProfiler) and
# plots the top terms of both groups in one diverging bar chart. Bar length is
# the gene count (negative = down, positive = up), bar colour is -log10(adjusted
# p-value).
#
# Input : tab-delimited text file, e.g. exported from Perseus, with columns
#           Gene.names  - gene symbols (one symbol per row)
#           Difference  - log2 fold change (> 0 up, < 0 down)
# Output: PNG and PDF of the plot and a CSV table of the plotted terms in 'output/'
#
# Requires: clusterProfiler, org.Hs.eg.db, ggplot2, dplyr
# ------------------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)   # human annotation; use another OrgDb for other species
library(ggplot2)
library(dplyr)

# ---- Settings (edit here) ----------------------------------------------------
input_file  <- "path/to/differential_expression.txt"
ontology    <- "CC"      # "BP", "MF" or "CC"
top_n       <- 6         # number of terms shown per direction
p_cutoff    <- 0.05      # adjusted p-value cutoff for enrichGO
plot_title  <- "GO term enrichment (top 6 up and down)"
output_dir  <- "output"
output_name <- "go_enrichment_diverging_barplot"
# ------------------------------------------------------------------------------

# Read data -------------------------------------------------------------------
gene_data <- read.table(input_file, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)

required_cols <- c("Gene.names", "Difference")
missing_cols  <- setdiff(required_cols, names(gene_data))
if (length(missing_cols) > 0) {
  stop("Missing column(s) in input file: ", paste(missing_cols, collapse = ", "))
}

up_genes   <- gene_data$Gene.names[gene_data$Difference > 0]
down_genes <- gene_data$Gene.names[gene_data$Difference < 0]

# GO enrichment ---------------------------------------------------------------
# Note: no custom background (universe) is set, so all annotated genes of the
# OrgDb are used. For proteomics data, the detected proteins can be passed to
# enrichGO(universe = ...) as a more suitable background.
run_go <- function(genes, direction) {
  ego <- enrichGO(gene          = genes,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "SYMBOL",
                  ont           = ontology,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = p_cutoff)
  df <- as.data.frame(ego)
  if (nrow(df) == 0) {
    warning("No enriched GO terms found for: ", direction)
    return(df)
  }
  df$direction <- direction
  df %>% arrange(p.adjust) %>% slice_head(n = top_n)
}

top_up   <- run_go(up_genes,   "Up")
top_down <- run_go(down_genes, "Down")

go_df <- bind_rows(top_up, top_down)
if (nrow(go_df) == 0) stop("No enriched GO terms in either direction.")

# Prepare plotting data -------------------------------------------------------
go_df <- go_df %>%
  mutate(CountSigned      = ifelse(direction == "Down", -Count, Count),
         logPadj          = -log10(p.adjust),
         Description_label = paste0(Description, " (", direction, ")"))

# Order the y-axis: Up terms by count, Down terms by (negative) count
up_df   <- go_df %>% filter(direction == "Up")   %>% arrange(desc(Count), p.adjust)
down_df <- go_df %>% filter(direction == "Down") %>% arrange(Count)

go_df_sorted <- bind_rows(up_df, down_df)
go_df_sorted$Description_label <- factor(go_df_sorted$Description_label,
                                         levels = go_df_sorted$Description_label)

x_max <- max(abs(go_df_sorted$CountSigned))   # symmetric x-axis

# Plot ------------------------------------------------------------------------
p <- ggplot(go_df_sorted,
            aes(x = CountSigned, y = Description_label, fill = logPadj)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 0, color = "black") +
  scale_x_continuous(limits = c(-x_max, x_max)) +
  scale_fill_gradient(low = "skyblue", high = "steelblue") +
  labs(title = plot_title,
       x = "Gene count (negative = down, positive = up)",
       y = "GO term",
       fill = "-log10(p.adjust)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p)

# Save ------------------------------------------------------------------------
dir.create(output_dir, showWarnings = FALSE)
ggsave(file.path(output_dir, paste0(output_name, ".png")), p,
       width = 9, height = 6, dpi = 300)
ggsave(file.path(output_dir, paste0(output_name, ".pdf")), p,
       width = 9, height = 6)
write.csv(go_df_sorted %>% select(-Description_label),
          file.path(output_dir, paste0(output_name, "_terms.csv")),
          row.names = FALSE)
