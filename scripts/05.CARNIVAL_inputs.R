# CARNIVAL inputs
# CARNIVAL renames PKN columns positionally, so order must be source, interaction, target.

ko_genes  <- c("PRF1", "RAB27A", "STX11", "STXBP2", "SH2D1A", "XIAP")
celltypes <- c("NK cells", "gdT cells", "CD8+ T cells")
tags      <- c("NK", "gdT", "CD8T")

# ---- Pruned, cell-type-specific ----
for (i in seq_along(celltypes)) {
  
  ct  <- celltypes[i]
  tag <- tags[i]
  
  pkn_ct <- pkn_list[[ct]]
  nodes  <- c(pkn_ct$source, pkn_ct$target)
  
  saveRDS(pkn_ct, sprintf("pkn_%s.rds", tag))
  
  tf_meas <- df %>%
    dplyr::filter(celltype == ct, source %in% nodes) %>%
    dplyr::select(source, mean) %>%
    tibble::deframe()
  
  saveRDS(tf_meas, sprintf("%s_baseline_measurements.rds", tag))
  
  ko_val <- min(tf_meas)
  
  for (gene in ko_genes) {
    if (!gene %in% nodes) { cat(tag, ":", gene, "skipped\n"); next }
    meas_ko <- c(tf_meas[names(tf_meas) != gene], setNames(ko_val, gene))
    saveRDS(meas_ko, sprintf("%s_ko%s_measurements.rds", tag, gene))
  }
  
  cat(tag, ":", nrow(pkn_ct), "edges,", length(tf_meas), "TFs\n")
}

# ---- Unpruned, gdT / PRF1 (pruning test) ----
pkn_nodes <- unique(c(pkn$source, pkn$target))

tf_meas <- df %>%
  dplyr::filter(celltype == "gdT cells", source %in% pkn_nodes) %>%
  dplyr::select(source, mean) %>%
  tibble::deframe()

ko_val  <- min(tf_meas)
meas_ko <- c(tf_meas[names(tf_meas) != "PRF1"], setNames(ko_val, "PRF1"))

saveRDS(pkn,     "pkn_full.rds")
saveRDS(tf_meas, "gdTfull_baseline_measurements.rds")
saveRDS(meas_ko, "gdTfull_koPRF1_measurements.rds")

# ---- Validation ----
sapply(list.files(pattern = "^pkn_.*\\.rds$"),
       function(f) paste(colnames(readRDS(f)), collapse = ", "))

for (tag in c("NK", "gdT", "CD8T")) {
  p    <- readRDS(sprintf("pkn_%s.rds", tag))
  base <- readRDS(sprintf("%s_baseline_measurements.rds", tag))
  cat(tag, ":", sum(names(base) %in% unique(c(p$source, p$target))),
      "of", length(base), "in PKN\n")
}