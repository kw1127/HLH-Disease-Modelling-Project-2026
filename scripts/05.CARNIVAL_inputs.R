# BUild CARNIVAL inputs
# NOTE: CARNIVAL's preprocessPriorKnowledgeNetwork() renames PKN columns
# positionally -- colnames(pkn) <- c("Node1", "Sign", "Node2") -- so the
# column order must be source, interaction, target. Wrong column order scrambles
# the network silently, with no error being thrown.

ko_genes <- c("PRF1", "RAB27A", "STX11", "STXBP2", "SH2D1A", "XIAP")
# UNC13D and LYST have no regulators in CollecTRI, so they cannot enter the network

celltypes <- c("NK cells", "gdT cells", "CD8+ T cells")
tags <- c("NK", "gdT", "CD8T")

for (i in seq_along(celltypes)) {
  
  ct  <- celltypes[i]
  tag <- tags[i]
  
  pkn_ct <- pkn_list[[ct]]
  nodes  <- c(pkn_ct$source, pkn_ct$target)
  
  saveRDS(pkn_ct, sprintf("pkn_%s.rds", tag))
  
  # Baseline TF activities
  tf_meas <- df %>%
    dplyr::filter(celltype == ct, source %in% nodes) %>%
    dplyr::select(source, mean) %>%
    tibble::deframe()
  
  saveRDS(tf_meas, sprintf("%s_baseline_measurements.rds", tag))
  
  # Knockout value scaled to the TF range
  ko_val <- min(tf_meas)
  
  # One knockout per HLH gene present in this cell type
  for (gene in ko_genes) {
    
    if (!gene %in% nodes) {
      cat(tag, ":", gene, "not in network, skipped\n")
      next
    }
    
    meas_ko <- c(tf_meas[names(tf_meas) != gene],
                 setNames(ko_val, gene))
    
    saveRDS(meas_ko, sprintf("%s_ko%s_measurements.rds", tag, gene))
  }
  
  cat(tag, ":", nrow(pkn_ct), "edges,", length(tf_meas), "TFs\n")
}

# Confirm column order before uploading
sapply(list.files(pattern = "^pkn_.*\\.rds$"), 
       function(f) paste(colnames(readRDS(f)), collapse = ", "))

list.files(pattern = "measurements\\.rds$")