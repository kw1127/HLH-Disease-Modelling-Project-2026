# Build CARNIVAL inputs for HPC

# Design: inverse CARNIVAL (no known perturbation -- healthy resting cells).
#   Baseline = measured TF activities
#   Knockout = measured TF activities + HLH gene fixed at loss-of-function

# HLH genes present in the coupled network (UNC13D absent)
ko_genes <- hlh[hlh %in% unique(unlist(lapply(pkn_list, function(p) 
  c(p$source, p$target))))]

# --- PILOT: NK cells, PRF1 knockout ---
ct <- "NK cells"

tf_meas <- df %>%
  dplyr::filter(celltype == ct,
                source %in% c(pkn_list[[ct]]$source, pkn_list[[ct]]$target)) %>%
  dplyr::select(source, mean) %>%
  tibble::deframe()

ko_val  <- min(tf_meas)
meas_ko <- c(tf_meas, setNames(ko_val, "PRF1"))

saveRDS(pkn_list[[ct]], "pkn_NK.rds")           
saveRDS(tf_meas, "NK_baseline_measurements.rds")
saveRDS(meas_ko, "NK_koPRF1_measurements.rds")

sinks_full <- setdiff(pkn$target, pkn$source)                    # coupled, 14,989 edges
sinks_nk   <- setdiff(pkn_list[["NK cells"]]$target,
                      pkn_list[["NK cells"]]$source)             # NK-pruned, 4,253 edges

length(sinks_full)
length(sinks_nk)
length(setdiff(sinks_nk, sinks_full))    # sinks created BY pruning

# the key question: were these already sinks before pruning?
c("PRF1","TBX21","GATA3","RUNX1","RELB","TCF7") %in% sinks_full