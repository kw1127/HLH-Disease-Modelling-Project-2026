# Build CARNIVAL inputs for HPC
# In: 04_pkn_list.rds, 03_tf_by_celltype.rds
# Out: outputs/carnival/*.rds
#
# Design: inverse CARNIVAL (no known perturbation -- healthy resting cells).
#   Baseline = measured TF activities
#   Knockout = measured TF activities + HLH gene fixed at loss-of-function

source("scripts/00_setup.R")

pkn_list <- readRDS(file.path(data_proc, "04_pkn_list.rds"))
df <- readRDS(file.path(data_proc, "03_tf_by_celltype.rds"))

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

length(tf_meas)

# Scale the KO to the TF range
ko_val <- min(tf_meas)
meas_ko <- c(tf_meas, setNames(ko_val, "PRF1"))

saveRDS(list(pkn = pkn_list[[ct]], meas = tf_meas),
        file.path(carn_dir, "baseline_NK.rds"))

saveRDS(list(pkn = pkn_list[[ct]], meas = meas_ko),
        file.path(carn_dir, "ko_PRF1_NK.rds"))