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

# confirm before uploading
head(readRDS("pkn_NK.rds"))   # source | interaction | target


# ---- Why the KO run failed ----
# CARNIVAL silently dropped PRF1 from the measurements,
# along with ~78 of the 192 TFs. So the baseline and KO runs solved the identical
# problem and returned the same objective.
#
# From carnival_pilot.err:
#   Warning: These measurement nodes are not in prior knowledge network and
#   will be ignored: ARID3A | ATF6 | ... | TBX21 | ... | ZNF148 | PRF1
#
# The reason is that PRF1 has no outgoing edges, so CARNIVAL cannot place it
# on a path from a perturbation to a measurement.

pkn_list[[ct]] %>% dplyr::filter(target == "PRF1")   # 11 incoming edges
pkn_list[[ct]] %>% dplyr::filter(source == "PRF1")   # 0 outgoing edges

# The dropped TFs are sinks for the same reason
sinks_nk <- setdiff(pkn_list[[ct]]$target, pkn_list[[ct]]$source)
c("PRF1", "TBX21", "GATA3", "RUNX1", "RELB", "TCF7") %in% sinks_nk