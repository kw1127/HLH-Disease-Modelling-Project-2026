# =============================================================================
# 14. CARNIVAL measurement input
#
# Two measurement sets are written, differing only in which TFs they contain.
# Both are solved over the same PKN, so the difference between the resulting
# networks isolates what the HLH anchoring contributed.
#
#   baseline - the most confident TFs by activity score. CARNIVAL is an ILP
#                  and runtime scales with the number of measurements, so this
#                  is capped at the top n by |score| among those clearing FDR.
#   hlh_anchored - the same TFs plus every CollecTRI regulator of an
#                  HLH-associated gene, included by design regardless of score
#                  or FDR. Without these the solver has no reason to reach the
#                  transcriptional layer, and the TRN join would be decorative:
#                  appended after solving rather than constrained during it.
#
# Forced TFs keep their observed scores. Nothing is inflated to make the solver
# attend to them; the origin column records which are which so the weakness of
# each can be reported rather than hidden.
# =============================================================================

# Build the measurement objects
hlh_tfs <- unique(trn$source)   # CollecTRI regulators of the 8 HLH genes

make_measobj <- function(contrast, n_top = 50, force = hlh_tfs) {
  d <- tf_contrast %>% dplyr::filter(condition == contrast)
  
  top <- d %>%
    dplyr::filter(p_adj < 0.05) %>%
    dplyr::slice_max(abs(score), n = n_top) %>%
    dplyr::mutate(origin = "top_significant")
  
  # a TF that is both significant and an HLH regulator appears once, as
  # top_significant, so the two sets are a union rather than a sum
  forced <- d %>%
    dplyr::filter(source %in% force, !source %in% top$source) %>%
    dplyr::mutate(origin = "forced_hlh_regulator")
  
  dplyr::bind_rows(top, forced) %>%
    dplyr::select(source, score, p_adj, origin)
}

meas_long <- dplyr::bind_rows(
  make_measobj("NK")        %>% 
    dplyr::mutate(contrast = "nk"),
  make_measobj("CD8 TEMRA") %>% 
    dplyr::mutate(contrast = "temra"))

# how much did forcing change the measurement set?
meas_long %>% 
  dplyr::count(contrast, origin)


# Audit of the forced TFs
#
# One row per forced TF per contrast. has_in_edge is the critical column: a TF
# with no incoming edge in the PKN cannot be explained by any solution, so it
# sits in the measurement set unreachable and is simply left unset.
pkn_targets <- unique(pkn$target)

tf_aud <- meas_long %>%
  dplyr::filter(origin == "forced_hlh_regulator") %>%
  dplyr::mutate(
    in_pkn      = source %in% pkn_nodes,
    has_in_edge = source %in% pkn_targets,
    hlh_targets = vapply(source, function(s)
      paste(sort(trn$target[trn$source == s]), collapse = ", "), character(1))) %>%
  dplyr::arrange(contrast, dplyr::desc(abs(score)))

tf_aud

# how many forced TFs are structurally reachable at all?
tf_aud %>%
  dplyr::count(contrast, in_pkn, has_in_edge)

saveRDS(tf_aud, "tf_aud.rds")


# Write the CARNIVAL inputs
#
# Restricting to pkn_nodes is the expression filter doing its work: the PKN was
# built on expressed_lax, so a TF not detected in these cell types cannot be
# asserted to be signalling in them. Forcing overrides the significance bar,
# not the expression bar.
#
# check.names = FALSE: the default would rewrite hyphenated gene symbols
# (NKX3-1 -> NKX3.1) and silently break the match against the PKN.
write_meas <- function(tag, variant) {
  m <- meas_long %>% dplyr::filter(contrast == tag)
  if (variant == "baseline")
    m <- m %>% dplyr::filter(origin == "top_significant")
  
  n_forced_before <- sum(m$origin == "forced_hlh_regulator")
  m <- m %>% dplyr::filter(source %in% pkn_nodes)
  n_forced_after  <- sum(m$origin == "forced_hlh_regulator")
  
  tf_v <- setNames(m$score, m$source)
  saveRDS(as.data.frame(as.list(tf_v), check.names = FALSE),
          sprintf("meas_%s_%s.rds", tag, variant))
  
  cat(sprintf("%-6s %-13s %2d measurements, range [%.1f, %.1f]; forced %d -> %d after PKN filter\n",
              tag, variant, length(tf_v), min(tf_v), max(tf_v),
              n_forced_before, n_forced_after))
}

for (tag in tags)
  for (variant in c("baseline", "hlh_anchored"))
    write_meas(tag, variant)


# Validation
for (tag in tags) {
  for (variant in c("baseline", "hlh_anchored")) {
    p <- readRDS("pkn_carnival.rds")
    m <- readRDS(sprintf("meas_%s_%s.rds", tag, variant))
    stopifnot(identical(colnames(p), c("source", "interaction", "target")),
              all(p$interaction %in% c(-1L, 1L)),
              nrow(m) == 1L,
              all(colnames(m) %in% c(p$source, p$target)))
    # a measured TF with no incoming edge cannot be explained by any solution
    cat(sprintf("%-6s %-13s OK: %d edges, %d measurements, %d with an in-edge\n",
                tag, variant, nrow(p), ncol(m),
                sum(colnames(m) %in% p$target)))
  }
}
