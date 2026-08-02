# ============================================================
# 13. Signalling network from OmniPath as prior knowledge network
# 
#     TF activities are the sole measurements
#     The transcriptional layer is joined after solving
#  
# ============================================================

# CARNIVAL is an ILP; runtime scales with the number of
# measurements, so take the most confident TFs rather than everything
# that clears FDR.
make_measobj <- function(contrast, n_top = 50) {
  tf_contrast %>%
    dplyr::filter(condition == contrast, p_adj < 0.05) %>%
    dplyr::slice_max(abs(score), n = n_top) %>%
    dplyr::select(source, score) %>%
    tidyr::pivot_wider(names_from = source, values_from = score) %>%
    as.data.frame()
}

meas_temra <- make_measobj("CD8 TEMRA")
meas_nk <- make_measobj("NK")

tf_all <- union(colnames(meas_temra), colnames(meas_nk))

# Signalling layer from OmniPath
ppi <- omnipath_interactions()      

sig_all <- ppi %>%
  dplyr::filter(consensus_direction == 1, # agreed direction source -> target
                consensus_stimulation + consensus_inhibition == 1) %>% # unambiguously signed
  dplyr::mutate(interaction = ifelse(consensus_stimulation == 1, 1L, -1L)) %>%
  dplyr::select(source = source_genesymbol, interaction,
                target = target_genesymbol, curation_effort) %>%
  dplyr::filter(source != "", # unmapped symbols would collapse into one phantom hub
                target != "",
                source != target) %>% # self-loops: one state per node, so vacuous or contradictory
  dplyr::group_by(source, interaction, target) %>%
  dplyr::summarise(curation_effort = max(curation_effort), .groups = "drop")

# nrow(sig_all) = 70,565 unique signed directed edges


# Choosing the curation and expression thresholds
# Two candidate gene universes:
#   strict = the DESeq2 filter (>=5 counts in >= min_grp samples). Designed
#   for stable dispersion estimates, so it is stricter than a
#   presence filter needs to be and drops low-abundance
#   signalling proteins.
#   lax = presence filter (>=1 count in >=5 lymphoid pseudobulk samples).
expressed_ct <- rownames(mat_ct)
expressed_lax <- rownames(pseudo_sub)[rowSums(pseudo_sub >= 1) >= 5]
c(strict = length(expressed_ct), lax = length(expressed_lax))

check <- function(ce, universe, label) {
  
  s <- sig_all %>%
    dplyr::filter(curation_effort >= ce,
                  source %in% universe, target %in% universe)
  nodes  <- unique(c(s$source, s$target))
  in_deg <- table(s$target)
  reach  <- intersect(tf_all, nodes)
  data.frame(universe = label, curation = ce,
             edges = nrow(s), nodes = length(nodes),
             tf_present = length(reach),
             # a TF with no incoming edge can never be explained by any solution
             tf_with_in = sum(reach %in% names(in_deg)))
}

grid <- do.call(rbind, c(
  lapply(1:4, check, universe = expressed_ct,  label = "strict"),
  lapply(1:4, check, universe = expressed_lax, label = "lax")))

grid

# Decision: lax universe, curation_effort >= 3.
#   The expression filter dominates: strict excludes 9-10 reachable TFs at
#   every curation level (32 vs 42 at curation effort >=1).
#   Within the lax universe, curation effort >=3 is the elbow: relaxing to >=2 buys one
#   extra TF for 2,331 extra edges; tightening to >=4 saves 1,283 edges but
#   costs two TFs.
sig <- sig_all %>%
  dplyr::filter(curation_effort >= 3,
                source %in% expressed_lax, target %in% expressed_lax) %>%
  dplyr::select(source, interaction, target)

nrow(sig)

stopifnot(nrow(dplyr::filter(dplyr::count(sig, source, target), n > 1)) == 0)

pkn <- as.data.frame(sig)
pkn_nodes <- unique(c(pkn$source, pkn$target))
saveRDS(pkn, "pkn_carnival.rds")

# Downstream transcriptional layer separate from the PPI pkn.
# Every CollecTRI regulator of an HLH gene is retained, including ones absent
# from the signalling layer, so the attrition can be counted rather than
# silently applied by a filter.
trn <- collectri %>%
  dplyr::filter(target %in% hlh_chr) %>%
  dplyr::mutate(interaction = as.integer(mor)) %>%
  dplyr::select(source, interaction, target) %>%
  dplyr::distinct()

saveRDS(trn, "trn.rds")

# Differential expression of the HLH genes.
# The genome-wide DESeq2 results from section 11 subset
# to the HLH genes of interest. This defines which genes are eligible to
# appear in the downstream layer.
res_list <- list(nk = res_nk, temra = res_temra)
tags <- c("nk", "temra")

hlh_de <- lapply(setNames(tags, tags), function(tag)
  as.data.frame(res_list[[tag]]) %>%
    tibble::rownames_to_column("target") %>%
    dplyr::filter(target %in% hlh_chr, !is.na(padj), padj < 0.05) %>%
    dplyr::select(target, log2FoldChange, padj))

sapply(hlh_de, nrow)   # how many of the 8 survive per contrast

saveRDS(hlh_de, "hlh_de.rds")
saveRDS(tf_contrast, "tf_contrast.rds")

# Audit: how can each HLH gene enter the model?
# It can only enter downstream, via a CollecTRI regulator that CARNIVAL
# recovers. in_signalling is reported because STXBP2 and STX11 carry
# protein-level OmniPath edges that this design deliberately does not use.
g_sig <- igraph::graph_from_data_frame(
  sig %>% dplyr::select(source, target), directed = TRUE)

hlh_audit <- tibble::tibble(gene = hlh_chr) %>%
  dplyr::mutate(
    n_tf_regulators = vapply(gene, function(g) sum(trn$target == g), integer(1)),
    n_regs_in_pkn = vapply(gene, function(g)
      sum(trn$target == g & trn$source %in% pkn_nodes), integer(1)),
    in_signalling = gene %in% pkn_nodes,
    de_nk = gene %in% hlh_de$nk$target,
    de_temra = gene %in% hlh_de$temra$target,
    can_enter = dplyr::case_when(
      n_regs_in_pkn > 0 ~ "yes, via TRN join",
      n_tf_regulators > 0 ~ "regulators exist but none in PKN",
      TRUE ~ "no CollecTRI regulator"))

hlh_audit

# Evidence that the exclusions are structural, not artefacts
# PRF1 has no usable outgoing edges anywhere in OmniPath, at any curation
# level: two of its four are unsigned, two have zero references.
all_int <- OmnipathR::import_all_interactions()

all_int %>%
  dplyr::filter(source_genesymbol == "PRF1") %>%
  dplyr::select(target_genesymbol, is_directed, is_stimulation,
                is_inhibition, curation_effort, n_references)

# usable outgoing edges per HLH gene across ALL OmniPath layers
all_int %>%
  dplyr::filter(source_genesymbol %in% hlh_chr,
                is_directed == 1,
                is_stimulation + is_inhibition == 1) %>%
  dplyr::count(source_genesymbol, name = "usable_outgoing") %>%
  dplyr::right_join(tibble::tibble(source_genesymbol = hlh_chr),
                    by = "source_genesymbol") %>%
  dplyr::mutate(usable_outgoing = tidyr::replace_na(usable_outgoing, 0L)) %>%
  dplyr::arrange(dplyr::desc(usable_outgoing))

# Genes with no regulators in CollecTRI cannot enter at all
collectri %>%
  dplyr::filter(target %in% hlh_chr) %>%
  dplyr::count(target)

setdiff(hlh_chr, collectri$target)

# PRF1's regulators, and whether they sit in the signalling layer
trn %>%
  dplyr::filter(target == "PRF1") %>%
  dplyr::mutate(regulator_in_pkn = source %in% pkn_nodes)

# End-to-end reachability now stops at the TF, not the gene
prf1_tfs <- trn %>%
  dplyr::filter(target == "PRF1", source %in% pkn_nodes) %>%
  dplyr::pull(source)

receptors <- intersect(c("IL12RB1","IL12RB2","IFNGR1","IFNGR2","IL2RB","IL18R1"),
                       igraph::V(g_sig)$name)

igraph::distances(g_sig, v = receptors,
                  to = intersect(prf1_tfs, igraph::V(g_sig)$name), mode = "out")

# canonical JAK-STAT route, receptor -> TF
# (substitute whichever regulator prf1_tfs actually returns)
igraph::V(g_sig)$name[
  igraph::shortest_paths(g_sig, from = "IL12RB1", to = "STAT4",
                         mode = "out")$vpath[[1]]]


# cell-type-specific regulation of PRF1
# Used ONLY to ask which regulators differ between cell types, not to build the
# networks that get solved. min_pct is a marker-gene threshold; with scRNA-seq
# dropout it is far too strict to use as a presence filter. The expression
# filter for the PKN was already applied at the pseudobulk level.
expressed_in <- function(ct, min_pct = 0.05) {
  cells <- colnames(pbmc.clean)[pbmc.clean$celltype == ct]
  cnt <- GetAssayData(pbmc.clean, assay = "RNA", layer = "counts")[, cells]
  rownames(cnt)[Matrix::rowMeans(cnt > 0) >= min_pct]
}

genes_nk <- expressed_in("NK")
genes_temra <- expressed_in("CD8 TEMRA")

prf1_regs <- trn %>% dplyr::filter(target == "PRF1") %>%
  dplyr::mutate(in_nk = source %in% genes_nk, in_temra = source %in% genes_temra)

prf1_regs

# regulators detected in NK but not CD8 TEMRA
prf1_regs %>%
  dplyr::filter(in_nk, !in_temra) %>%
  dplyr::pull(source)


# CARNIVAL inputs: one shared PKN, per-contrast measurements.
# Both cell types are solved over identical topology, so any difference between
# the solved networks reflects the measurements, not the network.
tf_list <- list(nk = meas_nk, temra = meas_temra)

for (tag in tags) {
  tf_v <- setNames(as.numeric(tf_list[[tag]][1, ]), colnames(tf_list[[tag]]))
  tf_v <- tf_v[names(tf_v) %in% pkn_nodes]
  
  # check.names = FALSE: the default would rewrite hyphenated gene symbols
  # (NKX3-1 -> NKX3.1) and silently break the match against the PKN.
  saveRDS(as.data.frame(as.list(tf_v), check.names = FALSE),
          sprintf("meas_%s_baseline.rds", tag))
  
  cat(sprintf("%-6s %d TF measurements, range [%.1f, %.1f]\n",
              tag, length(tf_v), min(tf_v), max(tf_v)))
}

# validation
for (tag in tags) {
  p <- readRDS("pkn_carnival.rds")
  m <- readRDS(sprintf("meas_%s_baseline.rds", tag))
  stopifnot(identical(colnames(p), c("source", "interaction", "target")),
            all(p$interaction %in% c(-1L, 1L)),
            nrow(m) == 1L,
            all(colnames(m) %in% c(p$source, p$target)))
  # a measured TF with no incoming edge cannot be explained by any solution
  cat(sprintf("%-6s OK: %d edges, %d measurements, %d with an in-edge\n",
              tag, nrow(p), ncol(m), sum(colnames(m) %in% p$target)))
}