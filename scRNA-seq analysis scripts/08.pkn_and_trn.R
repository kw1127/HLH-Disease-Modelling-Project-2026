# =============================================================================
# 13. Prior knowledge network and transcriptional layer
#
# Builds the two networks carnival needs.
#
#   pkn - a signed, directed protein signalling network from OmniPath. 
#
#   trn - CollecTRI edges from transcription factors onto the HLH-associated
#         genes. This is joined afterwards.
#
# Also records why each HLH gene can or cannot appear in the final model.
# =============================================================================
tags <- c("nk", "temra")   # the two contrasts, used as a suffix throughout


# Signalling layer from OmniPath
ppi <- omnipath_interactions()

sig_all <- ppi %>%
  filter(consensus_direction == 1, # agreed direction source -> target
         consensus_stimulation + consensus_inhibition == 1) %>% # unambiguously signed
  mutate(interaction = ifelse(consensus_stimulation == 1, 1L, -1L)) %>%
  select(source = source_genesymbol, interaction,
         target = target_genesymbol, curation_effort) %>%
  filter(source != "", # unmapped symbols would collapse into one phantom hub
         target != "",
         source != target) %>% # self-loops: one state per node, so vacuous or contradictory
  group_by(source, interaction, target) %>%
  summarise(curation_effort = max(curation_effort), .groups = "drop")

nrow(sig_all)   # 70,565 unique signed directed edges


# Choosing the curation and expression thresholds
#
# Two candidate gene universes:
#   strict - the DESeq2 filter (>=5 counts in >= min_grp samples). Designed for
#            stable dispersion estimates, so it is stricter than a presence
#            filter needs to be and drops low-abundance signalling proteins.
#   lax - a presence filter (>=1 count in >=5 lymphoid pseudobulk samples).
#
# tf_all below is the set of TFs used to score each candidate threshold. It is
# the top-50 significant TFs per contrast, which is the significance-only
# measurement set. The final measurement set (section 14) is a superset of this,
# since it also forces in the CollecTRI regulators of the HLH genes.
expressed_ct  <- rownames(mat_ct)
expressed_lax <- rownames(pseudo_sub)[rowSums(pseudo_sub >= 1) >= 5]
c(strict = length(expressed_ct), lax = length(expressed_lax))

tf_all <- tf_contrast %>%
  filter(p_adj < 0.05) %>%
  group_by(condition) %>%
  slice_max(abs(score), n = 50) %>%
  ungroup() %>%
  pull(source) %>%
  unique()

check <- function(ce, universe, label) {
  s <- sig_all %>%
    filter(curation_effort >= ce,
           source %in% universe, target %in% universe)
  nodes  <- unique(c(s$source, s$target))
  in_deg <- table(s$target)
  reach  <- intersect(tf_all, nodes)
  data.frame(universe = label, curation = ce,
             edges = nrow(s), nodes = length(nodes),
             tf_present = length(reach), 
             tf_with_in = sum(reach %in% names(in_deg))) # a TF with no incoming edge can never be explained by any solution
}

grid <- do.call(rbind, c(
  lapply(1:4, check, universe = expressed_ct,  label = "strict"),
  lapply(1:4, check, universe = expressed_lax, label = "lax")))

grid

# Export as image for thesis
grid_df <- as.data.frame(grid, stringsAsFactors = FALSE)

# if rbind produced a character matrix, restore the numeric columns
grid_df <- type.convert(grid_df, as.is = TRUE)

# rownames from rbind aren't a column; make them one if they carry information
grid_df <- tibble::rownames_to_column(grid_df, "row")   # skip if rownames are 1..8

grid_tbl <- grid_df |>
  gt(groupname_col = "label") |>        
  fmt_number(columns = where(is.numeric), decimals = 3) |>
  cols_align(align = "center", columns = everything()) |>
  tab_options(table.font.size = 11, data_row.padding = px(4))

gtsave(grid_tbl, "table_grid.png", vwidth = 700, expand = 5)

# Decision: lax universe, curation_effort >= 3.
#   The expression filter dominates: strict excludes 6-8 reachable TFs at every
#   curation level (32 vs 40 at curation effort >=1).
#   Within the lax universe, >=3 is the elbow: relaxing to >=2 buys one extra TF
#   for 2,331 extra edges; tightening to >=4 saves 1,283 edges but costs two TFs.

sig <- sig_all %>%
  filter(curation_effort >= 3,
         source %in% expressed_lax, target %in% expressed_lax) %>%
  select(source, interaction, target)

nrow(sig)

stopifnot(nrow(dplyr::filter(dplyr::count(sig, source, target), n > 1)) == 0)

pkn <- as.data.frame(sig)
pkn_nodes <- unique(c(pkn$source, pkn$target))
saveRDS(pkn, "pkn_carnival.rds")


# Transcriptional layer and HLH differential expression
#
# Every CollecTRI regulator of an HLH gene is retained, including those absent
# from the signalling layer, so the attrition can be counted in 13.4 rather than
# disappearing into a filter.
trn <- collectri %>%
  filter(target %in% hlh_chr) %>%
  mutate(interaction = as.integer(mor)) %>%
  select(source, interaction, target) %>%
  distinct()

saveRDS(trn, "trn.rds")

# The genome-wide DESeq2 results from section 11, subset to the HLH genes. This
# defines which genes are eligible to appear in the downstream layer at all.
res_list <- list(nk = res_nk, temra = res_temra)

hlh_de <- lapply(setNames(tags, tags), function(tag)
  as.data.frame(res_list[[tag]]) %>%
    rownames_to_column("target") %>%
    filter(target %in% hlh_chr, !is.na(padj), padj < 0.05) %>%
    select(target, log2FoldChange, padj))

sapply(hlh_de, nrow)   # how many of the 8 survive per contrast

saveRDS(hlh_de, "hlh_de.rds")
saveRDS(tf_contrast, "tf_contrast.rds")


# Audit: how can each HLH gene enter the model?
#
# Only downstream, via a CollecTRI regulator that CARNIVAL recovers.
# in_signalling is reported because STXBP2 and STX11 carry protein-level
# OmniPath edges that this design deliberately does not use.
g_sig <- graph_from_data_frame(select(sig, source, target), directed = TRUE)

hlh_audit <- tibble(gene = hlh_chr) %>%
  mutate(
    n_tf_regulators = vapply(gene, function(g) sum(trn$target == g), integer(1)),
    n_regs_in_pkn   = vapply(gene, function(g)
      sum(trn$target == g & trn$source %in% pkn_nodes), integer(1)),
    in_signalling = gene %in% pkn_nodes,
    de_nk = gene %in% hlh_de$nk$target, # hlh_de is filtered to padj < 0.05, so non-DE genes are absent by design
    de_temra = gene %in% hlh_de$temra$target,
    can_enter = case_when(
      n_regs_in_pkn > 0 ~ "yes, via TRN join",
      n_tf_regulators > 0 ~ "regulators exist but none in PKN",
      TRUE ~ "no CollecTRI regulator"))

hlh_audit

audit_tbl <- hlh_audit %>%
  select(gene, n_tf_regulators, n_regs_in_pkn, in_signalling, de_nk, de_temra) %>%
  gt() %>%
  cols_label(gene = "Gene",
             n_tf_regulators = "CollecTRI regulators",
             n_regs_in_pkn = "Regulators in PKN",
             in_signalling = "In signalling layer",
             de_nk = "DE in NK",
             de_temra = "DE in CD8 TEMRA") %>%
  fmt(columns = c(in_signalling, de_nk, de_temra),
      fns = function(x) ifelse(x, "Yes", "No")) %>%
  tab_footnote("Genes with no CollecTRI regulator cannot enter the model at any
                stage. STX11 and STXBP2 carry protein-level OmniPath edges, so
                appear in the signalling layer, but this design does not use
                them as signalling nodes.") %>%
  tab_options(table.font.size = 11)

gtsave(audit_tbl, "Table_hlh_coverage.png", vwidth = 900, expand = 10)


# Evidence that the exclusions are structural, not artefacts
#
# PRF1 has no usable outgoing edges anywhere in OmniPath at any curation level:
# two of its four are unsigned, two have zero references.
all_int <- import_all_interactions()

all_int %>%
  filter(source_genesymbol == "PRF1") %>%
  select(target_genesymbol, is_directed, is_stimulation,
         is_inhibition, curation_effort, n_references)

# usable outgoing edges per HLH gene, across all OmniPath layers
all_int %>%
  dplyr::filter(source_genesymbol %in% hlh_chr,
                is_directed == 1,
                is_stimulation + is_inhibition == 1) %>%
  dplyr::count(source_genesymbol, name = "usable_outgoing") %>%
  dplyr::right_join(tibble::tibble(source_genesymbol = hlh_chr),
                    by = "source_genesymbol") %>%
  dplyr::mutate(usable_outgoing = tidyr::replace_na(usable_outgoing, 0L)) %>%
  dplyr::arrange(dplyr::desc(usable_outgoing))

# genes with no regulators in CollecTRI cannot enter at all
collectri %>% 
  dplyr::filter(target %in% hlh_chr) %>% 
  dplyr::count(target)

setdiff(hlh_chr, collectri$target)

# PRF1's regulators, and whether they sit in the signalling layer
trn %>%
  dplyr::filter(target == "PRF1") %>%
  dplyr::mutate(regulator_in_pkn = source %in% pkn_nodes)

# End-to-end reachability stops at the TF, not the gene.
prf1_tfs <- trn %>%
  dplyr::filter(target == "PRF1", source %in% pkn_nodes) %>%
  dplyr::pull(source)

receptors <- intersect(c("IL12RB1", "IL12RB2", "IFNGR1", "IFNGR2", "IL2RB", "IL18R1"),
                       V(g_sig)$name)

distances(g_sig, v = receptors,
          to = intersect(prf1_tfs, V(g_sig)$name), mode = "out")

# canonical JAK-STAT route, receptor -> TF
# (substitute whichever regulator prf1_tfs actually returns)
V(g_sig)$name[shortest_paths(g_sig, from = "IL12RB1", to = "STAT4",
                             mode = "out")$vpath[[1]]]


# Cell-type-specific regulation of PRF1
#
# Used ONLY to ask which regulators differ between cell types, not to build the
# networks that get solved. min_pct is a marker-gene threshold; with scRNA-seq
# dropout it is far too strict to use as a presence filter, and the expression
# filter for the PKN was already applied at the pseudobulk level.
expressed_in <- function(ct, min_pct = 0.05) {
  cells <- colnames(pbmc.clean)[pbmc.clean$celltype == ct]
  cnt <- GetAssayData(pbmc.clean, assay = "RNA", layer = "counts")[, cells]
  rownames(cnt)[Matrix::rowMeans(cnt > 0) >= min_pct]
}

genes_nk <- expressed_in("NK")
genes_temra <- expressed_in("CD8 TEMRA")

prf1_regs <- trn %>%
  dplyr::filter(target == "PRF1") %>%
  dplyr::mutate(in_nk = source %in% genes_nk, in_temra = source %in% genes_temra)

prf1_regs

# regulators detected in NK but not CD8 TEMRA
prf1_regs %>% 
  dplyr::filter(in_nk, !in_temra) %>% 
  dplyr::pull(source)

# shared vs specific, both layers — these are your own helpers
sig_of <- function(x) x$nodes$id[x$nodes$layer %in% c("signalling","perturbation")]
list(shared = intersect(sig_of(pruned$anchored$nk), sig_of(pruned$anchored$temra)),
     nk_only = setdiff(sig_of(pruned$anchored$nk), sig_of(pruned$anchored$temra)),
     temra_only = setdiff(sig_of(pruned$anchored$temra), sig_of(pruned$anchored$nk)))

# TF -> pHLH edges
trn_of <- function(x) unique(paste(x$edges$source[x$edges$layer=="trn"],
                                   x$edges$target[x$edges$layer=="trn"], sep=" -> "))
list(shared = intersect(trn_of(pruned$anchored$nk), trn_of(pruned$anchored$temra)),
     nk_only = setdiff(trn_of(pruned$anchored$nk), trn_of(pruned$anchored$temra)),
     temra_only = setdiff(trn_of(pruned$anchored$temra), trn_of(pruned$anchored$nk)))
