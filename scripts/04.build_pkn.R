# Prior knowledge network: signalling + transcriptional layers
Idents(seurat_filt) <- "celltype_final"

# Import prior knowledge network from OmniPath
ppi <- omnipath_interactions() # 85,217 interactions

# Signalling layer includes directed and unambiguously signed interactions.
# Curation_effort >= 2 was chosen.
# >=3 loses key interactions involved in HLH gene regulation.
# >=1 gives ~70k edges which is too big for solving a network.

# Build a signed, directed prior knowledge network (PKN) for CARNIVAL from OmniPath PPIs
sig <- ppi %>%
  dplyr::filter(consensus_direction == 1, # keep interactions with an agreed direction (source -> target)
                consensus_stimulation + consensus_inhibition == 1, # keep only signed edges
                curation_effort >= 2) %>% # keep only edges supported by >= 2 curation sources
  dplyr::mutate(interaction = ifelse(consensus_stimulation == 1, 1, -1)) %>% # encode sign as CARNIVAL expects: +1 activation, -1 inhibition
  dplyr::select(source = source_genesymbol, # regulator gene (edge start)
                target = target_genesymbol, # target gene (edge end)
                interaction) %>% # the signed edge weight
  dplyr::distinct() # remove duplicate edges

# nrow(sig) returns 14,947 interactions survived

# Can HLH genes be used as CARNIVAL inputs?
# Which are in the signalling layer at all?
hlh[hlh %in% c(sig$source, sig$target)]  # STX11, STXBP2, SH2D1A, XIAP

# Which are absent entirely?
setdiff(hlh, c(sig$source, sig$target)) # PRF1, UNC13D, RAB27A, LYST

# Which have outgoing edges? A source node needs these for CARNIVAL to propagate a perturbation
sig %>% 
  dplyr::filter(source %in% hlh) %>% 
  dplyr::count(source) 

# STXBP2 has 1 outgoing edges and XIAP has 5

# Which are sinks? Sinks have regulators (incoming edges) but regulate nothing (no outgoing edges)
intersect(hlh, setdiff(sig$target, sig$source)) # STX11 and SH2D1A

# Of the two with outgoing edges, do they reach any of the measured TFs?
# A source that reaches no target contributes nothing
g_sig <- igraph::graph_from_data_frame(sig %>% 
                                         dplyr::select(source, target),
                                       directed = TRUE)

# Measured TFs present in the network: NK-cell TFs that also exist as nodes in g_sig
tfs_measured_nk <- df %>%
  dplyr::filter(celltype == "NK cells") %>% # NK TF activity results
  dplyr::pull(source) %>% # the TF names
  intersect(igraph::V(g_sig)$name) # keep only those that are nodes in the PKN

# For each HLH gene that has outgoing edges (STXBP2, XIAP), count how many measured TFs it can reach
for (s in intersect(hlh, sig$source)) {
  d <- igraph::distances(g_sig, v = s, to = tfs_measured_nk, mode = "out") # shortest-path distance, following edge direction
  cat(s, "reaches", sum(is.finite(d)), "of", length(tfs_measured_nk), "TFs\n") # finite distance = a directed path exists; Inf = unreachable
}

# Summary of signalling layer
#   PRF1, UNC13D, RAB27A, LYST are absent from the signed signalling layer
#   STX11, SH2D1A are have no outgoing edges so can't propagate
#   STXBP2 has an outgoing edge but reaches no measured TF
#   XIAP has five outgoing edges and reach 387/486 TFs; the only usable perturbation node

# Is this a filtering artefact, or structural?

# PRF1 is absent from the signalling layer before any filter is applied
"PRF1" %in% ppi$source_genesymbol   # FALSE
"PRF1" %in% ppi$target_genesymbol   # FALSE

# Across ALL OmniPath layers, PRF1 has only 4 outgoing edges, and none are usable
# for CARNIVAL, which needs edges that are both directed and signed
all_int <- OmnipathR::import_all_interactions()

prf1_all <- all_int %>%
  dplyr::filter(source_genesymbol == "PRF1") %>%
  dplyr::select(target_genesymbol, is_directed, is_stimulation,
                is_inhibition, curation_effort, n_references)

#   PRF1 -> CTSB is signed and directed, but curation_effort 0 and 0 references
#   PRF1 -> CDK2 is signed and directed, but curation_effort 0 and 0 references
#   PRF1 -> MTOR is directed, but unsigned
#   PRF1 -> PRKCA is directed, unsigned

# Count usable outgoing edges (directed AND unambiguously signed) for each HLH gene,
# across ALL OmniPath interaction layers.
hlh_all <- all_int %>%
  dplyr::filter(source_genesymbol %in% hlh,
                is_directed == 1,
                is_stimulation + is_inhibition == 1) %>%
  dplyr::count(source_genesymbol, name = "usable_outgoing_edges") %>%
  dplyr::right_join(tibble::tibble(source_genesymbol = hlh), by = "source_genesymbol") %>%
  dplyr::mutate(usable_outgoing_edges = tidyr::replace_na(usable_outgoing_edges, 0)) %>%
  dplyr::arrange(desc(usable_outgoing_edges))

# Usable outgoing edges per HLH gene
edges_in_pkn <- sig %>%
  dplyr::filter(source %in% hlh) %>%
  dplyr::count(source, name = "usable_outgoing_edges") %>%
  dplyr::rename(Gene = source)

# How many measured TFs can each gene reach? 
# Use the full measured-TF panel present in the network, not specific cell types
tfs_panel <- df %>%
  dplyr::pull(source) %>%
  unique() %>%
  intersect(igraph::V(g_sig)$name)

hlh <- as.character(hlh)

reach_counts <- vapply(hlh, function(s) {
  if (!s %in% igraph::V(g_sig)$name) return(0L)
  d <- igraph::distances(g_sig, v = s, to = tfs_panel, mode = "out")
  as.integer(sum(is.finite(d)))
}, integer(1))

# Assemble the full table of all 8 HLH genes
hlh_table <- tibble::tibble(Gene = hlh) %>%
  dplyr::left_join(edges_in_pkn, by = "Gene") %>%
  dplyr::mutate(
    usable_outgoing_edges = tidyr::replace_na(usable_outgoing_edges, 0L),
    tfs_reached           = reach_counts[Gene],
    in_pkn = dplyr::case_when(
      !Gene %in% c(sig$source, sig$target) ~ "Absent",
      Gene %in% sig$target & !Gene %in% sig$source ~ "Sink",
      TRUE ~ "Present"
    ),
    usable_input = ifelse(tfs_reached > 0, "Yes", "No"),
    reason = dplyr::case_when(
      in_pkn == "Absent" ~ "No signed directed edges in the network",
      in_pkn == "Sink" ~ "No outgoing edges",
      tfs_reached == 0 ~ "Has outgoing edges but reaches none of the measured TFs",
      TRUE ~ "Upstream regulator with broad reach to measured TFs"
    )
  ) %>%
  dplyr::arrange(desc(tfs_reached), desc(usable_outgoing_edges)) %>%
  dplyr::select(
    Gene,
    `In PKN` = in_pkn,
    `Usable outgoing edges` = usable_outgoing_edges,
    `Measured TFs reached` = tfs_reached,
    `Usable CARNIVAL input` = usable_input,
    Reason = reason) %>%
  gt() %>%
  tab_header(
    title = "Suitability of HLH-associated genes as CARNIVAL perturbation nodes",
    subtitle = "Edges and reachability computed on the signed signalling PKN"
  ) %>%
  tab_source_note(
    source_note = "Reachability: number of measured transcription factors reachable via a directed path. Genes absent from the PKN or acting as sinks cannot serve as perturbation nodes."
  ) %>%
  gtsave("hlh_carnival_suitability.png")

# ---- Couple signalling with transcriptional layer ----
# Adding TF -> HLH gene edges lets the HLH genes enter the network as
# downstream endpoints instead of inputs
trn <- net %>%
  dplyr::filter(target %in% hlh) %>%
  dplyr::select(source, target, interaction = mor)

pkn <- dplyr::bind_rows(sig, trn) %>% 
  dplyr::distinct() %>%
  dplyr::select(source, interaction, target)

dim(pkn)                  # 14,989 edges

saveRDS(pkn, "pkn_full.rds")

hlh[hlh %in% c(pkn$source, pkn$target)]   # 6 of 8

# UNC13D and LYST have no regulators in CollecTRI either
net %>% 
  dplyr::filter(target %in% c("UNC13D", "LYST"))   # 0 rows

# But PRF1 is still a sink in the coupled network: 11 regulators, 0 targets
pkn %>% 
  dplyr::filter(target == "PRF1")

pkn %>% 
  dplyr::filter(source == "PRF1")   # 0 rows

# Confirm the layers join: PRF1's regulators must be in the signalling layer
prf1_tfs <- net %>% 
  dplyr::filter(target == "PRF1") %>% 
  dplyr::pull(source)

prf1_tfs[prf1_tfs %in% c(sig$source, sig$target)]

# End-to-end reachability: cytokine receptor -> signalling -> TF -> PRF1
g <- igraph::graph_from_data_frame(pkn %>% dplyr::select(source, target),
                                   directed = TRUE)

receptors <- c("IL12RB1", "IL12RB2", "IFNGR1", "IFNGR2", "IL2RB", "IL18R1")
receptors <- receptors[receptors %in% igraph::V(g)$name]

for (i in receptors) {
  d <- igraph::distances(g, v = i, to = "PRF1", mode = "out")
  cat(i, "-> PRF1:", ifelse(is.finite(d), paste(d, "hops"), "unreachable"), "\n")
}

# Canonical JAK-STAT path recovered
igraph::V(g)$name[igraph::shortest_paths(g, from = "IL12RB1",
                                         to = "PRF1", mode = "out")$vpath[[1]]]


# ---- Cell-type-specific PKN ----
# An edge cannot be active if the protein is not expressed in that cell type
expressed <- function(obj, ct, min_pct = 0.1) {
  cells <- WhichCells(obj, idents = ct)
  cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")[, cells]
  rownames(cnt)[Matrix::rowMeans(cnt > 0) >= min_pct]
}

pkn_list <- lapply(celltypes, function(ct) {
  genes <- expressed(seurat_filt, ct)
  pkn %>% 
    dplyr::filter(source %in% genes, target %in% genes) %>%
    dplyr::select(source, interaction, target)  
})

names(pkn_list) <- celltypes

sapply(pkn_list, nrow)
sapply(pkn_list, function(p) "PRF1" %in% c(p$source, p$target))

sapply(pkn_list, function(p) paste(colnames(p), collapse = ", "))

# PRF1 is a sink in the NK network 
pkn_list[["NK cells"]] %>% 
  dplyr::filter(target == "PRF1")   

pkn_list[["NK cells"]] %>% 
  dplyr::filter(source == "PRF1")   


# ---- EOMES -> PRF1 is NK-specific ----
# EOMES and ELF4 are not expressed above threshold in healthy CD8 T cells,
# so these edges get pruned out. 
cd8_edges <- paste(pkn_list[["CD8+ T cells"]]$source,
                   pkn_list[["CD8+ T cells"]]$target)

nk_only <- pkn_list[["NK cells"]] %>%
  dplyr::filter(!paste(source, target) %in% cd8_edges)

nk_only %>% 
  dplyr::filter(target == "PRF1")   # ELF4 -> PRF1, EOMES -> PRF1


# ---- Fig 3: TF activity of effector-gene regulators ----
tf_map <- net %>%
  dplyr::filter(target %in% c("PRF1", "RAB27A")) %>%
  dplyr::distinct(source, target)

tf_annot <- tf_map %>%
  group_by(source) %>%
  summarise(regulates = paste(sort(unique(target)), collapse = " + ")) %>%
  column_to_rownames("source")

hlh_mat <- df %>%
  filter(source %in% rownames(tf_annot)) %>%
  pivot_wider(id_cols = celltype, names_from = source, values_from = mean) %>%
  column_to_rownames("celltype") %>%
  as.matrix()

tf_annot <- tf_annot[colnames(hlh_mat), , drop = FALSE]

annot_colors <- list(regulates = c("PRF1" = "#880808",
                                   "RAB27A" = "#0425E0",
                                   "PRF1 + RAB27A" = "#31ED0C"))

pheatmap(hlh_mat,
         color = colors.use,
         border_color = "white",
         cellwidth = 15, cellheight = 15,
         annotation_col = tf_annot,
         annotation_colors = annot_colors,
         main = "TF activity: regulators of PRF1 and RAB27A",
         filename = "fig3_hlh_regulator_heatmap.pdf",
         width = 10, height = 5)