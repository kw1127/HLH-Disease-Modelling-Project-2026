# Prior knowledge network: signalling + transcriptional layers
Idents(seurat_filt) <- "celltype_final"

# Import prior knowledge network
ppi <- omnipath_interactions()

nrow(ppi)   # 85,217 interactions

# Signalling layer: directed and unambiguously signed
# curation_effort >= 2 chosen: >=3 loses the IL12RB1 -> PRF1 path,
# >=1 gives ~70k edges which is intractable for the ILP solver
sig <- ppi %>%
  dplyr::filter(consensus_direction == 1,
                consensus_stimulation + consensus_inhibition == 1,
                curation_effort >= 2) %>%
  dplyr::mutate(interaction = ifelse(consensus_stimulation == 1, 1, -1)) %>%
  dplyr::select(source = source_genesymbol,
                target = target_genesymbol,
                interaction) %>%
  dplyr::distinct()

nrow(sig)


# HLH genes cannot serve as CARNIVAL input nodes

# Which are in the signalling layer at all?
hlh[hlh %in% c(sig$source, sig$target)]  

# Which are absent entirely?
setdiff(hlh, c(sig$source, sig$target))     

# Which have outgoing edges? A source node needs these to propagate a perturbation
sig %>% 
  dplyr::filter(source %in% hlh) %>% 
  dplyr::count(source)                      

# Which are sinks? Sinks have regulators but regulate nothing
intersect(hlh, setdiff(sig$target, sig$source)) 

# Of the two with outgoing edges, do they reach any of the measured TFs?
# A source that reaches no measurements contributes nothing to the solution
g_sig <- igraph::graph_from_data_frame(sig %>% 
                                         dplyr::select(source, target),
                                       directed = TRUE)

tfs_measured <- df %>%
  dplyr::filter(celltype == "NK cells") %>%
  dplyr::pull(source) %>%
  intersect(igraph::V(g_sig)$name)

for (s in intersect(hlh, sig$source)) {
  d <- igraph::distances(g_sig, v = s, to = tfs_measured, mode = "out")
  cat(s, "reaches", sum(is.finite(d)), "of", length(tfs_measured), "TFs\n")
}

# Summary:
#   PRF1, UNC13D, LYST, RAB27A  absent from the signalling layer
#   STX11, SH2D1A are present but sinks, no outgoing edges
#   STXBP2 has 1 outgoing edge
#   XIAP has 5 outgoing edges, but monocyte and DC-enriched, not NK/T
#
# The HLH genes are effectors of cytotoxic granule exocytosis. They carry out the
# final step of the pathway and do not signal onward, so they have no outgoing edges.

# ---- Is this a filtering artefact, or structural? ----

# PRF1 is absent from the signalling layer before any filter is applied
"PRF1" %in% ppi$source_genesymbol   # FALSE
"PRF1" %in% ppi$target_genesymbol   # FALSE

# Across ALL OmniPath layers, PRF1 has only 4 outgoing edges, and none are usable
# for CARNIVAL, which needs edges that are both directed and signed
all_int <- OmnipathR::import_all_interactions()

all_int %>%
  dplyr::filter(source_genesymbol == "PRF1") %>%
  dplyr::select(target_genesymbol, is_directed, is_stimulation,
                is_inhibition, curation_effort, n_references)

#   PRF1 -> CTSB signed, but curation_effort 0 and 0 references
#   PRF1 -> CDK2 signed, but curation_effort 0 and 0 references
#   PRF1 -> MTOR unsigned
#   PRF1 -> PRKCA unsigned

# Same check across all 8 HLH genes
all_int %>%
  dplyr::filter(source_genesymbol %in% hlh,
                is_directed == 1,
                is_stimulation + is_inhibition == 1) %>%
  dplyr::count(source_genesymbol, name = "usable_outgoing_edges")


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