# Prior knowledge network: signalling + transcriptional layers, coupled
Idents(seurat_filt) <- "celltype_final"

# Import prior knowledge network
ppi <- omnipath_interactions()

# Signalling layer: directed and signed
# curation_effort >= 2 chosen: >=3, >= 2 removes important HLH paths, >=1 gives 70k edges 
sig <- ppi %>%
  dplyr::filter(consensus_direction == 1,
                consensus_stimulation + consensus_inhibition == 1,   
                curation_effort >= 2) %>%
  dplyr::mutate(interaction = ifelse(consensus_stimulation == 1, 1, -1)) %>%
  dplyr::select(source = source_genesymbol,
                target = target_genesymbol,
                interaction) %>%
  dplyr::distinct()

# --- HLH genes cannot be used as input nodes ---
hlh[hlh %in% c(sig$source, sig$target)] # 4 of 8 present

sig %>% 
  dplyr::filter(source %in% hlh) %>% 
  dplyr::count(source) # 2 have outgoing edges

#   PRF1, UNC13D, RAB27A, LYST -- absent from signalling layer (terminal effectors)
#   STX11, SH2D1A  -- sinks, no outgoing edges
#   STXBP2 -- 1 outgoing edge, reaches 0 TFs
#   LYST, XIAP -- only viable sources, and both monocyte-enriched

# Couple signalling with transcriptional layer (TF -> HLH gene),
# so HLH genes enter the network as downstream endpoints, not inputs
trn <- net %>%
  dplyr::filter(target %in% hlh) %>%
  dplyr::select(source, target, interaction = mor)

pkn <- dplyr::bind_rows(sig, trn) %>% dplyr::distinct()

dim(pkn)                                  # 14,989 edges
hlh[hlh %in% c(pkn$source, pkn$target)]   # 7 of 8; UNC13D has no regulators anywhere

# Confirm the layers join: PRF1's TFs must be in the signalling layer
prf1_tfs <- net %>% dplyr::filter(target == "PRF1") %>% dplyr::pull(source)
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

# Canonical JAK-STAT path recovered: IL12RB1 -> JAK2 -> STAT3 -> PRF1
igraph::V(g)$name[igraph::shortest_paths(g, from = "IL12RB1", 
                                         to = "PRF1", mode = "out")$vpath[[1]]]

# --- Cell-type-specific PKNs ---
expressed <- function(obj, ct, min_pct = 0.1) {
  cells <- WhichCells(obj, idents = ct)
  cnt   <- GetAssayData(obj, assay = "RNA", layer = "counts")[, cells]
  rownames(cnt)[Matrix::rowMeans(cnt > 0) >= min_pct]
}

pkn_list <- lapply(celltypes, function(ct) {
  g <- expressed(seurat_filt, ct)
  pkn %>% dplyr::filter(source %in% g, target %in% g)
})
names(pkn_list) <- celltypes

sapply(pkn_list, nrow)                                          
sapply(pkn_list, function(p) "PRF1" %in% c(p$source, p$target)) 

# --- FINDING: EOMES -> PRF1 is NK-specific ---
# EOMES and ELF4 are not expressed above threshold in healthy CD8 T cells,
# so these edges are pruned out. Explains why naive CD8 T cells lack perforin
# while NK cells express it.
cd8_edges <- paste(pkn_list[["CD8+ T cells"]]$source, 
                   pkn_list[["CD8+ T cells"]]$target)

nk_only <- pkn_list[["NK cells"]] %>%
  dplyr::filter(!paste(source, target) %in% cd8_edges)

nk_only %>% dplyr::filter(target == "PRF1")   # ELF4 -> PRF1, EOMES -> PRF1

# --- Fig 3: TF activity of effector-gene regulators ---
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
