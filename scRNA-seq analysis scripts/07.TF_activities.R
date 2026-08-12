# ============================================================
# 12. Transcription factor activity (decoupleR + CollecTRI)
#
# A TF's own mRNA level is a poor measure of its activity, because activity is
# controlled by phosphorylation, localisation, and cofactors. Instead, activity is
# inferred from the behaviour of the genes it regulates.
# ============================================================

# CollecTRI is a curated set of TF-target relationships. Each edge carries a "mor"
# (mode of regulation): +1 if the TF activates that target, -1 if it represses.
# split_complexes = FALSE keeps complexes such as NF-kB as one entry rather than
# splitting them into subunits.
collectri <- get_collectri(organism = "human", split_complexes = FALSE)

# Part 1
#
# TF Activity of top 40 TFs per cell type by absolute activity score
# Collapse the pseudobulk samples to one mean profile per cell type.
# tapply() groups the values of each gene by cell type and averages within group.
ct <- droplevels(factor(meta_sub[colnames(mat_ct), "celltype"]))
mat_ct_mean <- t(apply(mat_ct, 1, function(x) tapply(x, ct, mean)))

# Centre and scale each gene across the four cell types (z-score per row).
#
# Why this is necessary: the linear model below regresses expression on regulon
# membership. On raw VST values the fit is driven by absolute expression level,
# which is shared by all four cell types, so every TF would score high. After
# centering, each value says "how far above or below this gene's own average". 
#
# What this means for interpretation: scores are relative to the mean of the four
# cell types, not to a global average. A TF uniformly high in all
# lymphocytes but absent in monocytes now looks flat, not high.
mat_z <- t(scale(t(mat_ct_mean)))
mat_z <- mat_z[stats::complete.cases(mat_z), ]

# Univariate linear model (ULM).
# For each TF and each cell type, fit: expression ~ regulon membership, where
# membership is 0 for genes the TF does not regulate, +1 for activated targets
# and -1 for repressed ones. The t-value of the fitted slope is the activity score.
#
# A high positive score means activated targets are up and repressed targets are
# down, which is what an active TF produces.
#
# "Univariate" means each TF is tested on its own. That is fast and robust, but
# it does not account for TFs that share target genes, so overlapping regulons
# can give correlated scores.
#
# minsize = 5 excludes TFs with less than 5 downstream targets.
tf_acts <- run_ulm(
  mat = mat_z, 
  network = collectri,
  .source = "source", 
  .target = "target", 
  .mor = "mor",
  minsize = 5) %>%
  dplyr::filter(statistic == "ulm") %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH"))

# Reshape to a TF x cell type matrix for plotting as a heatmap
tf_mat <- tf_acts %>%
  pivot_wider(id_cols = source, names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()

# Ranked by largest absolute score, not by significance and not by variance.
#
# Not by FDR: with only four conditions, the z-scoring above compresses the
# t-values, so only a few tests clear BH correction. Those that do are the TFs with the
# largest regulons, because a bigger regulon gives a more precise
# slope. This isn't very meaningful
#
# Not by standard deviation: with n = 4, and two of those (CD4 and CD8 naive)
# being near-identical, an SD over four numbers is unstable
tf_ranked <- tf_acts %>%
  group_by(source) %>%
  summarise(max_abs = max(abs(score)), # strongest activity in any cell type
            best_padj = min(p_adj), .groups = "drop") %>%
  arrange(desc(max_abs))

top_tfs <- head(tf_ranked$source, 40)

# Transpose so cell types are rows and TFs are columns
top_acts_mat <- t(tf_mat[top_tfs, ])

# Fix the row order naive -> effector
row_order <- c("CD4 naive T", "CD8 naive T", "CD8 TEMRA", "NK")
top_acts_mat <- top_acts_mat[intersect(row_order, rownames(top_acts_mat)), ]

# Module annotation
# Assigned from known biology
modules <- list(
  "Effector / IFN" = c("TBX21","STAT1","IRF1","IRF3","IRF5",
                       "NFKB","RUNX1","RUNX2","PML"),
  "MHC class II" = c("CIITA","RFX5","RFXANK","RFXAP"),
  "Growth / metabolic" = c("MYC","MZF1","SP1","NFYA","NFYB",
                           "SREBF1","SREBF2","NR1H3","ATF6"),
  "Stress / AP-1" = c("JUN","AP1","ATF2","DDIT3","NFE2L2",
                      "TP53","BACH1","LITAF"))

# Everything starts as "Other" and is reassigned if it appears in a module.
Function <- setNames(rep("Other", ncol(top_acts_mat)), colnames(top_acts_mat))
for (m in names(modules))
  Function[names(Function) %in% modules[[m]]] <- m

# pheatmap needs a data frame whose rownames match the column names of the matrix.
ann_col <- data.frame(Function = factor(
  Function, levels = c(names(modules), "Other")),
  row.names = names(Function))

ann_colors <- list(Function = c(
  "Effector / IFN" = "#C0504D",
  "MHC class II" = "#9970AB",
  "Growth / metabolic" = "#4F81BD",
  "Stress / AP-1" = "#E8A33D",
  "Other" = "grey85"))

# significance stars 
# Recovers the FDR information
star_mat <- tf_acts %>%
  dplyr::filter(source %in% colnames(top_acts_mat)) %>%
  dplyr::mutate(lab = dplyr::case_when(
    p_adj < 0.0001 ~ "****",
    p_adj < 0.001 ~ "***",
    p_adj < 0.01  ~ "**",
    p_adj < 0.05  ~ "*",
    TRUE ~ "")) %>%
  pivot_wider(id_cols = condition, names_from = source, values_from = lab) %>%
  column_to_rownames("condition") %>%
  as.matrix()

# Reorder to exactly match the heatmap or the stars would land on wrong cells.
star_mat <- star_mat[rownames(top_acts_mat), colnames(top_acts_mat)]

# Heatmap 
colors.use <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
lim <- quantile(abs(top_acts_mat), 0.95)
my_breaks <- c(seq(-lim, 0, length.out = 51),
               seq(lim / 50, lim, length.out = 50))

pheatmap(
  mat = t(top_acts_mat),
  color = colors.use,
  breaks = my_breaks,
  border_color = "white",
  cellwidth = 32,
  cellheight = 9,
  cluster_cols = FALSE,
  gaps_col = 2,
  cutree_rows = 3,
  treeheight_row = 30,
  angle_col = 45,
  annotation_row = ann_col,
  annotation_colors = ann_colors,
  annotation_names_row = FALSE,
  display_numbers = t(star_mat),
  number_color = "black",
  fontsize_number = 7,
  fontsize_row = 7,
  fontsize_col = 11,
  filename = "tf_acts_lymphoid.png",
  res = 600,
  width  = 4*32/72 + 30/72 + 0.7 + 2.3 + 0.2,
  height = 40*9/72 + 1.3)

# Part two
#
# TF activity per contrast
# The heatmap above describes cell types. This describes the two comparisons,
# which is what the network modelling needs: each column is already a difference
# so TFs shared by NK and TEMRA show up in both rather than
# cancelling against each other.
#
# The two contrasts were filtered independently, so their gene sets differ.
# Intersecting first prevents misalignment when they are bound into one matrix.
g <- intersect(names(stat_temra), names(stat_nk))
stat_mat <- cbind("CD8 TEMRA" = stat_temra[g], "NK" = stat_nk[g])

# No centering here. Unlike the expression matrix used for the heatmap, these values are already
# differences, so the baseline is built into the contrast.
# Centering across two columns would remove the signal being measured.
tf_contrast <- run_ulm(
  mat = stat_mat,
  network = collectri,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5) %>%
  dplyr::filter(statistic == "ulm") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup()

# Gene set enrichment analysis (GSEA)
net_homo <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") |>
  dplyr::select(source = gs_name, target = gene_symbol) |>
  dplyr::distinct() |>
  dplyr::mutate(mor = 1) |>
  dplyr::group_by(source) |>
  dplyr::filter(dplyr::n() <= 500) |>
  dplyr::ungroup()

# Run GSEA
gsea_contrast <- decoupleR::run_fgsea(
  mat = stat_mat,
  network = net_homo,
  .source = "source",
  .target = "target",
  minsize = 15, # filters the network
  maxSize = 500, # drops vague parent terms
  eps = 0, # multilevel p-values, no permutation floor
  times = 1000) |>
  dplyr::filter(statistic == "norm_fgsea") |>
  dplyr::group_by(condition) |>
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) |>
  dplyr::ungroup()

# Flag GO terms whose signal is driven predominantly by ribosomal protein genes 
# Several GO:BP terms (myoblast fusion, syncytium formation, muscle contraction)
# carry ~45 RPL genes as erroneous members. In this contrast their entire leading
# edge is ribosomal, so they report the translational program under an incorrect muscle label.
# Identify them from the leading edge which generates the NES.
# Correctly-annotated translation terms to retain as the representative finding
keep_translation <- c("GOBP_CYTOPLASMIC_TRANSLATION",
                      "GOBP_RIBOSOME_BIOGENESIS",
                      "GOBP_RRNA_PROCESSING")

pathways <- split(net_homo$target, net_homo$source)
ribo_re <- "^(RP[LS]|UBA52|FAU|RACK1)"

ribo_drop <- function(col) {
  f <- fgsea::fgsea(pathways, stat_mat[, col],
                    minSize = 15, maxSize = 500, eps = 0)
  f[, .(pathway,
        n_le = lengths(leadingEdge),
        ribo_le = sapply(leadingEdge, \(g) mean(grepl(ribo_re, g))))
  ][ribo_le > 0.5 & n_le >= 10, pathway] # n_le floor ignores 1-gene edges
}

drop <- unique(unlist(lapply(colnames(stat_mat), ribo_drop)))
drop <- setdiff(drop, keep_translation)

# Union across conditions so both panels are filtered by the same rule
drop2 <- c("GOBP_NEUROINFLAMMATORY_RESPONSE",
           "GOBP_MATURATION_OF_SSU_RRNA",
           "GOBP_CYTOPLASMIC_TRANSLATIONAL_INITIATION",
           "GOBP_RIBOSOMAL_SMALL_SUBUNIT_ASSEMBLY",
           "GOBP_RRNA_METABOLIC_PROCESS",
           "GOBP_PROTEIN_RNA_COMPLEX_ORGANIZATION")

gsea_plotting <- gsea_contrast |>
  filter(!source %in% c(drop, drop2)) |>
  mutate(pathway = gsub("_", " ", sub("^GOBP_", "", source)),
         pathway = stringr::str_trunc(pathway, 55))

# Diverging bar plot
top_n_each <- function(d, n = 10) {
  bind_rows(d |> group_by(condition) |> slice_max(score, n = n) |> ungroup(),
            d |> group_by(condition) |> slice_min(score, n = n) |> ungroup())
}

gsea_barplot <- gsea_plotting |>
  filter(p_adj < 0.05) |>
  top_n_each(8) |>
  mutate(pathway = reorder_within(pathway, score, condition)) |>
  ggplot(aes(score, pathway, fill = score > 0)) +
  geom_col() +
  scale_y_reordered() +
  facet_wrap(~ condition, scales = "free_y") +
  scale_fill_manual(values = c("steelblue", "firebrick"), guide = "none") +
  labs(x = "NES", y = NULL) +
  theme_bw()

ggsave("diverging_barplot.png", gsea_barplot, width = 10, height = 8, dpi = 300, bg = "white")

# Scatter plot of TF activities
tf_wide <- tf_contrast %>%
  mutate(cond = recode(condition, "CD8 TEMRA" = "temra", "NK" = "nk")) %>%
  select(source, cond, score, p_adj) %>%
  pivot_wider(names_from = cond, values_from = c(score, p_adj)) %>%
  mutate(
    sig_nk = p_adj_nk < 0.05,
    sig_temra = p_adj_temra < 0.05,
    class = case_when(
      sig_nk & sig_temra & sign(score_nk) == sign(score_temra) ~ "Shared",
      sig_nk & sig_temra ~ "Opposite",
      sig_nk ~ "NK only",
      sig_temra ~ "TEMRA only",
      TRUE ~ "n.s."
    ),
    delta = score_temra - score_nk
  ) %>%
  mutate(class = factor(class, levels = c("Shared", "NK only", "TEMRA only", "Opposite", "n.s.")))

lim_contrast <- max(abs(c(tf_wide$score_nk, tf_wide$score_temra)), na.rm = TRUE) * 1.05

# label the strongest shared TFs and the biggest outliers from the diagonal
label_contrast <- tf_wide %>%
  filter(class != "n.s.") %>%
  filter(rank(-(score_nk^2 + score_temra^2)) <= 15 | rank(-abs(delta)) <= 10)

tf_scatter <- ggplot(tf_wide, aes(score_nk, score_temra)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(aes(colour = class), size = 2, alpha = 0.85) +
  geom_text_repel(data = label_contrast, aes(label = source), size = 3,
                  max.overlaps = Inf, box.padding = 0.35, segment.colour = "grey60") +
  scale_colour_manual(values = c(
    "Shared" = "#8DA0CB", "NK only" = "#FC8D62",
    "TEMRA only" = "#66C2A5", "n.s." = "grey80"),
    drop = TRUE) +
  coord_equal(xlim = c(-lim_contrast, lim_contrast), ylim = c(-lim_contrast, lim_contrast)) +
  labs(x = "NK vs naive T",
       y = "CD8 TEMRA vs naive CD8 T",
       colour = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave("tf_scatter.png", tf_scatter, width = 10, height = 10, dpi = 300, bg = "white")

# Shared: Wald statistics for the HLH gene
hlh_obs_wide <- tibble(
  gene = hlh_chr,
  nk = as.numeric(stat_nk[hlh_chr]),
  temra = as.numeric(stat_temra[hlh_chr]),
  nk_padj = res_nk[hlh_chr, "padj"],
  temra_padj = res_temra[hlh_chr, "padj"]) %>%
  mutate(delta = temra - nk)

hlh_obs_long <- hlh_obs_wide %>%
  select(gene, NK = nk, `CD8 TEMRA` = temra) %>%
  pivot_longer(-gene, names_to = "condition", values_to = "obs")

p_hlh_obs <- hlh_obs_wide %>%
  filter(!is.na(nk), !is.na(temra)) %>%
  mutate(gene = reorder(gene, (nk + temra) / 2)) %>%
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
  geom_segment(aes(x = nk, xend = temra, y = gene, yend = gene),
               colour = "grey65", linewidth = 0.7) +
  geom_point(aes(nk, gene, colour = "NK", shape = nk_padj < 0.05),
             size = 3, alpha = 0.75, stroke = 0.9) +
  geom_point(aes(temra, gene, colour = "CD8 TEMRA", shape = temra_padj < 0.05),
             size = 3, alpha = 0.75, stroke = 0.9) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = "padj < 0.05", `FALSE` = "n.s."),
                     name = NULL) +
  scale_colour_manual(values = c("NK" = "#2166AC", "CD8 TEMRA" = "#B2182B"),
                      breaks = c("NK", "CD8 TEMRA"), name = NULL) +
  labs(x = "Wald statistic", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey92"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave("hlh_observed.png", p_hlh_obs, width = 6.5, height = 4, dpi = 300, bg = "white")