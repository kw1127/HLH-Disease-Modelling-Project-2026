# ============================================================
# 11. Differential expression
#
# Two contrasts, each comparing a cytotoxic effector to a naive baseline:
#   NK vs pooled naive (CD4 + CD8 naive T)
#   CD8 TEMRA vs CD8 naive T
#
# The naive populations are the resting counterparts, so the contrast isolates
# the effector programme rather than lineage differences.
#
# NK has no resting counterpart of its own, so it is compared to the pooled
# naive mean, while CD8 TEMRA is compared to its own resting lineage. The NK
# contrast therefore spans a greater transcriptional distance and carries a
# lineage component on top of the effector programme. A pooled-naive version of
# the TEMRA contrast is also computed as a concordance check.
# ============================================================

DefaultAssay(pbmc.clean) <- "RNA"

hlh <- c(
  "PRF1",   
  "UNC13D",  
  "STX11",  
  "STXBP2", 
  "RAB27A", 
  "LYST",    
  "SH2D1A",  
  "XIAP"
)

hlh <- factor(hlh, levels = hlh)
hlh_chr <- as.character(hlh)

# Colour: per gene, the population means are
# z-scored across the 12 populations. It encodes each gene's
# relative rank across populations. 
# Size is the percent of cells with non-zero counts, which is sensitive to capture depth.
dot_hlh <- DotPlot(pbmc.clean, features = levels(hlh)) +
  RotatedAxis() +
  labs(title = NULL,
       x = "HLH-associated gene",
       y = "Cell type") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 9)) + 
  geom_vline(xintercept = c(4.5, 6.5), linetype = "dashed", colour = "grey70")

ggsave("22_hlh_dotplot.png", dot_hlh, width = 8, height = 5)


# The four cell types used for the contrasts
effectors <- c("NK", "CD8 TEMRA")
naive <- c("CD8 naive T", "CD4 naive T")
used <- c(effectors, naive)

# Subset the pseudobulk matrix and metadata to those cell types, keeping the two
# objects aligned by using the metadata rownames to index the matrix columns.
keep_cols <- meta_pb$celltype %in% used
meta_sub <- droplevels(meta_pb[keep_cols, , drop = FALSE])
pseudo_sub <- pseudo[, rownames(meta_sub), drop = FALSE]

# Set CD8 naive T as the reference level, so model coefficients are expressed as
# "cell type X versus CD8 naive T" by default.
meta_sub$celltype <- relevel(factor(meta_sub$celltype), ref = "CD8 naive T")

# Keep only donors contributing at least two cell types.
# The design includes pb_sample as a term, which means each donor's baseline is
# estimated from that donor's own samples. A donor with only one column gives no
# within-donor comparison and makes the model unidentifiable.
cols_per_sample <- table(meta_sub$pb_sample)
usable_samples <- names(cols_per_sample)[cols_per_sample >= 2]
dropped <- setdiff(levels(factor(meta_sub$pb_sample)), usable_samples)

meta_sub <- droplevels(meta_sub[meta_sub$pb_sample %in% usable_samples, , drop = FALSE])
pseudo_sub <- pseudo_sub[, rownames(meta_sub), drop = FALSE]

# Realised design: all 21 pb_samples contribute all four cell types, so no
# samples are dropped and every cell-type comparison is fully within-donor.
cat(sprintf("dropped %d pb_samples; %d pb_samples retained\n",
            length(dropped), length(unique(meta_sub$pb_sample))))
print(table(meta_sub$celltype))
print(table(rowSums(table(meta_sub$pb_sample, meta_sub$celltype) > 0)))

# design = ~ pb_sample + celltype
#   pb_sample: a paired/blocking term. Each donor gets their own intercept, so
#   donor-to-donor differences are absorbed and do not inflate the
#   variance used to test cell type.
#   celltype: the term of interest, tested after that adjustment.
#
# Batch is absent from the formula because each pb_sample belongs to exactly one
# batch. Batch is nested within pb_sample and is therefore already absorbed by
# the donor term, and is not separately estimable.
dds_ct <- DESeqDataSetFromMatrix(
  countData = pseudo_sub,
  colData = meta_sub,
  design = ~ pb_sample + celltype
)

# Fail if the design matrix is rank-deficient.
mm <- model.matrix(~ pb_sample + celltype, meta_sub)
stopifnot(qr(mm)$rank == ncol(mm))

# Gene filter: keep genes with >=5 counts in at least as many samples as the
# smallest cell type group. This ensures a gene can  be detected
# throughout one whole group, so a real group-specific gene is not removed,
# while genes seen in only a couple of samples are.
min_grp <- min(table(meta_sub$celltype))
keep <- rowSums(counts(dds_ct) >= 5) >= min_grp
cat(sprintf("gene filter: %d of %d retained (min_grp = %d)\n",
            sum(keep), nrow(dds_ct), min_grp))
dds_ct <- dds_ct[keep, ]

# Size factors correct for differences in total sequencing depth between samples.
# type = "poscounts" computes them using only non-zero counts per gene, which is
# needed here because pseudobulk matrices still contain many zeros and the
# default method requires a gene to be non-zero in every sample.
# DESeq() does not overwrite existing size factors, so these are the ones used.
dds_ct <- estimateSizeFactors(dds_ct, type = "poscounts")

# Variance stabilising transformation.
# Raw counts have variance that grows with the mean, so high-expressed genes
# dominate distance calculations. VST puts everything on a roughly log scale with
# constant variance, which is what PCA and clustering need.
# blind = TRUE ignores the design, so the transform cannot be accused of having
# been shaped by the comparison it is used to inspect.
vsd_ct <- vst(dds_ct, blind = TRUE)
mat_ct <- assay(vsd_ct)

# PCA of the pseudobulk samples: checks that cell types separate and that batch
# is not the dominant source of variation. plotPCA uses the 500 most variable
# genes by default.
pca_df <- plotPCA(vsd_ct, intgroup = c("celltype", "batch"), returnData = TRUE)
pv <- round(100 * attr(pca_df, "percentVar"))
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = celltype, shape = batch)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(x = sprintf("PC1: %d%%", pv[1]), y = sprintf("PC2: %d%%", pv[2]),
       title = NULL) +
  theme_bw()
ggsave("23_pca_lymphoid.png", p_pca, width = 8, height = 5, dpi = 300)

# Fit the model. DESeq() runs three steps: size factor estimation, dispersion
# estimation (how variable each gene is beyond what counting noise explains),
# and the negative binomial fit with Wald tests.
dds_ct <- DESeq(dds_ct, quiet = TRUE)

# Dispersion plot: Black points are per-gene estimates, the red
# line is the fitted trend, blue points are the final shrunken values. Genes
# should sit near the trend and shrink toward it; a cloud with no trend means
# the model has not fitted well.
dev.new()
dispersion_plot <- ggplotify::as.ggplot(~plotDispEsts(dds_ct, main = NULL))
dev.off()

# Sample-to-sample distances on the VST values. Samples of the same cell type
# should cluster together. Useful for spotting a mislabelled or outlier sample.
sampleDists <- dist(t(mat_ct))
ann <- as.data.frame(colData(dds_ct)[, c("celltype", "batch")])
sample_heatmap <- pheatmap(as.matrix(sampleDists),
                           clustering_distance_rows = sampleDists,
                           clustering_distance_cols = sampleDists,
                           clustering_method = "ward.D2",
                           annotation_col = ann,
                           show_rownames = FALSE,
                           show_colnames = FALSE,
                           color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                           main = "",
                           silent = TRUE)

s_s_heatmap <- as.ggplot(sample_heatmap$gtable)


# Contrast 1: CD8 TEMRA vs CD8 naive T 
# A comparison of two levels of the same factor. Written as
# c(factor, numerator, denominator), so a positive fold change means higher in
# CD8 TEMRA.
# alpha = 0.05 sets the FDR level used for independent filtering and the summary.
res_temra <- results(dds_ct,
                     contrast = c("celltype", "CD8 TEMRA", "CD8 naive T"),
                     alpha = 0.05)

# Extract the Wald statistics before lfc shrinkage
# stat = log2FoldChange / standard error. It combines effect size and precision in one value.
# A large fold change measured noisily gets a small statistic whereas the fold change alone would not distinguish the two.
# This is the right input for TF activity inference, which needs a ranking that
# already accounts for uncertainty.
# lfcShrink overwrites res_temra and drops the stat column, so the order matters.
stat_temra <- setNames(res_temra$stat, rownames(res_temra))
stat_temra <- stat_temra[!is.na(stat_temra)]

res_temra <- lfcShrink(dds_ct,
                       contrast = c("celltype", "CD8 TEMRA", "CD8 naive T"),
                       type = "ashr",
                       res = res_temra)

stopifnot(length(stat_temra) > 0,
          !"stat" %in% names(as.data.frame(res_temra)))

# Contrast 2: NK vs the average of two naive types
# With CD8 naive T as the reference, the group means are:
#   CD8 naive T = b0
#   CD4 naive T = b0 + b_CD4
#   NK = b0 + b_NK
# The pooled naive mean is the unweighted average of the two naive groups:
#   pooled = b0 + b_CD4 / 2
# So the difference of interest is:
#   NK - pooled = b_NK - 0.5 * b_CD4
# which is the weight vector: +1 on b_NK, -0.5 on b_CD4, 0 elsewhere.
rn <- resultsNames(dds_ct)

# DESeq2 mangles level names into coefficient names (spaces become dots), so the
# names are looked up rather than typed by hand. stopifnot fails loudly if the
# pattern matches zero or more than one coefficient.
coef_of <- function(lvl) {
  hit <- grep(paste0("^celltype_", make.names(lvl), "_vs_"), rn, value = TRUE)
  stopifnot(length(hit) == 1L)
  hit
}

con <- setNames(numeric(length(rn)), rn)
con[coef_of("NK")] <-  1
con[coef_of("CD4 naive T")] <- -0.5
# CD8 naive T is the reference, its coefficient is 0 by construction

# unname() because results() expects a plain numeric vector in coefficient order.
res_nk <- results(
  dds_ct, 
  contrast = unname(con), 
  alpha = 0.05)

stat_nk <- setNames(res_nk$stat, rownames(res_nk))
stat_nk <- stat_nk[!is.na(stat_nk)]

res_nk <- lfcShrink(dds_ct,
                    contrast = unname(con),
                    type = "ashr",
                    res = res_nk)

stopifnot(length(stat_nk) > 0,
          !"stat" %in% names(as.data.frame(res_nk)))

# Diagnostics on the HLH panel.
# padj is set to NA by independent filtering and Cook's outlier flagging; the
# volcanoes and hlh_de in section 13 drop those rows silently.
print(sapply(list(nk = res_nk, temra = res_temra), \(r)
             c(NA_padj = sum(is.na(r$padj)), tested = sum(!is.na(r$padj)))))

# XIAP is down in NK, opposite to the other seven, so the panel is not
# coordinately upregulated.
print(as.data.frame(res_nk[hlh_chr, c("baseMean", "log2FoldChange", "padj")]))
print(as.data.frame(res_temra[hlh_chr, c("baseMean", "log2FoldChange", "padj")]))

# STX11 is low-abundance. Its small padj comes from near-total absence in the
# CD8 naive baseline rather than from broad detection in NK.
print(tapply(counts(dds_ct, normalized = TRUE)["STX11", ], meta_sub$celltype,
             \(x) round(mean(x), 2)))
print(tapply(counts(dds_ct)["STX11", ], meta_sub$celltype, \(x) sum(x >= 5)))

# Volcano plots for the two contrasts
# Well-known genes included purely for orientation, so a reader can confirm the
# plot is the right way round: cytotoxic genes on the effector side, naive
# markers (CCR7, SELL, TCF7) on the other.
anchors <- c("GNLY", "NKG7", "GZMB", "CCL5", "CD3D", "CCR7", "SELL", "TCF7")
hlh_chr <- c("PRF1", "UNC13D", "STX11", "STXBP2",
             "RAB27A", "LYST", "SH2D1A", "XIAP")

cap <- 150   # y-axis ceiling; capped points drawn as triangles
lfc_lim <- 8    # shared x range, so effect sizes are comparable across panels

volcano_panel <- function(res, title) {
  
  df <- as.data.frame(res) |>
    rownames_to_column("gene") |>
    filter(!is.na(padj)) |>          # DESeq2 sets padj to NA for filtered genes
    mutate(
      y = pmin(-log10(padj), cap),
      capped = -log10(padj) > cap,
      sig = case_when(
        padj < 0.05 & log2FoldChange > 1 ~ "Up in effector",
        padj < 0.05 & log2FoldChange < -1 ~ "Up in naive",
        TRUE ~ "n.s."))
  
  lab <- df |>
    dplyr::filter(gene %in% c(anchors, hlh_chr)) |>
    dplyr::mutate(is_hlh = gene %in% hlh_chr,
                  col = ifelse(is_hlh, "#1B7837", "grey25"),
                  face = ifelse(is_hlh, "bold", "italic"),
                  sz = ifelse(is_hlh, 3.5, 3))
  
  ggplot(df, aes(log2FoldChange, y)) +
    geom_point(aes(colour = sig, shape = capped), size = 0.7, alpha = 0.4) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    scale_colour_manual(values = c("Up in effector" = "#C0504D",
                                   "Up in naive" = "#4F81BD",
                                   "n.s." = "grey80"), name = NULL) +
    geom_vline(xintercept = c(-1, 1), lty = 2, colour = "grey60") +
    geom_hline(yintercept = -log10(0.05), lty = 2, colour = "grey60") +
    geom_point(data = filter(df, gene %in% hlh_chr),
               colour = "#1B7837", size = 2.5, shape = 18) +
    geom_text_repel(data = lab, aes(label = gene),
                    colour = lab$col, fontface = lab$face, size = lab$sz,
                    max.overlaps = Inf, seed = 42,
                    ylim = c(NA, cap * 0.9), nudge_y = -8, box.padding = 0.5,
                    segment.colour = lab$col, segment.size = 0.3) +
    coord_cartesian(xlim = c(-lfc_lim, lfc_lim)) +
    scale_x_continuous(breaks = seq(-lfc_lim, lfc_lim, 4)) +
    labs(x = "log2 fold change", y = "-log10 (padj)", title = title) +
    theme_bw() +
    theme(plot.title = element_text(size = 12))
}

# assemble the plot
x_obs <- range(c(res_temra$log2FoldChange, res_nk$log2FoldChange), na.rm = TRUE)
if (x_obs[1] < -lfc_lim || x_obs[2] > lfc_lim)
  warning(sprintf("observed LFC range [%.1f, %.1f] exceeds the plotted [-%g, %g]",
                  x_obs[1], x_obs[2], lfc_lim, lfc_lim))

fig_volcanoes <-
  (volcano_panel(res_temra, "CD8 TEMRA vs CD8 naive T") /
     volcano_panel(res_nk,    "NK vs pooled naive")) +
  plot_layout(guides = "collect", axis_titles = "collect_x") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

ggsave("Figure_06_volcanoes.png", fig_volcanoes,
       width = 6.5, height = 8.5, dpi = 300, bg = "white")

saveRDS(stat_temra, "stat_temra.rds")
saveRDS(stat_nk, "stat_nk.rds")

# HLH expression landscape and pseudobulk PCA
dot_hlh_panel <- dot_hlh +
  labs(tag = "A", x = NULL) +
  theme(legend.position = "right", legend.box = "horizontal") 

pca_panel <- p_pca +
  labs(tag = "B") +
  theme(legend.position = "bottom", legend.box = "horizontal") 

fig05_design <- "
AAAA
BBBB
"

fig_hlh_pca <- dot_hlh_panel + pca_panel +
  plot_layout(design = fig05_design, heights = c(1, 0.8))

ggsave("Figure_05_hlh_pca.png", fig_hlh_pca,
       width = 10, height = 10, dpi = 300, bg = "white")

# Supplementary: model diagnostics, stacked so each gets the full width.
fig_supp <- dispersion_plot / s_s_heatmap +
  plot_layout(heights = c(1, 1.4)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave("Supp_Figure_05_deseq_qc.png", fig_supp,
       width = 8, height = 12, dpi = 300, bg = "white")


# Additional sanity checks for the thesis
#
# 1. Is pooling CD4 and CD8 naive T a good idea?
# Counted at the same thresholds used for the effector contrasts, so the three
# numbers are directly comparable. The argument is not that the CD4 vs CD8 naive
# count is small in absolute terms but that it is an order of magnitude below
# either effector contrast, which together with the PC1/PC2 split (70% vs 7%)
# supports treating the two naive populations as one resting reference.
res_naive <- results(dds_ct,
                     contrast = c("celltype", "CD4 naive T", "CD8 naive T"),
                     alpha = 0.05)

n_de <- function(r) sum(r$padj < 0.05 & abs(r$log2FoldChange) > 1, na.rm = TRUE)

cat(sprintf(
  "DE genes (padj<0.05, |LFC|>1) of %d tested:\n  CD4 vs CD8 naive: %d (%.1f%%)\n  NK vs pooled naive: %d (%.1f%%)\n  CD8 TEMRA vs CD8 naive: %d (%.1f%%)\n",
  nrow(dds_ct),
  n_de(res_naive), 100 * n_de(res_naive) / nrow(dds_ct),
  n_de(res_nk), 100 * n_de(res_nk) / nrow(dds_ct),
  n_de(res_temra), 100 * n_de(res_temra) / nrow(dds_ct)))

# 2. Does the TEMRA result depend on which naive reference is used?
# The main TEMRA contrast uses CD8 naive T alone while NK uses the pooled mean,
# so the two are not on identical footing. Refitting TEMRA against the same
# pooled reference lets the concordance be reported rather than assumed.
con_t <- setNames(numeric(length(rn)), rn)
con_t[coef_of("CD8 TEMRA")] <-  1
con_t[coef_of("CD4 naive T")] <- -0.5

res_temra_pooled <- results(dds_ct, contrast = unname(con_t), alpha = 0.05)
stat_temra_pooled <- setNames(res_temra_pooled$stat, rownames(res_temra_pooled))
stat_temra_pooled <- stat_temra_pooled[!is.na(stat_temra_pooled)]

g_shared <- intersect(names(stat_temra), names(stat_temra_pooled))
cat(sprintf("TEMRA Wald stat correlation, CD8-only vs pooled reference: %.3f\n",
            cor(stat_temra[g_shared], stat_temra_pooled[g_shared])))

# 0,985 correlation between CD8-only and pooled naive reference

saveRDS(stat_temra_pooled, "stat_temra_pooled.rds")