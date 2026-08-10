# =============================================================================
# 16. Specificity check: do the HLH-gene regulators show activity outside the
#     cytotoxic populations?
#
# The anchored design forces every CollecTRI regulator of an HLH gene into the
# measurement set. If those regulators are equally active in populations that do
# not carry the cytotoxic machinery, the anchoring is recovering generic
# signalling rather than anything specific to the effector programme.
#
# Two potential controls, chosen from the HLH gene expression dot plot:
#   B cells - near-zero across the whole panel, a clean negative control
#   monocytes - high for the trafficking genes (STX11, STXBP2, LYST), but
#                negative for PRF1, so a partial negative only.
#
# This asks whether the regulators carry activity signal in these populations.
# =============================================================================

# Same naive baseline as the effector contrasts, so the comparison is like for
# like: each control population is tested against the same resting reference.
naive <- c("CD8 naive T", "CD4 naive T")
controls <- c("Naive B", "classical monocytes")

# Check representation before trusting a paired design. The effector contrasts
# had all 21 pb_samples in every cell type; anything less weakens the pairing.
meta_pb %>%
  dplyr::filter(celltype %in% c(controls, naive)) %>%
  dplyr::count(celltype)


# Differential expression, control population vs pooled naive
run_contrast <- function(ct) {
  
  used <- c(ct, naive)
  keep_col <- meta_pb$celltype %in% used
  meta_c <- droplevels(meta_pb[keep_col, , drop = FALSE])
  pseudo_c <- pseudo[, rownames(meta_c), drop = FALSE]
  
  # donors contributing only one cell type give no within-donor comparison
  usable <- names(which(table(meta_c$pb_sample) >= 2))
  meta_c <- droplevels(meta_c[meta_c$pb_sample %in% usable, , drop = FALSE])
  pseudo_c <- pseudo_c[, rownames(meta_c), drop = FALSE]
  
  meta_c$celltype <- relevel(factor(meta_c$celltype), ref = "CD8 naive T")
  
  dds <- DESeqDataSetFromMatrix(pseudo_c, meta_c, ~ pb_sample + celltype)
  
  min_grp <- min(table(meta_c$celltype))
  dds <- dds[rowSums(counts(dds) >= 5) >= min_grp, ]
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, quiet = TRUE)
  
  # control vs the unweighted mean of the two naive populations, matching the
  # NK contrast: +1 on the control coefficient, -0.5 on CD4 naive
  rn <- resultsNames(dds)
  con <- setNames(numeric(length(rn)), rn)
  coef_of <- function(lvl) {
    hit <- grep(paste0("^celltype_", make.names(lvl), "_vs_"), rn, value = TRUE)
    stopifnot(length(hit) == 1L); hit
  }
  con[coef_of(ct)]            <-  1
  con[coef_of("CD4 naive T")] <- -0.5
  
  res  <- results(dds, contrast = unname(con), alpha = 0.05)
  stat <- setNames(res$stat, rownames(res))
  stat[!is.na(stat)]
}

stat_ctrl <- lapply(setNames(controls, controls), run_contrast)

saveRDS(stat_ctrl, "stat_controls.rds")


# TF activity, all four populations on one matrix
#
# Restricted to genes tested in every contrast so the ULM fits are comparable.
# No centring: these are already differences relative to naive.
stat_all <- c(list("NK" = stat_nk, "CD8 TEMRA" = stat_temra), stat_ctrl)

g_all <- Reduce(intersect, lapply(stat_all, names))
length(g_all)

mat_all <- do.call(cbind, lapply(stat_all, function(x) x[g_all]))
colnames(mat_all) <- names(stat_all)

colSums(is.na(mat_all)) 

tf_all_ct <- run_ulm(mat = mat_all, 
                     network = collectri,
                     .source = "source", 
                     .target = "target", 
                     .mor = "mor",
                     minsize = 5) %>%
  dplyr::filter(statistic == "ulm") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup()


# Of the CollecTRI regulators of HLH genes, how many carry significant activity
# in each population, and how strong is it? If the controls look like the
# effectors, the anchoring is not specific.
hlh_tfs <- unique(trn$source)

spec <- tf_all_ct %>%
  dplyr::filter(source %in% hlh_tfs) %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(n_tested = dplyr::n(),
                   n_sig = sum(p_adj < 0.05),
                   median_abs = round(median(abs(score)), 2),
                   max_abs = round(max(abs(score)), 2),
                   .groups = "drop")

spec

# Per-TF, for the regulators that were forced in the anchored runs
tf_all_ct %>%
  dplyr::filter(source %in% tf_aud$source) %>%
  dplyr::select(source, condition, score, p_adj) %>%
  tidyr::pivot_wider(names_from = condition, values_from = c(score, p_adj)) %>%
  dplyr::arrange(dplyr::desc(abs(score_NK)))

# Background comparison: the same summary over all TFs, not just HLH regulators.
# If the HLH regulators are no more distinctive than the background, that is the
# result rather than the controls being uninformative.
tf_all_ct %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(n_sig_all = sum(p_adj < 0.05),
                   median_abs_all = round(median(abs(score)), 2),
                   .groups = "drop")

saveRDS(tf_all_ct, "tf_all_celltypes.rds")

# Are HLH regulators over-represented among significant TFs, relative to that cell type's own background?
enrich <- tf_all_ct %>%
  dplyr::mutate(is_hlh = source %in% hlh_tfs, sig = p_adj < 0.05) %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(
    hlh_sig     = sum(is_hlh & sig),
    hlh_n       = sum(is_hlh),
    other_sig   = sum(!is_hlh & sig),
    other_n     = sum(!is_hlh),
    hlh_rate    = round(100 * hlh_sig / hlh_n, 1),
    other_rate  = round(100 * other_sig / other_n, 1),
    fold        = round((hlh_sig / hlh_n) / (other_sig / other_n), 2),
    p_fisher    = signif(fisher.test(
      matrix(c(hlh_sig, hlh_n - hlh_sig,
               other_sig, other_n - other_sig), nrow = 2))$p.value, 3),
    .groups = "drop")

enrich


# Figure: Specificity of the HLH regulatory layer      
#
# A1: how much more likely an HLH regulator is to be active than a random TF,
#     within each cell type. Fold of 1 means no enrichment.
# A2: which HLH gene each cell type's signal sits on.
ct_order <- c("Naive B", "classical monocytes", "NK", "CD8 TEMRA")

pA1 <- enrich %>%
  mutate(condition = factor(condition, levels = ct_order),
         lab = ifelse(p_fisher < 0.001,
                      sprintf("%.1fx\np = %.0e", fold, p_fisher),
                      sprintf("%.1fx\np = %.3f", fold, p_fisher))) %>%
  ggplot(aes(condition, fold, fill = condition)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = lab), vjust = -0.25, size = 3, lineheight = 0.9) +
  scale_fill_manual(values = c("Naive B" = "grey70",
                               "classical monocytes" = "#4F81BD",
                               "NK" = "#2166AC",
                               "CD8 TEMRA" = "#B2182B"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL,
       y = "Enrichment of pHLH gene regulators\namong active TFs (fold)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1))

# proportion rather than count, so genes with large regulons do not dominate
# which HLH genes each regulator targets
tf_targets <- trn %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(targets = paste(sort(unique(target)), collapse = ", "),
                   .groups = "drop")

by_gene <- trn %>%
  dplyr::inner_join(tf_all_ct, by = "source") %>%
  dplyr::group_by(condition, target) %>%
  dplyr::summarise(n_reg      = dplyr::n(),
                   n_sig      = sum(p_adj < 0.05),
                   median_abs = round(median(abs(score)), 2),
                   .groups = "drop") %>%
  tidyr::pivot_wider(names_from = condition,
                     values_from = c(n_sig, median_abs))

by_gene

pA2 <- by_gene %>%
  select(target, n_reg, starts_with("n_sig_")) %>%
  pivot_longer(starts_with("n_sig_"),
               names_to = "condition", values_to = "n_sig") %>%
  mutate(condition = factor(sub("^n_sig_", "", condition), levels = ct_order),
         prop = n_sig / n_reg,
         target = factor(target, levels = rev(sort(unique(target))))) %>%
  ggplot(aes(condition, target, fill = prop)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%d/%d", n_sig, n_reg)), size = 3) +
  scale_fill_gradient(low = "white", high = "#B2182B",
                      limits = c(0, 1), name = "Fraction of\nregulators\nactive") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1))

figA <- pA1 + pA2 + plot_annotation(tag_levels = "A")
ggsave("Figure_specificity.png", figA, width = 11, height = 4.5,
       dpi = 300, bg = "white")

