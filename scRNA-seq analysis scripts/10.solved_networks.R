# =============================================================================
# 15. Solved networks: inspection, pruning and downstream join
#
# Runs over both carnival runs (forced TFs or not). The measurement set is the only difference, 
# so anything present in one network and absent from the other is a result of anchoring.
#
#   baseline - measurements were the top-50 significant TFs only
#   anchored - the same TFs + the forced CollecTRI regulators of HLH genes
#
# weightedSIF returns the expanded network: the PKN plus one synthetic
# perturbation node with an edge to every root node. Weight is the percentage
# of optimal solutions containing that edge, so weight 0 means the edge was
# never used. The real network is the non-zero subset with the virtual node
# removed.
#
# Three nested views are produced per variant:
#   full - the whole solved network
#   anchored - restricted to nodes upstream of the TFs regulating a
#              differentially expressed HLH gene (supplementary figure)
#   pruned - the above, minus nodes rare in the solution pool (main figure)
# =============================================================================
tags <- c("nk", "temra")                       
variants <- c("baseline", "anchored")            
cond_of <- c(nk = "NK", temra = "CD8 TEMRA")       

trn <- readRDS("trn.rds")          
hlh_de <- readRDS("hlh_de.rds") # HLH genes at padj < 0.05, per contrast
tf_contrast <- readRDS("tf_contrast.rds")   # ULM activity scores, both contrasts
tf_aud <- readRDS("tf_aud.rds") # the forced HLH regulators

# file naming differs between the two runs
carn_file <- function(tag, variant) {
  if (variant == "baseline") sprintf("carnival_sig_%s.rds", tag)
  else                       sprintf("carnival_anch_%s.rds", tag)
}

meas_file <- function(tag, variant) {
  if (variant == "baseline") sprintf("meas_%s_baseline.rds", tag)
  else                       sprintf("meas_%s_hlh_anchored.rds", tag)
}

# Extract the solved network and join the transcriptional layer to all of them
build_layers <- function(tag, variant) {
  
  carn <- readRDS(carn_file(tag, variant))
  
  raw <- as.data.frame(carn$weightedSIF)   # expanded network, all candidate edges
  names(raw) <- c("source", "sign", "target", "weight")
  raw[c("sign", "weight")] <- lapply(raw[c("sign", "weight")], as.numeric)  # arrive as character
  
  sig_layer <- raw %>% 
    filter(weight != 0, source != "Perturbation") # the actual solution
  
  solved <- unique(c(sig_layer$source, sig_layer$target)) # nodes CARNIVAL recovered
  
  # Entry points the solver chose. A node can appear twice with opposite signs
  # where the pool disagrees on direction, hence unique().
  pert_nodes <- unique(raw$target[raw$source == "Perturbation" & raw$weight != 0])
  
  att <- as.data.frame(carn$nodesAttributes, stringsAsFactors = FALSE)
  att[] <- lapply(att, function(x) if (is.list(x)) unlist(x) else x)  # flatten list columns
  att$AvgAct <- as.numeric(att$AvgAct) # UpAct - DownAct
  
  # NodeType is "M" for measured TFs, blank for inferred signalling proteins.
  # Under anchoring this set also contains the forced HLH regulators.
  # Perturbation status is not in NodeType; it comes from the edges above.
  tf_nodes <- att$Node[att$NodeType == "M"]
  
  de  <- hlh_de[[tag]] # DE HLH genes eligible for the downstream layer
  tfa <- tf_contrast %>% # measured activity of every TF in this contrast
    filter(condition == cond_of[[tag]]) %>%
    select(source, tf_score = score)
  
  # TF -> HLH edges. Not inferred: CollecTRI edges from TFs the model recovered.
  # Sign consistency is needed, so an inferred-active repressor cannot be drawn
  # pointing at an upregulated gene. This is the only place the asserted edges
  # are tested against data.
  trn_layer <- trn %>%
    filter(source %in% solved, target %in% de$target) %>% # regulator recovered, target DE
    left_join(select(de, target, log2FoldChange), by = "target") %>%  # observed direction
    left_join(tfa, by = "source") %>% # inferred TF state
    filter(sign(tf_score) * sign(interaction) == sign(log2FoldChange)) %>%
    transmute(source, target, sign = interaction,
              layer = "trn", weight = NA_real_)  # no pool frequency: never in the ILP
  
  edges <- bind_rows(sig_layer %>% mutate(layer = "signalling"), trn_layer)
  
  # consistency (% of pool solutions) and logFC are different quantities and are
  # kept in separate columns; only sign is shared across layers.
  nodes <- bind_rows(
    att %>% filter(Node %in% solved) %>%
      transmute(id = Node,
                layer = case_when(Node %in% pert_nodes ~ "perturbation",  # entry point
                                  NodeType == "M" ~ "TF", # measured
                                  TRUE ~ "signalling"),   # inferred
                sign = sign(AvgAct), consistency = abs(AvgAct), logFC = NA_real_),
    de %>% filter(target %in% trn_layer$target) %>% # only genes that kept an edge
      transmute(id = target, layer = "HLH_gene",
                sign = sign(log2FoldChange), consistency = NA_real_,
                logFC = log2FoldChange)) %>%
    left_join(rename(tfa, id = source), by = "id") %>%  # TFs carry model and measurement
    distinct(id, .keep_all = TRUE) # a gene in both layers appears once
  
  # Did the solver fit the measurements it was given? A sign-flipped TF has its
  # HLH targets hanging off a node whose inferred state contradicts the data.
  # Under anchoring, "dropped" and "unexplained" are expected for forced TFs with
  # near-zero scores, which contribute almost nothing to the objective.
  meas <- readRDS(meas_file(tag, variant))
  qc <- tibble(id  = as.character(colnames(meas)),
               ulm = as.numeric(unlist(meas[1, ]))) %>%          # what went in
    left_join(select(att, id = Node, AvgAct), by = "id") %>%     # what came out
    mutate(status = case_when(is.na(AvgAct) ~ "dropped",      # absent from solution
                              AvgAct == 0 ~ "unexplained",  # present, inferred inactive
                              sign(AvgAct) == sign(ulm) ~ "fit", # direction agrees
                              TRUE ~ "sign_flipped"))
  
  # Where each DE HLH gene is lost: no CollecTRI regulator, regulator not
  # recovered by CARNIVAL, or recovered but sign-inconsistent.
  attrition <- de %>%
    mutate(n_regulators = vapply(target, function(g) sum(trn$target == g), integer(1)),
           n_recovered  = vapply(target, function(g)
             sum(trn$target == g & trn$source %in% solved), integer(1)),
           n_consistent = vapply(target, function(g)
             sum(trn_layer$target == g), integer(1)))   # survives into the figure
  
  list(edges = edges, nodes = nodes, qc = qc, attrition = attrition,
       summary = c(expanded_edges    = nrow(raw),                            # PKN + virtual
                   virtual_edges     = sum(raw$source == "Perturbation"),    # = number of roots
                   solved_edges      = nrow(sig_layer),                      # the real network
                   solved_nodes      = length(solved),
                   perturbations     = length(pert_nodes),                   # inferred entry points
                   pert_that_are_tfs = sum(pert_nodes %in% tf_nodes),        # 0 = no shallow entries
                   trn_edges         = nrow(trn_layer)))                     # downstream layer size
}

net <- lapply(setNames(variants, variants), function(v)
  lapply(setNames(tags, tags), build_layers, variant = v))

# Figure: What anchoring changed, what did it add to the network?                 
#
# The transcriptional layer is the only part of the network the anchoring was
# meant to affect, so it is the direct measure of what forcing achieved.
anchor_effect <- do.call(rbind, lapply(variants, function(v)
  do.call(rbind, lapply(tags, function(tag)
    data.frame(variant  = v,
               contrast = cond_of[[tag]],
               trn      = net[[v]][[tag]]$summary[["trn_edges"]],
               sig      = net[[v]][[tag]]$summary[["solved_edges"]])))))

pB <- anchor_effect %>%
  mutate(variant = factor(variant, levels = c("baseline", "anchored"),
                          labels = c("Unanchored", "HLH-anchored"))) %>%
  ggplot(aes(contrast, trn, fill = variant)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = trn), position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("Unanchored" = "grey70",
                               "HLH-anchored" = "#B2182B"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "TF -> HLH gene edges recovered") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "bottom")

ggsave("Figure_anchoring_effect.png", pB, width = 5, height = 4.5,
       dpi = 300, bg = "white")


for (v in variants) for (tag in tags) {
  cat("\n========", v, "/", tag, "========\n")
  print(net[[v]][[tag]]$summary)  # sizes at each stage
  cat("\nmeasurement fit:\n");     print(table(net[[v]][[tag]]$qc$status))
  cat("\nHLH gene attrition:\n");  print(as.data.frame(net[[v]][[tag]]$attrition))
  cat("\nnode layers:\n");         print(table(net[[v]][[tag]]$nodes$layer))
}

# Sign-consistency rate: of the CollecTRI edges eligible for testing (regulator
# recovered, target DE), what fraction predicted the observed direction?
conc <- do.call(rbind, lapply(variants, function(v)
  do.call(rbind, lapply(tags, function(tag) {
    a <- net[[v]][[tag]]$attrition
    data.frame(variant = v, contrast = tag,
               tested = sum(a$n_recovered),
               passed = sum(a$n_consistent),
               rate   = round(100 * sum(a$n_consistent) / sum(a$n_recovered), 1))
  }))))


# Which forced TFs did the solver actually set?
#
# Forcing puts a TF in the objective; it does not guarantee the solution reaches
# it. Only the TFs the solver set support a claim about the HLH layer.
for (tag in tags) {
  forced <- tf_aud$source[tf_aud$contrast == tag]
  cat("\n", tag, " - forced HLH regulators (anchored solve):\n", sep = "")
  print(net$anchored[[tag]]$qc %>%
          filter(id %in% forced) %>%
          arrange(desc(abs(ulm))) %>%
          as.data.frame())
}

# How far upstream are the inferred entry points?
#
# The concern with free perturbations is that the solver enters the network
# immediately above each measured TF, giving no cascade. Distance 1 means a
# perturbation sits directly on a TF; 2+ means at least one inferred signalling
# protein lies between them.
hops <- function(x) {
  g <- graph_from_data_frame(   # signalling layer only: TRN edges are not paths
    x$edges %>% filter(layer == "signalling") %>% select(source, target),
    directed = TRUE)
  p <- intersect(x$nodes$id[x$nodes$layer == "perturbation"], V(g)$name)
  t <- intersect(x$nodes$id[x$nodes$layer == "TF"], V(g)$name)
  d <- distances(g, v = p, to = t, mode = "out")   # each entry point to each TF
  apply(d, 1, function(r) if (all(is.infinite(r))) NA_real_ else min(r))
}

for (v in variants) for (tag in tags) {
  h <- hops(net[[v]][[tag]])
  cat("\n", v, "/", tag, " - hops from perturbation to nearest measured TF:\n", sep = "")
  print(table(h, useNA = "ifany"))   # mode of 2 = a genuine intermediate layer
  cat("direct (1 hop): ", paste(names(h)[!is.na(h) & h == 1], collapse = ", "), "\n", sep = "")
}


# Evidence for the consistency filter
#
# The DUSP / PPM1 / PPP2 / PTPR phosphatases form a large family in OmniPath,
# all acting on the same few MAPK nodes. Where the ILP needs to set a MAPK state
# any member will do, so the pool substitutes them freely. Two signals that this
# is solver behaviour rather than biology: low pool frequency, and disagreement
# on sign between members acting on the same targets.
for (v in variants) for (tag in tags) {
  cat("\n", v, "/", tag, " - phosphatase family, pool frequency and sign:\n", sep = "")
  print(net[[v]][[tag]]$nodes %>%
          filter(grepl("^DUSP|^PPM1|^PPP2|^PTPR|^PTPN", id)) %>%
          select(id, layer, sign, consistency) %>%
          arrange(desc(consistency)))   # low values plus mixed sign = interchangeable
}

for (v in variants) for (tag in tags) {
  cat("\n", v, "/", tag, " - consistency quantiles:\n", sep = "")
  print(quantile(net[[v]][[tag]]$nodes$consistency, na.rm = TRUE,
                 probs = c(0, .25, .5, .75, .9, 1)))
}


# HLH-anchored subnetwork, then consistency pruning
subset_to_hlh <- function(x, order = 3, keep_fragments = TRUE) {
  e   <- x$edges
  sig <- e %>% filter(layer == "signalling")
  g   <- graph_from_data_frame(select(sig, source, target), directed = TRUE)
  
  anchor_tfs <- intersect(unique(e$source[e$layer == "trn"]), V(g)$name)  # TFs with an HLH target
  keep <- unique(unlist(lapply(ego(g, order, anchor_tfs, mode = "in"),    # everything upstream
                               function(v) V(g)$name[v])))
  
  # Weakly connected components containing a measured TF but no anchor: short
  # local explanations for TFs the model could not integrate into the cascade.
  if (keep_fragments) {
    tf_all <- x$nodes$id[x$nodes$layer == "TF"]
    comp   <- components(g, mode = "weak")$membership
    frag   <- names(comp)[comp %in% setdiff(comp[intersect(tf_all, names(comp))],
                                            comp[anchor_tfs])]
    keep <- union(keep, frag)
  }
  
  hlh <- unique(e$target[e$layer == "trn"])
  list(edges = e %>% filter((layer == "signalling" &          # signalling edges within the set
                               source %in% keep & target %in% keep) |
                              layer == "trn"),                # TRN kept whole
       nodes = x$nodes %>% filter(id %in% c(keep, hlh)))
}

# HLH genes are exempt from the consistency filter: they are not CARNIVAL nodes
# and have no pool frequency.
prune <- function(x, min_cons = 15) {
  keep <- x$nodes$id[is.na(x$nodes$consistency) | x$nodes$consistency >= min_cons]
  e <- x$edges %>% filter(source %in% keep & target %in% keep)          # edges between survivors
  n <- x$nodes %>% filter(id %in% unique(c(e$source, e$target)))        # drop isolated nodes
  list(edges = e, nodes = n)
}

sub <- lapply(net, function(vv) lapply(vv, subset_to_hlh, order = 3))

# Threshold sweep. The elbow is the sharp drop between 0 and 10, after which the
# node count is nearly flat: everything below ~10% is interchangeable.
ks <- c(0, 5, 10, 15, 20, 25, 50)
sweep <- do.call(rbind, lapply(variants, function(v)
  do.call(rbind, lapply(tags, function(tag)
    data.frame(variant = v, contrast = tag,
               setNames(as.list(sapply(ks, function(k)
                 nrow(prune(sub[[v]][[tag]], k)$nodes))), paste0(ks, "%")),
               check.names = FALSE)))))
cat("\nnodes retained by consistency threshold:\n"); print(sweep)

MIN_CONS <- 15   # one threshold for all four networks, or the panels are not comparable
pruned <- lapply(sub, function(vv) lapply(vv, prune, min_cons = MIN_CONS))

for (v in variants) for (tag in tags) {
  cat("\n========", v, "/", tag, "========\n")
  cat("full solution:  ", nrow(net[[v]][[tag]]$nodes), " nodes, ",
      nrow(net[[v]][[tag]]$edges), " edges\n", sep = "")
  cat("HLH-anchored:   ", nrow(sub[[v]][[tag]]$nodes), " nodes, ",
      nrow(sub[[v]][[tag]]$edges), " edges\n", sep = "")
  cat("pruned >=", MIN_CONS, "%:   ", nrow(pruned[[v]][[tag]]$nodes), " nodes, ",
      nrow(pruned[[v]][[tag]]$edges), " edges\n", sep = "")
  print(table(pruned[[v]][[tag]]$nodes$layer))   # all four tiers should survive
  cat("\nTF -> HLH edges retained:\n")
  print(pruned[[v]][[tag]]$edges %>% filter(layer == "trn") %>% select(source, target, sign))
}

# Figure: node-consistency plot                      
#
# Every node carries the percentage of optimal solutions in which it was
# assigned a non-zero state. This shows the distribution before pruning, so the
# threshold can be seen
cons_df <- do.call(rbind, lapply(variants, function(v)
  do.call(rbind, lapply(tags, function(tag)
    sub[[v]][[tag]]$nodes %>%
      filter(!is.na(consistency)) %>%      # HLH genes have no pool frequency
      select(id, layer, consistency) %>%
      mutate(variant  = v,
             contrast = cond_of[[tag]])))))

cons_df <- cons_df %>%
  mutate(panel = factor(paste0(ifelse(variant == "baseline",
                                      "Unanchored", "Anchored"),
                               " - ", contrast),
                        levels = c("Unanchored - NK", "Anchored - NK",
                                   "Unanchored - CD8 TEMRA",
                                   "Anchored - CD8 TEMRA")),
         layer = factor(layer, levels = c("perturbation", "signalling", "TF")))

# nodes ranked within each network, so the shape of the drop-off is visible
cons_df <- cons_df %>%
  group_by(panel) %>%
  arrange(desc(consistency), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

node_consistency <- ggplot(cons_df, aes(rank, consistency, colour = layer)) +
  geom_hline(yintercept = MIN_CONS, linetype = "dashed", colour = "grey40") +
  geom_point(size = 1.4, alpha = 0.85) +
  facet_wrap(~ panel, nrow = 2, scales = "free_x") +
  scale_colour_manual(values = c(perturbation = "#E8A33D",
                                 signalling = "#4F81BD",
                                 TF = "#B2182B"),
                      labels = c("Inferred entry point", "Signalling protein",
                                 "Measured TF"),
                      name = NULL) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "Nodes, ranked within each network",
       y = "Solutions containing the node (%)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"))

ggsave("Figure_pool_consistency.png", node_consistency, width = 9, height = 6,
       dpi = 300, bg = "white")

# How many nodes sit at 100%?
cons_df %>%
  group_by(panel) %>%
  summarise(n = n(),
            at_100 = sum(consistency == 100),
            below_thresh = sum(consistency < MIN_CONS), .groups = "drop")

# Export for Cytoscape
#
# pruned = main figure, hlh and full = supplementary.
# Import the edge file as a network, then the node file keyed on `id`.
# Node column `layer` drives shape; `sign` and `consistency` drive colour and
# transparency; edge column `layer` distinguishes solved signalling edges from
# asserted CollecTRI edges and should be styled differently (solid vs dashed).
for (v in variants) for (tag in tags) {
  for (view in c("pruned", "hlh", "full")) {
    x <- switch(view,
                pruned = pruned[[v]][[tag]],
                hlh    = sub[[v]][[tag]],
                full   = net[[v]][[tag]])
    write.table(x$edges, sprintf("edges_%s_%s_%s.tsv", tag, view, v),
                sep = "\t", row.names = FALSE, quote = FALSE, na = "")  # na = "" so Cytoscape
    write.table(x$nodes, sprintf("nodes_%s_%s_%s.tsv", tag, view, v),   # reads blanks as missing
                sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  }
}

saveRDS(list(full = net, anchored = sub, pruned = pruned, min_cons = MIN_CONS),
        "network_layers.rds")


# Comparison between contrasts
#
# The transcriptional layer is where the two effector populations can be
# compared directly, since both were solved over identical topology.
trn_of <- function(x) unique(paste(x$edges$source[x$edges$layer == "trn"],
                                   x$edges$target[x$edges$layer == "trn"],
                                   sep = " -> "))   # regulator-target pairs, not TFs alone

for (v in variants) {
  cat("\n", v, " - TF -> HLH edges, pruned networks:\n", sep = "")
  print(list(shared     = intersect(trn_of(pruned[[v]]$nk), trn_of(pruned[[v]]$temra)),
             nk_only    = setdiff(trn_of(pruned[[v]]$nk), trn_of(pruned[[v]]$temra)),
             temra_only = setdiff(trn_of(pruned[[v]]$temra), trn_of(pruned[[v]]$nk))))
}

# A regulator present in one contrast and not the other has three possible
# explanations, and only the last is about signalling:
#   (i)   never measured in that contrast (did not clear FDR / top-n cut)
#   (ii)  measured but not recovered by the solver
#   (iii) recovered but below the consistency threshold
regulator_status <- function(v, from, to) {
  tf_extra <- setdiff(unique(sub[[v]][[from]]$edges$source[sub[[v]][[from]]$edges$layer == "trn"]),
                      unique(sub[[v]][[to]]$edges$source[sub[[v]][[to]]$edges$layer == "trn"]))
  meas_to <- colnames(readRDS(meas_file(to, v)))   # what the other contrast had
  tibble(tf = tf_extra, measured_in_other = tf_extra %in% meas_to) %>%   # FALSE = case (i)
    left_join(tf_contrast %>% filter(condition == cond_of[[to]]) %>%
                select(tf = source, score, p_adj), by = "tf")
}

for (v in variants) {
  cat("\n", v, " - regulators in TEMRA but not NK, were they measurable in NK?\n", sep = "")
  print(regulator_status(v, "temra", "nk"))
  cat("\n", v, " - regulators in NK but not TEMRA, were they measurable in TEMRA?\n", sep = "")
  print(regulator_status(v, "nk", "temra"))
}

# Upstream layer: which signalling proteins are shared vs contrast-specific
sig_of <- function(x) x$nodes$id[x$nodes$layer %in% c("signalling", "perturbation")]

for (v in variants) {
  cat("\n", v, " - signalling nodes (pruned):\n", sep = "")
  print(list(shared     = intersect(sig_of(pruned[[v]]$nk), sig_of(pruned[[v]]$temra)),
             nk_only    = setdiff(sig_of(pruned[[v]]$nk), sig_of(pruned[[v]]$temra)),
             temra_only = setdiff(sig_of(pruned[[v]]$temra), sig_of(pruned[[v]]$nk))))
}


# What did anchoring change?
#
# Same everything, but the measurement set is the only
# difference. Nodes present in the anchored solve and absent from the baseline
# are those the forced HLH regulators recruited.
for (tag in tags) {
  a <- net$anchored[[tag]]$nodes$id
  b <- net$baseline[[tag]]$nodes$id
  cat("\n", tag, ": ", length(a), " nodes anchored vs ", length(b), " baseline\n", sep = "")
  cat("gained by anchoring: ", paste(setdiff(a, b), collapse = ", "), "\n", sep = "")
  cat("lost: ",               paste(setdiff(b, a), collapse = ", "), "\n", sep = "")
  cat("TF -> HLH edges: anchored ", sum(net$anchored[[tag]]$edges$layer == "trn"),
      ", baseline ",               sum(net$baseline[[tag]]$edges$layer == "trn"), "\n", sep = "")
}