# ============================================================
# 14. Solved networks: inspection, pruning and downstream join
#
#     weightedSIF returns the expanded network: the PKN plus one virtual
#     "Perturbation" node with an edge to every root node. Edges never
#     selected carry weight 0. The real network is the subset
#     with weight != 0 and the virtual node removed.
#
#     Three nested views are produced:
#       full - the full solved network
#       anchored - restricted to nodes upstream of the TFs regulating a
#                  differentially expressed HLH gene (supplementary figure)
#       pruned - anchored, minus nodes rare in the solution pool (main figure)
# ============================================================

library(dplyr); library(tibble); library(igraph)

tags <- c("nk", "temra") # the two contrasts, used as a suffix throughout
cond_of <- c(nk = "NK", temra = "CD8 TEMRA") # maps a tag to its label in tf_contrast
trn <- readRDS("trn.rds") # CollecTRI edges onto the HLH genes
hlh_de <- readRDS("hlh_de.rds") # HLH genes passing padj < 0.05, per contrast
tf_contrast <- readRDS("tf_contrast.rds") # ULM activity scores, both contrasts

# one CARNIVAL result object per contrast
carn <- lapply(setNames(tags, tags),
               function(tag) readRDS(sprintf("carnival_sig_%s.rds", tag)))


# ------------------------------------------------------------
# 14a. Extract the solved network and join the transcriptional layer
# ------------------------------------------------------------

build_layers <- function(tag) {
  
  raw <- as.data.frame(carn[[tag]]$weightedSIF) # expanded network, all candidate edges
  names(raw) <- c("source", "sign", "target", "weight")
  raw[c("sign", "weight")] <- lapply(raw[c("sign", "weight")], as.numeric) # arrive as character
  
  sig_layer <- raw %>% filter(weight != 0, source != "Perturbation") # the actual solution
  solved <- unique(c(sig_layer$source, sig_layer$target)) # nodes CARNIVAL recovered
  
  # Entry points the solver chose. A node can appear twice with opposite signs
  # where the pool disagrees on direction, hence unique().
  pert <- raw %>% filter(source == "Perturbation", weight != 0)
  pert_nodes <- unique(pert$target)
  
  att <- as.data.frame(carn[[tag]]$nodesAttributes, stringsAsFactors = FALSE)
  att[] <- lapply(att, function(x) if (is.list(x)) unlist(x) else x) # flatten list columns
  att$AvgAct <- as.numeric(att$AvgAct) # UpAct - DownAct, i.e. signed pool frequency
  
  # NodeType is "M" for measured TFs, blank for inferred signalling proteins.
  # Perturbation status is not in NodeType; it comes from the edges above.
  tf_nodes <- att$Node[att$NodeType == "M"]
  
  de <- hlh_de[[tag]] # DE HLH genes eligible for the downstream layer
  tfa <- tf_contrast %>% # measured activity of every TF in this contrast
    filter(condition == cond_of[[tag]]) %>%
    select(source, tf_score = score)
  
  # TF -> HLH edges. Not inferred: CollecTRI edges from TFs the model recovered.
  # Sign consistency required, so an inferred-active repressor cannot be drawn
  # pointing at an upregulated gene.
  trn_layer <- trn %>%
    filter(source %in% solved, target %in% de$target) %>% # regulator recovered, target DE
    left_join(select(de, target, log2FoldChange), by = "target") %>% # observed direction
    left_join(tfa, by = "source") %>% # inferred TF state
    filter(sign(tf_score) * sign(interaction) == sign(log2FoldChange)) %>% # consistency test
    transmute(source, target, sign = interaction,
              layer = "trn", weight = NA_real_) # no pool frequency: never in the ILP
  
  edges <- bind_rows(sig_layer %>% mutate(layer = "signalling"), trn_layer)
  
  # consistency (% of pool solutions) and logFC are different quantities and are
  # kept in separate columns; only sign is shared across layers.
  nodes <- bind_rows(
    att %>% filter(Node %in% solved) %>%
      transmute(id = Node,
                layer = case_when(Node %in% pert_nodes ~ "perturbation", # entry point
                                  NodeType == "M" ~ "TF", # measured
                                  TRUE ~ "signalling"), # inferred
                sign = sign(AvgAct), consistency = abs(AvgAct), logFC = NA_real_),
    de %>% filter(target %in% trn_layer$target) %>% # only genes that kept an edge
      transmute(id = target, layer = "HLH_gene",
                sign = sign(log2FoldChange), consistency = NA_real_,
                logFC = log2FoldChange)) %>%
    left_join(rename(tfa, id = source), by = "id") %>% # TFs carry both model and measurement
    distinct(id, .keep_all = TRUE) # a gene present in both layers appears once
  
  # Did the solver fit the measurements it was given? A sign-flipped TF has its
  # HLH targets hanging off a node whose inferred state contradicts the data.
  meas <- readRDS(sprintf("meas_%s_baseline.rds", tag))
  qc <- tibble(id = as.character(colnames(meas)),
               ulm = as.numeric(unlist(meas[1, ]))) %>% # what went in
    left_join(select(att, id = Node, AvgAct), by = "id") %>% # what came out
    mutate(status = case_when(is.na(AvgAct) ~ "dropped", # absent from the solution
                              AvgAct == 0 ~ "unexplained", # present but inferred inactive
                              sign(AvgAct) == sign(ulm) ~ "fit", # direction agrees
                              TRUE ~ "sign_flipped")) # direction contradicts the input
  
  # Where each DE HLH gene is lost: no CollecTRI regulator, regulator not
  # recovered by CARNIVAL, or recovered but sign-inconsistent.
  attrition <- de %>%
    mutate(n_regulators = vapply(target, function(g) sum(trn$target == g), integer(1)),
           n_recovered = vapply(target, function(g)
             sum(trn$target == g & trn$source %in% solved), integer(1)),
           n_consistent = vapply(target, function(g)
             sum(trn_layer$target == g), integer(1))) # survives into the figure
  
  list(edges = edges, nodes = nodes, qc = qc, attrition = attrition,
       summary = c(expanded_edges = nrow(raw), # PKN + virtual edges
                   virtual_edges = sum(raw$source == "Perturbation"), # = number of PKN roots
                   solved_edges = nrow(sig_layer), # the real network
                   solved_nodes = length(solved),
                   perturbations = length(pert_nodes), # inferred entry points
                   pert_that_are_tfs = sum(pert_nodes %in% tf_nodes), # 0 = no shallow entries
                   trn_edges = nrow(trn_layer))) # size of the downstream layer
}

net <- lapply(setNames(tags, tags), build_layers)

for (tag in tags) {
  cat("\n========", tag, "========\n")
  print(net[[tag]]$summary) # network sizes at each stage
  cat("\nmeasurement fit:\n"); print(table(net[[tag]]$qc$status)) # solver honoured the input?
  cat("\nHLH gene attrition:\n"); print(as.data.frame(net[[tag]]$attrition)) # why genes drop out
  cat("\nnode layers:\n"); print(table(net[[tag]]$nodes$layer)) # tier composition
}


# ------------------------------------------------------------
# 14b. How far upstream are the inferred entry points?
#
# The concern with free perturbations is that the solver enters the network
# immediately above each measured TF, giving no cascade. Distance 1 means a
# perturbation sits directly on a TF; 2+ means at least one inferred
# signalling protein lies between them.
# ------------------------------------------------------------

hops <- function(tag) {
  g <- graph_from_data_frame( # signalling layer only: TRN edges are not paths
    net[[tag]]$edges %>% filter(layer == "signalling") %>% select(source, target),
    directed = TRUE)
  p <- intersect(net[[tag]]$nodes$id[net[[tag]]$nodes$layer == "perturbation"], V(g)$name)
  t <- intersect(net[[tag]]$nodes$id[net[[tag]]$nodes$layer == "TF"], V(g)$name)
  d <- distances(g, v = p, to = t, mode = "out") # shortest path, each entry to each TF
  apply(d, 1, function(r) if (all(is.infinite(r))) NA_real_ else min(r)) # nearest TF per entry
}

for (tag in tags) {
  h <- hops(tag)
  cat("\n", tag, "- hops from perturbation to nearest measured TF:\n", sep = "")
  print(table(h, useNA = "ifany")) # mode of 2 = a genuine intermediate layer
  cat("direct (1 hop): ", paste(names(h)[!is.na(h) & h == 1], collapse = ", "), "\n", sep = "")
}


# ------------------------------------------------------------
# 14c. Evidence for the consistency filter
#
# The DUSP / PPM1 / PPP2 / PTPR phosphatases form a large family in OmniPath,
# all acting on the same few MAPK nodes. Where the ILP needs to set a MAPK
# state, any member will do, so the pool substitutes them freely. Two signals
# that this is solver behaviour rather than biology: low pool frequency, and
# disagreement on SIGN between members acting on the same targets.
# ------------------------------------------------------------

for (tag in tags) {
  cat("\n", tag, "- phosphatase family, pool frequency and sign:\n", sep = "")
  print(net[[tag]]$nodes %>%
          filter(grepl("^DUSP|^PPM1|^PPP2|^PTPR|^PTPN", id)) %>%
          select(id, layer, sign, consistency) %>%
          arrange(desc(consistency))) # low values plus mixed sign = interchangeable
}

for (tag in tags) {
  cat("\n", tag, "- consistency quantiles:\n", sep = "")
  print(quantile(net[[tag]]$nodes$consistency, na.rm = TRUE, # where the bulk of nodes sit
                 probs = c(0, .25, .5, .75, .9, 1)))
}


# ------------------------------------------------------------
# 14d. HLH-anchored subnetwork, then consistency pruning
# ------------------------------------------------------------

subset_to_hlh <- function(tag, order = 3, keep_fragments = TRUE) {
  e <- net[[tag]]$edges
  sig <- e %>% filter(layer == "signalling")
  g <- graph_from_data_frame(select(sig, source, target), directed = TRUE)
  
  anchor_tfs <- intersect(unique(e$source[e$layer == "trn"]), V(g)$name) # TFs with an HLH target
  keep <- unique(unlist(lapply(ego(g, order, anchor_tfs, mode = "in"), # everything upstream
                               function(v) V(g)$name[v])))
  
  # Weakly connected components containing a measured TF but no anchor: short
  # local explanations for TFs the model could not integrate into the cascade.
  if (keep_fragments) {
    tf_all <- net[[tag]]$nodes$id[net[[tag]]$nodes$layer == "TF"]
    comp <- components(g, mode = "weak")$membership
    frag <- names(comp)[comp %in% setdiff(comp[intersect(tf_all, names(comp))],
                                          comp[anchor_tfs])]
    keep <- union(keep, frag)
  }
  
  hlh <- unique(e$target[e$layer == "trn"])
  list(edges = e %>% filter((layer == "signalling" & # signalling edges within the kept set
                               source %in% keep & target %in% keep) | layer == "trn"), # TRN kept whole
       nodes = net[[tag]]$nodes %>% filter(id %in% c(keep, hlh)))
}

# HLH genes are exempt from the consistency filter: they are not CARNIVAL
# nodes and have no pool frequency.
prune <- function(x, min_cons = 15) {
  keep <- x$nodes$id[is.na(x$nodes$consistency) | x$nodes$consistency >= min_cons]
  e <- x$edges %>% filter(source %in% keep & target %in% keep) # edges between surviving nodes
  n <- x$nodes %>% filter(id %in% unique(c(e$source, e$target))) # drop nodes left isolated
  list(edges = e, nodes = n)
}

sub <- lapply(setNames(tags, tags), subset_to_hlh, order = 3)

# Threshold sweep. The elbow is the sharp drop between 0 and 10, after which
# the node count is nearly flat: everything below ~10% is interchangeable.
sweep <- sapply(c(0, 5, 10, 15, 20, 25, 50),
                function(k) sapply(sub, function(x) nrow(prune(x, k)$nodes)))
colnames(sweep) <- paste0(c(0, 5, 10, 15, 20, 25, 50), "%")
cat("\nnodes retained by consistency threshold:\n"); print(sweep)

MIN_CONS <- 15 # same threshold for both contrasts, or the panels are not comparable
pruned <- lapply(sub, prune, min_cons = MIN_CONS)

for (tag in tags) {
  cat("\n========", tag, "========\n")
  cat("full solution:  ", nrow(net[[tag]]$nodes), " nodes, ", # everything recovered
      nrow(net[[tag]]$edges), " edges\n", sep = "")
  cat("HLH-anchored:   ", nrow(sub[[tag]]$nodes), " nodes, ", # after the topology filter
      nrow(sub[[tag]]$edges), " edges\n", sep = "")
  cat("pruned >=", MIN_CONS, "%:   ", nrow(pruned[[tag]]$nodes), " nodes, ", # main figure
      nrow(pruned[[tag]]$edges), " edges\n", sep = "")
  print(table(pruned[[tag]]$nodes$layer)) # all four tiers should survive
  cat("\nTF -> HLH edges retained:\n")
  print(pruned[[tag]]$edges %>% filter(layer == "trn") %>% select(source, target, sign))
}


# ------------------------------------------------------------
# 14e. Export for Cytoscape
#     *_pruned = main figure, *_hlh = supplementary, *_full = supplementary
#     Import the edge file as a network, then the node file keyed on `id`.
# ------------------------------------------------------------

for (tag in tags) {
  for (v in c("pruned", "hlh", "full")) {
    x <- switch(v, pruned = pruned[[tag]], hlh = sub[[tag]], full = net[[tag]])
    write.table(x$edges, sprintf("edges_%s_%s.tsv", tag, v),
                sep = "\t", row.names = FALSE, quote = FALSE, na = "") # na = "" so Cytoscape
    write.table(x$nodes, sprintf("nodes_%s_%s.tsv", tag, v),          # reads blanks as missing
                sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  }
}

saveRDS(list(full = net, anchored = sub, pruned = pruned, min_cons = MIN_CONS),
        "network_layers.rds")


# ------------------------------------------------------------
# 14f. Comparison between contrasts
#
# The transcriptional layer is where the two effector populations can be
# compared directly, since both were solved over identical topology.
# ------------------------------------------------------------

trn_of <- function(x) unique(paste(x$edges$source[x$edges$layer == "trn"],
                                   x$edges$target[x$edges$layer == "trn"],
                                   sep = " -> ")) # regulator-target pairs, not TFs alone

cat("\nTF -> HLH edges, pruned networks:\n")
print(list(shared = intersect(trn_of(pruned$nk), trn_of(pruned$temra)), # common core
           nk_only = setdiff(trn_of(pruned$nk), trn_of(pruned$temra)),
           temra_only = setdiff(trn_of(pruned$temra), trn_of(pruned$nk))))

# A regulator present in one contrast and not the other has three possible
# explanations, and only the last is about signalling:
#   (i)   never measured in that contrast (did not clear FDR / top-n cut)
#   (ii)  measured but not recovered by the solver
#   (iii) recovered but below the consistency threshold
regulator_status <- function(from, to) {
  tf_extra <- setdiff(unique(sub[[from]]$edges$source[sub[[from]]$edges$layer == "trn"]),
                      unique(sub[[to]]$edges$source[sub[[to]]$edges$layer == "trn"]))
  meas_to <- colnames(readRDS(sprintf("meas_%s_baseline.rds", to))) # what the other contrast had
  tibble(tf = tf_extra, measured_in_other = tf_extra %in% meas_to) %>% # FALSE = case (i)
    left_join(tf_contrast %>% filter(condition == cond_of[[to]]) %>%
                select(tf = source, score, p_adj), by = "tf") # its score there, if any
}

cat("\nregulators in TEMRA but not NK - were they measurable in NK?\n")
print(regulator_status("temra", "nk"))

cat("\nregulators in NK but not TEMRA - were they measurable in TEMRA?\n")
print(regulator_status("nk", "temra"))

# Upstream layer: which signalling proteins are shared vs contrast-specific
sig_of <- function(x) x$nodes$id[x$nodes$layer %in% c("signalling", "perturbation")]
cat("\nsignalling nodes (pruned):\n")
print(list(shared = intersect(sig_of(pruned$nk), sig_of(pruned$temra)), # common wiring
           nk_only = setdiff(sig_of(pruned$nk), sig_of(pruned$temra)),
           temra_only = setdiff(sig_of(pruned$temra), sig_of(pruned$nk)))) # the NF-kB module