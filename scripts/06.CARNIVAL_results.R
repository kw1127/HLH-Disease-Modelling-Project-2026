library(dplyr)
res_dir <- "carnival_results"

resolved <- function(z, u, d) pmax(z, u, d) >= 90

# ---- Connectivity of each HLH gene per cell type ----
for (ct in c("NK", "gdT", "CD8T")) {
  p <- readRDS(file.path(res_dir, sprintf("pkn_%s.rds", ct)))
  cat("\n---", ct, "---\n")
  for (g in c("PRF1", "RAB27A", "STX11", "STXBP2", "SH2D1A", "XIAP")) {
    cat(sprintf("%-8s in=%2d out=%2d\n", g, sum(p$target == g), sum(p$source == g)))
  }
}

# Resolved changes per knockout (pruned grid)
for (ct in c("NK", "gdT", "CD8T")) {
  
  b <- readRDS(file.path(res_dir, sprintf("%s_baseline_carnival_result.rds", ct)))$nodesAttributes
  
  for (kf in list.files(res_dir, pattern = sprintf("^%s_ko.*_carnival_result.rds$", ct))) {
    
    gene <- sub(sprintf("^%s_ko(.*)_carnival_result.rds$", ct), "\\1", kf)
    k <- readRDS(file.path(res_dir, kf))$nodesAttributes
    
    m <- full_join(
      b %>% select(Node, bz = ZeroAct, bu = UpAct, bd = DownAct, base = AvgAct),
      k %>% select(Node, kz = ZeroAct, ku = UpAct, kd = DownAct, ko = AvgAct),
      by = "Node") %>%
      mutate(across(c(base, ko), as.numeric)) %>%
      filter(resolved(bz, bu, bd), resolved(kz, ku, kd), abs(ko - base) >= 50, Node != gene)
    
    cat(sprintf("\n=== %s / %s -- %d resolved ===\n", ct, gene, nrow(m)))
    if (nrow(m)) print(m %>% select(Node, base, ko) %>% arrange(desc(abs(ko - base))), n = 30)
  }
}

# TF programme under PRF1 knockout 
tfs <- c("PRF1", "EOMES", "RUNX3", "STAT4", "STAT3", "STAT1", "ELF4", "TBX21")
for (ct in c("NK", "gdT", "CD8T")) {
  b <- readRDS(file.path(res_dir, sprintf("%s_baseline_carnival_result.rds", ct)))$nodesAttributes
  k <- readRDS(file.path(res_dir, sprintf("%s_koPRF1_carnival_result.rds", ct)))$nodesAttributes
  cat("\n---", ct, "---\n")
  for (t in tfs) {
    bv <- b$AvgAct[b$Node == t]; kv <- k$AvgAct[k$Node == t]
    if (length(bv) && length(kv)) cat(sprintf("%-8s %6s -> %6s\n", t, bv, kv))
  }
}

# Unpruned gdT / PRF1 comparison
b <- readRDS(file.path(res_dir, "gdTfull_baseline_carnival_result.rds"))$nodesAttributes
k <- readRDS(file.path(res_dir, "gdTfull_koPRF1_carnival_result.rds"))$nodesAttributes

m <- full_join(
  b %>% 
    select(Node, bz = ZeroAct, bu = UpAct, bd = DownAct, base = AvgAct),
  k %>% 
    select(Node, kz = ZeroAct, ku = UpAct, kd = DownAct, ko = AvgAct),
  by = "Node") %>% 
  mutate(across(c(base, ko), as.numeric)) %>%
  filter(resolved(bz, bu, bd), resolved(kz, ku, kd), abs(ko - base) >= 50, Node != "PRF1")

print(m %>% 
        select(Node, base, ko) %>% 
        arrange(desc(abs(ko - base))), n = 50)