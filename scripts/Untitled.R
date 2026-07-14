k <- readRDS("carnival_results/gdT_koPRF1_carnival_result.rds")
b <- readRDS("carnival_results/gdT_baseline_carnival_result.rds")

sif <- as.data.frame(k$weightedSIF)
att <- as.data.frame(k$nodesAttributes)
batt <- as.data.frame(b$nodesAttributes)

sif$Weight <- as.numeric(sif$Weight)
att$AvgAct <- as.numeric(att$AvgAct)
batt$AvgAct <- as.numeric(batt$AvgAct)

# Core network: edges present in all 100 solutions
core <- sif %>% 
  filter(Weight == 100) 
nrow(core)   # 134

nodes <- unique(c(core$Node1, core$Node2))
length(nodes)

# Node attributes, with the baseline value and the delta
node_att <- att %>%
  filter(Node %in% nodes) %>%
  left_join(batt %>% select(Node, base = AvgAct), by = "Node") %>%
  mutate(delta = AvgAct - base,
         changed = abs(delta) >= 50)

table(node_att$changed)   # how many of the core nodes actually shifted

"PRF1" %in% nodes
core %>% 
  dplyr::filter(Node2 == "PRF1")

# and confirm which four changed
node_att %>% 
  dplyr::filter(changed) %>% 
  dplyr::select(Node, base, AvgAct, delta)

write.table(core, "gdT_koPRF1_core_network.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(node_att, "gdT_koPRF1_core_nodes.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)