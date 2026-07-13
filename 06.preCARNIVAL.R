pkn_files  <- list.files(pattern = "^pkn_(NK|gdT|CD8T)\\.rds$")
meas_files <- list.files(pattern = "measurements\\.rds$")

pkn_files
meas_files

sapply(list.files(pattern = "^pkn_(NK|gdT|CD8T)\\.rds$"), 
       function(f) paste(colnames(readRDS(f)), collapse = ", "))

for (tag in c("NK", "gdT", "CD8T")) {
  p     <- readRDS(sprintf("pkn_%s.rds", tag))
  nodes <- unique(c(p$source, p$target))
  base  <- readRDS(sprintf("%s_baseline_measurements.rds", tag))
  cat(tag, ":", sum(names(base) %in% nodes), "of", length(base), "in PKN\n")
}