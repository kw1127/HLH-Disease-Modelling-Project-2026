# ============================================================================
# GOSH cohort data cleaning and analysis
# ============================================================================

library(tidyverse)
library(janitor)
library(readxl)
library(openxlsx)
library(ggpubr)
library(ComplexHeatmap)
library(grid)
library(patchwork)

# ============================================================================
# Original dataset: PRF1 panel data
# ============================================================================

data <- read_excel("Patient Data /HLH_data.xlsx", sheet = "TOTAL")

# Only keep patients
# Perforin is in text so convert it to numeric for analysis
data <- data %>%
  slice(1:112) %>%
  mutate(Perforin_expression_percent = as.numeric(`Perforin expression %`)) %>%
  mutate(sex = as.factor(Sex)) %>%
  clean_names()

# Drop columns not used
data_clean <- data %>%
  rename(`Perforin Expression %` = perforin_expression_percent_2) %>%
  rename(`Patient ID` = id) %>%
  select(-"parent_nationality", -"familial_disease", -"parental_consanguinity",
         -"dob", -"episode_number", -"age", -"coll_date", -"rec_date", -"source",
         -"location", -"clinician", -"test_set_name", -"mutation",
         -"perforin_expression_percent", -"perforin_code", -"peak_mfi", -"gra",
         -"s_cd25_pg_ml", -"perforin_expression",
         -('hemoglobin_g_l':'treatment_applied'),
         -"sex_2") %>%
  mutate(across(starts_with("c_"), as.numeric))

# Only keep patients with recorded genetics
data_clean <- data_clean %>%
  filter(if_all(starts_with("c_"), ~ !is.na(.)))

# Three patients had no perforin value so impute them by state average
data_clean <- data_clean %>%
  mutate(`Perforin Expression %` = case_when(
    `Patient ID` == 91 ~ 0.56,
    `Patient ID` == 103 ~ 0.52,
    `Patient ID` == 21 ~ 0.19,
    TRUE ~ `Perforin Expression %`)) %>%
  mutate(`Perforin Expression %` = `Perforin Expression %` * 100)

# Clean up the diagnoses into three groups
data_labelled <- data_clean %>%
  mutate(diagnosis_group = case_when(
    diagnosis == "HLH" ~ "HLH",
    is.na(diagnosis) ~ "Unknown",
    TRUE ~ "Non-HLH")) %>%
  mutate(diagnosis_group = factor(diagnosis_group,
                                  levels = c("HLH", "Non-HLH", "Unknown"))) %>%
  select(-"sex", -"diagnosis") %>%
  arrange(diagnosis_group)

# ============================================================================
# Clinical data collected from GOSH
# ============================================================================

patients <- read_excel("Clinical Data/combined_patients.xlsx", sheet = "Original") %>%
  filter(!is.na(`GOSH MRN`)) %>%
  select(-starts_with("..."))

sCD25 <- read_excel("Clinical Data/sCD25.xlsx") %>%
  mutate(`sCD25 (pg/ml)` = as.numeric(gsub("[<>]", "", Value))) %>%
  filter(!is.na(`GOSH MRN`), !is.na(Value)) %>%
  group_by(`GOSH MRN`) %>%
  slice_max(`Collection Date (Best Available)`, n = 1, with_ties = FALSE) %>%
  ungroup()

patients <- patients %>%
  select(-any_of(c("Collection Date (sCD25)", "sCD25", "sCD25 (pg/ml)"))) %>%
  left_join(
    sCD25 %>%
      transmute(`GOSH MRN`,
                `Collection Date (sCD25)` = `Collection Date (Best Available)`,
                sCD25 = Value,
                `sCD25 (pg/ml)`),
    by = "GOSH MRN", relationship = "many-to-one", na_matches = "never")

triglycerides <- read_excel("Clinical Data/Tri.xlsx") %>%
  mutate(`triglycerides (mmol/g)` =
           as.numeric(gsub("[<>]", "", `Triglycerides (mmol/g)`))) %>%
  filter(!is.na(`GOSH MRN`), !is.na(`Triglycerides (mmol/g)`)) %>%
  group_by(`GOSH MRN`) %>%
  slice_max(`Collection Date (Best Available)`, n = 1, with_ties = FALSE) %>%
  ungroup()

patients <- patients %>%
  select(-any_of(c("Collection Date (Triglycerides)", "Triglycerides (mmol/g)",
                   "triglycerides (mmol/g)"))) %>%
  left_join(
    triglycerides %>%
      transmute(`GOSH MRN`,
                `Collection Date (Triglycerides)` = `Collection Date (Best Available)`,
                `Triglycerides (mmol/g)`,
                `triglycerides (mmol/g)`),
    by = "GOSH MRN", relationship = "many-to-one", na_matches = "never")

fibrinogen <- read_excel("Clinical Data/fibr.xlsx") %>%
  mutate(`fibrinogen (g/L)` = as.numeric(gsub("[<>]", "", `Fibrinogen (g/L)`))) %>%
  filter(!is.na(`GOSH MRN`), !is.na(`Fibrinogen (g/L)`)) %>%
  group_by(`GOSH MRN`) %>%
  slice_max(`Collection Date (Fibrinogen)`, n = 1, with_ties = FALSE) %>%
  ungroup()

patients <- patients %>%
  select(-any_of(c("Collection Date (Fibrinogen)", "Fibrinogen (g/L)",
                   "fibrinogen (g/L)"))) %>%
  left_join(
    fibrinogen %>%
      transmute(`GOSH MRN`, `Collection Date (Fibrinogen)`,
                `Fibrinogen (g/L)`, `fibrinogen (g/L)`),
    by = "GOSH MRN", relationship = "many-to-one", na_matches = "never")

ferritin <- read_excel("Clinical Data/ferritin.xlsx") %>%
  mutate(`ferritin (ug/L)` = as.numeric(gsub("[<>]", "", `Ferritin (ug/L)`))) %>%
  filter(!is.na(`GOSH MRN`), !is.na(`Ferritin (ug/L)`)) %>%
  group_by(`GOSH MRN`) %>%
  slice_max(`Collection Date (Ferritin)`, n = 1, with_ties = FALSE) %>%
  ungroup()

patients <- patients %>%
  select(-any_of(c("Collection Date (Ferritin)", "Ferritin (ug/L)",
                   "ferritin (ug/L)"))) %>%
  left_join(
    ferritin %>%
      transmute(`GOSH MRN`, `Collection Date (Ferritin)`,
                `Ferritin (ug/L)`, `ferritin (ug/L)`),
    by = "GOSH MRN", relationship = "many-to-one", na_matches = "never")

haemoglobin <- read_excel("Clinical Data/Haemoglobin.xlsx") %>%
  filter(!is.na(`GOSH MRN`), !is.na(`Haemoglobin (g/L)`)) %>%
  group_by(`GOSH MRN`) %>%
  slice_max(`Collection Date (Haemoglobin)`, n = 1, with_ties = FALSE) %>%
  ungroup()

patients <- patients %>%
  select(-any_of(c("Haemoglobin (g/L)", "Collection Date (Haemoglobin)"))) %>%
  left_join(
    haemoglobin %>%
      transmute(`GOSH MRN`, `Collection Date (Haemoglobin)`, `Haemoglobin (g/L)`),
    by = "GOSH MRN", relationship = "many-to-one", na_matches = "never")

neutrophils <- read_excel("Clinical Data/neuts.xlsx") %>%
  filter(!is.na(`GOSH MRN`), !is.na(`Neutrophils (x10*9/L)`)) %>%
  group_by(`GOSH MRN`) %>%
  slice_max(`Collection Date (Neutrophils)`, n = 1, with_ties = FALSE) %>%
  ungroup()

patients <- patients %>%
  select(-any_of(c("Neutrophils (x10^9/L)", "Neutrophils (x10*9/L)",
                   "Collection Date (Neutrophils)"))) %>%
  left_join(
    neutrophils %>%
      transmute(`GOSH MRN`, `Collection Date (Neutrophils)`,
                `Neutrophils (x10^9/L)` = `Neutrophils (x10*9/L)`),
    by = "GOSH MRN", relationship = "many-to-one", na_matches = "never")

# ============================================================================
# Genetics data collected from GOSH
# ============================================================================

# Keep unique patients
genetics_new <- read_excel("Clinical Data/mutation_matrix.xlsx") %>%
  filter(!if_all(everything(), is.na)) %>%
  distinct()

patients <- patients %>%
  left_join(genetics_new, by = "GOSH MRN",
            relationship = "many-to-one", na_matches = "never")

# Gene_affected is filled in for everyone who was sequenced, including the
# patients where nothing was found.
with_genetics <- patients %>%
  filter(!is.na(Gene_affected)) %>%
  rename_with(~ make_clean_names(.x), starts_with("c."))

genetics_use <- with_genetics %>%
  select(`GOSH MRN`, `Perforin state`, `CD3 (CD8+CD107A+) (%)`, GRA,
         `Perforin Expression (CD56+ cells)`, Gene_affected, starts_with("c_"))

# ============================================================================
# Combine the two datasets into one table
# ============================================================================

# List which variant belongs to which gene
prf1_v <- c("c.50del", "c.116C>A", "c.386G>C", "c.493G>A", "c.635A>G",
            "c.725G>A", "c.841_843del", "c.916G>T", "c.1018G>A", "c.1040A>T",
            "c.1117C>T", "c.1153C>T", "c.1304C>T")

sh2d1a_v <- c("c.137+2T>C")

stxbp2_v <- c("c.116del", "c.194G>A", "c.1247-1G>C", "c.1621G>A")

xiap_v <- c("c.145C>T", "c.389_392del", "c.446_450del", "c.553G>A", "c.608G>A",
            "c.664C>T", "c.712C>T", "c.802C>T", "c.1037dup", "c.1141C>T",
            "c.1261_1262del", "c.1349G>A", "c.1358G>A")

gene_lookup <- c(
  setNames(rep("PRF1", length(prf1_v)), make_clean_names(prf1_v)),
  setNames(rep("SH2D1A", length(sh2d1a_v)), make_clean_names(sh2d1a_v)),
  setNames(rep("STXBP2", length(stxbp2_v)), make_clean_names(stxbp2_v)),
  setNames(rep("XIAP", length(xiap_v)), make_clean_names(xiap_v))
)

# HGVS labels for the plots
variant_hgvs <- c(
  c_11g_a = "c.11G>A",
  c_46c_t = "c.46C>T",
  c_50del = "c.50del",
  c_268_t = "c.268C>T",
  c_272c_t = "c.272C>T",
  c_386g_c = "c.386G>C",
  c_445g_a = "c.445G>A",
  c_493g_a = "c.493G>A",
  c_666c_a = "c.666C>A",
  c_694c_t = "c.694C>T",
  c_725g_a = "c.725G>A",
  c_731t_g = "c.731T>G",
  c_841_843del = "c.841_843del",
  c_916g_t = "c.916G>T",
  c_1040a_t = "c.1040A>T",
  c_1117c_t = "c.1117C>T",
  c_1122g_a = "c.1122G>A",
  c_1153c_t = "c.1153C>T",
  c_1229_1230delins_cc = "c.1229_1230delinsCC",
  c_1357g_a = "c.1357G>A",
  c_1621g_a = "c.1621G>A",
  c_194g_a = "c.194G>A",
  c_1349g_a = "c.1349G>A",
  c_900c_t = "c.900C>T",
  c_822_c_t = "c.822C>T",
  c_726c_t = "c.726C>T",
  c_434_t_c = "c.434T>C",
  c_539_22g_c = "c.539+22G>C",
  c_539_39g_a = "c.539+39G>A",
  c_539_61g_a = "c.539+61G>A",
  c_539_82c_t = "c.539+82C>T"
)

# Rename variant to match
data_labelled <- data_labelled %>%
  rename(c_50del = c_50del_t)

# Anchor on a digit so that any other column starting with c_ does not get
# picked up as a variant by mistake.
gc_dl <- names(data_labelled) %>%
  str_subset("^c_\\d")

gc_gx <- names(genetics_use) %>% 
  str_subset("^c_\\d")

# Anything in the panel data that is not in the lists above must be PRF1,
# since the first cohort only tested for PRF1
extra <- setdiff(gc_dl, names(gene_lookup))
gene_lookup[extra] <- "PRF1"
message("Assumed PRF1 (confirm): ", paste(extra, collapse = ", "))

geno_cols <- union(gc_dl, gc_gx)
n_dl <- nrow(data_labelled)

# The cleaning left underscores in inconsistent places (c_900c_t vs c_822_c_t).
# If the same variant got spelled in different way then union() would keep both
# and split carriers across two columns
norm_variant <- function(x) tolower(gsub("_", "", x))
dupe_groups <- split(geno_cols, norm_variant(geno_cols))
dupe_groups <- dupe_groups[lengths(dupe_groups) > 1]
if (length(dupe_groups)) {
  stop("Same variant under multiple column names: ",
       paste(sapply(dupe_groups, paste, collapse = " / "), collapse = "; "))
}

# Putting both cohorts into the same column names
table_dl <- data_labelled %>%
  transmute(
    patient_id = row_number(),
    cohort = "data_labelled",
    original_id = as.character(`Patient ID`),
    gene_tested = "PRF1",
    perforin_pct = `Perforin Expression %`,
    perforin_state_recorded = NA_character_,
    gene_affected_recorded = NA_character_,
    result_clean = as.character(result),
    across(all_of(gc_dl))
  )

table_gx <- genetics_use %>%
  transmute(
    patient_id = n_dl + row_number(),
    cohort = "genetics",
    original_id = as.character(`GOSH MRN`),
    gene_tested = "PRF1/SH2D1A/STXBP2/XIAP",
    perforin_pct = `Perforin Expression (CD56+ cells)`,
    perforin_state_recorded = as.character(`Perforin state`),
    gene_affected_recorded = as.character(Gene_affected),
    result_clean = NA_character_,
    across(all_of(gc_gx))
  )

# A variant column missing from one cohort means it was not tested there, so
# fill those with 0 rather than leaving them as NA
combined <- bind_rows(table_dl, table_gx) %>%
  mutate(across(all_of(geno_cols), ~ replace_na(as.numeric(.x), 0)))

# Which genes each patient carries a variant in, of any class
combined$gene_derived <- apply(as.matrix(combined[geno_cols]), 1, function(r) {
  hit <- geno_cols[r > 0]
  if (!length(hit)) "nm" else paste(sort(unique(gene_lookup[hit])), collapse = "; ")
})

# Perforin state uses the lab thresholds: over 50 Normal, 10 to 50 Abnormal,
# under 10 Absent
combined <- combined %>%
  mutate(
    perforin_state_calc = case_when(
      is.na(perforin_pct) ~ NA_character_,
      perforin_pct < 10 ~ "Absent",
      perforin_pct <= 50 ~ "Abnormal",
      TRUE ~ "Normal"),
    perforin_state = perforin_state_calc,
    analysis_group = case_when(
      cohort == "genetics" & gene_affected_recorded == "nm" ~ "No mutation",
      cohort == "genetics" ~ "Mutation",
      is.na(result_clean) ~ NA_character_,
      str_detect(result_clean, regex("mutation", ignore_case = TRUE)) &
        !str_detect(result_clean, regex("^no mutation", ignore_case = TRUE)) ~ "Mutation",
      str_detect(result_clean, regex("sequence var", ignore_case = TRUE)) ~ "Sequence variance",
      str_detect(result_clean, regex("polymorphism", ignore_case = TRUE)) ~ "Polymorphism",
      str_detect(result_clean, regex("^no mutation", ignore_case = TRUE)) ~ "No mutation",
      TRUE ~ NA_character_),
    # gene_affected_recorded is the referring lab's own column. I am reading it
    # as "the gene with a pathogenic finding", but that is my interpretation of
    # their field, so I check it in the QC section below.
    gene_affected = case_when(
      gene_affected_recorded == "PRF1" ~ "PRF1 (Pathogenic)",
      gene_affected_recorded == "nm" ~ "No variant detected",
      !is.na(gene_affected_recorded) ~ gene_affected_recorded,
      is.na(analysis_group) ~ NA_character_,
      analysis_group == "Mutation" ~ "PRF1 (Pathogenic)",
      rowSums(across(all_of(geno_cols))) > 0 ~ "PRF1 (Non-pathogenic)",
      TRUE ~ "No PRF1 variant detected"))

# The two cohorts give different kinds of negative control. The genetics
# cohort was sequenced across all four genes, the panel cohort only for PRF1
combined <- combined %>%
  mutate(control_type = case_when(
    cohort == "genetics" & gene_affected_recorded == "nm" ~ "Full negative",
    cohort == "data_labelled" & analysis_group == "No mutation" ~ "PRF1-panel negative",
    TRUE ~ NA_character_))

# ============================================================================
# Variant classification
# ============================================================================

# Taken from VarSome and ClinVar, accessed 20 August 2026. For anything with
# no ClinVar entry or with conflicting submissions are recorded VUS, and those
# collapse into "Uncertain" for the figure. Note in the caption that
# "no evidence yet" and "the evidence disagrees" are not the same
variant_class_raw <- c(
  c_11g_a = "Likely benign",
  c_46c_t = "VUS",
  c_50del = "Pathogenic",
  c_268_t = "Likely benign",
  c_272c_t = "Hypomorphic",
  c_386g_c = "Pathogenic",
  c_434_t_c = "Benign",
  c_445g_a = "Pathogenic",
  c_493g_a = "VUS",
  c_539_22g_c = "Benign",
  c_539_39g_a = "VUS",
  c_539_61g_a = "Benign",
  c_539_82c_t = "Benign",
  c_666c_a = "Pathogenic",
  c_694c_t = "Pathogenic",
  c_725g_a = "VUS",
  c_726c_t = "Likely benign",
  c_731t_g = "VUS",
  c_822_c_t = "Benign",
  c_841_843del = "VUS",
  c_900c_t = "Benign",
  c_916g_t = "Pathogenic",
  c_1040a_t = "VUS",
  c_1117c_t = "VUS",
  c_1122g_a = "Pathogenic",
  c_1153c_t = "Benign",
  c_1229_1230delins_cc = "Likely pathogenic",
  c_1357g_a = "Likely benign",
  c_194g_a = "Pathogenic",
  c_1621g_a = "Benign",
  c_1349g_a = "VUS"
)

# Pathogenic variants at the top, benign at the bottom
plot_levels <- c("Pathogenic/LP", "Hypomorphic", "Uncertain", "Benign/LB")

collapse_class <- function(x) {
  x <- str_squish(str_to_lower(x))
  case_when(
    x %in% c("pathogenic", "likely pathogenic") ~ "Pathogenic/LP",
    x == "hypomorphic" ~ "Hypomorphic",
    x %in% c("benign", "likely benign") ~ "Benign/LB",
    TRUE ~ "Uncertain"
  )
}

# ============================================================================
# Zygosity, PRF1 only
# ============================================================================

# These are the variants not counted when deciding whether someone has one
# or two pathogenic PRF1 alleles.
poly_cols <- c(
  # A91V is a hypomorphic risk allele
  "c_272c_t",
  
  # Benign, intronic or synonymous according to ClinVar and gnomAD.
  "c_900c_t", "c_434_t_c",
  "c_539_82c_t", "c_539_61g_a", "c_539_39g_a", "c_539_22g_c",
  "c_822_c_t", "c_726c_t", "c_11g_a", "c_268_t", "c_1357g_a",
  
  # ClinVar calls this Benign. Reclassifying it moved patient 101 out of the
  # biallelic group. That patient was homozygous with 56.1% perforin and the
  # referring lab had recorded them as PRF1 affected. Main comparison went
  # from p = 7.13e-07 to p = 5.49e-07
  "c_1153c_t",
  
  # VarSome calls this Likely pathogenic but the evidence is thin: phyloP is
  # 0.529, so no conservation at that position, and only 2 predictors call it
  # damaging. I exclude it on that basis. I still plot it as Pathogenic/LP so
  # the figure shows the database call rather than my own judgement.
  # Sensitivity analysis (c) tests including it
  "c_1229_1230delins_cc"
)

prf1_cols <- names(gene_lookup)[gene_lookup == "PRF1"]

# Named these prf1_ because they only count PRF1 alleles. Someone with a
# STXBP2 or XIAP variant will show 0 here
prf1_path_cols <- setdiff(intersect(geno_cols, prf1_cols), poly_cols)

combined <- combined %>%
  mutate(prf1_path_alleles = rowSums(across(all_of(prf1_path_cols))),
         prf1_zygosity = case_when(
           prf1_path_alleles == 0 ~ "none",
           prf1_path_alleles == 1 ~ "monoallelic",
           TRUE ~ "biallelic"))

# The gene check comes first so that STXBP2, XIAP and SH2D1A patients are not
# sent down the PRF1 zygosity branch
combined <- combined %>%
  mutate(group = case_when(
    gene_affected %in% c("STXBP2", "XIAP", "SH2D1A") ~ gene_affected,
    analysis_group == "Mutation" & prf1_zygosity == "biallelic" ~ "PRF1 biallelic",
    analysis_group == "Mutation" & prf1_zygosity == "monoallelic" ~ "PRF1 monoallelic",
    analysis_group == "Mutation" & prf1_zygosity == "none" ~ "Mutation (VUS only)",
    TRUE ~ analysis_group))

# ============================================================================
# Quality control
# ============================================================================

# Is perforin expression on the same scale in both cohorts?
combined %>%
  group_by(cohort) %>%
  summarise(n = n(),
            min = min(perforin_pct, na.rm = TRUE),
            max = max(perforin_pct, na.rm = TRUE),
            n_missing = sum(is.na(perforin_pct)))

# Do my thresholds reproduce the lab's own state calls? Returned 67/67
combined %>%
  filter(!is.na(perforin_state_recorded)) %>%
  count(perforin_state_recorded, perforin_state_calc)

# Does the gene I derive agree with the gene the lab recorded?
combined %>%
  filter(!is.na(gene_affected_recorded),
         gene_affected_recorded != gene_derived) %>%
  select(patient_id, original_id, gene_affected_recorded, gene_derived)

# The check above misses the case where both fields say PRF1 but the only
# variant present is benign, so here ask whether each patient actually
# carries a pathogenic variant in the gene the lab named. Returned patient 101
combined %>%
  filter(gene_affected_recorded %in% c("PRF1", "STXBP2", "XIAP", "SH2D1A")) %>%
  select(patient_id, gene_affected_recorded, perforin_pct, perforin_state,
         all_of(setdiff(geno_cols, poly_cols))) %>%
  pivot_longer(all_of(setdiff(geno_cols, poly_cols)),
               names_to = "variant", values_to = "alleles") %>%
  mutate(variant_gene = gene_lookup[variant]) %>%
  group_by(patient_id, gene_affected_recorded, perforin_pct, perforin_state) %>%
  summarise(own_gene_alleles = sum(alleles[variant_gene == gene_affected_recorded]),
            .groups = "drop") %>%
  filter(own_gene_alleles == 0)

# Does the result text agree with whether any variant was actually recorded? Yes
combined %>%
  filter(cohort == "data_labelled") %>%
  mutate(has_variant = rowSums(across(all_of(geno_cols))) > 0) %>%
  count(result_clean, has_variant)

# Is anyone falling outside every group? No
combined %>%
  filter(is.na(group)) %>%
  count(cohort, result_clean)

# Group sizes
combined %>%
  count(cohort, analysis_group)

combined %>%
  count(cohort, gene_affected)

combined %>%
  count(group)

# ============================================================================
# Final table
# ============================================================================

final_table <- combined %>%
  select(patient_id, cohort, original_id, gene_tested,
         perforin_pct, perforin_state,
         gene_affected, analysis_group, result_clean,
         prf1_zygosity, prf1_path_alleles, group, control_type,
         all_of(geno_cols)) %>%
  arrange(patient_id)

# original_id holds GOSH MRNs, so drop it before further analyses
final_table_clean <- final_table %>%
  select(-original_id)

# ============================================================================
# Sensitivity analyses
# ============================================================================

# Three classification calls could move patients between the main comparison
# groups, so check each one and report both versions. The *_cols objects are
# defined here because the Stats section below uses them.
zyg_from <- function(df, cols) {
  n <- rowSums(df[, cols, drop = FALSE])
  case_when(n == 0 ~ "none", n == 1 ~ "monoallelic", TRUE ~ "biallelic")
}

uncertain_vars <- names(variant_class_raw)[
  collapse_class(variant_class_raw) == "Uncertain"]

# (a) Uncertain variants count toward zygosity. What if they do not?
# 4 of the 13 biallelic patients depend on a VUS as their second allele, so
# the strict biallelic group is 9.
# TODO record where those 4 land from the count() below. Either all 4 become
# monoallelic, or 2 become monoallelic and 2 become none. Both give n = 9, so
# the count is the only thing that separates them, and the earlier comments
# in this script disagreed on it.
strict_cols <- setdiff(prf1_path_cols, uncertain_vars)
combined %>%
  mutate(strict = zyg_from(., strict_cols)) %>%
  count(prf1_zygosity, strict) %>%
  print(n = Inf)

# All of the uncertain variants are heterozygous only, so none of them can
# create a homozygous call
combined %>%
  select(patient_id, any_of(uncertain_vars)) %>%
  pivot_longer(-patient_id, names_to = "variant", values_to = "alleles") %>%
  filter(alleles > 0) %>%
  count(variant, alleles)

# The strict comparison itself is run in the Stats section below, p = 3.319e-06

# (b) What if A91V counted as a pathogenic allele?
incl_a91v_cols <- union(prf1_path_cols, "c_272c_t")

combined %>%
  mutate(with_a91v = zyg_from(., incl_a91v_cols)) %>%
  count(prf1_zygosity, with_a91v) %>%
  print(n = Inf)

final_table_clean %>%
  mutate(z_a91v = zyg_from(., incl_a91v_cols)) %>%
  filter(group == "No mutation" |
           (analysis_group == "Mutation" & z_a91v == "biallelic")) %>%
  mutate(cmp = if_else(group == "No mutation", "No mutation", "PRF1 biallelic")) %>%
  wilcox.test(perforin_pct ~ cmp, data = ., exact = FALSE)
# p = 2.978e-07. Slightly stronger than the main result, conclusion unchanged

# (c) What if c.1229_1230delinsCC counted as pathogenic?
incl_delins_cols <- union(prf1_path_cols, "c_1229_1230delins_cc")

combined %>%
  mutate(with_delins = zyg_from(., incl_delins_cols)) %>%
  count(prf1_zygosity, with_delins) %>%
  print(n = Inf)

final_table_clean %>%
  mutate(z_delins = zyg_from(., incl_delins_cols)) %>%
  filter(group == "No mutation" |
           (analysis_group == "Mutation" & z_delins == "biallelic")) %>%
  mutate(cmp = if_else(group == "No mutation", "No mutation", "PRF1 biallelic")) %>%
  wilcox.test(perforin_pct ~ cmp, data = ., exact = FALSE)
# p = 2.45e-07 when run without exact = FALSE. Reconfirm under the normal
# approximation and record here. Conclusion unchanged either way

# ============================================================================
# Stats
# ============================================================================

# exact = FALSE everywhere. Perforin has heavy ties at 0, so wilcox.test was
# already falling back to the normal approximation and only warning about it.
# Setting it explicitly makes the methods section match what was computed.
# The p-values below were recorded before this was applied uniformly, so
# reconfirm each one and update both the comment and the figure, which
# hardcodes 0.000399, 5.491e-07, 0.14572 and 0.001624.

# Biallelic versus monoallelic: does one pathogenic allele reduce perforin?
final_table_clean %>%
  filter(group %in% c("PRF1 biallelic", "PRF1 monoallelic")) %>%
  wilcox.test(perforin_pct ~ group, data = ., exact = FALSE)
# p = 0.000399

# Does perforin differ between biallelic PRF1 defects and no mutation?
final_table_clean %>%
  filter(group %in% c("PRF1 biallelic", "No mutation")) %>%
  wilcox.test(perforin_pct ~ group, data = ., exact = FALSE)
# p = 5.491e-07, was 7.13e-07 before I reclassified c.1153C>T.
# Biallelic PRF1 defects abolish perforin expression

# Same comparison but only the 9 patients with two confidently pathogenic
# alleles, so the result does not rest on any VUS. This is sensitivity (a)
final_table_clean %>%
  filter(group %in% c("PRF1 biallelic", "No mutation")) %>%
  mutate(strict = zyg_from(., strict_cols)) %>%
  filter(group == "No mutation" | strict == "biallelic") %>%
  wilcox.test(perforin_pct ~ group, data = ., exact = FALSE)
# p = 3.319e-06. Same order of magnitude, so the finding is consistent

# Monoallelic versus no mutation
final_table_clean %>%
  filter(group %in% c("PRF1 monoallelic", "No mutation")) %>%
  wilcox.test(perforin_pct ~ group, data = ., exact = FALSE)
# p = 0.14572
# n = 12, so this is underpowered for a modest heterozygote shift. Report as
# no significant difference, not as evidence of no difference

# Polymorphism versus no mutation
final_table_clean %>%
  filter(group %in% c("Polymorphism", "No mutation")) %>%
  wilcox.test(perforin_pct ~ group, data = ., exact = FALSE)
# p = 1.222e-08
# Confounded by cohort. All the polymorphism carriers come from the panel
# data and 48 of 54 no-mutation patients come from the genetics data, so
# redo it within one cohort below

final_table_clean %>%
  filter(cohort == "data_labelled", group %in% c("Polymorphism", "No mutation")) %>%
  wilcox.test(perforin_pct ~ group, data = ., exact = FALSE)
# p = 0.980
# Within cohort there is no difference, so the pooled result was just cohort
# effects

# Is the main finding a result of cohort effects?
final_table_clean %>%
  filter(cohort == "data_labelled", group %in% c("PRF1 biallelic", "No mutation")) %>%
  wilcox.test(perforin_pct ~ group, data = ., exact = FALSE)
# p = 0.00600, no

final_table_clean %>%
  filter(cohort == "genetics", group %in% c("PRF1 biallelic", "No mutation")) %>%
  wilcox.test(perforin_pct ~ group, data = ., exact = FALSE)
# p = 0.001325
# It replicates separately in both cohorts, so it is not a cohort effect

# Biallelic PRF1 versus patients negative across all four tested genes
final_table_clean %>%
  filter(group == "PRF1 biallelic" | control_type == "Full negative") %>%
  wilcox.test(perforin_pct ~ group, data = ., exact = FALSE)
# p = 8.39e-10

# ----------------------------------------------------------------------------
# Panel B: are the two mutation-negative groups comparable controls?
# ----------------------------------------------------------------------------

# The panel cohort was referred for perforin flow cytometry on clinical
# suspicion and sequenced for PRF1 only if that came back low, so it is
# selected on the outcome variable. This test quantifies the ascertainment
# effect, it does not discover it. It is also its own family, so it is not
# Holm-adjusted with the panel A comparisons
final_table_clean %>%
  filter(group == "No mutation") %>%
  wilcox.test(perforin_pct ~ cohort, data = ., exact = FALSE, conf.int = TRUE)
# W = 29, p = 0.001624
# Hodges-Lehmann shift -45.8%, 95% CI -59.8 to -25.0. cohort is still
# data_labelled/genetics here, so negative means the panel cohort sits lower.
# CI is 35 points wide at n = 6, so quote the shift as approximate

# Medians and IQRs. The manuscript quotes these, and n = 6 makes the
# panel-cohort box unstable
final_table_clean %>%
  filter(group == "No mutation") %>%
  group_by(cohort) %>%
  summarise(n = n(),
            median = median(perforin_pct),
            q1 = quantile(perforin_pct, 0.25),
            q3 = quantile(perforin_pct, 0.75),
            n_not_normal = sum(perforin_state != "Normal"),
            .groups = "drop")

# Same point as a proportion below the normal threshold. Fisher handles the
# small cells where a chi-square would not
final_table_clean %>%
  filter(group == "No mutation") %>%
  mutate(not_normal = perforin_state != "Normal") %>%
  with(fisher.test(table(cohort, not_normal)))

# Three perforin values were imputed by state average, which pulls toward the
# centre of a group defined by the variable being tested. At n = 6 one such
# value moves the median, so check whether any landed in the panel-cohort
# negatives. Uses final_table because original_id is dropped in _clean
imputed_ids <- c("91", "103", "21")

final_table %>%
  filter(group == "No mutation", cohort == "data_labelled") %>%
  mutate(imputed = original_id %in% imputed_ids) %>%
  select(patient_id, original_id, perforin_pct, imputed)

# Rerun without them. Exclusion is scoped to the panel cohort so a GOSH MRN
# of 21, 91 or 103 in the genetics cohort cannot be dropped by accident
final_table %>%
  filter(group == "No mutation",
         !(cohort == "data_labelled" & original_id %in% imputed_ids)) %>%
  wilcox.test(perforin_pct ~ cohort, data = ., exact = FALSE, conf.int = TRUE)
# TODO record p and shift. If unchanged, one line in the text. If it moves,
# panel B needs a caveat in the caption

# ============================================================================
# Plots
# ============================================================================

group_levels <- c("PRF1 biallelic",
                  "PRF1 monoallelic",
                  "Sequence variance",
                  "Polymorphism",
                  "No mutation",
                  "Mutation (VUS only)",
                  "STXBP2",
                  "XIAP",
                  "SH2D1A")

perforin_levels <- c("Absent", "Abnormal", "Normal")

# Perforin by group, coloured by cohort so to see any cohort separation.
# Dashed lines are the state thresholds
final_table_clean %>%
  mutate(group = factor(group, levels = group_levels)) %>%
  ggplot(aes(group, perforin_pct)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour = cohort), width = 0.15, alpha = 0.6) +
  geom_hline(yintercept = c(10, 50), linetype = "dashed", alpha = 0.4) +
  labs(x = NULL, y = "Perforin expression (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Stacked bar chart of patients of each perforin state by condition
final_table_clean %>%
  count(group, perforin_state) %>%
  mutate(group = factor(group, levels = group_levels)) %>%
  ggplot(aes(group, n, fill = factor(perforin_state, levels = perforin_levels))) +
  geom_col(position = "fill") +
  labs(x = NULL, y = "Proportion", fill = "Perforin state") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Mutation negative patients split by cohort, to show the ascertainment issue
final_table_clean %>%
  filter(group == "No mutation") %>%
  ggplot(aes(cohort, perforin_pct)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  labs(title = "Mutation-negative patients by cohort", y = "Perforin (%)") +
  theme_minimal()

# Which variants are carried most often?
final_table_clean %>%
  summarise(across(all_of(geno_cols), ~ sum(.x > 0))) %>%
  pivot_longer(everything(), names_to = "variant", values_to = "n") %>%
  filter(n > 0) %>%
  mutate(gene = gene_lookup[variant],
         label = coalesce(variant_hgvs[variant], variant)) %>%
  ggplot(aes(reorder(label, n), n, fill = gene)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Patients carrying variant") +
  theme_minimal()

# Pathogenic PRF1 variants split by zygosity
final_table_clean %>%
  select(patient_id, all_of(prf1_path_cols)) %>%
  pivot_longer(-patient_id, names_to = "variant", values_to = "alleles") %>%
  filter(alleles > 0) %>%
  count(variant, alleles) %>%
  mutate(state = if_else(alleles == 1, "Heterozygous", "Homozygous"),
         label = coalesce(variant_hgvs[variant], variant)) %>%
  ggplot(aes(reorder(label, n), n, fill = state)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Patients carrying variant", fill = NULL) +
  theme_minimal()

# How many of the biallelic patients are homozygous versus compound heterozygous?
# 4 are compound heterozygous, 9 are homozygous
final_table_clean %>%
  select(patient_id, all_of(prf1_path_cols)) %>%
  pivot_longer(-patient_id, names_to = "variant", values_to = "alleles") %>%
  filter(alleles > 0) %>%
  group_by(patient_id) %>%
  summarise(n_variants = n(), total_alleles = sum(alleles), .groups = "drop") %>%
  filter(total_alleles >= 2) %>%
  count(genotype = if_else(n_variants == 1, "Homozygous", "Compound het"))

# Perforin by state, just a sanity check to show how expression differs by state
final_table_clean %>%
  mutate(perforin_state = factor(perforin_state, levels = perforin_levels)) %>%
  ggplot(aes(perforin_state, perforin_pct, fill = perforin_state)) +
  geom_hline(yintercept = c(10, 50), linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_boxplot(width = 0.5, outlier.shape = NA, colour = "black") +
  geom_jitter(width = 0.15, size = 1, alpha = 0.8, show.legend = FALSE) +
  stat_summary(
    fun.data = function(x) data.frame(y = -5, label = paste0("n = ", length(x))),
    geom = "text", size = 3) +
  labs(y = "Perforin expression (%)", x = "Perforin state") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)))

# ============================================================================
# Oncoprint
# ============================================================================

variant_cols <- final_table_clean %>%
  select(matches("^c_\\d")) %>%
  colnames()

# Subset to carriers once and build both the matrix and the annotations from
# that same subset
carrier_ids <- final_table_clean$patient_id[
  rowSums(final_table_clean[, variant_cols] > 0, na.rm = TRUE) > 0]

ftc_plot <- final_table_clean %>%
  filter(patient_id %in% carrier_ids)

# The panel widths show the perforin distribution among carriers, not the
# whole cohort
m <- t(as.matrix(ftc_plot[, variant_cols]))
colnames(m) <- ftc_plot$patient_id

# Zero the NAs to prevent rowSums from inserting all NA rows
m[is.na(m)] <- 0
m <- m[rowSums(m > 0) > 0, , drop = FALSE]
kept <- rownames(m)

vclass <- droplevels(factor(collapse_class(unname(variant_class_raw[kept])),
                            levels = plot_levels))

gene <- factor(unname(gene_lookup[kept]))

# Ordering by cDNA position
cds_pos <- as.numeric(sub("^c[._]?(\\d+).*$", "\\1", kept))
ord <- order(gene, cds_pos)

m <- m[ord, , drop = FALSE]
kept <- kept[ord]
gene <- gene[ord]
vclass <- vclass[ord]
rownames(m) <- coalesce(variant_hgvs[kept], kept)

# oncoPrint calculates percentages over the carriers only, so calculate for
# the whole cohort
pct_full <- sprintf("%.0f%%", 100 * rowSums(m > 0) / nrow(final_table_clean))

mat_chr <- matrix(c("", "HET", "HOM")[m + 1],
                  nrow = nrow(m), dimnames = dimnames(m))

col_alt <- c(HET = "#7FB3D5", HOM = "#B2182B")

box <- function(fill) {
  function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.4, "mm"), h - unit(0.4, "mm"),
              gp = gpar(fill = fill, col = NA))
  }
}

alter_fun <- list(
  background = box("grey92"),
  HET = box(col_alt["HET"]),
  HOM = box(col_alt["HOM"])
)

group_cols <- c(
  "PRF1 biallelic" = "#1B7837",
  "PRF1 monoallelic" = "#A6DBA0",
  "STXBP2" = "#762A83",
  "XIAP" = "#C2A5CF",
  "SH2D1A" = "#E7D4E8",
  "Mutation (VUS only)" = "#2166AC",
  "Sequence variance" = "#92C5DE",
  "Polymorphism" = "#F4A582",
  "No mutation" = "grey80"
)

perforin_cols <- c(
  "Absent" = "#d95f02",
  "Abnormal" = "#7570b3",
  "Normal" = "#1b9e77"
)

gene_cols <- c(PRF1 = "#E7298A",
               STXBP2 = "#66A61E",
               XIAP = "#7570B3",
               SH2D1A = "#E6AB02")

top_ann <- HeatmapAnnotation(
  Group = droplevels(factor(ftc_plot$group, levels = group_levels)),
  Perforin_state = droplevels(factor(ftc_plot$perforin_state,
                                     levels = perforin_levels)),
  col = list(Group = group_cols, Perforin_state = perforin_cols),
  annotation_name_gp = gpar(fontsize = 8),
  simple_anno_size = unit(3.5, "mm"),
  na_col = "white"
)

# Split columns by perforin state to link genotype to phenotype
# Rows split by variant class
ht <- oncoPrint(
  mat_chr,
  alter_fun = alter_fun,
  col = col_alt,
  row_order = 1:nrow(mat_chr),
  column_split = droplevels(factor(ftc_plot$perforin_state,
                                   levels = perforin_levels)),
  row_split = vclass,
  top_annotation = top_ann,
  show_pct = FALSE,
  left_annotation = rowAnnotation(
    Gene = gene,
    pct = anno_text(pct_full, just = "right",
                    location = unit(1, "npc"), gp = gpar(fontsize = 7)),
    col = list(Gene = gene_cols),
    annotation_name_gp = gpar(fontsize = 8)
  ),
  right_annotation = rowAnnotation(
    "Carriers" = anno_oncoprint_barplot(border = FALSE)
  ),
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 7),
  row_title_gp = gpar(fontsize = 8, fontface = "bold"),
  row_title_rot = 0,
  heatmap_legend_param = list(
    title = "Zygosity",
    at = c("HET", "HOM"),
    labels = c("Heterozygous", "Homozygous"))
)

pdf("oncoprint.pdf", width = 11, height = 10)
draw(ht, merge_legends = TRUE, heatmap_legend_side = "right")
dev.off()

# ============================================================================
# Final figure
# ============================================================================

# One threshold definition so both panels use identical cutoffs
star_label <- function(p) {
  case_when(p < 1e-4 ~ "****", p < 1e-3 ~ "***",
            p < 1e-2 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
}

n_lab <- function(x) data.frame(y = -6, label = paste0("n = ", length(x)))

# Match group names to the x positions ggplot actually draws, then stack
# brackets shortest-span lowest so they never cross
brackets <- function(cmp, x_levels, base, step = 9, adjust = "holm") {
  cmp %>%
    mutate(p.adj = p.adjust(p, method = adjust),
           label = star_label(p.adj),
           xmin = match(group1, x_levels),
           xmax = match(group2, x_levels)) %>%
    arrange(abs(xmax - xmin)) %>%
    mutate(y.position = base + step * (row_number() - 1))
}

stars <- function(df) {
  list(stat_pvalue_manual(df, label = "label", tip.length = 0.01,
                          bracket.size = 0.3, size = 3.2),
       expand_limits(y = max(df$y.position) + 4))
}

cohort_labs <- c(data_labelled = "PRF1 panel", genetics = "Four-gene cohort")

plot_df <- final_table_clean %>%
  mutate(group = factor(group, levels = group_levels),
         cohort = recode(cohort, !!!cohort_labs))

# Empty factor levels are dropped by the discrete scale, so read bracket
# positions off what is actually plotted. No SH2D1A carriers were identified
group_x <- levels(droplevels(plot_df$group))
cohort_x <- sort(unique(plot_df$cohort))

# Holm across the three comparisons shown in panel A
stat_group <- tibble::tribble(
  ~group1, ~group2, ~p,
  "PRF1 biallelic", "PRF1 monoallelic", 0.000399,
  "No mutation", "PRF1 biallelic", 5.491e-07,
  "No mutation", "PRF1 monoallelic", 0.14572
) %>%
  brackets(group_x, base = 94)

# Panel B is a separate question, so its own family and unadjusted
stat_neg <- tibble::tibble(
  group1 = "Four-gene cohort", group2 = "PRF1 panel", p = 0.001624
) %>%
  brackets(cohort_x, base = 108, adjust = "none")

# Shared layers. Threshold lines come after the points so they draw on top
common <- list(
  geom_hline(yintercept = c(10, 50), linetype = "dashed", alpha = 0.4),
  stat_summary(fun.data = n_lab, geom = "text", size = 2.8),
  labs(x = NULL, y = "Perforin expression (%)"),
  theme_minimal(base_size = 11)
)

# Lock the jitter so the figure regenerates identically
set.seed(1)

plot_group <- plot_df %>%
  add_count(group) %>%
  ggplot(aes(group, perforin_pct)) +
  geom_boxplot(data = ~ filter(.x, n >= 5), outlier.shape = NA) +
  geom_jitter(aes(colour = cohort), width = 0.15, alpha = 0.7, size = 1.6) +
  common +
  stars(stat_group) +
  labs(colour = "Cohort") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

plot_neg <- plot_df %>%
  filter(group == "No mutation") %>%
  ggplot(aes(cohort, perforin_pct)) +
  geom_boxplot(outlier.shape = NA, width = 0.45) +
  geom_jitter(width = 0.12, alpha = 0.7, size = 1.6) +
  common +
  stars(stat_neg)

clinical_figure <- plot_group / plot_neg +
  plot_layout(heights = c(1.8, 1.2)) +
  plot_annotation(tag_levels = "A")

ggsave("figure_perforin_by_group.png", clinical_figure, width = 8, height = 9)