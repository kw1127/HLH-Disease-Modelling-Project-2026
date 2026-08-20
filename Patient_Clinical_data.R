# ============================================================================
# GOSH cohort data cleaning for analysis
# ============================================================================

# Load libraries
library(tidyverse)
library(janitor)
library(readxl)
library(openxlsx)

# ============================================================================
# Original dataset: PRF1 panel data
# ============================================================================

# Load original patient dataset
data <- read_excel("Patient Data /HLH_data.xlsx", sheet = "TOTAL")

# Only keep patients
# Convert perforin expression to numerical factors
data <- data %>%
  slice(1:112) %>%
  mutate(Perforin_expression_percent = as.numeric(`Perforin expression %`)) %>%
  mutate(sex = as.factor(Sex)) %>%
  clean_names()

# Remove unimportant columns
# Convert all variants to numerical factors
# Filter for patients with variants recorded
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
  mutate(across(starts_with("c_"), as.numeric)) %>%
  filter(if_all(starts_with("c_"), ~ !is.na(.)))

# Impute missing perforin values by perforin state, then convert to percentage
data_clean <- data_clean %>%
  mutate(`Perforin Expression %` = case_when(
    `Patient ID` == 91 ~ 0.56,
    `Patient ID` == 103 ~ 0.52,
    `Patient ID` == 21 ~ 0.19,
    TRUE ~ `Perforin Expression %`)) %>%
  mutate(`Perforin Expression %` = `Perforin Expression %` * 100)

# Clean up clinician-derived diagnoses
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

# Load patient and clinical data collected from GOSH
patients <- read_excel("~/combined_patients.xlsx", sheet = "Original") %>%
  filter(!is.na(`GOSH MRN`)) %>%
  select(-starts_with("..."))

# Load sCD25 data
sCD25 <- read_excel("Clinical Data/sCD25.xlsx") %>%
  mutate(`sCD25 (pg/ml)` = as.numeric(gsub("[<>]", "", Value))) %>%
  filter(!is.na(`GOSH MRN`), !is.na(Value)) %>%
  group_by(`GOSH MRN`) %>%
  slice_max(`Collection Date (Best Available)`, n = 1, with_ties = FALSE) %>%
  ungroup()

# Join
patients <- patients %>%
  select(-any_of(c("Collection Date (sCD25)", "sCD25", "sCD25 (pg/ml)"))) %>%
  left_join(
    sCD25 %>%
      transmute(`GOSH MRN`,
                `Collection Date (sCD25)` = `Collection Date (Best Available)`,
                sCD25 = Value,
                `sCD25 (pg/ml)`),
    by = "GOSH MRN", relationship = "many-to-one", na_matches = "never")

# Load trigylceride data
triglycerides <- read_excel("Clinical Data/Tri.xlsx") %>%
  mutate(`triglycerides (mmol/g)` =
           as.numeric(gsub("[<>]", "", `Triglycerides (mmol/g)`))) %>%
  filter(!is.na(`GOSH MRN`), !is.na(`Triglycerides (mmol/g)`)) %>%
  group_by(`GOSH MRN`) %>%
  slice_max(`Collection Date (Best Available)`, n = 1, with_ties = FALSE) %>%
  ungroup()

# Join
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

# Load fibrinogen data
fibrinogen <- read_excel("Clinical Data/fibr.xlsx") %>%
  mutate(`fibrinogen (g/L)` = as.numeric(gsub("[<>]", "", `Fibrinogen (g/L)`))) %>%
  filter(!is.na(`GOSH MRN`), !is.na(`Fibrinogen (g/L)`)) %>%
  group_by(`GOSH MRN`) %>%
  slice_max(`Collection Date (Fibrinogen)`, n = 1, with_ties = FALSE) %>%
  ungroup()

# Join
patients <- patients %>%
  select(-any_of(c("Collection Date (Fibrinogen)", "Fibrinogen (g/L)",
                   "fibrinogen (g/L)"))) %>%
  left_join(
    fibrinogen %>%
      transmute(`GOSH MRN`, `Collection Date (Fibrinogen)`,
                `Fibrinogen (g/L)`, `fibrinogen (g/L)`),
    by = "GOSH MRN", relationship = "many-to-one", na_matches = "never")

# Load ferritin data
ferritin <- read_excel("Clinical Data/ferritin.xlsx") %>%
  mutate(`ferritin (ug/L)` = as.numeric(gsub("[<>]", "", `Ferritin (ug/L)`))) %>%
  filter(!is.na(`GOSH MRN`), !is.na(`Ferritin (ug/L)`)) %>%
  group_by(`GOSH MRN`) %>%
  slice_max(`Collection Date (Ferritin)`, n = 1, with_ties = FALSE) %>%
  ungroup()

# Join
patients <- patients %>%
  select(-any_of(c("Collection Date (Ferritin)", "Ferritin (ug/L)",
                   "ferritin (ug/L)"))) %>%
  left_join(
    ferritin %>%
      transmute(`GOSH MRN`, `Collection Date (Ferritin)`,
                `Ferritin (ug/L)`, `ferritin (ug/L)`),
    by = "GOSH MRN", relationship = "many-to-one", na_matches = "never")

# Load haemoglobin data
haemoglobin <- read_excel("Clinical Data/Haemoglobin.xlsx") %>%
  filter(!is.na(`GOSH MRN`), !is.na(`Haemoglobin (g/L)`)) %>%
  group_by(`GOSH MRN`) %>%
  slice_max(`Collection Date (Haemoglobin)`, n = 1, with_ties = FALSE) %>%
  ungroup()

# Join
patients <- patients %>%
  select(-any_of(c("Haemoglobin (g/L)", "Collection Date (Haemoglobin)"))) %>%
  left_join(
    haemoglobin %>%
      transmute(`GOSH MRN`, `Collection Date (Haemoglobin)`, `Haemoglobin (g/L)`),
    by = "GOSH MRN", relationship = "many-to-one", na_matches = "never")

# Load neutrophils data
neutrophils <- read_excel("~/Downloads/(No subject) (2) 2/neuts.xlsx") %>%
  filter(!is.na(`GOSH MRN`), !is.na(`Neutrophils (x10*9/L)`)) %>%
  group_by(`GOSH MRN`) %>%
  slice_max(`Collection Date (Neutrophils)`, n = 1, with_ties = FALSE) %>%
  ungroup()

# Join
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

# Load the genetics data
# Exclude dublicate patients
genetics_new <- read_excel("Clinical Data/mutation_matrix.xlsx") %>%
  filter(!if_all(everything(), is.na)) %>%
  distinct()

# Join to the rest of the patients by GOSH MRNs
# This keeps only the patients who are included in the
patients <- patients %>%
  left_join(genetics_new, by = "GOSH MRN",
            relationship = "many-to-one", na_matches = "never")

# Filter for patients with gene_affected filled in including no mutants
with_genetics <- patients %>%
  filter(!is.na(Gene_affected)) %>%
  rename_with(make_clean_names, starts_with("c."))

# Keep relevant columns
genetics_use <- with_genetics %>%
  select(`GOSH MRN`, `Perforin state`, `CD3 (CD8+CD107A+) (%)`, GRA,
         `Perforin Expression (CD56+ cells)`, Gene_affected, starts_with("c_"))

# ============================================================================
# Combine with original dataset into one table
# ============================================================================

# Classify each variant by affected gene
prf1_v <- c("c.50del","c.116C>A","c.386G>C","c.493G>A","c.635A>G","c.725G>A",
            "c.841_843del","c.916G>T","c.1018G>A","c.1040A>T","c.1117C>T",
            "c.1153C>T","c.1304C>T")
sh2d1a_v <- c("c.137+2T>C")
stxbp2_v <- c("c.116del","c.194G>A","c.1247-1G>C","c.1621G>A")
xiap_v <- c("c.145C>T","c.389_392del","c.446_450del","c.553G>A","c.608G>A",
            "c.664C>T","c.712C>T","c.802C>T","c.1037dup","c.1141C>T",
            "c.1261_1262del","c.1349G>A","c.1358G>A")

# Create table for each variant and affected gene
gene_lookup <- c(
  setNames(rep("PRF1", length(prf1_v)), make_clean_names(prf1_v)),
  setNames(rep("SH2D1A", length(sh2d1a_v)), make_clean_names(sh2d1a_v)),
  setNames(rep("STXBP2", length(stxbp2_v)), make_clean_names(stxbp2_v)),
  setNames(rep("XIAP", length(xiap_v)), make_clean_names(xiap_v))
)

data_labelled <- data_labelled %>%
  rename(c_50del = c_50del_t)

gc_dl <- names(data_labelled) %>%
  str_subset("^c_")

gc_gx <- names(genetics_use) %>%
  str_subset("^c_")

# Add PRF1 variants from original patients panel
extra <- setdiff(gc_dl, names(gene_lookup))
gene_lookup[extra] <- "PRF1"
message("Assumed PRF1 (confirm): ", paste(extra, collapse = ", "))

geno_cols <- union(gc_dl, gc_gx)
n_dl <- nrow(data_labelled)

# Create two tables with standardised column names for joining
table_dl <- data_labelled %>%
  transmute(
    patient_id = row_number(),
    cohort = "data_labelled",
    original_id = as.character(`Patient ID`),
    gene_tested = "PRF1",
    perforin_pct = `Perforin Expression %`,
    perforin_state_recorded = NA_character_,
    gene_affected_recorded = NA_character_,
    result_clean = as.character(result),   # Result column from Excel
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

# Join the two tables
combined <- bind_rows(table_dl, table_gx) %>%
  mutate(across(all_of(geno_cols), ~ replace_na(as.numeric(.x), 0)))

# Flags anyone carrying a non-pathogenic variant labelled as
# PRF1, so it is not the analysis variable.
combined$gene_derived <- apply(as.matrix(combined[geno_cols]), 1, function(r) {
  hit <- geno_cols[r > 0]
  if (!length(hit)) "nm" else paste(sort(unique(gene_lookup[hit])), collapse = "; ")
})

# Perforin state:  >50 Normal | 10-50 Abnormal | <10 Absent
# Gene_affected: gene with an actual pathogenic finding, not merely tested
combined <- combined %>%
  mutate(
    perforin_state_calc = case_when(
      is.na(perforin_pct) ~ NA_character_,
      perforin_pct < 10 ~ "Absent",
      perforin_pct <= 50 ~ "Abnormal",
      TRUE ~ "Normal"),
    perforin_state = coalesce(perforin_state_recorded, perforin_state_calc),
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
    gene_affected = case_when(
      gene_affected_recorded == "PRF1" ~ "PRF1 (Pathogenic)",
      gene_affected_recorded == "nm" ~ "No variant detected",
      !is.na(gene_affected_recorded) ~ gene_affected_recorded,
      is.na(analysis_group) ~ NA_character_,
      analysis_group == "Mutation" ~ "PRF1 (Pathogenic)",
      rowSums(across(all_of(geno_cols))) > 0 ~ "PRF1 (Non-pathogenic)",
      TRUE ~ "No PRF1 variant detected"))

combined <- combined %>%
  mutate(control_type = case_when(
    cohort == "genetics" & gene_affected_recorded == "nm" ~ "Full negative",
    cohort == "data_labelled" & analysis_group == "No mutation" ~ "PRF1-panel negative",
    TRUE ~ NA_character_))

# Benign/intronic/synonymous — not counted toward zygosity
poly_cols <- c("c_272c_t", "c_900c_t", "c_434_t_c",
               "c_539_82c_t", "c_539_61g_a", "c_539_39g_a", "c_539_22g_c",
               "c_822_c_t", "c_726c_t")

# Pathogenic PRF1 alleles only
prf1_cols <- names(gene_lookup)[gene_lookup == "PRF1"]
path_cols <- setdiff(intersect(geno_cols, prf1_cols), poly_cols)

combined <- combined %>%
  mutate(path_alleles = rowSums(across(all_of(path_cols))),
         zygosity = case_when(
           path_alleles == 0 ~ "none",
           path_alleles == 1 ~ "monoallelic",
           TRUE ~ "biallelic"))

combined <- combined %>%
  mutate(group = case_when(
    gene_affected %in% c("STXBP2", "XIAP") ~ gene_affected,
    analysis_group == "Mutation" & zygosity == "biallelic"   ~ "PRF1 biallelic",
    analysis_group == "Mutation" & zygosity == "monoallelic" ~ "PRF1 monoallelic",
    TRUE ~ analysis_group))

# ============================================================================
# QC
# ============================================================================

# Is the perforin expression on the same scale in both cohorts?
combined %>%
  group_by(cohort) %>%
  summarise(n = n(),
            min = min(perforin_pct, na.rm = TRUE),
            max = max(perforin_pct, na.rm = TRUE),
            n_missing = sum(is.na(perforin_pct)))

# Do my thresholds reproduce the lab's own state calls?
combined %>%
  filter(!is.na(perforin_state_recorded)) %>%
  count(perforin_state_recorded, perforin_state_calc)

# Does the derived gene agree with the recorded one?
combined %>%
  filter(!is.na(gene_affected_recorded),
         gene_affected_recorded != gene_derived) %>%
  select(patient_id, original_id, gene_affected_recorded, gene_derived)

# Does result_clean agree with whether any variant was recorded?
combined %>%
  filter(cohort == "data_labelled") %>%
  mutate(has_variant = rowSums(across(all_of(geno_cols))) > 0) %>%
  count(result_clean, has_variant)

# Group sizes
combined %>%
  count(cohort, analysis_group)

combined %>%
  count(cohort, gene_affected)

# ============================================================================
# Final table
# ============================================================================

final_table <- combined %>%
  select(patient_id, cohort, original_id, gene_tested,
         perforin_pct, perforin_state,
         gene_affected, analysis_group, result_clean,
         zygosity, group,
         all_of(geno_cols)) %>%
  arrange(patient_id)

stopifnot(nrow(final_table) == n_dl + nrow(genetics_use))
stopifnot(!any(is.na(final_table[geno_cols])))

view(final_table)

# Original_id holds GOSH MRNs for the genetics cohort
# Get rid
final_table_clean <- final_table %>%
  select(-original_id)

write.xlsx(final_table_clean, "final_table.xlsx")

# ============================================================================
# Stats
# ============================================================================

# Does PRF1 expression differ in having biallelic PRF1 defects from no mutation?
final_table_clean %>%
  filter(group %in% c("PRF1 biallelic", "No mutation")) %>%
  wilcox.test(perforin_pct ~ group, data = .)

# p=3.05e-07
# Biallelic defects in PRF1 completely impairs protein expression.

# Does perforin expression differ between monoallelic PRF1 mutations from having no mutation?
final_table_clean %>%
  filter(group %in% c("PRF1 monoallelic", "No mutation")) %>%
  wilcox.test(perforin_pct ~ group, data = .)

# p=0.1457

# Any difference in perforin expression from polymorphism to no mutation?
final_table_clean %>%
  filter(group %in% c("Polymorphism", "No mutation")) %>%
  wilcox.test(perforin_pct ~ group, data = .)

# p=1.222e-08

final_table_clean %>%
  filter(cohort == "data_labelled", group %in% c("Polymorphism", "No mutation")) %>%
  wilcox.test(perforin_pct ~ group, data = .)

# p=0.980

# Cohort effects in main finding
final_table_clean %>%
  filter(cohort == "data_labelled", group %in% c("PRF1 biallelic", "No mutation")) %>%
  wilcox.test(perforin_pct ~ group, data = .)

# p=0.0234
# Replicated in data_labelled.

# Cohort effects in main finding
final_table_clean %>%
  filter(cohort == "genetics", group %in% c("PRF1 biallelic", "No mutation")) %>%
  wilcox.test(perforin_pct ~ group, data = .)

# p=

# ============================================================================
# Plots
# ============================================================================

# How does perforin expression differ between conditions?
final_table_clean %>%
  mutate(group = factor(group, levels = c("PRF1 biallelic", "PRF1 monoallelic",
                                          "Sequence variance", "Polymorphism",
                                          "No mutation", "STXBP2", "XIAP"))) %>%
  ggplot(aes(group, perforin_pct)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour = cohort), width = 0.15, alpha = 0.6) +
  geom_hline(yintercept = c(10, 50), linetype = "dashed", alpha = 0.4) +
  labs(x = NULL, y = "Perforin expression (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# How many patients in each condition are by each perforin state?
final_table_clean %>%
  count(group, perforin_state) %>%
  ggplot(aes(group, n, fill = factor(perforin_state,
                                     levels = c("Absent", "Abnormal", "Normal")))) +
  geom_col(position = "fill") +
  labs(x = NULL, y = "Proportion", fill = "Perforin state")

# Mutation-negative patients by cohort
final_table_clean %>%
  filter(group == "No mutation") %>%
  ggplot(aes(cohort, perforin_pct)) +
  geom_boxplot() + geom_jitter(width = 0.1) +
  labs(title = "Mutation-negative patients by cohort", y = "Perforin (%)")

# What are the variants carried most?
final_table_clean %>%
  summarise(across(all_of(geno_cols), ~ sum(.x > 0))) %>%
  pivot_longer(everything(), names_to = "variant", values_to = "n") %>%
  filter(n > 0) %>%
  mutate(gene = gene_lookup[variant]) %>%
  ggplot(aes(reorder(variant, n), n, fill = gene)) +
  geom_col() + coord_flip() +
  labs(x = NULL, y = "Patients carrying variant")

# Pathogenic PRF1 variants 
final_table_clean %>%
  select(patient_id, all_of(path_cols)) %>%
  pivot_longer(-patient_id, names_to = "variant", values_to = "alleles") %>%
  filter(alleles > 0) %>%
  mutate(gene = gene_lookup[variant]) %>%
  filter(gene == "PRF1") %>%
  count(variant, alleles) %>%
  mutate(state = if_else(alleles == 1, "Heterozygous", "Homozygous")) %>%
  ggplot(aes(reorder(variant, n), n, fill = state)) +
  geom_col() + coord_flip() +
  labs(x = NULL, y = "Patients carrying variant", fill = NULL)

final_table_clean %>%
  select(patient_id, all_of(path_cols)) %>%
  pivot_longer(-patient_id, names_to = "variant", values_to = "alleles") %>%
  filter(alleles > 0) %>%
  group_by(patient_id) %>%
  summarise(n_variants = n(), total_alleles = sum(alleles)) %>%
  filter(total_alleles >= 2) %>%
  count(genotype = if_else(n_variants == 1, "Homozygous", "Compound het"))
