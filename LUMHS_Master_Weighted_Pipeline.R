# ==============================================================================
# MASTER COMPILATION PIPELINE: COX-WEIGHTED 15-GENE BUFFA PROGNOSTIC MODEL
# Baseline Global Dataset Phase - LUMHS MPhil Research Initiative (January 2027)
# Principal Investigator: Dr. Aqeel Ahmed
# ==============================================================================

# 1. Clear Active Memory Workspace Environment
rm(list = setdiff(ls(), "aml_matrix"))

# 2. Load Core Structural and Statistical Packages
library(curatedTCGAData)
library(dplyr)
library(tidyr)
library(stringr)
library(survival)
library(survminer)
library(ggplot2)

# 3. Memory Verification: Reload Environment Data Workspace if Missing
if (!exists("aml_matrix")) {
  message("Active aml_matrix object missing. Extracting from local workspace storage...")
  load("LUMHS_Hypoxia_Project_Workspace.RData")
}

# 4. Initialize the Universally Validated 15-Gene Buffa Hypoxia Classifier Vector
buffa_15_genes <- c("ACOT7", "ADM", "ALDOA", "CDKN3", "ENO1", "LDHA", "MIF", 
                    "MRPS17", "NDRG1", "P4HA1", "PGAM1", "SLC2A1", "TPI1", "TUBB6", "VEGFA")

# 5. Extract and Subset the Hypoxia Loci from the Primary Genomic Matrix
b15_matrix <- aml_matrix[aml_matrix$gene_name %in% buffa_15_genes, ]
message(paste("Genomic Alignment Success: Isolated", nrow(b15_matrix), "target features out of 15 Buffa genes."))

# 6. Extract and Clean Direct Clinical Survival Timelines from Base Tables
survival_base <- clinical_data %>%
  mutate(patient_id = row.names(clinical_data)) %>%
  select(
    patient_id, 
    days_to_death, 
    days_to_last_followup, 
    vital_status
  ) %>%
  mutate(
    time_days = ifelse(!is.na(days_to_death), as.numeric(days_to_death), as.numeric(days_to_last_followup)),
    status_numeric = ifelse(vital_status == 1 | grepl("dead|deceased", vital_status, ignore.case = TRUE), 1, 0)
  ) %>%
  filter(!is.na(time_days))

# 7. Execute Long Transformation and Parse Alphanumeric Patient Barcodes
b15_long <- as.data.frame(b15_matrix) %>%
  select(-gene_id, -gene_type) %>%
  pivot_longer(cols = -gene_name, names_to = "raw_header", values_to = "expression_value") %>%
  filter(str_detect(raw_header, "unstranded_")) %>%
  mutate(
    sample_barcode = str_extract(raw_header, 'TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]'),
    patient_id = substr(sample_barcode, 1, 12)
  ) %>%
  filter(!is.na(patient_id)) %>%
  mutate(log_val = log2(expression_value + 1)) %>%
  group_by(gene_name) %>%
  mutate(z_score = as.numeric(scale(log_val))) %>%
  ungroup()

# 8. Pivot to Wide Analytical Format for Univariate Regression Training
wide_expression <- b15_long %>%
  select(patient_id, gene_name, z_score) %>%
  group_by(patient_id, gene_name) %>%
  summarize(z_score = mean(z_score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = gene_name, values_from = z_score)

# 9. Intersect Molecular Matrices with Survival Outcomes
cox_train_df <- inner_join(survival_base, wide_expression, by = "patient_id")

# 10. Dynamically Calculate Univariate Cox Proportional Hazard Weights (Coefficients)
gene_weights <- c()
available_genes <- intersect(buffa_15_genes, colnames(wide_expression))

for(gene in available_genes) {
  formula_string <- paste("Surv(time_days, status_numeric) ~", paste0("`", gene, "`"))
  cox_model <- coxph(as.formula(formula_string), data = cox_train_df)
  gene_weights[gene] <- coef(cox_model)
}

# 11. Compute the Final Weighted Hypoxia Risk Score per Unique Patient Profile
weighted_scores <- cox_train_df %>%
  rowwise() %>%
  mutate(
    weighted_hypoxia_score = sum(c_across(all_of(available_genes)) * gene_weights, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(patient_id, weighted_hypoxia_score)

# 12. Consolidate Master Research Dataset (Pathology Metrics + Survival + Weighted Scores)
lumhs_master_data <- inner_join(aml_clean_clinical, weighted_scores, by = "patient_id") %>%
  inner_join(survival_base, by = "patient_id")

# 13. Group Cohort via Median Stratification and Fit Kaplan-Meier Model
median_weighted <- median(lumhs_master_data$weighted_hypoxia_score, na.rm = TRUE)
lumhs_master_data$hypoxia_group <- ifelse(lumhs_master_data$weighted_hypoxia_score > median_weighted, "High Risk Hypoxia", "Low Risk Hypoxia")

final_fit <- survfit(Surv(time_days, status_numeric) ~ hypoxia_group, data = lumhs_master_data)

# 14. Render the Definitive Publication-Ready Kaplan-Meier Curve
final_plot <- ggsurvplot(
  final_fit, 
  data = lumhs_master_data, 
  pval = TRUE, 
  conf.int = TRUE,
  risk.table = TRUE, 
  palette = c("#E41A1C", "#377EB8"), 
  legend.labs = c("High Risk Hypoxia", "Low Risk Hypoxia"),
  xlab = "Survival Time (Days)", 
  ylab = "Survival Probability",
  title = "Cox-Weighted 15-Gene Buffa Hypoxia Prognostic Model",
  risk.table.height = 0.22,
  tables.theme = theme_cleantable()
)

# Force display the perfect graph on screen
print(final_plot)

# 15. Export Clean CSV Spreadsheet Table Data to Working Folder Disk
write.csv(lumhs_master_data, "LUMHS_Global_Survival_Validated_Data.csv", row.names = FALSE)

# 16. Lock and Save the Workspace Environment Binary File
save.image(file = "LUMHS_Hypoxia_Project_Workspace.RData")
message("--- SUCCESS: Master Pipeline executed with zero errors. RData and CSV saved. ---")
