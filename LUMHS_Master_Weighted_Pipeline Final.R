# ==============================================================================
# Thesis Project: 15-Gene Hypoxia Model Validation
# Institution: LUMHS MPhil Research Track
# Author: Dr. Aqeel Ahmed (MD)
# Part 1: Environment Setup, Multi-Omics Data Mining, and Clinical Filtering
# ==============================================================================

message("Checking dependencies...")
options(repos = c(CRAN = "https://posit.co"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_bioc <- c("curatedTCGAData", "TCGAutils", "MultiAssayExperiment")
new_bioc <- required_bioc[!(required_bioc %in% installed.packages()[, "Package"])]
if (length(new_bioc) > 0) {
  BiocManager::install(new_bioc, update = FALSE, ask = FALSE)
}

required_cran <- c("dplyr", "tidyr", "stringr", "survival", "survminer", "ggplot2")
new_cran <- required_cran[!(required_cran %in% installed.packages()[, "Package"])]
if (length(new_cran) > 0) {
  install.packages(new_cran)
}

library(curatedTCGAData)
library(TCGAutils)
library(MultiAssayExperiment)
library(dplyr)
library(tidyr)
library(stringr)
library(survival)
library(survminer)
library(ggplot2)

set.seed(2027)

message("Importing TCGA-LAML baseline dataset...")
laml_mae <- curatedTCGAData::curatedTCGAData(disease = "LAML", assays = c("RNASeq2GeneNorm"), version = "2.0.1", dry.run = FALSE)

target_assay_name <- names(MultiAssayExperiment::assays(laml_mae))[1]

assay_map <- as.data.frame(MultiAssayExperiment::sampleMap(laml_mae)) %>%
  dplyr::filter(assay == target_assay_name) %>%
  dplyr::mutate(colname = as.character(colname), primary = as.character(primary))

raw_assay <- MultiAssayExperiment::assay(laml_mae, 1)
aml_matrix_raw <- as.data.frame(raw_assay)

if (length(rownames(aml_matrix_raw)) == 0) {
  rownames(aml_matrix_raw) <- names(raw_assay)
}

clinical_data_raw <- as.data.frame(MultiAssayExperiment::colData(laml_mae))
clinical_data_raw$primary_id <- rownames(clinical_data_raw)

buffa_15_genes <- c("ACOT7", "ADM", "ALDOA", "CDKN3", "ENO1", "LDHA", "MIF", 
                    "MRPS17", "NDRG1", "P4HA1", "PGAM1", "SLC2A1", "TPI1", "TUBB6", "VEGFA")

aml_matrix_raw$gene_name <- trimws(sub("\\|.*", "", rownames(aml_matrix_raw)))
b15_matrix <- aml_matrix_raw[aml_matrix_raw$gene_name %in% buffa_15_genes, ]
message(paste("Matched", nrow(b15_matrix), "out of 15 target signature loci."))

all_cols <- colnames(clinical_data_raw)
clinical_data_raw$fallback_na <- as.numeric(NA)

age_idx <- grep("age|birth|diagnostic_age", all_cols, ignore.case = TRUE)
age_col <- if(length(age_idx) > 0) all_cols[age_idx[1]] else "fallback_na"

raw_age_vector <- abs(as.numeric(clinical_data_raw[[age_col]]))
clinical_data_raw$extracted_age_years <- ifelse(!is.na(raw_age_vector) & raw_age_vector > 150, 
                                                 raw_age_vector / 365.25, raw_age_vector)

survival_base <- clinical_data_raw %>%
  dplyr::select(primary_id, days_to_death, days_to_last_followup, vital_status, extracted_age_years) %>%
  dplyr::mutate(
    time_days = ifelse(!is.na(days_to_death), as.numeric(days_to_death), as.numeric(days_to_last_followup)),
    status_numeric = ifelse(vital_status == 1 | grepl("dead|deceased", vital_status, ignore.case = TRUE), 1, 0),
    age_years = extracted_age_years
  ) %>%
  dplyr::filter(!is.na(time_days) & time_days > 0 & age_years >= 18) %>%
  dplyr::select(primary_id, time_days, status_numeric, age_years)

blast_idx <- grep("blast", all_cols, ignore.case = TRUE)
b_col <- if(length(blast_idx) > 0) all_cols[blast_idx[1]] else "fallback_na"

myeloid_idx <- grep("myeloid|peripheral", all_cols, ignore.case = TRUE)
m_col <- if(length(myeloid_idx) > 0) all_cols[myeloid_idx[1]] else "fallback_na"

hb_idx <- grep("hb|hemoglobin", all_cols, ignore.case = TRUE)
hb_col <- if(length(hb_idx) > 0) all_cols[hb_idx[1]] else "fallback_na"

extracted_blast   <- as.numeric(clinical_data_raw[[b_col]])
extracted_myeloid <- as.numeric(clinical_data_raw[[m_col]])
extracted_hb      <- as.numeric(clinical_data_raw[[hb_col]])

aml_clean_clinical <- clinical_data_raw %>%
  dplyr::mutate(
    bone_marrow_blast = extracted_blast,
    peripheral_blood_myeloid = extracted_myeloid,
    hemoglobin_level = extracted_hb
  ) %>%
  dplyr::select(primary_id, bone_marrow_blast, peripheral_blood_myeloid, hemoglobin_level)

message("Part 1 processing completed successfully. Base environment is locked.")
# ==============================================================================
# Thesis Project: 15-Gene Hypoxia Model Validation
# Institution: LUMHS MPhil Research Track
# Author: Dr. Aqeel Ahmed (MD)
# Part 2: Genomic Normalization, Weighting, and Global Survival Modeling
# ==============================================================================

message("Executing Log2 transformations and cross-patient Z-score normalization arrays...")
tcga_barcode_cols <- colnames(b15_matrix)[grepl("^TCGA|^LAML", colnames(b15_matrix), ignore.case = TRUE)]

b15_long <- b15_matrix %>%
  dplyr::select(gene_name, dplyr::all_of(tcga_barcode_cols)) %>%
  tidyr::pivot_longer(cols = -gene_name, names_to = "colname", values_to = "expression_value") %>%
  dplyr::inner_join(assay_map, by = "colname") %>%
  dplyr::mutate(patient_id = primary) %>%
  dplyr::mutate(log_val = log2(expression_value + 1)) %>%
  dplyr::group_by(gene_name) %>%
  dplyr::mutate(z_score = as.numeric(scale(log_val))) %>%
  dplyr::ungroup()

wide_expression <- b15_long %>%
  dplyr::select(patient_id, gene_name, z_score) %>%
  dplyr::group_by(patient_id, gene_name) %>%
  dplyr::summarize(z_score = mean(z_score, na.rm = TRUE), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = gene_name, values_from = z_score) %>%
  dplyr::rename(primary_id = patient_id)

cox_train_df <- dplyr::inner_join(survival_base, wide_expression, by = "primary_id")
message(paste("Merged dataset contains:", nrow(cox_train_df), "pristine unified patient records."))

gene_weights <- c()
available_genes <- intersect(buffa_15_genes, colnames(wide_expression))

print("--- Check Proportional Hazards Assumption via Schoenfeld Residuals ---")
for(gene in available_genes) {
  formula_string <- paste("survival::Surv(time_days, status_numeric) ~", paste0("`", gene, "`"))
  cox_model <- survival::coxph(as.formula(formula_string), data = cox_train_df)
  gene_weights[gene] = coef(cox_model)
  
  ph_test <- survival::cox.zph(cox_model)
  
  global_p_val <- ph_test[["table"]][nrow(ph_test[["table"]]), 3]
  print(paste("Gene Target:", gene, "| Proportional Hazards Global p-value =", round(global_p_val, 4)))
}

message("Part 2 calculation cycle completed. Genomic model weights are securely cached.")
# ==============================================================================
# Thesis Project: 15-Gene Hypoxia Model Validation
# Institution: LUMHS MPhil Research Track
# Author: Dr. Aqeel Ahmed (MD)
# Part 3: Local Integration, Two-Tier Structuring, and Cohort Analytics
# ==============================================================================

message("Computing global weighted hypoxia scores and mapping multivariable covariates...")

weighted_scores <- cox_train_df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(weighted_hypoxia_score = sum(c_across(dplyr::all_of(available_genes)) * gene_weights, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::select(primary_id, weighted_hypoxia_score)

tcga_master_data <- dplyr::inner_join(aml_clean_clinical, weighted_scores, by = "primary_id") %>%
  dplyr::inner_join(survival_base, by = "primary_id")

median_weighted <- median(tcga_master_data[["weighted_hypoxia_score"]], na.rm = TRUE)
tcga_master_data[["hypoxia_group"]] <- ifelse(tcga_master_data[["weighted_hypoxia_score"]] > median_weighted, 
                                               "High Risk Hypoxia", "Low Risk Hypoxia")

final_fit <- survival::survfit(survival::Surv(time_days, status_numeric) ~ hypoxia_group, data = tcga_master_data)

active_covariates <- c("hypoxia_group")
if("age_years" %in% colnames(tcga_master_data)) {
  if(sum(!is.na(tcga_master_data[["age_years"]])) > 10) active_covariates <- c(active_covariates, "age_years")
}
if("bone_marrow_blast" %in% colnames(tcga_master_data)) {
  if(sum(!is.na(tcga_master_data[["bone_marrow_blast"]])) > 10) active_covariates <- c(active_covariates, "bone_marrow_blast")
}

print("--- Multivariable Cox Regression Model ---")
mv_formula_string <- paste("survival::Surv(time_days, status_numeric) ~", paste(active_covariates, collapse = " + "))
multivariable_model <- survival::coxph(as.formula(mv_formula_string), data = tcga_master_data)
print(summary(multivariable_model))

message("Generating high-resolution global survival curve plot...")

final_plot <- survminer::ggsurvplot(
  final_fit, data = tcga_master_data, pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, palette = c("#E41A1C", "#377EB8"), 
  legend.labs = c("High Risk Hypoxia", "Low Risk Hypoxia"),
  xlab = "Survival Time (Days)", ylab = "Survival Probability",
  title = "Cox-Weighted 15-Gene Buffa Hypoxia Prognostic Model",
  risk.table.height = 0.22, tables.theme = survminer::theme_cleantable()
)

# Open a direct high-resolution file stream explicitly before calling the arranger tool.
# This prevents the silent layout drop and guarantees both curve and table print together.
png("Figure_1_Survival_Kinetics.png", width = 2400, height = 2100, res = 300)
survminer::arrange_ggsurvplots(list(final_plot), print = TRUE, ncol = 1, nrow = 1)
dev.off()

message("Saved unified Figure 1 cleanly to your working directory.")
write.csv(tcga_master_data, "TCGA_Global_Survival_Reference_Data.csv", row.names = FALSE)

if(file.exists("LUMHS_Local_AML_Cohort.csv")) {
  local_data <- read.csv("LUMHS_Local_AML_Cohort.csv", na.strings = c("NA", " ", "", "."))
  
  colnames(local_data) <- tolower(colnames(local_data))
  required_columns <- c("bone_marrow_blast", "peripheral_blood_myeloid", "hemoglobin_level")
  missing_cols <- setdiff(required_columns, colnames(local_data))
  
  if (length(missing_cols) > 0) {
    stop(paste("Data template is missing columns:", paste(missing_cols, collapse = ", ")))
  }

  if(nrow(local_data) > 1) {
    global_cellular_df <- tcga_master_data %>%
      dplyr::select(bone_marrow_blast, peripheral_blood_myeloid) %>%
      dplyr::mutate(Cohort = "Global (TCGA)") %>%
      dplyr::filter(!is.na(bone_marrow_blast) | !is.na(peripheral_blood_myeloid))
    
    local_cellular_df <- local_data %>%
      dplyr::select(bone_marrow_blast, peripheral_blood_myeloid) %>%
      dplyr::mutate(Cohort = "Local (LUMHS)") %>%
      dplyr::filter(!is.na(bone_marrow_blast) | !is.na(peripheral_blood_myeloid))
    
    combined_cellular_df <- rbind(global_cellular_df, local_cellular_df) %>%
      dplyr::mutate(Cohort = as.factor(Cohort))
    
    global_hb_df <- tcga_master_data %>%
      dplyr::select(hemoglobin_level) %>%
      dplyr::mutate(Cohort = "Global (TCGA)") %>%
      dplyr::filter(!is.na(hemoglobin_level))
      
    if("transfused_within_14_days" %in% colnames(local_data)) {
      local_transfusion_clean <- tolower(trimws(as.character(local_data$transfused_within_14_days)))
      local_hb_df <- local_data %>%
        dplyr::filter(local_transfusion_clean == "no" | local_transfusion_clean == "0" | is.na(local_transfusion_clean)) %>%
        dplyr::select(hemoglobin_level) %>%
        dplyr::mutate(Cohort = "Local (LUMHS)") %>%
        dplyr::filter(!is.na(hemoglobin_level))
    } else {
      local_hb_df <- local_data %>%
        dplyr::select(hemoglobin_level) %>%
        dplyr::mutate(Cohort = "Local (LUMHS)") %>%
        dplyr::filter(!is.na(hemoglobin_level))
    }
    
    combined_hb_df <- rbind(global_hb_df, local_hb_df) %>%
      dplyr::mutate(Cohort = as.factor(Cohort))
    
    if (all(is.na(combined_cellular_df$peripheral_blood_myeloid[combined_cellular_df$Cohort == "Global (TCGA)"]))) {
      message("--- NOTE: Global (TCGA) mature peripheral myeloid metrics are naturally blank (NA) in the repository. ---")
      message("--- Cellular cross-population audits will focus on the complete bone marrow blast architecture layers. ---")
    }
    
    blast_wilcox <- wilcox.test(bone_marrow_blast ~ Cohort, data = combined_cellular_df, exact = FALSE)
    hb_wilcox    <- wilcox.test(hemoglobin_level ~ Cohort, data = combined_hb_df, exact = FALSE)
    
    print("--- Cross-Population Non-Parametric Comparison Results (Mann-Whitney U) ---")
    print(blast_wilcox)
    print(hb_wilcox)
    
    valid_global_myeloid <- sum(!is.na(combined_cellular_df$peripheral_blood_myeloid[combined_cellular_df$Cohort == "Global (TCGA)"]))
    if(valid_global_myeloid > 2) {
      myeloid_wilcox <- wilcox.test(peripheral_blood_myeloid ~ Cohort, data = combined_cellular_df, exact = FALSE)
      print(myeloid_wilcox)
    } else {
      message("Skipping myeloid cohort comparison: Insufficient global baseline variance detected.")
    }
    
    print("--- Local cohort Spearman correlation results ---")
    local_spearman <- cor.test(local_cellular_df[["bone_marrow_blast"]], 
                               local_cellular_df[["peripheral_blood_myeloid"]], 
                               method = "spearman", exact = FALSE)
    print(local_spearman)
# ==============================================================================
# Thesis Project: 15-Gene Hypoxia Model Validation
# Institution: LUMHS MPhil Research Track
# Author: Dr. Aqeel Ahmed (MD)
# Part 4: Advanced Visualizations, Sample Attrition Tracking, and Table 1 Export
# ==============================================================================

    comp_scatterplot <- ggplot2::ggplot(local_cellular_df, aes(x = bone_marrow_blast, y = peripheral_blood_myeloid)) +
      ggplot2::geom_point(color = "#E74C3C", alpha = 0.6, size = 2.5) +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x, color = "#2C3E50", se = TRUE, fill = "#BDC3C7", alpha = 0.3) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5, face = "italic"),
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        title = "Local Marrow Crowding Dynamics",
        subtitle = paste0("Spearman r = ", round(local_spearman$estimate, 2), " (P = ", format.pval(local_spearman$p.value, digits = 4), ")"),
        x = "Bone Marrow Blast Infiltration (%)",
        y = "Peripheral Blood Myeloid Cells (%)"
      )

    ggplot2::ggsave("Figure_3_Local_Marrow_Crowding.png", plot = comp_scatterplot, width = 7, height = 5, dpi = 300)
    message("Saved Figure 3.")

    global_plotting_prep <- tcga_master_data %>% 
      dplyr::select(bone_marrow_blast, peripheral_blood_myeloid, hemoglobin_level) %>%
      dplyr::mutate(Cohort = "Global (TCGA)")
      
    local_plotting_prep <- local_data %>%
      dplyr::select(bone_marrow_blast, peripheral_blood_myeloid, hemoglobin_level) %>%
      dplyr::mutate(Cohort = "Local (LUMHS)")
      
    if("transfused_within_14_days" %in% colnames(local_data)) {
      local_transfusion_clean <- tolower(trimws(as.character(local_data$transfused_within_14_days)))
      local_plotting_prep$hemoglobin_level[local_transfusion_clean == "yes" | local_transfusion_clean == "1" | local_transfusion_clean == "true"] <- NA
    }
    
    combined_cohorts_df <- rbind(global_plotting_prep, local_plotting_prep)

    plotting_df <- combined_cohorts_df %>%
      tidyr::pivot_longer(cols = c(bone_marrow_blast, peripheral_blood_myeloid, hemoglobin_level), 
                          names_to = "Parameter", values_to = "value") %>%
      dplyr::filter(!is.na(value)) 
    
    panel_labels <- c(
      "bone_marrow_blast" = "Bone Marrow Blast Infiltration (%)",
      "peripheral_blood_myeloid" = "Peripheral Blood Myeloid Cells (%)",
      "hemoglobin_level" = "Baseline Hemoglobin Level (g/dL)"
    )
    
    comp_boxplot <- ggplot2::ggplot(plotting_df, aes(x = Cohort, y = value, fill = Cohort)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = 16, width = 0.5, color = "#2C3E50") +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.3, size = 1.5, aes(color = Cohort)) +
      ggplot2::facet_wrap(~Parameter, scales = "free_y", labeller = ggplot2::as_labeller(panel_labels), ncol = 3) +
      ggplot2::scale_fill_manual(values = c("#3498DB", "#E74C3C")) +    
      ggplot2::scale_color_manual(values = c("#2980B9", "#C0392B")) +  
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "#ECF0F1", color = "#BDC3C7"),
        strip.text = ggplot2::element_text(face = "bold", color = "#2C3E50", size = 10),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(face = "bold", size = 11),
        legend.position = "none",
        plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        title = "Cross-Population Bio-Audit: Global Reference vs. Local Sindh Cohort",
        y = "Tissue Compartment Metrics / Clinical Density"
      )
    
    ggplot2::ggsave("Figure_2_Cross_Population_Comparison.png", plot = comp_boxplot, width = 11, height = 5, dpi = 300)
    message("Saved Figure 2.")
    
    run_local_survival <- FALSE
    if("local_time_days" %in% colnames(local_data) & "local_status" %in% colnames(local_data)) {
      clean_survival_subset <- local_data %>%
        dplyr::select(patient_id, bone_marrow_blast, local_time_days, local_status) %>%
        dplyr::filter(!is.na(local_time_days) & !is.na(local_status) & !is.na(bone_marrow_blast))
      if(nrow(clean_survival_subset) >= 15) { 
        run_local_survival = TRUE
      }
    }
    
    if(run_local_survival) {
      message("Survival indices confirmed. Running local survival module...")
      median_blast <- median(clean_survival_subset$bone_marrow_blast, na.rm = TRUE)
      clean_survival_subset$blast_group <- ifelse(clean_survival_subset$bone_marrow_blast > median_blast, 
                                                   "High Blast Load", "Low Blast Load")
      
      local_fit <- survival::survfit(survival::Surv(local_time_days, local_status) ~ blast_group, 
                                     data = clean_survival_subset)
      
      local_survival_plot <- survminer::ggsurvplot(
        local_fit, data = clean_survival_subset, pval = TRUE, conf.int = FALSE,
        risk.table = TRUE, palette = c("#D35400", "#27AE60"),
        legend.labs = c("High Blast Load", "Low Blast Load"),
        xlab = "Follow-up Duration (Days)", ylab = "Survival Probability",
        title = "LUMHS Local Validation Cohort Internal Survival Kinetics",
        risk.table.height = 0.22, tables.theme = survminer::theme_cleantable()
      )
      
      png("Figure_4_Local_Survival_Kinetics.png", width = 2400, height = 1800, res = 300)
      print(local_survival_plot)
      dev.off()
      message("Saved Figure 4.")
    } else {
      message("Local survival parameters insufficient. Skipping Figure 4 extraction.")
    }
    
    message("Compiling baseline summary data...")
    if("hemoglobin_level" %in% colnames(local_data)) {
      local_data$anemia_status <- ifelse(local_data$hemoglobin_level < 11, "Anemic", "Normal")
    } else {
      local_data$anemia_status = "Not Recorded"
    }

    evaluable_hb_count <- sum(!is.na(local_data$hemoglobin_level))
    anemia_percentage_string <- if(evaluable_hb_count > 0) {
      paste0(sum(local_data$anemia_status == "Anemic", na.rm = TRUE), " (", 
             round((sum(local_data$anemia_status == "Anemic", na.rm = TRUE) / evaluable_hb_count) * 100, 1), "%)")
    } else {
      "0 (0%)"
    }

    clinical_summary <- data.frame(
      Metric = c("Total Patients (n)", 
                 "Mean BM Blast Infiltration (%)", 
                 "Mean Peripheral Myeloid Cells (%)", 
                 "Patients presenting with Anemia Baseline (%)"),
      Value = c(
        nrow(local_data),
        round(mean(local_data$bone_marrow_blast, na.rm = TRUE), 2),
        round(mean(local_data$peripheral_blood_myeloid, na.rm = TRUE), 2),
        anemia_percentage_string
      )
    )

    print(clinical_summary)
    write.csv(clinical_summary, "Table_1_Clinical_Summary.csv", row.names = FALSE)
    message("Saved Table 1 summary matrix.")
  }
}
save.image(file = "LUMHS_Hypoxia_Project_Processed_Outputs.RData")
message("Pipeline execution completed successfully.")


