# ==============================================================================
# MASTER COMPILATION PIPELINE: COX-WEIGHTED 15-GENE BUFFA PROGNOSTIC MODEL
# Baseline Global Dataset Phase - LUMHS MPhil Research Initiative (2027)
# Principal Investigator: Dr. Aqeel Ahmed (MD)
# Target Journal Submission: PLOS ONE (Technical Soundness Track)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. AUTOMATED BOOTSTRAPPING FOR CRAN AND BIOCONDUCTOR INFRASTRUCTURES
# ------------------------------------------------------------------------------
message("Initializing package dependency audit toolsets...")
options(repos = c(CRAN = "https://r-project.org"))

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

# ------------------------------------------------------------------------------
# 2. DYNAMIC PROGRAMMATIC FETCH & OFFICIAL MAP LINKING
# ------------------------------------------------------------------------------
message("Downloading raw TCGA-LAML multi-assay experiment data directly from API...")
laml_mae <- curatedTCGAData::curatedTCGAData(disease = "LAML", assays = c("RNASeq2GeneNorm"), version = "2.0.1", dry.run = FALSE)

target_assay_name <- names(MultiAssayExperiment::assays(laml_mae))

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
# ------------------------------------------------------------------------------
# 3. INITIALIZE THE BUFFA 15-GENE HYPOXIA CLASSIFIER VECTOR
# ------------------------------------------------------------------------------
buffa_15_genes <- c("ACOT7", "ADM", "ALDOA", "CDKN3", "ENO1", "LDHA", "MIF", 
                    "MRPS17", "NDRG1", "P4HA1", "PGAM1", "SLC2A1", "TPI1", "TUBB6", "VEGFA")

aml_matrix_raw$gene_name <- sub("\\|.*", "", rownames(aml_matrix_raw))
b15_matrix <- aml_matrix_raw[aml_matrix_raw$gene_name %in% buffa_15_genes, ]
message(paste("Genomic Alignment Success: Isolated", nrow(b15_matrix), "target features out of 15 Buffa genes."))

# ------------------------------------------------------------------------------
# 4. EXTRACT AND CLEAN DIRECT CLINICAL SURVIVAL TIMELINES
# ------------------------------------------------------------------------------
survival_base <- clinical_data_raw %>%
  dplyr::select(primary_id, days_to_death, days_to_last_followup, vital_status, age) %>%
  dplyr::mutate(
    time_days = ifelse(!is.na(days_to_death), as.numeric(days_to_death), as.numeric(days_to_last_followup)),
    status_numeric = ifelse(vital_status == 1 | grepl("dead|deceased", vital_status, ignore.case = TRUE), 1, 0),
    age_years = as.numeric(age)
  ) %>%
  dplyr::filter(!is.na(time_days) & time_days > 0)

all_cols <- colnames(clinical_data_raw)
blast_candidates <- c("percent_bone_marrow_blasts", "percent_blasts", 
                      "leukemia_blast_cell_cellularity_percentage", "bone_marrow_blasts")
blast_col_match <- intersect(blast_candidates, all_cols)

myeloid_candidates <- c("percent_myeloid_cells_peripheral_blood", "percent_myeloid_cells",
                        "peripheral_blood_myeloid_percentage", "myeloid_cells")
myeloid_col_match <- intersect(myeloid_candidates, all_cols)

clinical_data_raw$fallback_na <- as.numeric(NA)

b_col <- if(length(blast_col_match) > 0) blast_col_match else "fallback_na"
m_col <- if(length(myeloid_col_match) > 0) myeloid_col_match else "fallback_na"

aml_clean_clinical <- clinical_data_raw %>%
  dplyr::mutate(
    bone_marrow_blast = as.numeric(clinical_data_raw[[b_col]]),
    peripheral_blood_myeloid = as.numeric(clinical_data_raw[[m_col]])
  ) %>%
  dplyr::select(primary_id, bone_marrow_blast, peripheral_blood_myeloid)

# ------------------------------------------------------------------------------
# 5. LONG TRANSFORMATION USING MULTI-ASSAY EXPERIMENT CONNECTIONS
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# 6. WIDE ANALYTICAL FORMATTING & COHORT INTERSECTION
# ------------------------------------------------------------------------------
wide_expression <- b15_long %>%
  dplyr::select(patient_id, gene_name, z_score) %>%
  dplyr::group_by(patient_id, gene_name) %>%
  dplyr::summarize(z_score = mean(z_score, na.rm = TRUE), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = gene_name, values_from = z_score) %>%
  dplyr::rename(primary_id = patient_id)

cox_train_df <- dplyr::inner_join(survival_base, wide_expression, by = "primary_id")
message(paste("CRITICAL MERGE VERIFICATION: Successfully synchronized", nrow(cox_train_df), "matching patient records."))

# ------------------------------------------------------------------------------
# 7. DYNAMIC UNIVARIATE COX MODEL GENERATION & PROACTIVE PROPORTIONAL HAZARDS AUDIT
# ------------------------------------------------------------------------------
gene_weights <- c()
available_genes <- intersect(buffa_15_genes, colnames(wide_expression))

print("--- ADVANCED COMPUTATIONAL AUDIT: PROPORTIONAL HAZARDS ASSUMPTION LOG ---")
for(gene in available_genes) {
  formula_string <- paste("Surv(time_days, status_numeric) ~", paste0("`", gene, "`"))
  cox_model <- survival::coxph(as.formula(formula_string), data = cox_train_df)
  gene_weights[gene] <- coef(cox_model)
  
  # PROACTIVE FIX: Executes Schoenfeld residual verification automatically for each weight locus
  ph_test <- survival::cox.zph(cox_model)
  print(paste("Locus:", gene, "| Schoenfeld Residual p-value =", round(ph_test$table[1, 3], 4)))
}

# ------------------------------------------------------------------------------
# 8. RISK CALCULATIONS, STRATIFICATION & MULTIVARIABLE COVARYING ADJUSTMENT
# ------------------------------------------------------------------------------
weighted_scores <- cox_train_df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(weighted_hypoxia_score = sum(c_across(dplyr::all_of(available_genes)) * gene_weights, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::select(primary_id, weighted_hypoxia_score)

lumhs_master_data <- dplyr::inner_join(aml_clean_clinical, weighted_scores, by = "primary_id") %>%
  dplyr::inner_join(survival_base, by = "primary_id")

median_weighted <- median(lumhs_master_data$weighted_hypoxia_score, na.rm = TRUE)
lumhs_master_data$hypoxia_group <- ifelse(lumhs_master_data$weighted_hypoxia_score > median_weighted, 
                                          "High Risk Hypoxia", "Low Risk Hypoxia")

final_fit <- survival::survfit(survival::Surv(time_days, status_numeric) ~ hypoxia_group, data = lumhs_master_data)

# PROACTIVE FIX: Run multivariable model to demonstrate score independence against clinical confounders
print("--- PROACTIVE INDEPENDENT COVARYING AUDIT: MULTIVARIABLE SURVIVAL MODEL ---")
multivariable_model <- survival::coxph(
  survival::Surv(time_days, status_numeric) ~ hypoxia_group + age_years + bone_marrow_blast, 
  data = lumhs_master_data
)
print(summary(multivariable_model))
# ------------------------------------------------------------------------------
# 9. GRAPHICAL RENDERING & TABULAR OUTPUT STORAGE (FIGURE 1)
# ------------------------------------------------------------------------------
final_plot <- survminer::ggsurvplot(
  final_fit, data = lumhs_master_data, pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, palette = c("#E41A1C", "#377EB8"), 
  legend.labs = c("High Risk Hypoxia", "Low Risk Hypoxia"),
  xlab = "Survival Time (Days)", ylab = "Survival Probability",
  title = "Cox-Weighted 15-Gene Buffa Hypoxia Prognostic Model",
  risk.table.height = 0.22, tables.theme = survminer::theme_cleantable()
)

png("Figure_1_Survival_Kinetics.png", width = 2400, height = 1800, res = 300)
survminer:::print.ggsurvplot(final_plot, newpage = FALSE)
dev.off()
message("Figure 1 successfully generated and saved to your project directory!")

write.csv(lumhs_master_data, "LUMHS_Global_Survival_Validated_Data.csv", row.names = FALSE)

# ==============================================================================
# PHASE 3: LOCAL TO GLOBAL CROSS-POPULATION COMPARISON (LUMHS VALIDATION)
# ==============================================================================
if(file.exists("LUMHS_Local_AML_Cohort.csv")) {
  local_data <- read.csv("LUMHS_Local_AML_Cohort.csv", na.strings = c("NA", " ", "", "."))
  
  colnames(local_data) <- tolower(colnames(local_data))
  required_columns <- c("bone_marrow_blast", "peripheral_blood_myeloid")
  missing_cols <- setdiff(required_columns, colnames(local_data))
  
  if (length(missing_cols) > 0) {
    stop(paste("CRITICAL ERROR: The template is missing columns:", paste(missing_cols, collapse = ", ")))
  }

  if(nrow(local_data) > 1) {
    global_comparison_df <- lumhs_master_data %>%
      dplyr::select(bone_marrow_blast, peripheral_blood_myeloid) %>%
      dplyr::mutate(Cohort = "Global (TCGA)") %>%
      dplyr::filter(!is.na(bone_marrow_blast) & !is.na(peripheral_blood_myeloid))
    
    local_comparison_df <- local_data %>%
      dplyr::select(bone_marrow_blast, peripheral_blood_myeloid) %>%
      dplyr::mutate(Cohort = "Local (LUMHS)") %>%
      dplyr::filter(!is.na(bone_marrow_blast) & !is.na(peripheral_blood_myeloid))
    
    combined_cohorts_df <- rbind(global_comparison_df, local_comparison_df)
    
    blast_t_test <- t.test(bone_marrow_blast ~ Cohort, data = combined_cohorts_df)
    myeloid_t_test <- t.test(peripheral_blood_myeloid ~ Cohort, data = combined_cohorts_df)
    
    print("--- CROSS-POPULATION BASELINE DIFFERENCES ---")
    print(blast_t_test)
    print(myeloid_t_test)
    
    local_spearman <- cor.test(local_comparison_df$bone_marrow_blast, 
                               local_comparison_df$peripheral_blood_myeloid, 
                               method = "spearman", exact = FALSE)
    
    print("--- LOCAL SINDH COHORT CORRELATION RESULTS ---")
    print(local_spearman)
    
    comp_scatterplot <- ggplot2::ggplot(local_comparison_df, aes(x = bone_marrow_blast, y = peripheral_blood_myeloid)) +
      ggplot2::geom_point(color = "#E74C3C", alpha = 0.6, size = 2.5) +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x, color = "#2C3E50", se = TRUE, fill = "#BDC3C7", alpha = 0.3) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"),
        panel.grid.minor = element_blank()
      ) +
      ggplot2::labs(
        title = "Local Marrow Crowding Dynamics",
        subtitle = paste0("Spearman r = ", round(local_spearman$estimate, 2), " (p = ", format.pval(local_spearman$p.value, digits = 4), ")"),
        x = "Bone Marrow Blast Infiltration (%)",
        y = "Peripheral Blood Myeloid Cells (%)"
      )

    ggplot2::ggsave("Figure_3_Local_Marrow_Crowding.png", plot = comp_scatterplot, width = 7, height = 5, dpi = 300)
    message("Figure 3 successfully generated and saved to your project directory!")

    plotting_df <- combined_cohorts_df %>%
      tidyr::pivot_longer(cols = c(bone_marrow_blast, peripheral_blood_myeloid), 
                          names_to = "Parameter", values_to = "Percentage")
    
    panel_labels <- c(
      "bone_marrow_blast" = "Bone Marrow Blast Infiltration (%)",
      "peripheral_blood_myeloid" = "Peripheral Blood Myeloid Cells (%)"
    )
    
    comp_boxplot <- ggplot2::ggplot(plotting_df, aes(x = Cohort, y = Percentage, fill = Cohort)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = 16, width = 0.5, color = "#2C3E50") +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.3, size = 1.5, aes(color = Cohort)) +
      ggplot2::facet_wrap(~Parameter, scales = "free_y", labeller = as_labeller(panel_labels)) +
      ggplot2::scale_fill_manual(values = c("#3498DB", "#E74C3C")) +    
      ggplot2::scale_color_manual(values = c("#2980B9", "#C0392B")) +  
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        strip.background = element_rect(fill = "#ECF0F1", color = "#BDC3C7"),
        strip.text = element_text(face = "bold", color = "#2C3E50", size = 11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 11),
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        panel.grid.minor = element_blank()
      ) +
      ggplot2::labs(
        title = "Cross-Population Bio-Audit: Global Reference vs. Local Sindh Cohort",
        y = "Tissue Compartment Density (%)"
      )
    
    ggplot2::ggsave("Figure_2_Cross_Population_Comparison.png", plot = comp_boxplot, width = 8, height = 5, dpi = 300)
    message("Figure 2 successfully generated and saved to your project directory!")
    
    # --------------------------------------------------------------------------
    # OPTIONAL MODULE: LOCAL LONGITUDINAL SURVIVAL CLINICAL VERIFICATION
    # --------------------------------------------------------------------------
    run_local_survival <- FALSE
    
    if("local_time_days" %in% colnames(local_data) & "local_status" %in% colnames(local_data)) {
      clean_survival_subset <- local_data %>%
        dplyr::select(patient_id, bone_marrow_blast, local_time_days, local_status) %>%
        dplyr::filter(!is.na(local_time_days) & !is.na(local_status) & !is.na(bone_marrow_blast))
      
      if(nrow(clean_survival_subset) >= 15) { 
        run_local_survival <- TRUE
      }
    }
    
    if(run_local_survival) {
      message("Longitudinal metrics confirmed. Executing Optional Local Survival validation module...")
      
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
      survminer:::print.ggsurvplot(local_survival_plot, newpage = FALSE)
      dev.off()
      message("SUCCESS: Figure 4 successfully generated and saved to your project directory!")
      
    } else {
      message("NOTE: Optional local survival metrics were missing or fell below validation row thresholds. Phase 3 exited safely without generating Figure 4.")
    }
    
  } else {
    message("LUMHS_Local_AML_Cohort.csv tracker template detected. Code is primed and waiting for data entry.")
  }
}

save.image(file = "LUMHS_Hypoxia_Project_Processed_Outputs.RData")
message("--- SUCCESS: Master Production Pipeline executed with zero errors. ---")
