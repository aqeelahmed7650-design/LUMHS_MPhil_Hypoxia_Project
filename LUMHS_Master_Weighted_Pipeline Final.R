# ==============================================================================
# MASTER COMPILATION PIPELINE: COX-WEIGHTED 15-GENE BUFFA PROGNOSTIC MODEL
# Baseline Global Dataset Phase - LUMHS MPhil Research Initiative
# Principal Investigator: Dr. Aqeel Ahmed
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. AUTOMATED BOOTSTRAPPING FOR CRAN AND BIOCONDUCTOR INFRASTRUCTURES
# ------------------------------------------------------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))

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

# ------------------------------------------------------------------------------
# 2. DYNAMIC PROGRAMMATIC FETCH OF RAW GLOBAL GENOMIC & CLINICAL COHORT
# ------------------------------------------------------------------------------
message("Downloading raw TCGA-LAML multi-assay experiment data directly from API...")
laml_mae <- curatedTCGAData(disease = "LAML", assays = c("RNASeq2GeneNorm"), version = "2.0.1", dry.run = FALSE)

raw_assay <- assays(laml_mae)[[1]]
aml_matrix_raw <- as.data.frame(raw_assay)
aml_matrix_raw$gene_name <- rownames(raw_assay)

clinical_data_raw <- as.data.frame(colData(laml_mae))

# ------------------------------------------------------------------------------
# 3. INITIALIZE THE BUFFA 15-GENE HYPOXIA CLASSIFIER VECTOR
# ------------------------------------------------------------------------------
buffa_15_genes <- c("ACOT7", "ADM", "ALDOA", "CDKN3", "ENO1", "LDHA", "MIF", 
                    "MRPS17", "NDRG1", "P4HA1", "PGAM1", "SLC2A1", "TPI1", "TUBB6", "VEGFA")

# FIXED: Strip out any Entrez ID vertical bar extensions (e.g., converting "ACOT7|641" to "ACOT7")
# to guarantee alignment with the 15-gene signature vector.
aml_matrix_raw$gene_name <- sub("\\|.*", "", rownames(raw_assay))

b15_matrix <- aml_matrix_raw[aml_matrix_raw$gene_name %in% buffa_15_genes, ]
message(paste("Genomic Alignment Success: Isolated", nrow(b15_matrix), "target features out of 15 Buffa genes."))

# ------------------------------------------------------------------------------
# 4. EXTRACT AND CLEAN DIRECT CLINICAL SURVIVAL TIMELINES
# ------------------------------------------------------------------------------
survival_base <- clinical_data_raw %>%
  # FIXED: Removed 'stringr::' namespace tag since toupper is a standard base R function
  dplyr::mutate(patient_id = toupper(patientID)) %>%
  dplyr::select(patient_id, days_to_death, days_to_last_followup, vital_status) %>%
  dplyr::mutate(
    time_days = ifelse(!is.na(days_to_death), as.numeric(days_to_death), as.numeric(days_to_last_followup)),
    status_numeric = ifelse(vital_status == 1 | grepl("dead|deceased", vital_status, ignore.case = TRUE), 1, 0)
  ) %>%
  dplyr::filter(!is.na(time_days) & time_days > 0)

# FIXED: Peer-Review Safeguard. Dynamically detect columns by text patterns 
# to protect against varying Bioconductor/Firehose API column renamings.
blast_col_match <- colnames(clinical_data_raw)[grep("blast", colnames(clinical_data_raw), ignore.case = TRUE)]
myeloid_col_match <- colnames(clinical_data_raw)[grep("myeloid|peripheral", colnames(clinical_data_raw), ignore.case = TRUE)]

# Fallback defaults if matching yields an unexpected empty index
if(is.na(blast_col_match[1])) blast_col_match <- "percent_bone_marrow_blasts"
if(is.na(myeloid_col_match[1])) myeloid_col_match <- "percent_myeloid_cells_peripheral_blood"

aml_clean_clinical <- clinical_data_raw %>%
  dplyr::mutate(
    patient_id = toupper(patientID),
    # Uses programmatic indexing (.[[]]) to pull columns flexibly by name string
    bone_marrow_blast = as.numeric(.[[blast_col_match[1]]]),
    peripheral_blood_myeloid = as.numeric(.[[myeloid_col_match[1]]]) 
  ) %>%
  dplyr::filter(!is.na(bone_marrow_blast) & !is.na(peripheral_blood_myeloid)) %>%
  dplyr::select(patient_id, bone_marrow_blast, peripheral_blood_myeloid)
# ------------------------------------------------------------------------------
# 5. LONG TRANSFORMATION AND CONVERSION VIA PRIMARY SOLID/BLOOD ALIGNMENT
# ------------------------------------------------------------------------------
b15_long <- b15_matrix %>%
  tidyr::pivot_longer(
    cols = -gene_name, 
    names_to = "raw_header", 
    values_to = "expression_value"
  ) %>%
  # FIXED: Removed structural anchor constraint to dynamically extract the TCGA barcode from anywhere in the text string
  dplyr::mutate(sample_barcode = stringr::str_extract(raw_header, "TCGA[\\.-][A-Z0-9]{2}[\\.-][A-Z0-9]{4}[\\.-][0-9]{2}[A-Z]")) %>%
  dplyr::filter(!is.na(sample_barcode)) %>%
  dplyr::mutate(sample_barcode = stringr::str_replace_all(sample_barcode, "\\.", "-")) %>%
  dplyr::mutate(patient_id = substr(sample_barcode, 1, 12)) %>%
  dplyr::mutate(log_val = log2(expression_value + 1)) %>%
  dplyr::group_by(gene_name) %>%
  dplyr::mutate(z_score = as.numeric(scale(log_val))) %>%
  dplyr::ungroup()

# ------------------------------------------------------------------------------
# 6. WIDE ANALYTICAL FORMATTING & COHORT INTERSECTION
# ------------------------------------------------------------------------------
wide_expression <- b15_long %>%
  select(patient_id, gene_name, z_score) %>%
  group_by(patient_id, gene_name) %>%
  summarize(z_score = mean(z_score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = gene_name, values_from = z_score)

cox_train_df <- inner_join(survival_base, wide_expression, by = "patient_id")

# ------------------------------------------------------------------------------
# 7. DYNAMIC UNIVARIATE COX MODEL GENERATION
# ------------------------------------------------------------------------------
gene_weights <- c()
available_genes <- intersect(buffa_15_genes, colnames(wide_expression))

for(gene in available_genes) {
  formula_string <- paste("Surv(time_days, status_numeric) ~", paste0("`", gene, "`"))
  cox_model <- coxph(as.formula(formula_string), data = cox_train_df)
  gene_weights[gene] <- coef(cox_model)
}

# ------------------------------------------------------------------------------
# 8. RISK CALCULATIONS AND STRATIFICATION
# ------------------------------------------------------------------------------
weighted_scores <- cox_train_df %>%
  rowwise() %>%
  mutate(weighted_hypoxia_score = sum(c_across(all_of(available_genes)) * gene_weights, na.rm = TRUE)) %>%
  ungroup() %>%
  select(patient_id, weighted_hypoxia_score)

lumhs_master_data <- inner_join(aml_clean_clinical, weighted_scores, by = "patient_id") %>%
  inner_join(survival_base, by = "patient_id")

median_weighted <- median(lumhs_master_data$weighted_hypoxia_score, na.rm = TRUE)
lumhs_master_data$hypoxia_group <- ifelse(lumhs_master_data$weighted_hypoxia_score > median_weighted, 
                                          "High Risk Hypoxia", "Low Risk Hypoxia")

final_fit <- survfit(Surv(time_days, status_numeric) ~ hypoxia_group, data = lumhs_master_data)

# ------------------------------------------------------------------------------
# 9. GRAPHICAL RENDERING & TABULAR OUTPUT STORAGE
# ------------------------------------------------------------------------------
final_plot <- ggsurvplot(
  final_fit, data = lumhs_master_data, pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, palette = c("#E41A1C", "#377EB8"), 
  legend.labs = c("High Risk Hypoxia", "Low Risk Hypoxia"),
  xlab = "Survival Time (Days)", ylab = "Survival Probability",
  title = "Cox-Weighted 15-Gene Buffa Hypoxia Prognostic Model",
  risk.table.height = 0.22, tables.theme = theme_cleantable()
)

png("Figure_1_Survival_Kinetics.png", width = 2400, height = 1800, res = 300)
survminer:::print.ggsurvplot(final_plot, newpage = FALSE)
dev.off()
message("Figure 1 successfully generated and saved to your project directory!")

write.csv(lumhs_master_data, "LUMHS_Global_Survival_Validated_Data.csv", row.names = FALSE)
save.image(file = "LUMHS_Hypoxia_Project_Processed_Outputs.RData")
message("--- SUCCESS: Master Pipeline executed with zero errors. ---")

# ==============================================================================
# PHASE 3: LOCAL TO GLOBAL CROSS-POPULATION COMPARISON (LUMHS VALIDATION)
# ==============================================================================

# 1. Load local clinical spreadsheet with robust missing value strings
if(file.exists("LUMHS_Local_AML_Cohort.csv")) {
  # Insulates against blank entries, dot markers, and standard NA text configurations
  local_data <- read.csv("LUMHS_Local_AML_Cohort.csv", na.strings = c("NA", " ", "", "."))
  
  # Force lowercase header mapping to prevent data typing and case mismatches
  colnames(local_data) <- tolower(colnames(local_data))
  
  # FIXED: Aligned target metrics with the true biological compartment definitions
  required_columns <- c("bone_marrow_blast", "peripheral_blood_myeloid")
  missing_cols <- setdiff(required_columns, colnames(local_data))
  
  if (length(missing_cols) > 0) {
    stop(paste("CRITICAL ERROR: The template is missing columns:", paste(missing_cols, collapse = ", ")))
  }

  # Ensure there is actual data logged past row 1 to prevent empty array crashes
  if(nrow(local_data) > 1) {
    
    # 2. Isolate and clean target metrics from the global cohort
    global_comparison_df <- lumhs_master_data %>%
      dplyr::select(bone_marrow_blast, peripheral_blood_myeloid) %>%
      dplyr::mutate(Cohort = "Global (TCGA)") %>%
      # Defensive Filter: Ensure comparison values are strictly numeric and complete
      dplyr::filter(!is.na(bone_marrow_blast) & !is.na(peripheral_blood_myeloid))
    
    # 3. Isolate and clean target metrics from your local cohort
    local_comparison_df <- local_data %>%
      dplyr::select(bone_marrow_blast, peripheral_blood_myeloid) %>%
      dplyr::mutate(Cohort = "Local (LUMHS)") %>%
      dplyr::filter(!is.na(bone_marrow_blast) & !is.na(peripheral_blood_myeloid))
    
    # 4. Bind both cohorts into a master comparative matrix
    combined_cohorts_df <- rbind(global_comparison_df, local_comparison_df)
    
    # 5. Run Independent T-Tests to check for regional baseline variations
    blast_t_test <- t.test(bone_marrow_blast ~ Cohort, data = combined_cohorts_df)
    myeloid_t_test <- t.test(peripheral_blood_myeloid ~ Cohort, data = combined_cohorts_df)
    
    # Print statistical differences to the console for your Results Chapter
    print("--- CROSS-POPULATION BASELINE DIFFERENCES ---")
    print(blast_t_test)
    print(myeloid_t_test)
    
    # 6. Run Local Non-Parametric Spearman Rank Correlations
    local_spearman <- cor.test(local_comparison_df$bone_marrow_blast, 
                               local_comparison_df$peripheral_blood_myeloid, 
                               method = "spearman", exact = FALSE)
    
    print("--- LOCAL SINDH COHORT CORRELATION RESULTS ---")
    print(local_spearman)
    
    comp_scatterplot <- ggplot(local_comparison_df, aes(x = bone_marrow_blast, y = peripheral_blood_myeloid)) +
      geom_point(color = "#E74C3C", alpha = 0.6, size = 2.5) +
      geom_smooth(method = "loess", formula = y ~ x, color = "#2C3E50", se = TRUE, fill = "#BDC3C7", alpha = 0.3) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"),
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = "Local Marrow Crowding Dynamics",
        subtitle = paste0("Spearman r = ", round(local_spearman$estimate, 2), " (p = ", format.pval(local_spearman$p.value, digits = 4), ")"),
        x = "Bone Marrow Blast Infiltration (%)",
        y = "Peripheral Blood Myeloid Cells (%)"
      )

    ggsave("Figure_3_Local_Marrow_Crowding.png", plot = comp_scatterplot, width = 7, height = 5, dpi = 300)
    message("Figure 3 successfully generated and saved to your project directory!")

    # 7. Reshape data for publication-grade ggplot2 dual-panel visualization
    plotting_df <- combined_cohorts_df %>%
      tidyr::pivot_longer(cols = c(bone_marrow_blast, peripheral_blood_myeloid), 
                          names_to = "Parameter", values_to = "Percentage")
    
    # FIXED: Updated clean panel text labels to accurately reflect peripheral blood metrics
    panel_labels <- c(
      "bone_marrow_blast" = "Bone Marrow Blast Infiltration (%)",
      "peripheral_blood_myeloid" = "Peripheral Blood Myeloid Cells (%)"
    )
    
    # Generate Figure 2 Boxplot
    comp_boxplot <- ggplot(plotting_df, aes(x = Cohort, y = Percentage, fill = Cohort)) +
      geom_boxplot(alpha = 0.7, outlier.shape = 16, width = 0.5, color = "#2C3E50") +
      geom_jitter(width = 0.15, alpha = 0.3, size = 1.5, aes(color = Cohort)) +
      facet_wrap(~Parameter, scales = "free_y", labeller = as_labeller(panel_labels)) +
      scale_fill_manual(values = c("#3498DB", "#E74C3C")) +    
      scale_color_manual(values = c("#2980B9", "#C0392B")) +  
      theme_bw(base_size = 12) +
      theme(
        strip.background = element_rect(fill = "#ECF0F1", color = "#BDC3C7"),
        strip.text = element_text(face = "bold", color = "#2C3E50", size = 11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 11),
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = "Cross-Population Bio-Audit: Global Reference vs. Local Sindh Cohort",
        y = "Tissue Compartment Density (%)"
      )
    
    ggsave("Figure_2_Cross_Population_Comparison.png", plot = comp_boxplot, 
           width = 8, height = 5, dpi = 300)
    message("Figure 2 successfully generated and saved to your project directory!")
    
  } else {
    message("LUMHS_Local_AML_Cohort.csv tracker template detected. Code is primed and waiting for data entry.")
  }
}
    # =========================================================================
    # OPTIONAL MODULE: LOCAL LONGITUDINAL SURVIVAL CLINICAL VERIFICATION
    # =========================================================================
    # System checks if optional tracking columns contain actual numerical events
    run_local_survival <- FALSE
    
    if("local_time_days" %in% colnames(local_data) & "local_status" %in% colnames(local_data)) {
      clean_survival_subset <- local_data %>%
        dplyr::select(patient_id, bone_marrow_blast, local_time_days, local_status) %>%
        dplyr::filter(!is.na(local_time_days) & !is.na(local_status) & !is.na(bone_marrow_blast))
      
      if(nrow(clean_survival_subset) >= 15) { # Requires a safe absolute minimum tracking threshold
        run_local_survival <- TRUE
      }
    }
    
    if(run_local_survival) {
      message("Longitudinal metrics confirmed. Executing Optional Local Survival validation module...")
      
      # Risk stratify your local cohort based on their marrow blast cell density median
      median_blast <- median(clean_survival_subset$bone_marrow_blast, na.rm = TRUE)
      clean_survival_subset$blast_group <- ifelse(clean_survival_subset$bone_marrow_blast > median_blast, 
                                                   "High Blast Load", "Low Blast Load")
      
      # Fit Kaplan-Meier parameters to local hospital cohort tracking metrics
      local_fit <- survival::survfit(survival::Surv(local_time_days, local_status) ~ blast_group, 
                                     data = clean_survival_subset)
      
      # Render the optional Figure 4 chart
      local_survival_plot <- survminer::ggsurvplot(
        local_fit, data = clean_survival_subset, pval = TRUE, conf.int = FALSE,
        risk.table = TRUE, palette = c("#D35400", "#27AE60"),
        legend.labs = c("High Blast Load", "Low Blast Load"),
        xlab = "Follow-up Duration (Days)", ylab = "Survival Probability",
        title = "LUMHS Local Cohort Internal Survival Kinetic Kinetics",
        risk.table.height = 0.22, tables.theme = survminer::theme_cleantable()
      )
      
      png("Figure_4_Local_Survival_Kinetics.png", width = 2400, height = 1800, res = 300)
      survminer:::print.ggsurvplot(local_survival_plot, newpage = FALSE)
      dev.off()
      message("SUCCESS: Figure 4 successfully generated and saved to your project directory!")
      
    } else {
      message("NOTE: Optional local survival metrics were missing or fell below validation row thresholds. Phase 3 exited safely without generating Figure 4.")
    }
