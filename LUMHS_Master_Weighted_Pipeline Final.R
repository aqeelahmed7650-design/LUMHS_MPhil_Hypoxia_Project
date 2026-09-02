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
# Render graph on screen and automatically write to disk
png("Figure_1_Survival_Kinetics.png", width = 2400, height = 1800, res = 300)
print(final_plot)
dev.off()
message("Figure 1 successfully generated and saved to your project directory!")


# 15. Export Clean CSV Spreadsheet Table Data to Working Folder Disk
write.csv(lumhs_master_data, "LUMHS_Global_Survival_Validated_Data.csv", row.names = FALSE)

# 16. Lock and Save the Workspace Environment Binary File
save.image(file = "LUMHS_Hypoxia_Project_Workspace.RData")
message("--- SUCCESS: Master Pipeline executed with zero errors. RData and CSV saved. ---")
## =============================================================================
## PHASE 3: LOCAL TO GLOBAL CROSS-POPULATION COMPARISON (LUMHS VALIDATION)
## =============================================================================

# 1. Load your local clinical data spreadsheet (Target Sample: 30-50 patients)
if(file.exists("LUMHS_Local_AML_Cohort.csv")) {
  local_data <- read.csv("LUMHS_Local_AML_Cohort.csv", na.strings = "NA")
  
  # Ensure there is actual data logged past row 1 to prevent empty array crashes
  if(nrow(local_data) > 1) {
    
    # 2. Extract and format the matching global parameters from your active workspace
    global_comparison_df <- lumhs_master_data %>%
      dplyr::select(bone_marrow_blast, bone_marrow_neutrophil) %>%
      dplyr::mutate(Cohort = "Global (TCGA)")
    
    # 3. Format and isolate your local parameters
    local_comparison_df <- local_data %>%
      dplyr::select(bone_marrow_blast, bone_marrow_neutrophil) %>%
      dplyr::mutate(Cohort = "Local (LUMHS)")
    
    # 4. Bind both cohorts into a master comparative matrix
    combined_cohorts_df <- rbind(global_comparison_df, local_comparison_df)
    
    # 5. Run Independent T-Tests to check for regional baseline variations
    blast_t_test <- t.test(bone_marrow_blast ~ Cohort, data = combined_cohorts_df)
    neutrophil_t_test <- t.test(bone_marrow_neutrophil ~ Cohort, data = combined_cohorts_df)
    
    # Print statistical differences to the console for your Results Chapter
    print("--- CROSS-POPULATION BASELINE DIFFERENCES ---")
    print(blast_t_test)
    print(neutrophil_t_test)
    
    # 6. Run Local Non-Parametric Spearman Rank Correlations
    local_spearman <- cor.test(local_data$bone_marrow_blast, 
                               local_data$bone_marrow_neutrophil, 
                               method = "spearman", exact = FALSE)
    
    print("--- LOCAL SINDH COHORT CORRELATION RESULTS ---")
    print(local_spearman)
    comp_scatterplot <- ggplot(local_data, aes(x = bone_marrow_blast, y = bone_marrow_neutrophil)) +
  geom_point(color = "#E74C3C", alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", color = "#2C3E50", se = TRUE, fill = "#BDC3C7", alpha = 0.3) +
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
    y = "Bone Marrow Mature Neutrophils (%)"
  )

ggsave("Figure_3_Local_Marrow_Crowding.png", plot = comp_scatterplot, width = 7, height = 5, dpi = 300)
print("Figure 3 successfully generated and saved to your project directory!")

    # 7. Reshape data for publication-grade ggplot2 dual-panel visualization
    library(tidyr)
    library(ggplot2)
    
    plotting_df <- combined_cohorts_df %>%
      tidyr::pivot_longer(cols = c(bone_marrow_blast, bone_marrow_neutrophil), 
                          names_to = "Parameter", values_to = "Percentage")
    
    # Create clean panel labels
    panel_labels <- c(
      "bone_marrow_blast" = "Bone Marrow Blast Infiltration (%)",
      "bone_marrow_neutrophil" = "Bone Marrow Mature Neutrophils (%)"
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
    
    # Save the chart automatically
    ggsave("Figure_2_Cross_Population_Comparison.png", plot = comp_boxplot, 
           width = 8, height = 5, dpi = 300)
    print("Figure 2 successfully generated and saved to your project directory!")
    
  } else {
    print("LUMHS_Local_AML_Cohort.csv tracker template detected. Code is primed and waiting for 2027 data entry.")
  }
}
