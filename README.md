# Cox-Weighted 15-Gene Buffa Prognostic Hypoxia Model in AML
### Baseline Global Dataset Phase & Local Validation Cohort Bio-Audit

**Principal Investigator:** Dr. Aqeel Ahmed  
**Institution:** Liaquat University of Medical & Health Sciences (LUMHS), MPhil Research Initiative  
**Target Publication:** *Scientific Reports* (Nature Portfolio)  
**Timeline Configuration:** January 2027 Baseline Deploy  

---

## 🔬 Project Overview
This repository contains the complete computational reproducibility pipeline for evaluating the prognostic strength of the universally validated **15-Gene Buffa Hypoxia Classifier Vector** ("ACOT7", "ADM", "ALDOA", "CDKN3", "ENO1", "LDHA", "MIF", "MRPS17", "NDRG1", "P4HA1", "PGAM1", "SLC2A1", "TPI1", "TUBB6", "VEGFA") inside Acute Myeloid Leukemia (AML) cohorts. 

The pipeline dynamically harvests global molecular-clinical matrices from the **TCGA-LAML** study via the Bioconductor API to train an optimized Cox Proportional Hazards regression score model, and bridges these metrics against a local validation cohort in Sindh to audit regional marrow crowding dynamics.

## Repository File Blueprint
* **`LUHMS_Master_Weighted_Pipeline.R`**: The unified, self-contained multi-phase production R script. Handles package bootstrapping, data ingestion, feature scaling, model training, plot rendering, and local cohort t-testing.
* **`LUMHS_Local_AML_Cohort.csv`**: The formatted clinical tracking template required to execute the cross-population comparison module.

## Execution Instructions
To execute this analysis pipeline and reproduce the publication figures, clone this repository to your local machine and run the core script within an established R Project environment:

```r
# Execute the entire master pipeline dynamically
source("LUHMS_Master_Weighted_Pipeline.R")
```

### Generated Publication Outputs
Upon flawless completion, the script automatically exports the following files to your working directory:
1. **`Figure_1_Survival_Kinetics.png`**: High-resolution Kaplan-Meier curves with cross-stratified log-rank risk boundaries and patient risk-table intervals.
2. **`Figure_2_Cross_Population_Comparison.png`**: Comparative dual-panel boxplots illustrating tissue compartment densities. *(Requires Local CSV Data)*
3. **`Figure_3_Local_Marrow_Crowding.png`**: Non-parametric LOESS smooth trend line correlation plot mapping bone marrow blast infiltration against peripheral blood myeloid distributions. *(Requires Local CSV Data)*
4. **`LUMHS_Global_Survival_Validated_Data.csv`**: Master tabular dataset combining computed patient-specific risk scores with matched pathobiology markers.

