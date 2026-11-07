# 15-Gene Hypoxia Model Validation & Cross-Population Bio-Audit Pipeline

### 🚀 Clinical Multi-Omics Data Mining & Translational Field Validation
**Author:** Dr. Aqeel Ahmed, MD (Postgraduate M.Phil Scholar)  
**Institution:** Liaquat University of Medical and Health Sciences (LUMHS), Sindh, Pakistan  
**Target Repository Status:** Production / Clinical Validation Track  
**Regulatory Baseline:** ICH-GCP E6(R3) Compliant Framework  

---

## 📋 Project Overview

This repository houses the complete, reproducible four-part R analytics pipeline developed to validate the prognostic stability of high-dimensional transcriptomic signatures and bridge them with low-cost clinical surrogates. 

The framework employs a **"Genomic Discovery + Pragmatic Clinical Translation"** model:
1. **Global Discovery:** It extracts raw bulk RNA-seq data from the global **TCGA-LAML** registry ($n=151$), normalizes gene vectors, and uses outcome-weighted Cox regression to map microenvironmental tissue hypoxia.
2. **Local Validation:** It cross-audits these high-level signatures against prospective, real-world field metrics ($n=45$) collected at a public tertiary center (**LUMHS, interior Sindh**), explicitly mapping how endemic systemic comorbidities (multifactorial nutritional iron-deficiency anemia) accelerate tissue ischemia.

---

## 🛠️ Pipeline Architecture & Execution Flow

The analytics engine is segmented into four sequential, error-checked execution blocks to protect script runtime integrity:

### 🔹 Part 1: Environment Setup, Multi-Omics Mining, & Clinical Filtering
*   Configures global repositories to the stable Cloud CRAN mirror to bypass institutional network access constraints.
*   Automates installation and background loading of core Bioconductor infrastructures (`curatedTCGAData`, `MultiAssayExperiment`, `TCGAutils`).
*   Extracts pre-normalized global sequencing expressions (`RNASeq2GeneNorm`) and harmonizes case barcodes.
*   Pivots to an accessible **Total White Blood Cell (WBC) Count ($x10^9/\text{L}$)** tracking array to optimize the pipeline for low-resource global health clinical deployment.

### 🔹 Part 2: Genomic Normalization & Weighting Arrays
*   Isolates expression tracks for the validated **15-gene Buffa Hypoxia Signature** (`ACOT7`, `ADM`, `ALDOA`, `CDKN3`, `ENO1`, `LDHA`, `MIF`, `MRPS17`, `NDRG1`, `P4HA1`, `PGAM1`, `SLC2A1`, `TPI1`, `TUBB6`, `VEGFA`).
*   Executes $\log_2$ transformations and row-wise, cross-patient Z-score normalization arrays.
*   Runs independent univariate Cox Proportional Hazards loops for each locus against clinical survival timelines to extract directional log-hazard weights ($\beta$ coefficients).
*   Evaluates mathematical model stability by checking the Proportional Hazards assumption via automated Schoenfeld residual diagnostics (`cox.zph`).

### 🔹 Part 3: Local Integration, Predictive Proxy Modeling, & Cohort Analytics
*   Aggregates patient-specific expression matrices into a final **Weighted Hypoxia Risk Score**:
    $$\text{Hypoxia Score}_j = \sum (\text{Z-score of Gene}_i \times \beta_i)$$
*   **Niche Proxy Core:** Builds an internal **Multivariable Logistic Regression Model** (`glm(..., family = binomial)`) within TCGA to mathematically evaluate how accurately routine clinical indicators (Hb, WBC, Blasts) reflect a deep genomic "High Hypoxia" state.
*   Merges global reference tables with the local `LUMHS_Local_AML_Cohort.csv` data.
*   Applies a strict, automated **14-day transfusion-naive filter masking block** to screen out recent donor blood artifacts from host baseline Hemoglobin ($Hb$) tracking.

### 🔹 Part 4: Advanced Visualizations & Algorithmic Safety Gates
*   Executes **non-parametric Mann-Whitney U tests** (`wilcox.test`) to evaluate population distribution shifts, safely managing skewed leukemic profiles.
*   Maps intramarrow vs. peripheral tumor density tracking using non-linear **Spearman Rank Correlation matrices** combined with local regression smoothers (**LOESS**).
*   Generates high-resolution, publication-ready data plots (`survminer::ggsurvplot`) with hardcoded time-split ticks (`break.time.by = 500`) to guarantee perfect risk-table axis alignment.
*   **Built-in Safety Gate:** If patient attrition or clinical tracking limitations drop active Tier 1 follow-up rows below a critical baseline ($n < 15$), the pipeline automatically bypasses long-term survival curves to avoid code crashes, outputting the study exclusively as a stable, cross-sectional descriptive bio-audit (`Table_1_Clinical_Summary.csv`).

---

## 📊 Data Schema Specification (`LUMHS_Local_AML_Cohort.csv`)

To execute the local clinical cross-population modules, the prospective field registry file must reside in the working directory utilizing this exact flat-file database schema:

| Column Name | Data Type | Field Constraints / Description |
| :--- | :--- | :--- |
| `patient_id` | Character | Primary key mapping from `LUMHS_AML_001` onwards. |
| `bone_marrow_blast` | Numeric | Manual bone marrow aspirate blast percentage (Must be $\geq 20.0\%$). |
| `peripheral_blood_myeloid` | Numeric | Absolute Total White Blood Cell (WBC) count ($x10^9/\text{L}$). |
| `hemoglobin_level` | Numeric | Diagnostic baseline Hemoglobin level ($g/\text{dL}$) from initial blood draw. |
| `transfused_within_14_days` | Character | Binary flag (`yes` or `no` in lowercase) for automated software transfusion masking. |
| `local_time_days` | Numeric | Longitudinal observational tracking timeline from Time-Zero ($T_0$). |
| `local_status` | Integer | Binary event flag ($1 = \text{deceased event}$, $0 = \text{censored / alive}$). |

---

## 🏃 Testing with Mock Data Simulation

A standalone, bug-fixed validation script (`LUMHS_Mock_Data_Generator.R`) is included in the project folder to test your local environment without live patient files. 
To validate your local environment, run:

```r
# 1. Run the generator script to create the local file framework
source("LUMHS_Mock_Data_Generator.R")

# 2. Execute pipeline components sequentially
source("Part1_Environment_Setup.R")
source("Part2_Genomic_Weighting.R")
source("Part3_Local_Integration.R")
source("Part4_Visualizations_Export.R")
```

---

## 📂 Expected Output Deliverables

Upon successful pipeline completion, the engine exports these five publication-ready files directly to your active working directory:
*   `Figure_1_Survival_Kinetics.png`: High-resolution Cox-weighted global survival timeline curves.
*   `Figure_2_Cross_Population_Comparison.png`: Facet-wrapped cross-population distribution boxplots.
*   `Figure_3_Local_Marrow_Crowding.png`: LOESS-smoothed scatter matrix tracking compartment crowding.

"Data Status Note: The file LUMHS_Local_AML_Cohort.csv currently contains the production-ready database schema headers. Prospective clinical field enrollment is actively running at the Liaquat University of Medical and Health Sciences (LUMHS) hospital counters. The file will dynamically populate with real-world patient matrices upon completion of the prospective validation track."

*   
*   `Figure_4_Local_Survival_Kinetics.png`: Local cohort internal Kaplan-Meier validation kinetics (if $n \geq 15$).
*   `Table_1_Clinical_Summary.csv`: Excel-readable aggregated descriptive clinical snapshot matrix.

---
*Disclaimer: This analytics pipeline is intended strictly for multi-omics transcriptomic validation, translational biostatistics research, and peer-reviewed publication frameworks. It does not constitute direct medical advice or clinical diagnosis tools.*
