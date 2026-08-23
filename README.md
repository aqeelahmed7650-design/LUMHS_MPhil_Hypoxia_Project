#  LUMHS MPhil Hypoxia Project: Computational Oncology Pipeline

This repository hosts the translational data engine and clinical tracking architecture designed for investigating bone marrow microenvironmental stress parameters in adult Acute Myeloid Leukemia (AML) cohorts at Liaquat University of Medical and Health Sciences (LUMHS).

##  Core Architecture Overview
* **`LUMHS_Master_Weighted_Pipeline.R`**: Multi-phase script executing univariate Cox Proportional Hazards regression scaling across 151 adult de novo datasets (TCGA-LAML) to optimize a 15-gene Buffa metagene hypoxia profile, yielding highly significant prognostic trajectories (**Log-rank p < 0.0001**).
* **`LUMHS_Local_AML_Cohort.csv`**: Pre-formatted, lowercase clinical template optimized to track local tissue densities (blasts, neutrophils, basophils, lymphocytes) across prospective regional samplings without data array mismatch risks.
* **`.gitignore`**: Operational configuration file shielding cloud sync runs from heavy data storage footprints.

##  Operational Focus
This computational framework bridges bulk multi-omic transcriptomic datasets with cost-effective, routine clinical bone marrow metrics, providing a validated predictive launchpad suitable for downstream integration with Digital Pathology and Artificial Intelligence whole-slide image analysis.
