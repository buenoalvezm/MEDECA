# Cancer Detection using Blood Proteome Profiling

![License](https://img.shields.io/badge/license-Apache2.0-blue)
![Version](https://img.shields.io/badge/version-1.0.0-green)

## Description

This project explores proximity extension assay-based proteomics for early cancer detection in patients with non-specific symptoms. We analyzed 1463 proteins in plasma from 456 patients prior to cancer diagnosis, identifying 29 proteins associated with cancer. A classification model was developed with an AUC of 0.80, achieving an AUC of 0.82 in an independent cohort. The model distinguishes cancer from autoimmune, inflammatory, and infectious diseases.

## Installation

Clone the repository:
   ```bash
   git clone https://github.com/buenoalvezm/MEDECA.git
   ```

# Usage
Run the scripts in the scripts folder for data processing and classification model generation for cancer detection.
- **01_data_preprocessing.Rmd**: quality control of Olink dataset
- **02_metadata_preprocessing.Rmd**: format MEDECA (discovery) and ALLVOS (replication) cohort metadata
- **03_differential_expression.Rmd**: perform differential expression analyses in the discovery and replication cohorts
- **04_cancer_classification.Rmd**: classification models to predict presence of cancer
- **05_cancer_type_analyses.Rmd**: investigate specific cancer types 

# License
This project is licensed under the Apache License 2.0. See the LICENSE file for more information.

# Version
Current Version: 1.0.0

# Contact
For any questions or inquiries, please contact María Bueno Álvez:
- Email: maria.bueno@scilifelab.se
- GitHub: buenoalvezm
