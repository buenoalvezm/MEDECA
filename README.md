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

Required hardware includes a computer with R, RStudio, and the necessary packages installed. Note: At least 4GB of RAM is recommended, especially for larger datasets.

Typical installation time on a standard desktop computer is less than 5 minutes.

# Usage

## Synthetic Data
The provided data is synthetic and randomly generated based on the original dataset's age, NPX, and group proportions. It retains some biological differences to ensure the code runs as intended. **No biological conclusions should be drawn from this data**, it is for code evaluation only and should not be used for research purposes.

## Running the Scripts
Scripts `01_data_preprocessing.Rmd` and `02_metadata_preprocessing.Rmd` do not need to be run when testing the code with the provided synthetic data. Instead, start with `03_differential_expression.Rmd` for differential expression analyses, followed by `04_cancer_classification.Rmd` for classification models, and `05_cancer_type_analyses.Rmd` to investigate specific cancer types.

For generating some figures, two external files are required: HPA data, which can be downloaded from proteinatlas.tsv, and pan-cancer markers, available as supplementary data from the publication [Next-generation pan-cancer blood proteome profiling using proximity extension assay](https://doi.org/10.1038/s41467-023-39765-y).


## Important Steps
Ensure the correct data paths are set in the scripts before running them. To test the code with the synthetic datasets, replace paths in scripts with:
   ```R
metadata <- read_csv("synthetic_data/combined_metadata_synthetic.csv")
meta_medeca <- read_csv("synthetic_data/medeca_metadata_synthetic.csv")
meta_allvos <- read_csv("synthetic_data/allvos_metadata_synthetic.csv")
data <- read_csv("synthetic_data/final_olink_data_synthetic.csv")
   ```
## Expected Outputs
The demo runtime on a standard desktop computer is approximately one hour. Outputs include plots and summary files of differential and machine learning results, stored in the results/ and data/processed/ directories.

# License
This project is licensed under the Apache License 2.0. See the LICENSE file for more information.

# Version
Current Version: 1.0.0

# Contact
For any questions or inquiries, please contact María Bueno Álvez:
- Email: maria.bueno@scilifelab.se
- GitHub: buenoalvezm
