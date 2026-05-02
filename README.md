# \# Metabolomics Data Analysis Pipeline 🧪

# 

# A modular, reproducible R pipeline for preprocessing, multivariate analysis, feature selection, pathway enrichment, and quality assessment of LC-MS metabolomics data.

# 

# \## 📁 Repository Structure

metabolomics-preprocessing/

├── README.md

├── .gitignore

├── analysis\_scripts

│ ├── 00\_setup\_environment.R # Install/load all required packages

│ ├── 01\_preprocessing\_data.R # Filtering, imputation, PQN, log2

│ ├── 02\_pca\_analysis.R # Pareto scaling \& PCA

│ ├── 03\_plsda\_analysis.R # PLS-DA, LOOCV, permutation test, VIP

│ ├── 04\_loocv\_validation.R # Leave-one-out CV on important metabolites

│ ├── 05\_bootstrap\_stability\_analysis.R # Bootstrap stability of VIP-selected metabolites

│ ├── 06\_random\_forest.R # Random Forest feature selection

│ ├── 07\_univariate\_analysis.R # t-test/Wilcoxon, effect size, fold change

│ ├── 08\_pathway\_analysis.R # Pathway enrichment bubble plot

│ ├── 09\_name\_mapping.R # Map safe R names → original metabolite names

│ ├── 10\_volcano\_plot.R # Publication‑ready volcano plots

│ └── 11\_quality\_control.R # QC sample PCA, outlier detection, RSD





\## 🚀 How to Run



1\. \*\*Clone or download\*\* this repository and open it in RStudio (or run from the terminal).

2\. \*\*Place your raw data\*\* in `data/raw\\\\\\\_metabolomics\\\\\\\_data.csv`.

&#x20;  - The file should contain columns `Sample\\\\\\\_ID`, `Group` (levels: `Fit-Good`, `Fit-Poor`), and one column per metabolite.

3\. \*\*Run the scripts in numerical order\*\* (00 → 01 → … → 11).

&#x20;  - Each script saves its outputs in the `output/` folder, which are then read by subsequent scripts.

&#x20;  - Script `11` (QC) requires `data/QC\\\\\\\_metabolites.csv` and can be run independently.



> \\\\\\\*\\\\\\\*Tip:\\\\\\\*\\\\\\\* You can source each script with `source("R/00\\\\\\\_setup\\\\\\\_environment.R")`, etc.



\## 📦 Dependencies



All packages are centrally managed by the `pacman` package.

Script `01` installs (if missing) and loads:



`dplyr`, `tidyr`, `readr`, `ggplot2`, `ggrepel`, `mixOmics`, `randomForest`,

`effsize`, `pROC`, `patchwork`, `viridis`, `cowplot`, `scales`, `RColorBrewer`,

`imputeLCMD`, `nortest`, `MASS`, `gridExtra`, `stringr`, `magrittr`, `tibble`, `tidyverse`



\## 📈 Analysis Workflow



| Step | Script | Description |

|------|--------|-------------|

| Setup | `00` | Loads all required R packages |

| Preprocessing | `01` | Missing value filtering (≥50% presence), half‑minimum imputation, PQN normalisation, log2 transformation |

| Exploratory | `02` | Pareto scaling + PCA score plots (basic \& journal‑style) |

| Supervised | `03` | PLS‑DA with manual LOOCV, VIP scores, permutation test (n=1000) |

| Validation | `04` | LOOCV using only important metabolites (VIP>1.5), confusion matrix, 95% CI |

| Stability | `05` | Bootstrap analysis (×100) to assess selection frequency of VIP metabolites |

| Feature Selection | `06` | Random Forest (ntree=1000) on VIP metabolites, integrated RF‑Bootstrap ranking |

| Univariate | `07` | t‑test/Wilcoxon, Cliff’s delta, fold change, FDR correction |

| Visualisation | `08` | Pathway enrichment bubble plot |

| Names | `09` | Converts safe R column names to original, human‑readable metabolite names |

| Quality Control | `10` | Optimised volcano plots with smart label placement |

| Final Figure | `11` | PCA of QC samples, outlier detection, RSD distribution |



\## 📊 Key Outputs



All results are saved in `output/` and `output/figures/`.

Important files include:



\- `preprocessed\\\\\\\_data.csv` – cleaned, normalised, log2‑transformed data

\- `preprocessed\\\\\\\_data\\\\\\\_scaled.csv` – Pareto‑scaled data for multivariate analysis

\- `VIP\\\\\\\_\\\\\\\*.csv` – VIP scores (all, important, top 50)

\- `PLSDA\\\\\\\_performance\\\\\\\_metrics.csv` – LOOCV accuracy, optimal ncomp, permutation p‑value

\- `LOOCV\\\\\\\_final\\\\\\\_report.csv` – final model accuracy and 95% CI

\- `bootstrap\\\\\\\_stability\\\\\\\_results.csv` – metabolites with selection frequencies and confidence levels

\- `RandomForest\\\\\\\_importance.csv` – RF importance scores

\- `univariate\\\\\\\_all\\\\\\\_results.csv` – complete univariate test results (p‑values, effect sizes, fold changes)

\- `univariate\\\\\\\_significant\\\\\\\_metabolites.csv` – significantly differential metabolites

\- `\\\\\\\*\\\\\\\_Mapped.csv` – versions of the above with original metabolite names

\- `QC\\\\\\\_assessment\\\\\\\_summary.csv` – QC metrics (RSD, outliers)



\## ⚙️ Customization



\- \*\*Group names\*\*: Change `Fit-Good` / `Fit-Poor` in the scripts if your design differs.

\- \*\*Missing value threshold\*\*: Modify `threshold` in `01`.

\- \*\*VIP threshold\*\*: Modify `vip\\\\\\\_threshold` in `04` and downstream scripts.

\- \*\*Number of bootstrap iterations\*\*: Adjust `n\\\\\\\_bootstrap` in `05`.

\- \*\*Label count in volcano\*\*: Change `top\\\\\\\_n\\\\\\\_labels` in `11`.



\## 📄 License



This project is provided under the \[MIT License](LICENSE) (or your preferred license).



\## ✉️ Contact



For questions or suggestions, please open an issue or contact the repository maintainer.

