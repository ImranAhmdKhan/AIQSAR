# AIQSAR — Comprehensive Statistical Analysis and QSAR Modeling App

## Overview

This R Shiny application provides a comprehensive, interactive, and user-friendly interface for end-to-end statistical analysis and Quantitative Structure-Activity Relationship (QSAR) modeling. It is designed to guide users from initial data exploration and preprocessing to advanced model building, rigorous validation, and the generation of publication-quality reports and visualizations.

The app consolidates numerous complex statistical and machine learning workflows into accessible modules, making it suitable for researchers, students, and data scientists in fields like cheminformatics, bioinformatics, and computational chemistry.

---

## Screenshots

### Data Explorer & Preprocessing
![Data Explorer tab showing file upload and missing-value heatmap](screenshots/01_data_explorer.png)

### Correlation Heatmap
![Interactive correlation heatmap of molecular descriptors](screenshots/02_correlation_heatmap.png)

### Full QSAR Pipeline Results
![Full QSAR Pipeline showing model comparison table and validation metrics](screenshots/03_qsar_pipeline.png)

### Applicability Domain — Williams Plot
![Williams Plot for applicability domain analysis](screenshots/04_williams_plot.png)

### AI Model Grader
![AI-powered model grader output with letter grade and recommendation](screenshots/05_model_grader.png)

### Analysis Dashboard
![Consolidated analysis dashboard with all key diagnostic plots](screenshots/06_analysis_dashboard.png)

> **Want to contribute screenshots?** Run the app locally, capture each tab, and open a pull request adding `.png` files to the [`screenshots/`](screenshots/) folder. See [`screenshots/README.md`](screenshots/README.md) for naming conventions.

---

## Key Features 🧪

- **Flexible Data Input**: Load data from CSV/Excel files or by pasting directly into the application.
- **Data Preprocessing**: A dedicated module for handling missing values (mean, median, k-NN imputation) and applying data transformations (Log, Square Root, Box-Cox).
- **Exploratory Analysis**: Generate descriptive statistics and visualize inter-variable relationships with a correlation heatmap.
- **Dimensionality Reduction & Clustering**: Perform Principal Component Analysis (PCA) and explore data structure with K-Means and Hierarchical Clustering.

### Automated Modeling Pipelines

- **Full QSAR Pipeline**: A one-click module that performs data splitting (Kennard-Stone), feature selection, model training (MLR, Random Forest, XGBoost), and validation.
- **Automated MLR**: A specialized pipeline for building robust Multiple Linear Regression models with a rigorous multi-step feature selection protocol.

### Advanced Machine Learning

Train and compare a suite of modern models including LightGBM, CatBoost, Cubist, and various Neural Networks (MLP, DNN). Leverage H2O's AutoML to automatically discover the best-performing models for your dataset.

### Rigorous Model Validation

- **Applicability Domain (AD)**: Generate and analyze Williams Plots to assess the reliability of model predictions.
- **Y-Scrambling**: Perform permutation testing to ensure the model is not based on chance correlations.
- **Validation Loop**: Test the stability of your model formula by retraining it on multiple random data splits.

### AI-Powered Insights

- **Model Grader**: Get an AI-driven analysis of your model comparison results, with a final grade and recommendation based on the performance vs. interpretability trade-off.
- **Model Interpretation**: Use DALEX-based methods (Break Down and Partial Dependence plots) to understand why your models make certain predictions.

### Publication-Ready Outputs

- **Analysis Dashboard**: A consolidated view of all key plots from your final analysis.
- **Reporting**: Generate and download comprehensive analysis logs, plot data, and final summary reports.

---

## Prerequisites & Installation

This application is built in R using the Shiny framework. To run it, you need to have **R** and **RStudio** installed on your system.

Install all required R packages by running the following command in your R console:

```r
install.packages(c(
  "shiny", "shinythemes", "DT", "readxl", "tidyverse", "corrplot", "pls",
  "glmnet", "prospectr", "car", "caret", "ggpubr", "randomForest", "e1071",
  "gbm", "shinyjs", "lmtest", "colourpicker", "xgboost", "plotly", "naniar",
  "tidymodels", "rmarkdown", "ggcorrplot", "ggthemes", "Ckmeans.1d.dp",
  "rpart", "DALEX", "DALEXtra", "promises", "future", "lightgbm", "catboost",
  "Cubist", "h2o", "RSNNS", "nortest", "moments"
))
```

> **Note:** Some packages, like `catboost`, may have specific system dependencies. Please refer to the official documentation for each package if you encounter installation issues.

---

## How to Run the Application

1. Save the R script as `AIQSAR_1.R` (or clone this repository).
2. Open `AIQSAR_1.R` in RStudio.
3. Click the **Run App** button at the top of the script editor.

Alternatively, run the app from the R console:

```r
shiny::runApp("AIQSAR_1.R")
```

---

## Workflow Guide

A typical analysis workflow using this tool:

1. **Data Explorer** — Upload your dataset (CSV/Excel) or paste it directly. Select your compound name column for labeling.
2. **Data Preprocessing** — Visualize missing data, apply imputation or transformation methods. The processed data carries through to all subsequent analyses.
3. **Automated Pipelines** — For a fast, comprehensive analysis, use the **Full QSAR Pipeline**: select your response variable, model type, and run to get a fully validated model with diagnostic plots.
4. **Manual Analysis** *(optional)* — Navigate tabs sequentially for a customized approach:
   - **Descriptive Analysis** — Understand data distribution and correlations.
   - **Dimensionality Reduction** — Explore structure with PCA and clustering.
   - **Modeling** — Build and compare models. The *Publication Analysis* tab is recommended for a final, robust MLR model.
5. **Advanced Validation** — Use Applicability Domain and Y-Scrambling to rigorously validate model robustness.
6. **AI Insights** — Use the AI-Powered Model Grader and Model Interpretation tabs for deeper performance and behavior insights.
7. **Reporting** — Use the Analysis Dashboard for a consolidated view, or the Final Report tab to generate a text-based summary.

---

## License

This project is open source.

> **Contribute:** We invite interested contributors to help improve AI-based QSAR modeling within this application. Pull requests and issues are welcome!
