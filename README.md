Comprehensive Statistical Analysis and QSAR Modeling App
Overview
This R Shiny application provides a comprehensive, interactive, and user-friendly interface for end-to-end statistical analysis and Quantitative Structure-Activity Relationship (QSAR) modeling. It is designed to guide users from initial data exploration and preprocessing to advanced model building, rigorous validation, and the generation of publication-quality reports and visualizations.

The app consolidates numerous complex statistical and machine learning workflows into accessible modules, making it suitable for researchers, students, and data scientists in fields like cheminformatics, bioinformatics, and computational chemistry.

Key Features 🧪
Flexible Data Input: Load data from CSV/Excel files or by pasting directly into the application.

Data Preprocessing: A dedicated module for handling missing values (mean, median, k-NN imputation) and applying data transformations (Log, Square Root, Box-Cox).

Exploratory Analysis: Generate descriptive statistics and visualize inter-variable relationships with a correlation heatmap.

Dimensionality Reduction & Clustering: Perform Principal Component Analysis (PCA) and explore data structure with K-Means and Hierarchical Clustering.

Automated Modeling Pipelines:

Full QSAR Pipeline: A one-click module that performs data splitting (Kennard-Stone), feature selection, model training (MLR, Random Forest, XGBoost), and validation.

Automated MLR: A specialized pipeline for building robust Multiple Linear Regression models with a rigorous multi-step feature selection protocol.

Advanced Machine Learning:

Train and compare a suite of modern models including LightGBM, CatBoost, Cubist, and various Neural Networks (MLP, DNN).

Leverage H2O's AutoML to automatically discover the best-performing models for your dataset.

Rigorous Model Validation:

Applicability Domain (AD): Generate and analyze Williams Plots to assess the reliability of model predictions.

Y-Scrambling: Perform permutation testing to ensure the model is not based on chance correlations.

Validation Loop: Test the stability of your model formula by retraining it on multiple random data splits.

AI-Powered Insights:

Model Grader: Get an AI-driven analysis of your model comparison results, with a final grade and recommendation based on the performance vs. interpretability trade-off.

Model Interpretation: Use DALEX-based methods (Break Down and Partial Dependence plots) to understand why your models make certain predictions.

Publication-Ready Outputs:

Analysis Dashboard: A consolidated view of all key plots from your final analysis.

Reporting: Generate and download comprehensive analysis logs, plot data, and final summary reports.

Prerequisites & Installation
This application is built in R using the Shiny framework. To run it, you need to have R and RStudio installed on your system.

First, you must install all the required R packages. You can do this by running the following command in your R console:

R

install.packages(c(
  "shiny", "shinythemes", "DT", "readxl", "tidyverse", "corrplot", "pls", 
  "glmnet", "prospectr", "car", "caret", "ggpubr", "randomForest", "e1071", 
  "gbm", "shinyjs", "lmtest", "colourpicker", "xgboost", "plotly", "naniar", 
  "tidymodels", "rmarkdown", "ggcorrplot", "ggthemes", "Ckmeans.1d.dp", 
  "rpart", "DALEX", "DALEXtra", "promises", "future", "lightgbm", "catboost", 
  "Cubist", "h2o", "RSNNS", "nortest", "moments"
))
Note: Some packages, like catboost, might have specific system dependencies. Please refer to the official documentation for each package if you encounter installation issues.

How to Run the Application
Save the entire R script provided as a file named app.R.

Open the app.R file in RStudio.

Click the "Run App" button that appears at the top of the script editor.

Alternatively, you can run the app from the R console by navigating to the directory where you saved the file and running:

R

shiny::runApp("app.R")
Workflow Guide
A typical analysis workflow using this tool would be:

Data Explorer: Start by uploading your dataset using the file input or by pasting it into the text area. Select your compound name column if you wish to retain it for labeling.

Data Preprocessing: Navigate to this tab to visualize missing data. Apply imputation or transformation methods as needed. The processed data will be used for all subsequent analyses.

Automated Pipelines: For a fast and comprehensive analysis, go to the Full QSAR Pipeline. Select your response variable, model type, and run the pipeline to get a fully validated model and a suite of diagnostic plots.

Manual Analysis (Optional): For a more customized approach, navigate through the tabs sequentially:

Descriptive Analysis: To understand your data's distribution and correlations.

Dimensionality Reduction: To explore the underlying structure with PCA and clustering.

Modeling: To build and compare various models. The Publication Analysis tab is highly recommended for building a final, robust MLR model.

Advanced Validation: After building a model in the Publication Analysis tab, use the tools in this menu (like Applicability Domain and Y-Scrambling) to rigorously validate its robustness and reliability.

AI Insights: Use the AI-Powered Model Grader and Model Interpretation tabs to gain deeper insights into your model's performance and behavior.

Reporting: Use the Analysis Dashboard for a consolidated view of your main results or the Final Report tab to generate a text-based summary of all analyses performed.

License
This project is open source.
Note:
We invite interested one to contribute in the improvisation towards AI-based QSAR
