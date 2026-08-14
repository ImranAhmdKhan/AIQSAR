---
title: 'AIQSAR: An Interactive R Shiny Application for AI-Powered Quantitative Structure-Activity Relationship Modeling and Statistical Analysis'
tags:
  - R
  - QSAR
  - cheminformatics
  - machine learning
  - Shiny
  - drug discovery
  - computational chemistry
authors:
  - name: Huma Basheer
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Annet
    orcid: 0000-0000-0000-0000
    affiliation: 2
  - name: Shaurya Raghuvanshi
    orcid: 0000-0000-0000-0000
    affiliation: 3
  - name: Imran A. Khan
    orcid: 0000-0000-0000-0000
    corresponding: true
    affiliation: 1
affiliations:
  - name: Department of Biosciences, COMSATS University Islamabad, Islamabad, Pakistan
    index: 1
  - name: Affiliation to be provided
    index: 2
  - name: Department of Computational Biology and Bioinformatics, Affiliation to be provided
    index: 3
date: 14 August 2026
bibliography: paper.bib
---

# Summary

`AIQSAR` is an open-source, browser-based R Shiny application that provides a
comprehensive and interactive platform for Quantitative Structure-Activity
Relationship (QSAR) modeling and statistical analysis. Designed for researchers,
students, and data scientists in computational chemistry, cheminformatics, and
drug discovery, `AIQSAR` consolidates the entire QSAR workflow—from raw data
ingestion and preprocessing through model training, rigorous validation, and
publication-quality reporting—into a single accessible interface that requires no
programming expertise. The application integrates classical statistical methods
(Multiple Linear Regression, Principal Component Analysis) with modern ensemble
and deep learning algorithms (Random Forest, XGBoost, LightGBM, CatBoost,
H2O AutoML, and neural networks), and couples them with explainability tools
powered by DALEX [@biecek2021dalex], enabling transparent, interpretable, and
reproducible computational chemistry research.

# Statement of Need

QSAR modeling is a cornerstone of rational drug design and toxicity prediction,
translating molecular descriptors into quantitative predictions of biological
activity [@chung2020machine; @tropsha2010best]. Despite its importance, the
practical application of QSAR methods remains a significant challenge for many
researchers due to fragmented toolchains, steep programming requirements, and the
lack of integrated validation frameworks. Existing solutions are either
commercially licensed (e.g., Schrödinger, MOE) or require advanced scripting
skills (e.g., RDKit, DeepChem). There is a clear and unmet need for a
fully open-source, GUI-driven platform that adheres to best-practice guidelines
for QSAR model development and validation [@gramatica2013principles], while
incorporating state-of-the-art machine learning algorithms.

`AIQSAR` fills this gap by providing:

1. A zero-code graphical interface accessible to domain experts without
   computational backgrounds.
2. A fully integrated validation suite covering internal cross-validation,
   external test-set evaluation, Y-scrambling, and applicability domain
   (Williams Plot) analysis [@becker2004applicability].
3. Native support for modern gradient boosting libraries
   [@chen2016xgboost; @ke2017lightgbm; @prokhorenkova2018catboost] and AutoML
   [@ledell2020h2o] alongside classical QSAR methods.
4. AI-powered model grading and DALEX-based model interpretation
   [@biecek2021dalex] to guide researchers toward the most informative and
   reliable models.

# Functionality

Figure 1 provides an overview of the `AIQSAR` graphical interface, illustrating
the tab-based navigation that guides users through the full QSAR workflow.

![Figure 1: Overview of the AIQSAR Shiny application interface showing the tab-based navigation panel (left) and the Full QSAR Pipeline results panel (right), including model comparison metrics and diagnostic plots.](screenshots/03_qsar_pipeline.png)

## Data Input and Preprocessing

`AIQSAR` accepts datasets in CSV or Microsoft Excel format, or allows users to
paste data directly into the application. A dedicated preprocessing module
provides multiple missing-value imputation strategies (mean, median, k-nearest
neighbours) and data transformation options (logarithmic, square root, Box-Cox
[@box1964analysis]). An interactive correlation heatmap and descriptive
statistics panel assist users in understanding data quality prior to modeling.

## Dimensionality Reduction and Clustering

The application integrates Principal Component Analysis (PCA) for dimensionality
reduction together with K-Means and Hierarchical Clustering to uncover latent
structure in chemical descriptor spaces, supporting informed feature selection and
dataset partitioning.

## Automated QSAR Pipelines

A one-click **Full QSAR Pipeline** automates the end-to-end workflow:

- **Dataset Splitting**: Kennard-Stone algorithm [@kennard1969computer] for
  representative training/test set selection.
- **Feature Selection**: Automated multi-step variable selection protocols
  tailored for Multiple Linear Regression (MLR).
- **Model Training**: Simultaneous training and comparison of MLR, Random Forest
  [@liaw2002classification], and XGBoost [@chen2016xgboost].
- **Validation Metrics**: R², Q²LOO, RMSE, MAE, and external validation
  statistics calculated automatically.

A dedicated **Automated MLR** module provides a rigorous stepwise feature
selection pipeline optimized for interpretable QSAR models conforming to OECD
QSAR validation principles [@gramatica2013principles].

## Advanced Machine Learning

Beyond the automated pipelines, `AIQSAR` provides manual access to an extended
model library:

- Gradient Boosting: LightGBM [@ke2017lightgbm], CatBoost
  [@prokhorenkova2018catboost], XGBoost [@chen2016xgboost]
- Neural Networks: MLP and DNN architectures via RSNNS
- AutoML: H2O's automated model search [@ledell2020h2o]
- Regularized Regression: Ridge, Lasso, and Elastic Net via glmnet
- Rule-Based Regression: Cubist

All models are trained and compared within the `caret` framework [@kuhn2008caret]
using unified cross-validation protocols.

## Model Validation

`AIQSAR` implements a comprehensive validation suite aligned with accepted QSAR
best practices [@tropsha2010best; @gramatica2013principles]:

- **Applicability Domain (AD)**: Williams Plot analysis based on leverage
  statistics and standardized residuals [@williams1986regression;
  @becker2004applicability], identifying compounds outside the model's reliable
  prediction space.
- **Y-Scrambling**: Permutation testing to confirm that model performance is not
  attributable to chance correlations between descriptor and response variables.
- **Validation Loop**: Repeated retraining of the final model formula on multiple
  random data splits to assess predictive stability.

## Explainability and AI Insights

A dedicated **AI Insights** module integrates DALEX [@biecek2021dalex] to
provide:

- **Break-Down Plots**: Instance-level feature attribution explaining individual
  predictions.
- **Partial Dependence Plots (PDP)**: Global feature-effect profiles across the
  descriptor space.
- **Model Grader**: An AI-driven scoring system that evaluates the trade-off
  between predictive performance and model interpretability, providing a final
  letter grade and actionable recommendation.

## Reporting

`AIQSAR` generates downloadable, publication-ready outputs including analysis
logs, processed data tables, plot data exports, and consolidated summary reports.
An **Analysis Dashboard** tab provides a unified view of all key diagnostic
figures from the final analysis session.

# Implementation

`AIQSAR` is implemented entirely in R [@rcoreteam2023r] using the Shiny framework
[@chang2015shiny] with `shinythemes` and `shinyjs` for enhanced user-interface
components, `DT` for interactive data tables, and `plotly` for interactive
visualizations. The application is organized into modular UI panels and reactive
server logic, ensuring responsiveness and scalability with moderate-sized
cheminformatics datasets typical in QSAR research (hundreds to low thousands of
compounds).

# Availability and Installation

`AIQSAR` is freely available under an open-source license at
[https://github.com/ImranAhmdKhan/AIQSAR](https://github.com/ImranAhmdKhan/AIQSAR).
The application can be run locally in RStudio by installing the required
dependencies and executing:

```r
shiny::runApp("AIQSAR_1.R")
```

A full list of required packages and installation instructions is provided in the
repository `README.md`.

# Acknowledgements

The authors acknowledge the developers of all open-source R packages integrated
within `AIQSAR`, without which this work would not have been possible.

# References
