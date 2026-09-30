
# ImputationPackage

<!-- badges: start -->
<!-- badges: end -->

The goal of ImputationPackage is to simulate missingness in PatientLevelPrediction data, apply imputation techniques, and evaluate downstream model performance.

## Overview
This package is for simulation studies that want to evaluate how different missing-data structures and imputation methods affect predictive modeling results, using PatientLevelPrediction data objects.

## Features
- Build PLP-ready data objects
- Simulate missingness under different mechanisms and allow different missingness structures in one data object
- Run multiple imputation methods
- Fit prediction models
- Aggregate and summarise evaluation outputs

## Prerequisites
- R >= 4.2.0
- Core dependencies used by the package include:
  - Andromeda
  - CirceR
  - DatabaseConnector
  - DeepPatientLevelPrediction
  - doParallel
  - dplyr
  - foreach
  - magrittr
  - PatientLevelPrediction
  - readr
  - rlang
  - SqlRender
  - tibble


## Installation

From the package root, install the package with:

```r
devtools::install(upgrade = "never", dependencies = TRUE)
```

If you only want to load the package for development work without installing it, use:

```r
devtools::load_all()
```

## Workflow
1. Build a (complete-case) PLP data object
2. Simulate missingness
3. Apply imputation methods
4. Fit prediction models
5. Evaluate and summarise prediction results
6. Visualize the results and run Friedman and Nemenyi tests

## Important Functions
- `buildPopulationPLPData()`: to build a (complete-case) PLP data object.
- `runMissingnessSimulation()`: to run a full simulation study.
- `createPLPImputationMethods()`: to generate the PLP imputation methods from the PatientLevelPrediction package.
- `runPlpImputers()`: to run imputation methods on train- and test data.
- `runPlpModels()`: to run prediction models on train- and test data and store the evaluation metrics.
- `collectEvaluationTables()`: to collect model evaluation outputs into tables for comparison and evaluation summaries.

## Layout
- `R/`: main package functions
- `man/`: generated documentation files
- `extras/`: example workflows, plotting scripts, and significance testing

## Example

Example scripts for running simulations are in `extras/univariate/` and `extras/multivariate/`.
