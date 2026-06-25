
# ImputationPackage

<!-- badges: start -->
<!-- badges: end -->

The goal of ImputationPackage is to simulate missingness in PatientLevelPrediction data, apply imputation techniques and evaluate downstream model performance.

## Overview
This package is designed for studies that want to evaluate how different missing-data structures and imputation methods affect predicitive modeling results.

## Features
- Build PLP-ready data objects
- Simulate missingness under different mechanisms and allow different missingness structures in one data object
- Run multiple imputation methods
- Fit prediction models
- Aggregate and summarise evaluation outputs

## Prerequisites
- R >= 4.2.0
- Required packages:
  - Andromeda
  - CirceR
  - DatabaseConnector
  - dplyr
  - FeatureExtraction
  - margrittr
  - PatientLevelPrediction
  - readr
  - rlang
  - sqlRender
  - tibble
  


## Installation

You can install the development version of ImputationPackage like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Workflow
1. Build a (complete-case) PLP data object
2. Simulate missingness
3. Apply imputation methods
4. Fit prediction models
5. Evaluate and summarise prediction results

## Important Functions
- `buildPopulationPLPData()`: to build a (complete-case) PLP data object.
- `runMissingnessSimulation()`: to run a full simulation study.
- `createPLPImputationMethods()`: to generate the PLP imputation methods from the PatientLevelPrediction package.
- `runPlpImputers()`: to run imputation methods on train- and test data.
- `runPlpModels()`: to run prediction models on train- and test data and store the evaluation metrics.
- `collectEvaluationTables()`: to collect model evaluation outputs into tables for comparison and evaluation summaries.

## Layout
- `R/ `: main package functions
- `extras/`: data exploration scripts and example scripts
- `tests/`: test folder

## Example

An example script to run a simulation is in `extras/run_workflow.R`

