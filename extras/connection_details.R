################################################################################
# Study connection and configuration
#
# Purpose:
# Centralize project paths, cohort definitions, database connection settings,
# and core PLP configuration used by downstream scripts.
#
# This file should only define configuration and reusable settings.
################################################################################

################################################################################
# Libraries
################################################################################
library(DatabaseConnector)
library(CirceR)
library(readr)
library(SqlRender)
library(PatientLevelPrediction)
library(RSQLite)
library(FeatureExtraction)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(VIM)
library(corrplot)
library(stringr)
library(paletteer)
library(RPostgres)
library(robustbase)


################################################################################
# Specify paths
################################################################################
projectPath <- " "

cohortPath <- file.path(projectPath,"cohort")

figuresPath <- file.path(projectPath,"plots")

resultsFolder <- file.path(projectPath,"Results","Missingness Simulation")

dir.create(
  figuresPath,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  resultsFolder,
  recursive = TRUE,
  showWarnings = FALSE
)

################################################################################
# Cohort definitions
################################################################################

targetJson <- file.path(
  cohortPath,
  "target_cohort_mace_age40_79_strict.json"
)

outcomeJson <- file.path(
  cohortPath,
  "outcome_cohort_mace 1.json"
)

################################################################################
# Connect to database
################################################################################

connectionDetails <- DatabaseConnector::createConnectionDetails()

################################################################################
# Database schemas
################################################################################
cdmDbSchema <- "cdm"

vocabularyDbSchema <- "cdm"

cohortDbSchema <- " "

cohortDbTable <- "cohort"

################################################################################
# Cohort identifiers
################################################################################
targetCohortId <- 1

outcomeCohortId <- 2

################################################################################
# Measurement concepts
################################################################################
measurementConceptIds <- c(
  3004249, # Systolic Blood Pressure
  3038553, # BMI
  3019900, # Total cholesterol, 3027114 3019900
  3023602, # HDL cholesterol, 3023602 3007070
  42870529   # LDL cholesterol,42870529  3009966
)

measurementConceptLabels <- c(
  "3004249" = "Systolic BP",
  "3038553" = "BMI",
  "3019900" = "Total Cholesterol",
  "3023602" = "HDL",
  "42870529" = "LDL"
)

measurementLabels <- c(
  "Systolic BP",
  "LDL",
  "HDL",
  "Total Cholesterol",
  "BMI"
)

################################################################################
# Demographic covariates
################################################################################
demographicCovariateIds <- c(
  1002, # Age
  8507001, 8532001 # Gender
)

################################################################################
# PLP settings
################################################################################
covariateSettings <- FeatureExtraction::createCovariateSettings(
  useDemographicsGender = TRUE,
  useDemographicsAge = TRUE,
  useDemographicsAgeGroup = FALSE, # Change to true for data exploration
  useConditionOccurrenceLongTerm = TRUE,
  useDrugExposureLongTerm = TRUE,
  useMeasurementValueLongTerm = TRUE,
  longTermStartDays = -365,
  endDays = 0
)

populationSettings <- PatientLevelPrediction::createStudyPopulationSettings(
  binary = TRUE,
  includeAllOutcomes = TRUE,
  firstExposureOnly = TRUE,
  washoutPeriod = 365,
  removeSubjectsWithPriorOutcome = TRUE,
  priorOutcomeLookback = 99999,
  requireTimeAtRisk = FALSE,
  minTimeAtRisk = 0,
  riskWindowStart = 1,
  startAnchor = "cohort start",
  riskWindowEnd = 1095,
  endAnchor = "cohort start"
)


