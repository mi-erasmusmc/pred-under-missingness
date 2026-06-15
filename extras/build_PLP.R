################################################################################
# Libraries
################################################################################
library(DatabaseConnector)
library(CirceR)
library(SqlRender)
library(readr)
library(dplyr)
library(PatientLevelPrediction)
library(FeatureExtraction)

################################################################################
# Build PLP data objects
#
# Purpose:
# Create the analysis-ready population, covariates, and outcome objects used
# for the simulation study and for visualizations.
#
# Inputs:
# - Database connection details
# - Cohort / outcome definitions
# - Feature extraction settings
#
# Outputs:
# - fullPopulation
# - fullPopulationPLPData
# - completePopulationAll
# - completePopulationAllPLPData
#
# Notes:
# This script is intended to be sourced by downstream analysis scripts.
################################################################################

# set path
source("/code/connection_details.R")

################################################################################
# Helper functions
################################################################################
buildAndExecuteCohort <- function(connection,
                                  jsonPath,
                                  cohortId,
                                  cdmSchema = cdmDbSchema,
                                  vocabularySchema = vocabularyDbSchema,
                                  cohortSchema = cohortDbSchema,
                                  cohortTable = cohortDbTable) {
  cohortJSON <- readr::read_file(jsonPath)
  cohortExpression <- CirceR::cohortExpressionFromJson(cohortJSON)

  cohortSQL <- CirceR::buildCohortQuery(
    expression = cohortExpression,
    options = CirceR::createGenerateOptions()
  )

  renderedSQL <- SqlRender::render(
    cohortSQL,
    cdm_database_schema = cdmSchema,
    vocabulary_database_schema = vocabularySchema,
    target_database_schema = cohortSchema,
    target_cohort_table = cohortTable,
    target_cohort_id = cohortId
  )

  translatedSQL <- SqlRender::translate(
    renderedSQL,
    targetDialect = "postgresql"
  )

  DatabaseConnector::executeSql(connection, translatedSQL)
}

plpDataHelper <- function(labels,
                          folds = NULL,
                          covariates,
                          covariateRef,
                          analysisRef,
                          templatePLPData = plpData) {
  covariateData <- Andromeda::copyAndromeda(templatePLPData$covariateData)

  covariateData$covariates <- covariates
  covariateData$covariateRef <- covariateRef
  covariateData$analysisRef <- analysisRef

  if (!is.null(templatePLPData$covariateData$timeRef)) {
    covariateData$timeRef <- templatePLPData$covariateData$timeRef
  }

  class(covariateData) <- class(templatePLPData$covariateData)
  attr(covariateData, "metaData") <- attr(templatePLPData$covariateData, "metaData")

  out <- list(labels = as.data.frame(labels))

  if (!is.null(folds)) {
    out$folds <- folds
  }

  out$covariateData <- covariateData

  class(out) <- "plpData"
  attr(out, "metaData") <- attr(templatePLPData, "metaData")
  out
}

populationSubset <- function(population, selectedRowIds) {
  filteredPopulation <- population %>%
    dplyr::filter(.data$rowId %in% selectedRowIds)

  attr(filteredPopulation, "metaData") <- attr(population, "metaData")
  filteredPopulation
}

buildPopulationPLPData <- function(selectedRowIds, includeFolds = TRUE) {
  selectedRowIds <- sort(unique(selectedRowIds))

  finalPopulation <- populationSubset(population, selectedRowIds)

  analysisCovariates <- bind_rows(
    demographicsCovariateValues %>% filter(.data$rowId %in% selectedRowIds),
    measurementCovariateValues %>% filter(.data$rowId %in% selectedRowIds)
  )

  folds <- NULL
  if (includeFolds) {
    folds <- data.frame(
      rowId = selectedRowIds,
      index=1L
    )
  }

  list(
    plpData = plpDataHelper(
      labels = finalPopulation,
      folds = folds,
      covariates = analysisCovariates,
      covariateRef = analysisCovariateRef,
      analysisRef = analysisCovariateAnalysis
    ), population = finalPopulation
  )
}

################################################################################
# Connect to database
################################################################################
connection <- DatabaseConnector::connect(connectionDetails)

################################################################################
# Build cohorts
################################################################################
buildAndExecuteCohort(connection = connection, jsonPath = targetJson, cohortId = targetCohortId)

buildAndExecuteCohort(connection = connection, jsonPath = outcomeJson, cohortId = outcomeCohortId)

################################################################################
# Extract PLP data
################################################################################
databaseDetails <- PatientLevelPrediction::createDatabaseDetails(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDbSchema,
  cdmDatabaseName = "",
  cohortDatabaseSchema = cohortDbSchema,
  cohortTable = cohortDbTable,
  targetId = 1,
  outcomeDatabaseSchema = cohortDbSchema,
  outcomeTable = cohortDbTable,
  outcomeIds = 2,
  cdmVersion = 5
)

restrictPlpDataSettings <- PatientLevelPrediction::createRestrictPlpDataSettings()

plpData <- PatientLevelPrediction::getPlpData(
  databaseDetails = databaseDetails,
  covariateSettings = covariateSettings,
  restrictPlpDataSettings = restrictPlpDataSettings
)

################################################################################
# Create population
################################################################################
population <- PatientLevelPrediction::createStudyPopulation(
  plpData = plpData,
  outcomeId = 2,
  populationSettings = populationSettings
)

################################################################################
# Create filtered datasets
################################################################################
populationRowIds <- population$rowId

covariateAnalysis <- plpData$covariateData$analysisRef %>% collect()
covariateRef <- plpData$covariateData$covariateRef %>% collect()

covariatePopulation <- plpData$covariateData$covariates %>%
  filter(.data$rowId %in% populationRowIds) %>%
  collect()

demographicsRef <- covariateRef %>%
  filter(.data$covariateId %in% demographicCovariateIds)

demographicsAnalysis <- covariateAnalysis %>%
  semi_join(demographicsRef, by = "analysisId")

demographicsCovariateValues <- covariatePopulation %>%
  filter(.data$covariateId %in% demographicCovariateIds) %>%
  select(rowId, covariateId, covariateValue)

measurementRef <- covariateRef %>%
  filter(
    !is.na(.data$conceptId),
    .data$conceptId %in% measurementConceptIds
  ) %>%
  distinct(covariateId, .keep_all = TRUE)

measurementAnalysis <- covariateAnalysis %>%
  semi_join(measurementRef, by = "analysisId")

measurementCovariateValues <- covariatePopulation %>%
  filter(.data$covariateId %in% measurementRef$covariateId) %>%
  select(rowId, covariateId, covariateValue)

analysisCovariateRef <- bind_rows(demographicsRef, measurementRef)
analysisCovariateAnalysis <- bind_rows(demographicsAnalysis, measurementAnalysis)


################################################################################
# Build full and complete-case datasets
################################################################################
buildPopulation <- buildPopulationPLPData(populationRowIds)
fullPopulationPLPData <- buildPopulation$plpData
fullPopulation <- buildPopulation$population

# Complete on all 5 measurements
allCompleteCaseRowIds <- measurementCovariateValues %>%
  filter(.data$covariateId %in% measurementRef$covariateId) %>%
  distinct(rowId, covariateId) %>%
  count(rowId, name = "nrMeasurementsObserved") %>%
  filter(.data$nrMeasurementsObserved == length(unique(measurementRef$covariateId))) %>%
  pull(rowId) %>%
  sort()

buildPopulationCompleteAll <- buildPopulationPLPData(allCompleteCaseRowIds)
completePopulationAllPLPData <- buildPopulationCompleteAll$plpData
completePopulationAll <- buildPopulationCompleteAll$population

populationCovData <- fullPopulationPLPData$covariateData
completePopulationAllCovData <- completePopulationAllPLPData$covariateData














