# End-to-end workflow example for imputationPackage
#
# This script shows one way to:
# 1. connect to the database
# 2. build target/outcome cohorts
# 3. extract PLP data
# 4. create a complete-case population
# 5. run the missingness simulation study
#
# Update the paths, credentials, schemas, and cohort JSON filenames before use.


library(devtools)
library(DatabaseConnector)
library(FeatureExtraction)
library(PatientLevelPrediction)
library(dplyr)
library(ImputationPackage)


################################################################################
# Project-specific settings:
# - Fill in the required paths
# - Create connection details
################################################################################

#######################
# Specify project path
projectPath <- ""
#######################

cohortPath <- file.path(projectPath, "cohort")
resultsFolder <- file.path(projectPath, "Results", "Missingness Simulation")

dir.create(resultsFolder, recursive = TRUE, showWarnings = FALSE)

targetJson <- file.path(
  cohortPath,
  "target_cohort_mace_age40_79_strict.json"
)

outcomeJson <- file.path(
  cohortPath,
  "outcome_cohort_mace 1.json"
)

targetCohortId <- 1
outcomeCohortId <- 2

cdmDbSchema <- "cdm"
vocabularyDbSchema <- "cdm"
# Specify cohort schema
cohortDbSchema <- " "
cohortDbTable <- "cohort"


connectionDetails <- DatabaseConnector::createConnectionDetails()

measurementConceptIds <- c(
  3004249,  # Systolic Blood Pressure
  3019900,  # Total cholesterol
  3023602,  # HDL cholesterol
  42870529  # LDL cholesterol
)

covariateSettings <- FeatureExtraction::createCovariateSettings(
  useDemographicsGender = TRUE,
  useDemographicsAge = TRUE,
  useDemographicsAgeGroup = FALSE,
  useConditionGroupEraLongTerm = FALSE,
  useConditionOccurrenceLongTerm = TRUE,
  useDrugExposureLongTerm = TRUE,
  useDrugGroupEraLongTerm = FALSE,
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

################################################################################
# Build cohorts and extract PLP data
################################################################################

connection <- DatabaseConnector::connect(connectionDetails)

buildAndExecuteCohort(
  connection = connection,
  jsonPath = targetJson,
  cohortId = targetCohortId,
  cdmSchema = cdmDbSchema,
  vocabularySchema = vocabularyDbSchema,
  cohortSchema = cohortDbSchema,
  cohortTable = cohortDbTable
)

buildAndExecuteCohort(
  connection = connection,
  jsonPath = outcomeJson,
  cohortId = outcomeCohortId,
  cdmSchema = cdmDbSchema,
  vocabularySchema = vocabularyDbSchema,
  cohortSchema = cohortDbSchema,
  cohortTable = cohortDbTable
)

databaseDetails <- PatientLevelPrediction::createDatabaseDetails(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDbSchema,
  cdmDatabaseName = "",
  cohortDatabaseSchema = cohortDbSchema,
  cohortTable = cohortDbTable,
  targetId = targetCohortId,
  outcomeDatabaseSchema = cohortDbSchema,
  outcomeTable = cohortDbTable,
  outcomeIds = outcomeCohortId,
  cdmVersion = 5
)

restrictPlpDataSettings <- PatientLevelPrediction::createRestrictPlpDataSettings()

plpData <- PatientLevelPrediction::getPlpData(
  databaseDetails = databaseDetails,
  covariateSettings = covariateSettings,
  restrictPlpDataSettings = restrictPlpDataSettings
)

population <- PatientLevelPrediction::createStudyPopulation(
  plpData = plpData,
  outcomeId = outcomeCohortId,
  populationSettings = populationSettings
)

################################################################################
# Prepare full and complete-case populations
################################################################################

populationRowIds <- population$rowId

covariateAnalysis <- plpData$covariateData$analysisRef %>% collect()
covariateRef <- plpData$covariateData$covariateRef %>% collect()

covariatePopulation <- plpData$covariateData$covariates %>%
  dplyr::filter(.data$rowId %in% populationRowIds) %>%
  collect()

allCovariateAnalysis <- covariateAnalysis
allCovariateRef <- covariateRef
allCovariateValues <- covariatePopulation %>%
  dplyr::select(.data$rowId, .data$covariateId, .data$covariateValue)

fullBuild <- buildPopulationPLPData(
  selectedRowIds = populationRowIds,
  population = population,
  covariateValues = allCovariateValues,
  covariateRef = allCovariateRef,
  analysisRef = allCovariateAnalysis,
  measurementConceptIds = measurementConceptIds,
  templatePLPData = plpData,
  includeFolds = TRUE
)

fullPopulationPLPData <- fullBuild$plpData
fullPopulation <- fullBuild$population

measurementRef <- covariateRef %>%
  dplyr::filter(
    !is.na(.data$conceptId),
    .data$conceptId %in% measurementConceptIds
  ) %>%
  dplyr::distinct(.data$covariateId, .keep_all = TRUE)

measurementCovariateValues <- covariatePopulation %>%
  dplyr::filter(.data$covariateId %in% measurementRef$covariateId) %>%
  dplyr::select(.data$rowId, .data$covariateId, .data$covariateValue)

completeCaseRowIds <- measurementCovariateValues %>%
  dplyr::distinct(.data$rowId, .data$covariateId) %>%
  dplyr::count(.data$rowId, name = "nrMeasurementsObserved") %>%
  dplyr::filter(.data$nrMeasurementsObserved == length(unique(measurementRef$covariateId))) %>%
  dplyr::pull(.data$rowId) %>%
  sort()

completeBuild <- buildPopulationPLPData(
  selectedRowIds = completeCaseRowIds,
  population = population,
  covariateValues = allCovariateValues,
  covariateRef = allCovariateRef,
  analysisRef = allCovariateAnalysis,
  measurementConceptIds = measurementConceptIds,
  templatePLPData = plpData,
  includeFolds = TRUE
)

completePopulationPLPData <- completeBuild$plpData
completePopulation <- completeBuild$population

################################################################################
# Select variables for the simulation study
################################################################################

bpConceptId <- 3004249
cholConceptId <- 3019900
hdlConceptId <- 3023602
ldlConceptId <- 42870529

bpCovariateId <- measurementRef %>%
  dplyr::filter(.data$conceptId == bpConceptId) %>%
  dplyr::pull(.data$covariateId) %>%
  unique()

cholCovariateId <- measurementRef %>%
  dplyr::filter(.data$conceptId == cholConceptId) %>%
  dplyr::pull(.data$covariateId) %>%
  unique()

hdlCovariateId <- measurementRef %>%
  dplyr::filter(.data$conceptId == hdlConceptId) %>%
  dplyr::pull(.data$covariateId) %>%
  unique()

ldlCovariateId <- measurementRef %>%
  dplyr::filter(.data$conceptId == ldlConceptId) %>%
  dplyr::pull(.data$covariateId) %>%
  unique()

simulationInputCovariateIds <- c(
  bpCovariateId,
  cholCovariateId,
  hdlCovariateId,
  ldlCovariateId
)

ageCovariateId <- 1002

################################################################################
# Single run
################################################################################

system.time({
  singleRun <- runMissingnessSimulation(
  data = completePopulationPLPData,
  population = completePopulation,
  targetCovariateId = bpCovariateId,
  causeCovariateIds = ageCovariateId,
  completeInputCovariateIds = simulationInputCovariateIds,
  completeCase = FALSE,
  mechanisms = "MAR",
  missingnessRatios = 0.8,
  imputationMethods =
    "simpleMedian_withIndicator",
  predictionModels = "lasso",
  runs = 1,
  outputFolder = resultsFolder,
  seed = 123
)})



# ################################################################################
# # Run the simulation study
# ################################################################################
#
#
# simulationResults <- runMissingnessSimulation(
#   data = completePopulationPLPData,
#   population = completePopulation,
#   targetCovariateId = bpCovariateId,
#   causeCovariateIds = ageCovariateId,
#   completeInputCovariateIds = simulationInputCovariateIds,
#   completeCase = TRUE,
#   mechanisms = c("MCAR", "MAR", "MNAR"),
#   missingnessRatios = seq(0, 0.8, by = 0.2),
#   imputationMethods = c(
#     "simpleMean_noIndicator",
#     "simpleMean_withIndicator",
#     "simpleMedian_noIndicator",
#     "simpleMedian_withIndicator",
#     "iterativePMM_noIndicator",
#     "iterativePMM_withIndicator",
#     "sklearnIterative_noIndicator",
#     "sklearnIterative_withIndicator"
#   ),
#   predictionModels = c("lasso", "xgboost"),
#   runs = 150,
#   outputFolder = resultsFolder,
#   seed = 123
# )
#

