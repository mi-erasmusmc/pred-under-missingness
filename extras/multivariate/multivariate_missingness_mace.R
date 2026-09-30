library(DatabaseConnector)
library(FeatureExtraction)
library(PatientLevelPrediction)
library(dplyr)
#library(ImputationPackage)

packageRoot <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
while (!file.exists(file.path(packageRoot, "DESCRIPTION")) &&
       dirname(packageRoot) != packageRoot) {
  packageRoot <- dirname(packageRoot)
}
if (!file.exists(file.path(packageRoot, "DESCRIPTION"))) {
  stop("Could not locate the package root from the current working directory.", call. = FALSE)
}
devtools::load_all(path = packageRoot, quiet = TRUE, export_all = FALSE)


################################################################################
# Project-specific settings:
# - Fill in the required paths
# - Create connection details
################################################################################

#######################
# Specify project path
projectPath <- "path/to/project"
#######################

cohortPath <- file.path(projectPath, "cohort")
resultsFolder <- file.path(projectPath, "results", "Missingness Simulation")

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
cohortDbSchema <- "your-cohort-schema"
cohortDbTable <- "cohort"

jdbc_driver <- "path/to/jdbc/driver"
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "your-dbms",
  server = "your-server",
  user = "your-user",
  password = "your-password",
  port = NULL,
  pathToDriver = jdbc_driver
)

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

# simulate missingness
patterns <- list(
  list(targetCovariateIds = integer(0)
  ),
  list(
    targetCovariateIds = c(ldlCovariateId, hdlCovariateId, cholCovariateId, bpCovariateId),
    causeCovariateIds = ageCovariateId,
    mechanism = "MAR",
    type = "LEFT"
  ),
  list(
    targetCovariateIds = c(ldlCovariateId, hdlCovariateId, cholCovariateId),
    causeCovariateIds = ageCovariateId,
    mechanism = "MAR",
    type = "LEFT"
  ),
  list(
    targetCovariateIds = c(bpCovariateId),
    causeCovariateIds = ageCovariateId,
    mechanism = "MAR",
    type = "LEFT"
  )
)

#mild
# run mild again
mild <- c(0.4, 0.35, 0.20, 0.05)
realistic <- c(0.18, 0.61, 0.17, 0.04)
extreme <- c(0.05, 0.80, 0.10, 0.05)
veryMild <- c(0.5, 0.25, 0.2, 0.05)
moderate <- c(0.30, 0.45, 0.2, 0.05)
minimal <- c(0.6, 0.15, 0.2, 0.05)
nearComplete <- c(0.8, 0.05, 0.1, 0.05)

imputationOverview <- createPLPImputationMethods(includeMissingIndicator = c(FALSE, TRUE))
imputationOverview$xgboost_noIndicator <- createXgboostImputer(
  missingThreshold = 0.95,
  nrounds = 20,
  maxDepth = 3L,
  eta = 0.3,
  subsample = 0.7,
  colsampleBytree = 1,
  minChildWeight = 1,
  maxiter = 5,
  numThreads = 1,
  treeMethod = "hist",
  maxBin = 256,
  addMissingIndicator = FALSE,
  verbose = TRUE
)
imputationOverview$xgboost_withIndicator <- createXgboostImputer(
  missingThreshold = 0.95,
  nrounds = 20,
  maxDepth = 3L,
  eta = 0.3,
  subsample = 0.7,
  colsampleBytree = 1,
  minChildWeight = 1,
  maxiter = 5,
  numThreads = 1,
  treeMethod = "hist",
  maxBin = 256,
  addMissingIndicator = TRUE,
  verbose = TRUE
)
system.time({
  testDynamicRun <- runDynamicMissingnessSimulation(
    data = completePopulationPLPData,
    population = completePopulation,
    patterns = patterns,
    freq = nearComplete,
    completeInputCovariateIds = c(ldlCovariateId, hdlCovariateId, cholCovariateId, bpCovariateId),
    completeCase = FALSE,
    savePlpModelResults = TRUE,
    imputationOverview = imputationOverview,
    imputationMethods = c("xgboost_noIndicator", "xgboost_withIndicator"
    ),
    predictionModels = c("lasso","xgboost","transformer"),
    runs = 50,
    outputFolder = resultsFolder,
    #startSimulation = 36,
    transformerDevice = "cuda:3",
    seed = 123,
  )
})

