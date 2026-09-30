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

cohortPath <- file.path(projectPath, "cohort_ckd")
resultsFolder <- file.path(projectPath, "results", "Missingness Simulation")

dir.create(resultsFolder, recursive = TRUE, showWarnings = FALSE)

targetJson <- file.path(
  cohortPath,
  "target_ckd.json"
)

outcomeJson <- file.path(
  cohortPath,
  "outcome_ckd.json"
)

targetCohortId <- 3
outcomeCohortId <- 4

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
  40762887,  # Creatinine
  3018251,  # Fasting Glucose
  3025315  # Body Weight
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
  riskWindowEnd = 1825,
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
creatinineConceptId <- 40762887
fastingGlucoseConceptId <- 3018251
bodyWeightConceptId <- 3025315

bpCovariateId <- measurementRef %>%
  dplyr::filter(.data$conceptId == bpConceptId) %>%
  dplyr::pull(.data$covariateId) %>%
  unique()

creatinineCovariateId <- measurementRef %>%
  dplyr::filter(.data$conceptId == creatinineConceptId) %>%
  dplyr::pull(.data$covariateId) %>%
  unique()

fastingGlucoseCovariateId <- measurementRef %>%
  dplyr::filter(.data$conceptId == fastingGlucoseConceptId) %>%
  dplyr::pull(.data$covariateId) %>%
  unique()

bodyWeightCovariateId <- measurementRef %>%
  dplyr::filter(.data$conceptId == bodyWeightConceptId) %>%
  dplyr::pull(.data$covariateId) %>%
  unique()

simulationInputCovariateIds <- c(
  bpCovariateId,
  creatinineCovariateId,
  fastingGlucoseCovariateId,
  bodyWeightCovariateId
)

ageCovariateId <- 1002

################################################################################
# Single run
################################################################################
#
# system.time({
#   singleRun <- runMissingnessSimulation(
#   data = completePopulationPLPData,
#   population = completePopulation,
#   targetCovariateId = creatinineCovariateId,
#   causeCovariateIds = ageCovariateId,
#   completeInputCovariateIds = simulationInputCovariateIds,
#   savePlpModelResults = TRUE,
#   completeCase = TRUE,
#   mechanisms = c("MCAR","MAR"),
#   typePerMech = list(MAR = "LEFT"),
#   missingnessRatios = seq(0,0.8, by = 0.2),
#   #imputationMethods = c("xgb"),
#   imputationMethods = "noImputation",
#   predictionModels = c("lasso","xgboost"),
#   runs = 50,
#   outputFolder = resultsFolder,
#   startSimulation = 36,
#   startScenarioId = 2,
#   seed = 123,
#   #transformerDevice = "cuda:3"
#
# )})

################################################################################
# XGBoost imputation run
################################################################################
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
  singleRun <- runMissingnessSimulation(
    data = completePopulationPLPData,
    population = completePopulation,
    targetCovariateId = creatinineCovariateId,
    causeCovariateIds = ageCovariateId,
    completeInputCovariateIds = simulationInputCovariateIds,
    savePlpModelResults = TRUE,
    completeCase = TRUE,
    mechanisms = c("MCAR","MAR"),
    typePerMech = list(MAR = "LEFT"),
    missingnessRatios = seq(0, 0.8, by = 0.2),
    imputationOverview = imputationOverview,
    imputationMethods = c("xgboost_noIndicator","xgboost_withIndicator",
                              "simpleMean_noIndicator",
                              "simpleMean_withIndicator",
                              "simpleMedian_noIndicator",
                              "simpleMedian_withIndicator",
                              "iterativePMM_noIndicator",
                              "iterativePMM_withIndicator",
                          "noImputation"),
    # xgboostTuningGrid = expand.grid(
    #   nrounds = c(100,300),
    #   maxDepth = c(4,6,8),
    #   eta = c(0.05,0.1,0.3),
    #   stringsAsFactors = FALSE
    # ),
    xgboostTuningGrid = NULL,
    xgboostTuningWorkers = 18L,
    xgboostValidationFraction = 0.2,
    xgboostMinHoldoutPerColumn = 1,
    xgboostTuningMetric = "rmse",
    xgboostVerbose = TRUE,
    predictionModels = "transformer",
    transformerDevice = "cuda:3",
    runs = 50,
    outputFolder = resultsFolder,
    # startSimulation = 4,
    # startScenarioId = 10,
    seed = 123
  )})





