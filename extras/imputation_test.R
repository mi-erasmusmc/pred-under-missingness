# Libraries
library(DatabaseConnector)
library(CirceR)
library(readr)
library(SqlRender)
library(PatientLevelPrediction)
library(FeatureExtraction)
library(dplyr)
library(tidyr)

################################################################################
# Connect to database
################################################################################
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "duckdb",
  server = "C:\\Users\\maud5\\Documents\\database-1M_filtered.duckdb"
)

connection <- DatabaseConnector::connect(connectionDetails)

################################################################################
# Generate target and outcome cohorts
################################################################################
targetJSON <- readr::read_file(
  "C:\\Users\\maud5\\OneDrive\\MSc 2025-2026\\Thesis applications\\Code\\target_cohort_mace_age40_79_strict.json"
)

cohortExpressionTarget <- CirceR::cohortExpressionFromJson(targetJSON)

sqlTarget <- CirceR::buildCohortQuery(
  expression = cohortExpressionTarget,
  options = CirceR::createGenerateOptions()
)

sqlRenderTarget <- SqlRender::render(
  sqlTarget,
  cdm_database_schema = "main",
  vocabulary_database_schema = "main",
  target_database_schema = "main",
  target_cohort_table = "cohort",
  target_cohort_id = 1
)

sqlDuckDBTarget <- SqlRender::translate(
  sqlRenderTarget,
  targetDialect = "duckdb"
)

DatabaseConnector::executeSql(connection, sqlDuckDBTarget)

outcomeJSON <- readr::read_file(
  "C:\\Users\\maud5\\OneDrive\\MSc 2025-2026\\Thesis applications\\Code\\outcome_cohort_mace 1.json"
)

cohortExpressionOutcome <- CirceR::cohortExpressionFromJson(outcomeJSON)

sqlOutcome <- CirceR::buildCohortQuery(
  expression = cohortExpressionOutcome,
  options = CirceR::createGenerateOptions()
)

sqlRenderOutcome <- SqlRender::render(
  sqlOutcome,
  cdm_database_schema = "main",
  vocabulary_database_schema = "main",
  target_database_schema = "main",
  target_cohort_table = "cohort",
  target_cohort_id = 2
)

sqlDuckDBOutcome <- SqlRender::translate(
  sqlRenderOutcome,
  targetDialect = "duckdb"
)

DatabaseConnector::executeSql(connection, sqlDuckDBOutcome)

################################################################################
# Extract PLP data
################################################################################
cdmDbSchema <- "main"
cohortsDbSchema <- "main"
cohortsDbTable <- "cohort"

measurementConceptIds <- c(
  3004249, # Systolic Blood Pressure
  3038553, # BMI
  3027114, # Total cholesterol
  3007070, # HDL cholesterol
  3009966  # LDL cholesterol
)

baseCovariateSettings <- FeatureExtraction::createCovariateSettings(
  useDemographicsGender = TRUE,
  useDemographicsAge = TRUE,
  useDemographicsAgeGroup = TRUE,
  useDemographicsEthnicity = TRUE,
  useDemographicsRace = TRUE,
  useConditionOccurrenceAnyTimePrior = TRUE,
  useDrugExposureAnyTimePrior = TRUE,
  useDeviceExposureAnyTimePrior = TRUE,
  useProcedureOccurrenceAnyTimePrior = TRUE,
  useMeasurementAnyTimePrior = FALSE,
  useMeasurementValueAnyTimePrior = FALSE,
  useObservationAnyTimePrior = TRUE
)

measurementValueSettings <- FeatureExtraction::createCovariateSettings(
  useMeasurementValueAnyTimePrior = TRUE,
  includedCovariateConceptIds = measurementConceptIds,
  addDescendantsToInclude = FALSE
)


covariateSettings <- list(
  baseCovariateSettings,
  measurementValueSettings
)

databaseDetails <- PatientLevelPrediction::createDatabaseDetails(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDbSchema,
  cdmDatabaseName = "",
  cohortDatabaseSchema = cohortsDbSchema,
  cohortTable = cohortsDbTable,
  targetId = 1,
  outcomeDatabaseSchema = cohortsDbSchema,
  outcomeTable = cohortsDbTable,
  outcomeIds = 2,
  cdmVersion = 5
)

# Can use this to restrict the sample size
restrictPlpDataSettings <- PatientLevelPrediction::createRestrictPlpDataSettings()

plpData <- PatientLevelPrediction::getPlpData(
  databaseDetails = databaseDetails,
  covariateSettings = covariateSettings,
  restrictPlpDataSettings = restrictPlpDataSettings
)

################################################################################
# Create population
################################################################################
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

population <- PatientLevelPrediction::createStudyPopulation(
  plpData = plpData,
  outcomeId = 2,
  populationSettings = populationSettings
)

populationRowIds <- population$rowId
covariateAnalysis <- plpData$covariateData$analysisRef %>% collect()
covariateRef <- plpData$covariateData$covariateRef %>% collect()

covariatePopulation <- plpData$covariateData$covariates %>%
  filter(.data$rowId %in% !!populationRowIds) %>%
  collect()

# PLP demographic covariates
demographicCovariateIds <- c(
  8003, 9003, 10003, 11003, 12003, 13003, 14003, 15003, # Age group
  1002, # Age
  8507001, 8532001, # Gender
  38003564005, 38003563005, # Ethnicity
  8516004, 8515004, 8527004 # Race
)

# Helper functions for PLP long-format data
plpDataHelper <- function(labels,
                                folds,
                                covariates,
                                covariateRef,
                                analysisRef,
                                templatePLPData = plpData) {
  covariateData <- Andromeda::andromeda(
    covariates = covariates,
    covariateRef = covariateRef,
    analysisRef = analysisRef
  )
  
  class(covariateData) <- "CovariateData"
  attr(covariateData, "metaData") <- attr(templatePLPData$covariateData, "metaData")
  
  plpData <- list(
    labels = labels,
    folds = folds,
    covariateData = covariateData
  )
  
  class(plpData) <- "plpData"
  attr(plpData, "metaData") <- attr(templatePLPData, "metaData")
  
  plpData
}

# To get population subset from PLPData 
populationSubset <- function(population, selectedRowIds) {
  filteredPopulation <- population %>%
    filter(.data$rowId %in% selectedRowIds)
  
  attr(filteredPopulation, "metaData") <- attr(population, "metaData")
  filteredPopulation
}

# Extract values of one covariate from population 
extractCovariateValues <- function(plpData,
                                     covariateId,
                                     missingMeansZero = FALSE) {
  values <- plpData$labels %>%
    select("rowId") %>%
    left_join(
      plpData$covariateData$covariates %>%
        filter(.data$covariateId == covariateId) %>% # get selected covariate
        collect() %>%
        select("rowId", "covariateValue"), # keep only the necessary columns
      by = "rowId"
    )
  
  # If missing means zero, replace NA by 0
  if (isTRUE(missingMeansZero)) {
    values <- values %>%
      mutate(covariateValue = tidyr::replace_na(.data$covariateValue, 0))
  }
  
  values
}

################################################################################
# Build clean analysis PLP data
################################################################################
demographicsCovariateValues <- covariatePopulation %>%
  filter(.data$covariateId %in% demographicCovariateIds) %>%
  select("rowId", "covariateId", "covariateValue")

demographicsRef <- covariateRef %>%
  filter(.data$covariateId %in% demographicCovariateIds)

demographicsAnalysis <- covariateAnalysis %>%
  semi_join(demographicsRef, by = "analysisId")

measurementRef <- covariateRef %>%
  filter(
    .data$analysisId == 705,
    .data$conceptId %in% measurementConceptIds
  )

measurementAnalysis <- covariateAnalysis %>%
  semi_join(measurementRef, by = "analysisId")

measurementCovariateValues <- covariatePopulation %>%
  filter(.data$covariateId %in% measurementRef$covariateId) %>%
  select("rowId", "covariateId", "covariateValue")

populationRef <- bind_rows(demographicsRef, measurementRef)
populationAnalysis <- bind_rows(demographicsAnalysis, measurementAnalysis)

buildPopulationPLPData <- function(selectedRowIds) {
  selectedRowIds <- sort(unique(selectedRowIds))
  
  studyLabels <- data.frame(rowId = selectedRowIds)
  studyFolds <- data.frame(rowId = selectedRowIds, index = 1L)
  
  populationCovariates <- bind_rows(
    demographicsCovariateValues %>%
      filter(.data$rowId %in% selectedRowIds),
    measurementCovariateValues %>%
      filter(.data$rowId %in% selectedRowIds)
  )
  
  finalPopulation <- populationSubset(population, selectedRowIds)
  
  list(
    plpData = plpDataHelper(
      labels = studyLabels,
      folds = studyFolds,
      covariates = populationCovariates,
      covariateRef = populationRef,
      analysisRef = populationAnalysis
    ),
    population = finalPopulation
  )
}

# Full analysis data: all subjects in the PLP population
buildPopulation <- buildPopulationPLPData(populationRowIds)
fullPopulationPLPData <- buildPopulation$plpData
fullPopulation <- buildPopulation$population

# Find systolic blood pressure covariate extracted by FeatureExtraction
bpCovariateIds <- measurementRef %>%
  filter(.data$conceptId == 3004249) %>%
  pull(.data$covariateId) %>%
  unique()

bpCovariateId <- bpCovariateIds[[1]]

# Complete-case analysis data: patients with observed systolic blood pressure
bpCompleteCaseRowIds <- measurementCovariateValues %>%
  filter(.data$covariateId == bpCovariateId) %>%
  distinct(.data$rowId) %>%
  pull(.data$rowId) %>%
  sort()

buildPopulationComplete <- buildPopulationPLPData(bpCompleteCaseRowIds)
completePopulationPLPData <- buildPopulationComplete$plpData
completePopulation <- buildPopulationComplete$population

################################################################################
# Missing data mechanism functions
#
# Approach follows from:
# Zhang, X. (2023). How to generate missing data for simulation studies. 
# The Quantitative Methods for Psychology, 19(2), 100-122.
#
# The current approach is similar to a sub case of the Ampute package within the MICE package
#
# For MAR and MNAR: Logit approach.
# For MCAR: use binomial
# 
# TO DO:
# - Research possible values of gamma2 (now it is simply set to 1)
# - Create function that allows missing data simulation of multiple variables at the same time
#   with different missing data mechanisms (use Ampute as inspiration, but apply to long format data)
# 
################################################################################
# Use inverse logit
logitInverse <- function(x) {
  1 / (1 + exp(-x))
}

# Function that calculates gamma1 depending on the value of gamma2
computeGamma1 <- function(z, targetProb, gamma2 = 1) {
  uniroot(
    function(gamma1) mean(logitInverse(gamma1 + gamma2 * z)) - targetProb,
    c(-50, 50),
    extendInt = "yes"
  )$root
}

simMCAR <- function(plpData, targetCovariateId, p) {
  # Get all rows where the target covariate exists
  targetRows <- plpData$covariateData$covariates %>%
    filter(.data$covariateId == targetCovariateId) %>%
    collect() %>%
    pull(.data$rowId)
  
  # MCAR: each observed value has probability p of becoming missing
  missing <- as.logical(rbinom(length(targetRows), size = 1, prob = p))
  dropRowIds <- targetRows[missing]
  
  newCovariates <- plpData$covariateData$covariates %>%
    collect() %>%
    filter(!(
      .data$covariateId == targetCovariateId &
        .data$rowId %in% dropRowIds
    ))
  
  plpDataHelper(
    labels = plpData$labels,
    folds = plpData$folds,
    covariates = newCovariates,
    covariateRef = plpData$covariateData$covariateRef %>% collect(),
    analysisRef = plpData$covariateData$analysisRef %>% collect(),
    templatePLPData = plpData
  )
}

simMAR <- function(plpData,
                             targetCovariateId,
                             causeCovariateId,
                             ratio,
                             gamma2 = 1) {
  targetRows <- plpData$covariateData$covariates %>%
    filter(.data$covariateId == targetCovariateId) %>%
    collect() %>%
    select("rowId")
  
  causeValues <- plpData$covariateData$covariates %>%
    filter(.data$covariateId == causeCovariateId) %>%
    collect() %>%
    inner_join(targetRows, by = "rowId") %>%
    filter(!is.na(.data$covariateValue)) %>%
    group_by(.data$rowId) %>%
    summarise(
      covariateValue = first(.data$covariateValue),
      .groups = "drop"
    ) %>%
    mutate(z = as.numeric(scale(.data$covariateValue)))
  
  gamma1 <- computeGamma1(causeValues$z, targetProb = ratio, gamma2 = gamma2)
  prob <- logitInverse(gamma1 + gamma2 * causeValues$z)
  dropRowIds <- causeValues %>%
    mutate(drop = rbinom(n(), size = 1, prob = prob) == 1) %>%
    filter(drop) %>%
    pull(.data$rowId)
  
  newCovariates <- plpData$covariateData$covariates %>%
    collect() %>%
    filter(!(.data$covariateId == targetCovariateId & .data$rowId %in% dropRowIds))
  
  plpDataHelper(
    labels = plpData$labels,
    folds = plpData$folds,
    covariates = newCovariates,
    covariateRef = plpData$covariateData$covariateRef %>% collect(),
    analysisRef = plpData$covariateData$analysisRef %>% collect(),
    templatePLPData = plpData
  )
}

simMNAR <- function(plpData,
                              targetCovariateId,
                              ratio,
                              gamma2 = 1) {
  targetValues <- plpData$covariateData$covariates %>%
    filter(.data$covariateId == targetCovariateId) %>%
    collect() %>%
    filter(!is.na(.data$covariateValue)) %>%
    group_by(.data$rowId) %>%
    summarise(
      covariateValue = first(.data$covariateValue),
      .groups = "drop"
    ) %>%
    mutate(y = as.numeric(scale(.data$covariateValue)))
  
  gamma1 <- computeGamma1(targetValues$y, targetProb = ratio, gamma2 = gamma2)
  prob <- logitInverse(gamma1 + gamma2 * targetValues$y)
  dropRowIds <- targetValues %>%
    mutate(drop = rbinom(n(), size = 1, prob = prob) == 1) %>%
    filter(drop) %>%
    pull(.data$rowId)
  
  newCovariates <- plpData$covariateData$covariates %>%
    collect() %>%
    filter(!(.data$covariateId == targetCovariateId & .data$rowId %in% dropRowIds))
  
  plpDataHelper(
    labels = plpData$labels,
    folds = plpData$folds,
    covariates = newCovariates,
    covariateRef = plpData$covariateData$covariateRef %>% collect(),
    analysisRef = plpData$covariateData$analysisRef %>% collect(),
    templatePLPData = plpData
  )
}

# Example: simulate one incomplete dataset per mechanism with 40% missingness in blood pressure
set.seed(123)

ageCovariateId <- 1002
missingnessRatio <- 0.4
gamma2 <- 1

bpMCAR40 <- simMCAR(
  plpData = completePopulationPLPData,
  targetCovariateId = bpCovariateId,
  p = missingnessRatio
)

# Let missingness be related to age
bpMAR40 <- simMAR(
  plpData = completePopulationPLPData,
  targetCovariateId = bpCovariateId,
  causeCovariateId = ageCovariateId,
  ratio = missingnessRatio,
  gamma2 = gamma2
)

bpMNAR40 <- simMNAR(
  plpData = completePopulationPLPData,
  targetCovariateId = bpCovariateId,
  ratio = missingnessRatio,
  gamma2 = gamma2
)

# Test whether right proportion of data is missing
# Check how much BP is now observed/missing.
# targetCovariateId is covariate id of measurement with missingness
checkMissingness <- function(originalPLPData, incompletePLPData, targetCovariateId) {
  originalNr <- originalPLPData$covariateData$covariates %>%
    filter(.data$covariateId == targetCovariateId) %>%
    collect() %>%
    nrow()
  
  incompleteNr <- incompletePLPData$covariateData$covariates %>%
    filter(.data$covariateId == targetCovariateId) %>%
    collect() %>%
    nrow()
  
  data.frame(
    originalObserved = originalNr,
    incompleteObserved = incompleteNr,
    removed = originalNr - incompleteNr,
    removedFraction = (originalNr - incompleteNr) / originalNr
  )
}

checkMissingness(completePopulationPLPData, bpMCAR40, bpCovariateId)
checkMissingness(completePopulationPLPData, bpMAR40, bpCovariateId)
checkMissingness(completePopulationPLPData, bpMNAR40, bpCovariateId)

################################################################################
# Test calling PLP imputation methods
# TO DO:
# Zoek uit missing threshold
# ZOek uit hoe je kan aangeven dat hij maar 1 variable impute (blood pressure)
################################################################################
plpImputationMethods <- list(
  simpleMean = PatientLevelPrediction::createSimpleImputer(
    method = "mean",
    missingThreshold = 0.95,
    addMissingIndicator = FALSE
  ),
  simpleMedian = PatientLevelPrediction::createSimpleImputer(
    method = "median",
    missingThreshold = 0.95,
    addMissingIndicator = FALSE
  ),
  iterativePMM = PatientLevelPrediction::createIterativeImputer(
    missingThreshold = 0.95,
    method = "pmm",
    methodSettings = list(
      pmm = list(
        k = 5,
        iterations = 5,
        alpha = 1
      )
    ),
    addMissingIndicator = FALSE
  )
)

if ("createSklearnIterativeImputer" %in% getNamespaceExports("PatientLevelPrediction")) {
  plpImputationMethods$sklearn_iterative <- PatientLevelPrediction::createSklearnIterativeImputer(
    missingThreshold = 0.95,
    methodSettings = list(
      maxIter = 10,
      initialStrategy = "mean",
      randomState = 432
    ),
    addMissingIndicator = FALSE
  )
}

# This is to work around the error of imputation functions using `||`, which doesn't work with dplyr::filter??
getPLPImputer <- function(funName) {
  imputerFun <- getFromNamespace(funName, "PatientLevelPrediction")
  
  if (funName %in% c("simpleImpute", "iterativeImpute")) {
    patchedBody <- gsub(
      "dplyr::filter(is.na(.data$missing) || ",
      "dplyr::filter(is.na(.data$missing) | ",
      paste(deparse(body(imputerFun)), collapse = "\n"),
      fixed = TRUE
    )
    body(imputerFun) <- parse(text = patchedBody)[[1]]
  }
  
  imputerFun
}

runPLPImputer <- function(plpData, imputerSettings) {
  imputerFun <- getPLPImputer(attr(imputerSettings, "fun"))
  
  imputerFun(
    trainData = plpData,
    featureEngineeringSettings = imputerSettings,
    done = FALSE
  )
}

# Example: impute the MCAR dataset using mean imputation.
bpMCARMeanImputer <- runPLPImputer(
  plpData = bpMCAR40,
  imputerSettings = plpImputationMethods$simpleMean
)

# Example: impute the MAR dataset using iterative PMM.
bpMARPMMImputer <- runPLPImputer(
  plpData = bpMAR40,
  imputerSettings = plpImputationMethods$iterativePMM
)
