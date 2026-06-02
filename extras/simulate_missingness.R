################################################################################
# Libraries
################################################################################
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
buildAndExecuteCohort <- function(connection,
                                  jsonPath,
                                  cohortId,
                                  cdmSchema = "main",
                                  vocabularySchema = "main",
                                  cohortSchema = "main",
                                  cohortTable = "cohort") {
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
    targetDialect = "duckdb"
  )
  
  DatabaseConnector::executeSql(connection, translatedSQL)
}

buildAndExecuteCohort(
  connection = connection,
  jsonPath = "C:\\Users\\maud5\\OneDrive\\MSc 2025-2026\\Thesis applications\\Code\\target_cohort_mace_age40_79_strict.json",
  cohortId = 1
)

buildAndExecuteCohort(
  connection = connection,
  jsonPath = "C:\\Users\\maud5\\OneDrive\\MSc 2025-2026\\Thesis applications\\Code\\outcome_cohort_mace 1.json",
  cohortId = 2
)

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

bpConceptId <- 3004249
ageCovariateId <- 1002

covariateSettings <- FeatureExtraction::createCovariateSettings(
  useDemographicsGender = TRUE,
  useDemographicsAge = TRUE,
  useDemographicsAgeGroup = TRUE,
  useMeasurementValueLongTerm = TRUE,
  longTermStartDays = -365,
  endDays = 0
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
  filter(.data$rowId %in% populationRowIds) %>%
  collect()

################################################################################
# Keep only the covariates needed for simulation/imputation
################################################################################
demographicCovariateIds <- c(
  8003, 9003, 10003, 11003, 12003, 13003, 14003, 15003, # Age group
  1002, # Age
  8507001, 8532001 # Gender
)

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

missingMeasurementConceptIds <- setdiff(
  measurementConceptIds,
  unique(measurementRef$conceptId)
)

if (length(missingMeasurementConceptIds) > 0) {
  warning(
    paste(
      "No PLP measurement covariate found for conceptId(s):",
      paste(missingMeasurementConceptIds, collapse = ", ")
    )
  )
}

bpCovariateIds <- measurementRef %>%
  filter(.data$conceptId == bpConceptId) %>%
  pull(.data$covariateId) %>%
  unique()

if (length(bpCovariateIds) != 1) {
  stop("Expected exactly one systolic blood pressure covariate.")
}

bpCovariateId <- bpCovariateIds[[1]]

analysisCovariateRef <- bind_rows(demographicsRef, measurementRef)
analysisCovariateAnalysis <- bind_rows(demographicsAnalysis, measurementAnalysis)

################################################################################
# Helper functions to rebuild minimal PLP data
################################################################################
# plpDataHelper <- function(labels,
#                           folds = NULL,
#                           covariates,
#                           covariateRef,
#                           analysisRef,
#                           templatePLPData = plpData) {
#   covariateData <- Andromeda::andromeda(
#     covariates = covariates,
#     covariateRef = covariateRef,
#     analysisRef = analysisRef
#   )
#   
#   class(covariateData) <- "CovariateData"
#   attr(covariateData, "metaData") <- attr(templatePLPData$covariateData, "metaData")
#   
#   plpDataOut <- list(
#     labels = labels,
#     covariateData = covariateData
#   )
#   
#   if (!is.null(folds)) {
#     plpDataOut$folds <- folds
#   }
#   
#   class(plpDataOut) <- "plpData"
#   attr(plpDataOut, "metaData") <- attr(templatePLPData, "metaData")
#   
#   plpDataOut
# }

plpDataHelper <- function(labels,
                          folds = NULL,
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
  attr(class(covariateData), "package") <- "FeatureExtraction"
  attr(covariateData, "metaData") <- attr(templatePLPData$covariateData, "metaData")
  
  plpDataOut <- list(
    labels = labels,
    covariateData = covariateData,
    metaData = templatePLPData$metaData
  )
  
  if (!is.null(folds)) {
    plpDataOut$folds <- folds
  }
  
  class(plpDataOut) <- "plpData"
  plpDataOut
}

populationSubset <- function(population, selectedRowIds) {
  filteredPopulation <- population %>%
    dplyr::filter(.data$rowId %in% selectedRowIds)
  
  attr(filteredPopulation, "metaData") <- attr(population, "metaData")
  filteredPopulation
}

buildPopulationPLPData <- function(selectedRowIds, includeFolds = FALSE) {
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
# Build full and BP-complete analysis datasets
################################################################################
buildPopulation <- buildPopulationPLPData(populationRowIds)
fullPopulationPLPData <- buildPopulation$plpData
fullPopulation <- buildPopulation$population

bpCompleteCaseRowIds <- measurementCovariateValues %>%
  filter(.data$covariateId == bpCovariateId) %>%
  distinct(rowId) %>%
  pull(rowId) %>%
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
# The current approach implements the multivariate amputation structure of the Ampute package within the MICE package 
# + combined with a function to calibrate gamma1
#
# For MAR and MNAR: Logit approach.
# For MCAR: use binomial
# 
# TO DO:
# - Research possible values of gamma2 (now it is simply set to 1)
# - [DONE] Create function that allows missing data simulation of multiple variables at the same time
#   with different missing data mechanisms (use Ampute as inspiration, but apply to long format data)
# 
################################################################################








################################################################################
# Helper functions
################################################################################
# Use inverse logit
logitInverse <- function(x) {
  1 / (1 + exp(-x))
}

# Function that calculates gamma1 depending on the value of gamma2
# targetprob = the target probability of missingness
computeGamma1 <- function(z, targetProb, gamma2 = 1) {
  uniroot(
    function(gamma1) mean(logitInverse(gamma1 + gamma2 * z)) - targetProb,
    c(-50, 50),
    extendInt = "yes"
  )$root
}

# Function that calculates the score of a patient
# causeCovariateIds = covariate ids of the variables that determine missingness
calcScore <- function(covariates, rowIds, causeCovariateIds, weights=NULL, standardize=TRUE) {
  
  if (is.null(weights)) {
    weights <- rep(1, length(causeCovariateIds))
  }
  
  if (length(causeCovariateIds) != length(weights)) {
    stop("Weights and number of variabels causing missingness (causeCovariateIds) must have the same length!")
  }
  
  weightTable <- tibble(covariateId = causeCovariateIds,
                        weight = weights)
  causeData <- covariates %>%
    filter(.data$rowId %in% rowIds,
           .data$covariateId %in% causeCovariateIds) %>%
    collect() %>%
    inner_join(weightTable, by = "covariateId") %>%
    # Multiply value of covariate by its weight
    mutate(weightedVal = .data$covariateValue * .data$weight) %>%
    group_by(.data$rowId) %>% 
    # sum up the weighted values of each covariate for a patient, this is the score of one patient
    summarise(
      score = sum(.data$weightedVal),
      .groups = "drop"
    ) %>% collect()
  
  # include patients that do not have the covariates causing missingness, then set their score to zero
  # Otherwise these patients can never become missing
  causeData <- tibble(rowId = rowIds) %>% left_join(causeData, by="rowId") %>%
    mutate(score = dplyr::coalesce(.data$score, 0))
  
  # Check whether there are unique scores, otherwise, don't standardize
  if (standardize && length(unique(causeData$score)) >1 ) {
    causeData$score <- as.numeric(scale(causeData$score))
  }
  
  causeData
  
}

transformScore <- function(score, type="RIGHT") {
  # Make sure it is upper case
  type <- toupper(type)
  
  meanScore <- mean(score, na.rm=TRUE)
  
  transformed <- switch(
    type,
    "RIGHT" = score - meanScore, # Larger scores (= further away from the mean) more likely to be missing (positive score), below the average gives a negative score --> less likely to be missing
    "LEFT" = meanScore - score, # Smaller scores more likely to be missing
    "MID" = -abs(score - meanScore), # People with average values become most likely to be missing. After logistic transformation, average patients more missing and extreme patients lower probability
    "TAIL"= abs(score - meanScore), # People  far away from the average become most likely to be missing. After logistic transformation, extreme patients more missing and average patients less missing
    stop("Type must be either RIGHT, LEFT, MID or TAIL")
  )
  
  if (length(unique(transformed)) > 1) {
    transformed <- as.numeric(scale(transformed))
  }
  
  transformed
}

# Check if the target covariate is complete before the start of the simulation
observedCovariateCheck <- function(covariatesDf, targetCovariateIds) {
  observedRows <- covariatesDf %>% 
    filter(.data$covariateId %in% targetCovariateIds) %>%
    distinct(.data$rowId, .data$covariateId)
  nrRows <- length(unique(covariatesDf$rowId))
  observedCheck <- observedRows %>% group_by(.data$covariateId) %>%
    summarise(
      nrObserved = n_distinct(.data$rowId),
      .groups = "drop"
    ) %>% mutate(nrMissing = nrRows - nrObserved)
  
  if (any(observedCheck$nrMissing >0)) {
    stop("At least one target covariate contains missing values before simulating missingness.")
  } 
  TRUE
  
}


################################################################################
# Simulate multivariate missingness using different patterns
# This implementation allows variables to follow different missingness mechanisms
# instead of assuming a single mechanism for all variables
################################################################################
simulateMultivariateMissingness <- function(data, 
                                            patterns,
                                             ratio, 
                                            freq = NULL,
                                            seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # If one of the targetCovariateIds accidentally contains missing values; stop procedure
  

  
  # If no patterns are given, stop amputation procedure
  if(!is.list(patterns) || length(patterns) == 0) {
    # Change to running the default pattern function?
    stop("List of patterns cannot be empty.", call.=FALSE)
  }
  
  nrOfPatterns <- length(patterns)
  
  # If no frequency is given, assume all patterns have equal probability
  if (is.null(freq)) {
    freq <- rep(1/nrOfPatterns, nrOfPatterns)
  }
  
  # Frequency length should be equal to the number of patterns (each pattern should be assigned a probability)
  if (length(freq) != nrOfPatterns) {
    stop("Frequency list must have the same length as the number of patterns.", call.=FALSE)
  }
  
  # Frequencies must be nonnegative
  if (any(freq <0)) {
    stop("Frequencies must be nonnegative.", call.=FALSE)
  }
  
  # Ensure frequencies add up to a positive value
  if(sum(freq) == 0) {
    stop("Frequencies cannot sum to 0.", call.=FALSE)
  }
  
  # Ensure frequencies are proportions and recalculate to a proper scale if the frequencies do not add up to 1 or 100
  freq <- freq / sum(freq)
  
  # Extract the different rowIds 
  rowIds <- data$labels %>% dplyr::pull(.data$rowId)
  
  # Select the covariateIds and values of each rowId
  covariatesTable <- data$covariateData$covariates
  covariatesDf <- covariatesTable %>%
    dplyr::collect() %>% 
    dplyr::select("rowId", "covariateId", "covariateValue")
  
  # Assign each row to a missingness pattern
  assignRows <- tibble(
    rowId = rowIds,
    # Randomly assign a pattern to each row
    patternId = sample.int(
      n = nrOfPatterns,
      size=length(rowIds),
      replace=TRUE,
      prob=freq
    )
  )
  
  # This becomes a list of rowIds, covariateIds and patternIds
  droppedPairs <- vector("list", nrOfPatterns)
  # This becomes a list of patterns and assigned rows to the patterns
  patternList <- vector("list", nrOfPatterns)
  
  for (i in seq_along(patterns)) {
    pattern <- patterns[[i]]
    
    if (is.null(pattern$targetCovariateIds)) {
      stop(paste("Pattern", i, "must have targetCovariateIds defined."), call.=FALSE)
      
    }
    
    # Extract the target covariate ids within a pattern & check if they are complete
    targetCovariateIds <- unique(pattern$targetCovariateIds)
    observedCovariateCheck(covariatesDf, targetCovariateIds)

    
    # Extract the covariate ids causing missingness in the target
    causeCovariateIds <- pattern$causeCovariateIds
    # Extract the weights assigned to the variables causing missingness within a pattern
    causeWeights <- pattern$causeWeights
    # Default is "RIGHT"
    type <- if (is.null(pattern$type)) "RIGHT" else toupper(pattern$type)
    gamma2 <- if (is.null(pattern$gamma2)) 1 else pattern$gamma2
    standardizeScore <- if (is.null(pattern$standardizeScore)) TRUE else pattern$standardizeScore
    mechanism <- if (is.null(pattern$mechanism)) {
      if (is.null(causeCovariateIds)) "MCAR" else "CUSTOM"
    } else {
      toupper(pattern$mechanism)
    }
    
    # Filter on the row Ids assigned to the current pattern
    assignedRowIds <- assignRows %>% 
      dplyr::filter(.data$patternId == i) %>% 
      dplyr::pull(.data$rowId)
    
    if (length(assignedRowIds) == 0) {
      droppedPairs[[i]] <- tibble(
        rowId = integer(),
        covariateId = integer(),
        patternId = integer()
      )
      patternList[[i]] <- tibble(
        patternId = i,
        mechanism = mechanism,
        assignedRows = 0L,
        candidateRows = 0L,
        droppedRows = 0L,
        droppedCells = 0L
      )
      next
    }
    
    # rowIds of the rows assigned to a pattern and the target covariateId
    # If assigned rowId does not have the targetCovariateId observed, 
    targetPairs <- covariatesDf %>% 
      filter(.data$rowId %in% assignedRowIds,
             .data$covariateId %in% targetCovariateIds) %>% 
      distinct(.data$rowId, .data$covariateId)
    
    # MCAR if no cause covariates supplied
    if (is.null(causeCovariateIds)) {
      probData <- tibble(
        rowId = assignedRowIds,
        prob = ratio
      )
    # If not MCAR: check if weights are set manually, 
    # else assign equal weights to each covariate causing missingness
    } else {
      if (is.null(causeWeights)) {
        causeWeights <- rep(1, length(causeCovariateIds))
      }
      # Calculate the score in case of MAR and MNAR
      causeData <- calcScore(
        covariates = covariatesTable,
        rowIds = assignedRowIds,
        causeCovariateIds = causeCovariateIds,
        weights = causeWeights,
        standardize = standardizeScore
      ) %>% 
        dplyr::mutate(z=transformScore(.data$score, type))
      
      gamma1 <- computeGamma1(
        z = causeData$z,
        targetProb = ratio,
        gamma2 = gamma2
      )
      
      # Calculate the probabilities based on the scores
      probData <- causeData %>% 
        dplyr::mutate(prob = logitInverse(gamma1 + gamma2 * .data$z)) %>%
        dplyr::select("rowId", "prob")
    }
    
    # Randomize which rows actually become missing
    dropRowIds <- probData %>%
      dplyr::mutate(drop = rbinom(dplyr::n(), size = 1, prob = .data$prob) == 1) %>% 
      dplyr::filter(.data$drop) %>% 
      dplyr::pull(.data$rowId)
    
    # Select the rows and the corresponding pattern id
    droppedPairsSelection <- targetPairs %>% 
      dplyr::filter(.data$rowId %in% dropRowIds) %>%
      dplyr::mutate(patternId = i)
    
    droppedPairs[[i]] <- droppedPairsSelection
    
    # Add the number of rows assign
    patternList[[i]] <- tibble(
      patternId = i,
      mechanism = mechanism,
      assignedRows = length(assignedRowIds), # Number of rows assigned to the pattern
      droppedRows = length(unique(dropRowIds)), # Number of rows that become missing
      droppedCells = nrow(droppedPairsSelection) # Number of cells that become missing
    )
  } 
    
    droppedPairsSelection <- bind_rows(droppedPairs)
    patternSummary <- bind_rows(patternList)
    
    if (nrow(droppedPairsSelection) == 0) {
      newCovariates <- covariatesDf
    } else {
      newCovariates <- covariatesDf %>%
        dplyr::anti_join(
          droppedPairsSelection %>% dplyr::select("rowId", "covariateId"),
          by = c("rowId", "covariateId")
        )
    }
    
    newPlpData <- plpDataHelper(
      labels = data$labels, 
      folds = data$folds, 
      covariates = newCovariates,
      covariateRef = data$covariateData$covariateRef %>% dplyr::collect(), 
      analysisRef = data$covariateData$analysisRef %>% dplyr::collect(),
      templatePLPData = data
    )
    
    list(
      plpData = newPlpData,
      droppedPairs = droppedPairsSelection, 
      assignments = assignRows,
      patternSummary = patternSummary,
      patterns = patterns,
      freq = freq,
      ratio=ratio
    )
}

# To check whether it simulated the correct ratio of missingness
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

################################################################################
# Example: simulate one incomplete dataset per mechanism
# with 40% missingness in blood pressure
################################################################################

set.seed(123)

ageCovariateId <- 1002 # For MAR
missingnessRatio <- 0.4
# Gamma2 determines how strongly the score affects the missingness probability
# Gamma2 = 1 means moderate dependence --> high scores increase missingness, but not deterministically
# If gamma2 is larger, patients with a high score become more likely to be missing
gamma2 <- 1

patterns <- list(
  list(
    mechanism = "MCAR",
    targetCovariateIds = bpCovariateId
  ),
  list(
    mechanism = "MAR",
    targetCovariateIds = bpCovariateId,
    causeCovariateIds = ageCovariateId, # If multiple variables: add causeWeights below
    type = "RIGHT", # people with higher age are more likely to have bp missing
    gamma2 = 1
  ),
  list(
    mechanism = "MNAR",
    targetCovariateIds = c(bpCovariateId),
    causeCovariateIds = c(bpCovariateId), # larger values of bp are more likely to become missing
    type = "RIGHT",
    gamma2 = 1
  )
)

# TO DO:
# Write a default pattern function: one variable missing at a time

multiMissing <- simulateMultivariateMissingness(
  data = completePopulationPLPData,
  patterns = patterns,
  ratio = 0.4,
  freq = c(0.2, 0.4, 0.4),
  seed = 123
)

multiMissing$patternSummary

checkMissingness(
  originalPLPData = completePopulationPLPData,
  incompletePLPData = multiMissing$plpData,
  targetCovariateId = bpCovariateId
)

multiMissing$droppedPairs %>%
  dplyr::count(.data$patternId, .data$covariateId)


## Test whether univariate missingness case works
singlePattern <- list(
  list(mechanism = "MCAR", 
       targetCovariateIds = bpCovariateId)
)

univariateMissingness <- simulateMultivariateMissingness(
  data = completePopulationPLPData,
  patterns = singlePattern,
  ratio = 0.4,
  seed = 123
)
univariateMissingness$patternSummary

checkMissingness(
  originalPLPData = completePopulationPLPData,
  incompletePLPData = univariateMissingness$plpData,
  targetCovariateId = bpCovariateId
)

univariateMissingness$droppedPairs %>%
  dplyr::count(.data$patternId, .data$covariateId)



# Define nr simulations
# In each simulation: runMissingnessSimulation for all 3 mechanisms and all missingness ratios
# This results in multiple datasets per simulation run
# Then take each dataset with simulated missingness: 
# First split into test and training set
# On the training set
# - Complete case analysis
# - mean/median imputation (with and without missingness indicator)
# - ICE + PMM (with and without missingness indicator)
# - PMM (with and without missingness indicator)
# - RF imputation (with and without indicator)








################################################################################
# Test calling PLP imputation methods
# TO DO:
# Investigate missingness threshold (default = 0.95)
# Implement Random Forest imputation technique
# Build loop/function over all possible imputation techniques (including with and without
# missingness indicator)
################################################################################
createPLPImputationMethods <- function(includeMissingIndicator = c(FALSE, TRUE)) {
  methods <- list()
  
  for (addIndicator in includeMissingIndicator) {
    suffix <- if (addIndicator) "_withIndicator" else "_noIndicator"
    
    methods[[paste0("simpleMean", suffix)]] <-
      PatientLevelPrediction::createSimpleImputer(
        method = "mean",
        missingThreshold = 0.95,
        addMissingIndicator = addIndicator
      )
    
    methods[[paste0("simpleMedian", suffix)]] <-
      PatientLevelPrediction::createSimpleImputer(
        method = "median",
        missingThreshold = 0.95,
        addMissingIndicator = addIndicator
      )
    
    methods[[paste0("iterativePMM", suffix)]] <-
      PatientLevelPrediction::createIterativeImputer(
        missingThreshold = 0.95,
        method = "pmm",
        methodSettings = list(
          pmm = list(
            k = 5,
            iterations = 5,
            alpha = 1
          )
        ),
        addMissingIndicator = addIndicator
      )
    
    if ("createSklearnIterativeImputer" %in% getNamespaceExports("PatientLevelPrediction")) {
      methods[[paste0("sklearnIterative", suffix)]] <-
        PatientLevelPrediction::createSklearnIterativeImputer(
          missingThreshold = 0.95,
          methodSettings = list(
            maxIter = 10,
            initialStrategy = "mean",
            randomState = 432
          ),
          addMissingIndicator = addIndicator
        )
    }
  }
  
  methods
}

plpImputationMethods <- createPLPImputationMethods()

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

# runPLPImputer <- function(plpData, imputerSettings) {
#   imputerFun <- getPLPImputer(attr(imputerSettings, "fun"))
#   
#   imputerFun(
#     trainData = plpData,
#     featureEngineeringSettings = imputerSettings,
#     done = FALSE
#   )
# }

runPLPImputers <- function(trainData,
                           testData = NULL,
                           imputationMethods,
                           imputationOverview = plpImputationMethods) {
  if (missing(imputationMethods) || length(imputationMethods) == 0) {
    stop("Provide at least one imputation name.")
  }
  
  unknownMethods <- setdiff(imputationMethods, names(imputationOverview))
  if (length(unknownMethods) > 0) {
    stop(
      "Unknown imputation method(s): ",
      paste(unknownMethods, collapse = ", ")
    )
  }
  
  results <- lapply(imputationMethods, function(methodName) {
    imputerSettings <- imputationOverview[[methodName]]
    funName <- attr(imputerSettings, "fun")
    imputerFun <- getPLPImputer(funName)
    
    # fit on train
    trainImputed <- imputerFun(
      trainData = trainData,
      featureEngineeringSettings = imputerSettings,
      done = FALSE
    )
    
    fittedSettings <- getFittedPLPImputerSettings(
      imputedTrainData = trainImputed,
      imputerSettings = imputerSettings
    )
    
    # apply to test if supplied
    testImputed <- NULL
    if (!is.null(testData)) {
      testImputed <- imputerFun(
        trainData = testData,
        featureEngineeringSettings = fittedSettings,
        done = TRUE
      )
    }
    
    list(
      methodName = methodName,
      functionName = funName,
      originalSettings = imputerSettings,
      fittedSettings = fittedSettings,
      trainImputed = trainImputed,
      testImputed = testImputed
    )
  })
  
  names(results) <- imputationMethods
  results
}

getPLPImputationMetaKey <- function(funName) {
  switch(
    funName,
    simpleImpute = "simpleImputer",
    iterativeImpute = "iterativeImputer",
    sklearnIterativeImpute = "sklearnIterativeImputer",
    stop("Unknown PLP imputation function: ", funName)
  )
}

getFittedPLPImputerSettings <- function(imputedTrainData, imputerSettings) {
  funName <- attr(imputerSettings, "fun")
  metaKey <- getPLPImputationMetaKey(funName)
  
  attr(imputedTrainData$covariateData, "metaData")$
    featureEngineering[[metaKey]]$settings$featureEngineeringSettings
}

################################################################################
# Prediction models


################################################################################
plpPredictionModels <- list(
  lasso = PatientLevelPrediction::setLassoLogisticRegression(seed = 12),
  xgboost = PatientLevelPrediction::setGradientBoostingMachine(seed = 12)
  # ... add DL model
)

runPlpModels <- function(
    imputationResults, 
    modelOverview = plpPredictionModels, 
    preprocessSettings  = PatientLevelPrediction::createPreprocessSettings(
      minFraction = 0.001, 
      normalize = TRUE, 
      removeRedundancy = TRUE
    ), 
    hyperparameterSettings = PatientLevelPrediction::createHyperparameterSettings(),
    analysisPrefix = "sim", 
    analysisPath = tempdir()) {
  results <- list()
  index <- 1
  
  for (imputation in names(imputationResults)) {
    imputationResult <- imputationResults[[imputation]]
    
    for (predModel in names(modelOverview)) {
      trainData <- imputationResult$trainImputed
      testData <- imputationResult$testImputed
      
      trainData$covariateData <- PatientLevelPrediction::preprocessData(
        covariateData = trainData$covariateData, 
        preprocessSettings = preprocessSettings
      )
      
      model <- PatientLevelPrediction::fitPlp(
        trainData = trainData,
        modelSettings = modelOverview[[predModel]],
        hyperparameterSettings = hyperparameterSettings, 
        analysisId = paste0(analysisPrefix, "_", imputation, "_", predModel),
        analysisPath = analysisPath
      )
      
      predictionTrain <- model$prediction
      predictionTest <- PatientLevelPrediction::predictPlp(
        plpModel = model, 
        plpData = testData, 
        population = testData$labels
      )
       # Continue here
      
    }
  }
}








################################################################################
# Run the simulation


################################################################################
runMissingnessSimulation <- function(data,
                                     population, 
                                     targetCovariateId, 
                                     causeCovariateIds = NULL, 
                                     mechanisms = c("MCAR", "MAR", "MNAR"), 
                                     missingnessRatios = seq(0, 0.9, by = 0.1), 
                                     typePerMech = list(MAR = "RIGHT", MNAR = "RIGHT"), 
                                     imputationOverview = plpImputationMethods,
                                     runs = 10,
                                     
                                     seed = 123) {
  results <- list()
  count <- 1
  
  for (simulation in seq_len(runs)) {
    cat("Simulation: ", simulation, "\n")
    
    splitSettings <- createDefaultSplitSetting(
      testFraction = 0.25, 
      trainFraction = 0.75, 
      nfold = 3,
      splitSeed = seed + count,
      type = "stratified"
    )
    
    splitPlpData <- splitData(
      plpData = data, 
      population = population, 
      splitSettings = splitSettings
    )
    
    trainData <- splitPlpData$Train
    testData <- splitPlpData$Test
    
    
    for (mech in mechanisms) {
      for (ratio in missingnessRatios) {
        type <- typePerMech[[mech]]
        if (is.null(type)) {
          type <- "RIGHT"
        }
        
        
        
        pattern <- switch(mech, 
                          "MCAR" = list(targetCovariateIds = targetCovariateId, mechanism = "MCAR"), 
                          "MAR" = list(targetCovariateIds = targetCovariateId, causeCovariateIds = causeCovariateIds, 
                                       mechanism = "MAR", type = type), 
                          "MNAR" = list(targetCovariateIds = targetCovariateId, causeCovariateIds = targetCovariateId, 
                                        mechanism = "MNAR", type = type))
        
        simulateTrain <- simulateMultivariateMissingness(
          data = trainData, 
          patterns = list(pattern), 
          ratio = ratio, 
          seed = seed+count
        )
        
        simulateTest <- simulateMultivariateMissingness(
          data = testData, 
          patterns = list(pattern), 
          ratio = ratio,
          seed = seed + count + 10000
        )
        
        trainMissingData <- simulateTrain$plpData
        testMissingData <- simulateTest$plpData
        
        imputationResults <- runPLPImputers(
          trainData = trainMissingData,
          testData = testMissingData,
          imputationMethods = c("simpleMean_noIndicator","simpleMean_withIndicator", 
                                "simpleMedian_noIndicator", "simpleMedian_withIndicator", 
                                "iterativePMM_noIndicator", "iterativePMM_withIndicator"),
          imputationOverview
        )
        
        
        results[[count]] <- list(
          simulation = simulation,
          mechanism = mech,
          ratio = ratio,
          type = type,
          trainMissingSummary = simulateTrain$patternSummary,
          testMissingSummary = simulateTest$patternSummary,
          imputationResults = imputationResults,
          modelResults = modelResults
        )
        
        count <- count + 1
      }
    }
    results
    
  }
  
  
}

