################################################################################
# Simulation framework
#
# Purpose:
# Simulate missing data under predefined mechanisms (for example MCAR, MAR,
# or MNAR), apply the chosen imputation method(s), predict outcomes using PLP
# and store evaluationn metrics.
#
# Example run is in line 1424
################################################################################


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
library(glmnet)
library(Cyclops)
library(reticulate)


reticulate::py_config()

# Specify path
source("/build_PLP.R")

################################################################################
# Keep only the covariates needed for simulation/imputation
################################################################################
bpConceptId <- 3004249
bmiConceptId <- 3038553
ageCovariateId <- 1002
ldlConceptId <- 42870529
cholConceptId <- 3019900
hdlConceptId <- 3023602


bpCovariateId <- measurementRef %>%
  filter(.data$conceptId == bpConceptId) %>%
  pull(.data$covariateId) %>%
  unique()

bmiCovariateId <- measurementRef %>%
  filter(.data$conceptId == bmiConceptId) %>%
  pull(.data$covariateId) %>%
  unique()


ldlCovariateId <- measurementRef %>%
  filter(.data$conceptId == ldlConceptId) %>%
  pull(.data$covariateId) %>%
  unique()

hdlCovariateId <- measurementRef %>%
  filter(.data$conceptId == hdlConceptId) %>%
  pull(.data$covariateId) %>%
  unique()

cholCovariateId <- measurementRef %>%
  filter(.data$conceptId == cholConceptId) %>%
  pull(.data$covariateId) %>%
  unique()

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

# Derive stable, named seeds from one master seed.
makeSeed <- function(masterSeed, simulationId = 0, scenarioId = 0, seedLabel = "split") {
  seedOffsets <- c(
    split = 1,
    trainMissing = 2,
    testMissing = 3,
    sklearnImputer = 4,
    lasso = 5,
    xgboost = 6
  )

  if (!seedLabel %in% names(seedOffsets)) {
    stop("Unknown seed label: ", seedLabel)
  }

  as.numeric(masterSeed + simulationId * 10000 + scenarioId * 100 + seedOffsets[[seedLabel]])
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

# Check if input data for simulation study is complete
checkPlpDataComplete <- function(plpData, requiredCovariateIds) {
  requiredCovariateIds <- sort(unique(requiredCovariateIds))
    totalRows <- n_distinct(plpData$labels$rowId)

    observed <- plpData$covariateData$covariates %>%
      filter(.data$covariateId %in% requiredCovariateIds) %>%
      distinct(.data$rowId, .data$covariateId) %>%
      collect() %>% count(.data$covariateId, name = "nObserved")

    check <- tibble(covariateId = requiredCovariateIds) %>% left_join(observed, by="covariateId") %>% mutate(nObserved = coalesce(.data$nObserved, 0L),
                                                                                                             nMissing = totalRows - .data$nObserved)
    if(any(check$nMissing >0)) {
      stop("Simulation input is not complete for the required covariates: ",
           paste(check$covariateId[check$nMissing > 0], collapse = ", "))
    }
    invisible(check)
}

reportMissingness <- function(plpData, covariateIds = NULL) {
  covariates <- plpData$covariateData$covariates %>%
    dplyr::collect() %>%
    dplyr::select(.data$rowId, .data$covariateId)

  # If selected covariateIds are passed: use only those covariates
  if (!is.null(covariateIds)) {
    covariates <- covariates %>%
      dplyr::filter(.data$covariateId %in% covariateIds)
  # If no selected covariateIds are passed: use all covariates
  } else {
    covariateIds <- covariates %>%
      dplyr::distinct(.data$covariateId) %>%
      dplyr::pull(.data$covariateId)
  }

  # Calculate total number of patients
  totalRows <- dplyr::n_distinct(plpData$labels$rowId)

  # For each covariate, count how many unique rows have that covariate observed
  tibble(covariateId = sort(unique(covariateIds))) %>%
    dplyr::left_join(
      covariates %>%
        dplyr::distinct(.data$rowId, .data$covariateId) %>%
        dplyr::count(.data$covariateId, name = "nObserved"),
      by = "covariateId"
    ) %>%
    dplyr::mutate(
      nObserved = dplyr::coalesce(.data$nObserved, 0L),
      nMissing = totalRows - .data$nObserved,
      missingFraction = .data$nMissing / totalRows # Actual missingness (should match the target missingness ratio)
    ) %>%
    dplyr::arrange(dplyr::desc(.data$nMissing), .data$covariateId)
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

################################################################################
# Test calling PLP imputation methods
# TO DO:
# Investigate missingness threshold (default = 0.95)
# Implement Random Forest imputation technique
# Build loop/function over all possible imputation techniques (including with and without
# missingness indicator)
################################################################################
createPLPImputationMethods <- function(includeMissingIndicator = c(FALSE, TRUE),
                                       sklearnRandomState = 432) {
  methods <- list()

  # Create suffix for adding a missingness indicator or not
  for (addIndicator in includeMissingIndicator) {
    suffix <- if (addIndicator) "_withIndicator" else "_noIndicator"

    # Label the imputation methods
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
            randomState = sklearnRandomState
          ),
          addMissingIndicator = addIndicator
        )
    }
  }

  methods
}

# Create the list of imputation methods
plpImputationMethods <- createPLPImputationMethods()

# Create a safe copy of covariate data object, such that it is not affected by imputation and can be reused again for another imputation method
safeCopyCovariateData <- function(covariateData) {
  out <- Andromeda::andromeda(
    covariates = covariateData$covariates %>% dplyr::collect(),
    covariateRef = covariateData$covariateRef %>% dplyr::collect(),
    analysisRef = covariateData$analysisRef %>% dplyr::collect()
  )

  if (!is.null(covariateData$timeRef)) {
    out$timeRef <- covariateData$timeRef %>% dplyr::collect()
  }

  class(out) <- class(covariateData)
  attr(out, "metaData") <- attr(covariateData, "metaData")
  out
}

# Create a safe copy of PLP data object, such that it is not affected by imputation and can be reused again for another imputation method
safeCopyPlpData <- function(x) {
  out <- list(
    labels = as.data.frame(x$labels)
  )

  if (!is.null(x$folds)) {
    out$folds <- as.data.frame(x$folds)
  }

  out$covariateData <- safeCopyCovariateData(x$covariateData)

  class(out) <- class(x)
  attr(out, "metaData") <- attr(x, "metaData")
  out
}

# Call the imputation method from the PLP package
getPLPImputer <- function(funName) {
  baseFun <- getFromNamespace(funName, "PatientLevelPrediction")

  function(trainData, featureEngineeringSettings, done = FALSE) {
    baseFun(
      trainData = safeCopyPlpData(trainData),
      featureEngineeringSettings = featureEngineeringSettings,
      done = done
    )
  }
}

# Run the imputation methods
runPlpImputers <- function(trainData,
                           testData = NULL,
                           imputationMethods,
                           imputationOverview = plpImputationMethods) {
  # At least one imputation method should be supplied
  # TO DO: CHANGE THIS S.T. DEFAULT IS COMPLETE CASE ANALYSIS
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
    # Get the imputation settings from the trained method
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

# Get the imputation settings
getFittedPLPImputerSettings <- function(imputedTrainData, imputerSettings) {
  funName <- attr(imputerSettings, "fun")
  metaKey <- getPLPImputationMetaKey(funName)

  attr(imputedTrainData$covariateData, "metaData")$
    featureEngineering[[metaKey]]$settings$featureEngineeringSettings
}

################################################################################
# Prediction models
################################################################################
# Define PLP prediction models
createPLPPredictionModels <- function(lassoSeed = 12L, xgboostSeed = 12L) {
  list(
    lasso = PatientLevelPrediction::setLassoLogisticRegression(seed = lassoSeed),
    # XGBoost was trained on the complete case to get a fixed set of hyperparameters
    xgboost = PatientLevelPrediction::setGradientBoostingMachine(
      ntrees = 300,
      learnRate = 0.01,
      minChildWeight = 1,
      maxDepth = 4,
      seed = xgboostSeed
    )
    # ... add DL model
  )
}

# Create a list of prediction models
plpPredictionModels <- createPLPPredictionModels()

getPlpMetrics <- function(evaluation, imputation, predictionModel) {

  cat("\n--- Evaluation object structure ---\n")
  print(names(evaluation))
  print(str(evaluation, max.level = 2))

  performance <- evaluation$evaluationStatistics

  cat("\n--- Extracted performance object ---\n")
  print(performance)
  print(str(performance))

  performance %>%
    dplyr::mutate(
      imputation = imputation,
      predictionModel = predictionModel
    )
}

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
  metrics <- list()
  index <- 1

  dir.create(analysisPath, recursive = TRUE, showWarnings = FALSE)

  for (imputation in names(imputationResults)) {
    imputationResult <- imputationResults[[imputation]]

    for (predModel in names(modelOverview)) {
      trainData <- imputationResult$trainImputed
      testData <- imputationResult$testImputed

      # Preprocess the data
      trainData$covariateData <- PatientLevelPrediction::preprocessData(
        covariateData = trainData$covariateData,
        preprocessSettings = preprocessSettings
      )

      cat("\nModel fitting check:\n")
      cat("Imputation:", imputation, "\n")
      cat("Prediction model:", predModel, "\n")
      cat("N train:", nrow(trainData$labels), "\n")
      cat("N outcomes:", sum(trainData$labels$outcomeCount >0), "\n")
      cat("Outcome prevalence:", mean(trainData$labels$outcomeCount > 0), "\n")

      fitResult <- tryCatch(
        {
          model <- PatientLevelPrediction::fitPlp(
            trainData = trainData,
            modelSettings = modelOverview[[predModel]],
            hyperparameterSettings = hyperparameterSettings,
            analysisId = paste0(analysisPrefix, "_", imputation, "_", predModel),
            analysisPath = analysisPath
          )

          predictionTest <- PatientLevelPrediction::predictPlp(
            plpModel = model,
            plpData = testData,
            population = testData$labels
          )

          predictionTrain <- as.data.frame(model$prediction)
          predictionTrain$evaluationType <- "Train"

          predictionTestEval <- as.data.frame(predictionTest)
          predictionTestEval$evaluationType <- "Test"

          combinedPrediction <- dplyr::bind_rows(
            predictionTrain,
            predictionTestEval
          )

          evaluation <- PatientLevelPrediction::evaluatePlp(
            prediction = combinedPrediction,
            typeColumn = "evaluationType"
          )

          list(
            success = TRUE,
            model = model,
            predictionTest = predictionTest,
            evaluation = evaluation,
            errorMessage = NA_character_
          )
        },
        error = function(e) {
          message(
            "Model fit failed for ",
            analysisPrefix,
            " | imputation=", imputation,
            " | predictionModel=", predModel,
            " | error=", conditionMessage(e)
          )

          list(
            success = FALSE,
            model = NULL,
            predictionTest = NULL,
            evaluation = NULL,
            errorMessage = conditionMessage(e)
          )
        }
      )

      results[[index]] <- list(
        imputation = imputation,
        predictionModel = predModel,
        plpModel = fitResult$model,
        predictionTest = fitResult$predictionTest,
        evaluation = fitResult$evaluation,
        success = fitResult$success,
        errorMessage = fitResult$errorMessage
      )

      if (isTRUE(fitResult$success)) {
        metrics[[index]] <- getPlpMetrics(
          evaluation = fitResult$evaluation,
          imputation = imputation,
          predictionModel = predModel
        )
      }

      index <- index + 1
    }
  }

  list(modelResults = results, metrics = dplyr::bind_rows(metrics))
}


################################################################################
# Complete Case Helpers
################################################################################
subsetPlpDataRows <- function(plpData, rowIds) {
  rowIds <- sort(unique(rowIds))

  labels <- as.data.frame(plpData$labels) %>%
    dplyr::filter(.data$rowId %in% rowIds)

  covariates <- plpData$covariateData$covariates %>%
    dplyr::filter(.data$rowId %in% rowIds) %>%
    dplyr::collect()

  covariateRef <- plpData$covariateData$covariateRef %>%
    dplyr::collect()

  analysisRef <- plpData$covariateData$analysisRef %>%
    dplyr::collect()

  folds <- NULL
  if (!is.null(plpData$folds)) {
    folds <- as.data.frame(plpData$folds) %>%
      dplyr::filter(.data$rowId %in% rowIds)
  }

  plpDataHelper(
    labels = labels,
    folds = folds,
    covariates = covariates,
    covariateRef = covariateRef,
    analysisRef = analysisRef,
    templatePLPData = plpData
  )
}

completeCasePlp <- function(trainData, testData, targetCovariateIds) {
  targetCovariateIds <- sort(unique(targetCovariateIds))

  # Collect rows that are observed in the target covariate Ids (other covariates should be complete by definition)
  getObservedRows <- function(plpData, covariateIds) {
    plpData$covariateData$covariates %>%
      dplyr::filter(.data$covariateId %in% covariateIds) %>%
      dplyr::distinct(.data$rowId, .data$covariateId) %>%
      dplyr::collect() %>%
      count(.data$rowId, name = "nObservedTargets") %>%
      filter(.data$nObservedTargets == length(covariateIds)) %>%
      pull(.data$rowId)
      # group_by(.data$rowId) %>%
      # summarise(
      #   nObservedTargets = n_distinct(.data$covariateId),
      #   .groups = "drop"
      # ) %>% filter(.data$nObservedTargets == length(unique(targetCovariateIds))) %>%
      # dplyr::pull(.data$rowId)
  }

  trainRows <- getObservedRows(trainData, targetCovariateIds)
  testRows <- getObservedRows(testData, targetCovariateIds)

  list(
    methodName = "completeCase",
    trainImputed = subsetPlpDataRows(trainData, trainRows),
    testImputed = subsetPlpDataRows(testData, testRows)
  )
}


################################################################################
# Evaluation helpers
################################################################################
extractEvaluationObjects <- function(evaluation, simulation, mechanism, ratio, type, targetVariables, causeVariables, imputation, predictionModel) {
  addTable <- function(x) {
    if (is.null(x)) return(NULL)

    as.data.frame(x) %>% mutate(
      simulation = simulation,
      mechanism = mechanism,
      ratio = ratio,
      type = type,
      targetVariables = targetVariables,
      causeVariables = causeVariables,
      imputation = imputation,
      predictionModel = predictionModel
    )
  }

  list(
    evaluationStatistics = addTable(evaluation$evaluationStatistics),
    calibrationSummary = addTable(evaluation$calibrationSummary),
    thresholdSummary = addTable(evaluation$thresholdSummary),
    demographicSummary = addTable(evaluation$demographicSummary),
    predictionDistribution = addTable(evaluation$predictionDistribution)
  )
}

collectEvaluationTables <- function(modelResults, simulation, mechanism, ratio, type, targetVariables, causeVariables) {
  allTables <- list(
    evaluationStatistics = list(),
    calibrationSummary = list(),
    thresholdSummary = list(),
    demographicSummary = list(),
    predictionDistribution = list())

    for (i in seq_along(modelResults$modelResults)) {
      result <- modelResults$modelResults[[i]]
      tables <- extractEvaluationObjects(
        evaluation = result$evaluation,
        simulation = simulation,
        mechanism = mechanism,
        ratio = ratio,
        type = type,
        targetVariables = targetVariables,
        causeVariables = causeVariables,
        imputation = result$imputation,
        predictionModel = result$predictionModel
      )

      for (tableName in names(tables)) {
        allTables[[tableName]][[i]] <- tables[[tableName]]
      }
    }

  lapply(allTables, bind_rows)

}

# Average results over simulations
summariseOverSims <- function(table) {
  if (is.null(table) || nrow(table) == 0) {
    return(table)
  }

  # Set value to numeric type, such that the average can be calculated
  if ("value" %in% names(table)) {
    table <- table %>% dplyr::mutate(value = as.numeric(value))
  }
  # Define the variables that define a scenario and are available across the scenarios
  groupVars <- intersect(c(
    "simulation","mechanism", "ratio","type","imputation", "targetVariables","causeVariables","predictionModel","evaluation","evaluationType","metric","category","threshold"
  ), names(table))

  # Define the numeric variables that should be summarized
  # Remove columns that are identifiers/group labels
  numericV <- table %>% select(where(is.numeric)) %>% names()
  numericV <- setdiff(numericV, c("simulation",groupVars))

  # Do not group over simulation, because this over what should be summarized
  table %>% group_by(dplyr::across(dplyr::all_of(setdiff(groupVars, "simulation")))) %>%
    summarise(across(all_of(numericV),
                     list(
                       mean = ~ mean(.x, na.rm=TRUE),
                       sd = ~ sd(.x, na.rm= TRUE)
                     ), .names = "{.col}_{.fn}"),
              nSimulation = n_distinct(.data$simulation),
              .groups = "drop")
}

# Add the table to the csv if it already exists, otherwise create a new CSV
appendCsvTable <- function(table, filePath) {
  # If table is empty or null; do nothing
  if (is.null(table) || nrow(table) == 0) {
    return(invisible(NULL))
  }

  # If the folder for the CSV does not exist yet, create a new folder
  dir.create(dirname(filePath), recursive = TRUE, showWarnings = FALSE)

  readr::write_csv(
    table,
    file = filePath,
    append = file.exists(filePath), # If CSV exists, append new rows and do not repeat the column names
    col_names = !file.exists(filePath) # If CSV does not exist, create a new CSV with column names
  )

  invisible(filePath)
}

# Save all outputs from one completed scenario and append to csv
appendScenarioResults <- function(evaluationTables, metrics = NULL, progress = NULL, folder) {
  # Define the raw output folder
  rawFolder <- file.path(folder, "raw")

  # Loop over all evaluation tables
  for (tableName in names(evaluationTables)) {
    table <- evaluationTables[[tableName]]

    if (tableName == "evaluationStatistics" && !is.null(table) && nrow(table) > 0) {
      table <- table %>% dplyr::mutate(
        metric = as.character(metric),
        evaluation = as.character(evaluation)
      )
    }
    # Append each table to its matchin CSV in the raw folder
    appendCsvTable(table, file.path(rawFolder, paste0(tableName, ".csv")))
  }

  # Store evaluation metrics in a CSV
  appendCsvTable(metrics, file.path(rawFolder, "metrics.csv"))
  # Store progress in a CSV
  appendCsvTable(progress, file.path(rawFolder, "scenarioProgress.csv"))
}

# Read raw CSVs back into R
loadRawEvaluationTables <- function(folder) {
  rawFolder <- file.path(folder, "raw")
  tableNames <- c(
    "evaluationStatistics",
    "calibrationSummary",
    "thresholdSummary",
    "demographicSummary",
    "predictionDistribution"
  )

  tables <- lapply(tableNames, function(tableName) {
    filePath <- file.path(rawFolder, paste0(tableName, ".csv"))

    if (!file.exists(filePath)) {
      return(NULL)
    }

    readr::read_csv(filePath, show_col_types = FALSE)
  })

  names(tables) <- tableNames
  tables
}

# Write final summarized output files after raw scenario results have been appended to CSVs
# and results across simulations have been summarized with summariseOverSims()
saveSummarisedResults <- function(aggregatedTables, folder) {
  mainFolder <- file.path(folder, "master")
  specificFolder <- file.path(folder, "specifications")

  dir.create(mainFolder, recursive = TRUE, showWarnings = FALSE)
  dir.create(specificFolder, recursive = TRUE, showWarnings = FALSE)

  # Loop over each summarized table and write it to a CSV
  for (tableName in names(aggregatedTables)) {
    table <- aggregatedTables[[tableName]]

    if (is.null(table) || nrow(table) == 0) {
      next
    }

    if (tableName == "evaluationStatistics") {
      table <- table %>% dplyr::mutate(
        metric = as.character(metric),
        evaluation = as.character(evaluation)
      )

      readr::write_csv(table, file.path(mainFolder, "evaluationStatistics.csv"))
    }

    readr::write_csv(table, file.path(mainFolder, paste0(tableName, ".csv")))

    configurations <- table %>% distinct(mechanism, ratio, type, imputation, predictionModel)

    for (i in seq_len(nrow(configurations))) {
      config <- configurations[i, ]

      configTab <- table %>% filter(
        .data$mechanism == config$mechanism,
        .data$ratio == config$ratio,
        .data$type == config$type,
        .data$imputation == config$imputation,
        .data$predictionModel == config$predictionModel
      )

      outputDirectory <- file.path(
        specificFolder, config$mechanism, paste0("ratio_", config$ratio),
        config$imputation,
        config$predictionModel
      )

      dir.create(outputDirectory, recursive = TRUE, showWarnings = FALSE)

      readr::write_csv(configTab, file.path(outputDirectory, paste0(tableName, ".csv")))

    }
  }
}

# This function is useful for the dynamic case
# If there are multiple target or cause variables, it joins them by using an underscore
# This serves as a label for the simulation output
collapseIds <- function(x) {
  if (is.null(x)) return(NA_character_)
  paste(sort(unique(x)), collapse = "_")
}

# Collect all patternIds
# Extract the covariate ids as separate objects, such that they can be passed to collapseIds
extractPatternIds <- function(patterns, field) {
  ids <- unlist(lapply(patterns, function(pattern) pattern[[field]]), use.names = FALSE)

  if (length(ids) == 0) {
    return(NULL)
  }

  ids
}
################################################################################
# Run the simulation
################################################################################
fiveIds <- c(bpCovariateId, bmiCovariateId, ldlCovariateId, hdlCovariateId, cholCovariateId)
runMissingnessSimulation <- function(data,
                                     population,
                                     targetCovariateId,
                                     causeCovariateIds = NULL,
                                     completeInputCovariateIds = measurementRef$covariateId,
                                     completeCase = FALSE,
                                     mechanisms = c("MCAR", "MAR", "MNAR"),
                                     missingnessRatios = seq(0, 0.8, by = 0.2),
                                     patterns = NULL,
                                     freq = NULL,
                                     typePerMech = list(MAR = "RIGHT", MNAR = "RIGHT"),
                                     imputationOverview = plpImputationMethods,
                                     imputationMethods = c("simpleMean_noIndicator","simpleMean_withIndicator",
                                                           "simpleMedian_noIndicator", "simpleMedian_withIndicator",
                                                           "iterativePMM_noIndicator", "iterativePMM_withIndicator",
                                                           "sklearnIterative_noIndicator", "sklearnIterative_withIndicator"),
                                     predictionModels = c("lasso","xgboost"),
                                     runs = 10,
                                     outputFolder = resultsFolder,
                                     appendIntermediateCsv = TRUE,
                                     startSimulation = 1L,
                                     startScenarioId = 1L,
                                     seed = 123L) {
  results <- list()
  count <- 1

  if (startSimulation < 1L) {
    stop("startSimulation must be at least 1.")
  }

  if (startSimulation > runs) {
    stop("startSimulation cannot be larger than runs.")
  }

  if (startScenarioId < 1L) {
    stop("startScenarioId must be at least 1.")
  }

  for (simulation in seq_len(runs)) {
    if (simulation < startSimulation) {
      next
    }

    cat("Simulation: ", simulation, "\n")

    splitSeed <- makeSeed(seed, simulationId = simulation, seedLabel = "split")
    splitSettings <- createDefaultSplitSetting(
      testFraction = 0.25,
      trainFraction = 0.75,
      nfold = 3,
      splitSeed = splitSeed,
      type = "stratified"
    )

    splitPlpData <- splitData(
      plpData = data,
      population = population,
      splitSettings = splitSettings
    )

    trainData <- splitPlpData$Train
    testData <- splitPlpData$Test

    checkPlpDataComplete(trainData, completeInputCovariateIds)
    checkPlpDataComplete(testData, completeInputCovariateIds)

    scenarioId <- 0L
    firstScenarioThisSimulation <- if (simulation == startSimulation) startScenarioId else 1L

    for (mech in mechanisms) {
      for (ratio in missingnessRatios) {
        scenarioId <- scenarioId + 1L

        if (scenarioId < firstScenarioThisSimulation) {
          next
        }

        type <- typePerMech[[mech]]
        if (is.null(type)) {
          type <- "RIGHT"
        }

        if (!is.null(patterns)) {
          patternsSettings <- patterns
        } else {
          patternsSettings <- list(
            switch(mech,
                   "MCAR" = list(targetCovariateIds = targetCovariateId, mechanism = "MCAR"),
                   "MAR" = list(targetCovariateIds = targetCovariateId, causeCovariateIds = causeCovariateIds, mechanism = "MAR", type = type),
                   "MNAR" = list(targetCovariateIds = targetCovariateId, causeCovariateIds = targetCovariateId, mechanism = "MNAR", type = type))
          )
        }

        targetVariables <- collapseIds(
          extractPatternIds(patternsSettings, "targetCovariateIds")
        )
        causeVariables <- collapseIds(
          extractPatternIds(patternsSettings, "causeCovariateIds")
        )

        trainMissingSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "trainMissing")
        testMissingSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "testMissing")
        sklearnSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "sklearnImputer")
        lassoSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "lasso")
        xgboostSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "xgboost")

        modelOverview <- createPLPPredictionModels(
          lassoSeed = lassoSeed,
          xgboostSeed = xgboostSeed
        )[predictionModels]

        if (identical(imputationOverview, plpImputationMethods)) {
          imputationOverviewCurrent <- createPLPImputationMethods(
            includeMissingIndicator = c(FALSE, TRUE),
            sklearnRandomState = sklearnSeed
          )
        } else {
          imputationOverviewCurrent <- imputationOverview
        }

        cat("\nPatterns used:\n")
        print(patternsSettings)
        simulateTrain <- simulateMultivariateMissingness(
          data = trainData,
          patterns = patternsSettings,
          ratio = ratio,
          freq = freq,
          seed = trainMissingSeed
        )

        simulateTest <- simulateMultivariateMissingness(
          data = testData,
          patterns = patternsSettings,
          ratio = ratio,
          freq = freq,
          seed = testMissingSeed
        )

        trainMissingData <- simulateTrain$plpData
        testMissingData <- simulateTest$plpData

       imputationResults <- runPlpImputers(
          trainData = trainMissingData,
          testData = testMissingData,
          imputationMethods = imputationMethods,
          imputationOverview = imputationOverviewCurrent
        )


       if (completeCase) {
         completeCaseResult <- completeCasePlp(
           trainData = trainMissingData,
           testData = testMissingData,
           targetCovariateIds = targetCovariateId)

         imputationResults <- c(list(completeCase = completeCaseResult),
                                imputationResults)
       }

        modelResults <- runPlpModels(
          imputationResults = imputationResults,
          modelOverview = modelOverview,
          analysisPrefix = paste0("sim", simulation, "_", mech, "_ratio", ratio),
          analysisPath = file.path(tempdir(), "missingnessSimulation")
        )

        evaluationTables <- collectEvaluationTables(modelResults = modelResults,
                                                    simulation = simulation,
                                                    mechanism = mech,
                                                    ratio = ratio,
                                                    type = type,
                                                    targetVariables,
                                                    causeVariables)
        # Store simulation progress
        scenarioProgress <- tibble(
          simulation = simulation,
          scenarioId = scenarioId,
          resultIndex = count,
          mechanism = mech,
          ratio = ratio,
          type = type,
          targetVariables = targetVariables,
          causeVariables = causeVariables,
          splitSeed = splitSeed,
          trainMissingSeed = trainMissingSeed,
          testMissingSeed = testMissingSeed,
          sklearnSeed = sklearnSeed,
          lassoSeed = lassoSeed,
          xgboostSeed = xgboostSeed,
          status = "completed"
        )

        # Append intermediate simulation results
        if (appendIntermediateCsv) {
          appendScenarioResults(
            evaluationTables = evaluationTables,
            metrics = modelResults$metrics %>%
              dplyr::mutate(
                simulation = simulation,
                scenarioId = scenarioId,
                mechanism = mech,
                ratio = ratio,
                type = type,
                targetVariables = targetVariables,
                causeVariables = causeVariables
              ),
            progress = scenarioProgress,
            folder = outputFolder
          )
        }

        # Store results
        results[[count]] <- list(
          simulation = simulation,
          mechanism = mech,
          ratio = ratio,
          type = type,
          targetVariables = targetVariables,
          causeVariables = causeVariables,
          seeds = list(
            split = splitSeed,
            trainMissing = trainMissingSeed,
            testMissing = testMissingSeed,
            sklearnImputer = sklearnSeed,
            lasso = lassoSeed,
            xgboost = xgboostSeed
          ),
          trainMissingSummary = simulateTrain$patternSummary,
          testMissingSummary = simulateTest$patternSummary,
          imputationResults = imputationResults,
          modelResults = modelResults,
          evaluationTables = evaluationTables,
          metrics = modelResults$metrics %>%
            dplyr::mutate(
              simulation = simulation,
              scenarioId = scenarioId,
              mechanism = mech,
              ratio = ratio,
              type = type
            )
        )
        # Increase count for the next simulation run
        count <- count + 1
      }
    }

  }

  # Collect the results from all runs and aggregate the results over the simulation runs
  if (appendIntermediateCsv) {
    allEvaluationTables <- loadRawEvaluationTables(outputFolder)
  } else {
    allEvaluationTables <- list(
      evaluationStatistics = dplyr::bind_rows(lapply(results, function(x) x$evaluationTables$evaluationStatistics)),
      calibrationSummary = dplyr::bind_rows(lapply(results, function(x) x$evaluationTables$calibrationSummary)),
      thresholdSummary = dplyr::bind_rows(lapply(results, function(x) x$evaluationTables$thresholdSummary)),
      demographicSummary = dplyr::bind_rows(lapply(results, function(x) x$evaluationTables$demographicSummary)),
      predictionDistribution = dplyr::bind_rows(lapply(results, function(x) x$evaluationTables$predictionDistribution))
    )
  }

  aggregatedTables <- lapply(allEvaluationTables, summariseOverSims)

  saveSummarisedResults(
    aggregatedTables = aggregatedTables,
    folder = outputFolder
  )

  list(
    results = results,
    allEvaluationTables = allEvaluationTables,
    aggregatedTables = aggregatedTables
  )

}

# Univariate missingness simulation
system.time({

  testSingleRun <- runMissingnessSimulation(
    data =completePopulationAllPLPData,
    population = completePopulationAll,
    targetCovariateId = ldlCovariateId,
    causeCovariateIds = cholCovariateId,
    completeCase = TRUE,

    mechanisms = "MAR",
    missingnessRatios = 0.8,

    imputationMethods = "sklearnIterative_withIndicator",
    predictionModels = "lasso",

    runs = 1,
    seed = 123,
  )

})

