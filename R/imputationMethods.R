#' @title Create default PLP imputation methods
#' @description
#' This function builds a list of default imputation method settings used in the
#' simulation workflow. It creates variants with and without missingness indicators and checks
#' whether the required packages are installed
#' @param includeMissingIndicator Logical vector indicating whether to create imputationmethods
#' with or without missingness  indicators or both
#' @param sklearnRandomState The random state used fot the sklearn iterative imputer when available
#' @param missForestSeed The seed used for missForest
#' @param iterativeXGBoostSeed The seed used for iterative XGBoost imputation
#' @param pmmTrainSeed The seed used when fitting PMM on the train data
#' @param pmmTestSeed The seed used when applying PMM to the test data
#' @return A list of `featureEngineeringSettings` objects
#' @export
createPLPImputationMethods <- function(includeMissingIndicator = c(FALSE, TRUE),
                                       sklearnRandomState = 432,
                                       missForestSeed = 432L,
                                       iterativeXGBoostSeed = 432L,
                                       pmmTrainSeed = 432L,
                                       pmmTestSeed = 432L) {
  # Initialize the list of imputation method settings
  methods <- list()

  # Create one method set for each requested missing indicator configuration
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
            alpha = 1,
            trainSeed = pmmTrainSeed,
            testSeed = pmmTestSeed
          )
        ),
        addMissingIndicator = addIndicator
      )

    # Include the sklearn iterative imputer only when the PLP export is available
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

    # Include missForest imputer if ranger package is installed
    if (requireNamespace("ranger", quietly = TRUE)) {
      methods[[paste0("missForest", suffix)]] <-
        createMissForestImputer(
          missingThreshold = 0.95,
          addMissingIndicator = addIndicator,
          verbose = TRUE,
          seed = missForestSeed
        )
    }

    # Include xgboost imputer if required packages are installed
    if (requireNamespace("xgboost", quietly = TRUE) && requireNamespace("Matrix", quietly = TRUE)) {
      methods[[paste0("xgboost", suffix)]] <-
        createXgboostImputer(
          missingThreshold = 0.95,
          addMissingIndicator = addIndicator,
          seed = iterativeXGBoostSeed
        )
    }
  }

  methods
}

#' @title Check whether settings correspond to PMM imputer
#' @description
#' This function identifies whether an imputer settings object uses PMM
#' @param imputerSettings An imputer settings object
#' @return Logical indicating whether the imputer method is `PMM`
#'
isPmmImputerSettings <- function(imputerSettings) {
  identical(imputerSettings$method, "pmm")
}

#' @title Get PMM execution seed
#' @description
#'  This function returns the PMM seed used for training or test set imputation. Non-PMM imputers return NULL
#'  @param imputerSettings An imputer settings object
#' @param done Logical parameter indicating whether the imputer is being applied to new data using fitted settings
#' @return A seed for PMM execution (or NULL for any other imputation method)
#'
getPmmExecutionSeed <- function(imputerSettings, done = FALSE) {
  if (!isPmmImputerSettings(imputerSettings)) {
    return(NULL)
  }

  pmmSettings <- imputerSettings$methodSettings$pmm
  if (is.null(pmmSettings)) {
    return(NULL)
  }

  if (isTRUE(done)) {
    pmmSettings$testSeed
  } else {
    pmmSettings$trainSeed
  }
}

#'@title Return a no imputation result
#'@description This function creates a results object for workflow where no
#'imputation is applied and the input PLP data are passed through unchanged
#'@param trainData The training PLP data object
#' @param testData An optional test PLP data object
#' @param methodName The method label used in the returned result
#' @return A list matching the standard imputation result structure
noImputationPlp <- function(trainData, testData = NULL, methodName = "noImputation") {
  # Return standard imputed data structure without changing the data
  list(
    methodName = methodName,
    functionName = NA_character_,
    originalSettings = NULL,
    fittedSettings = NULL,
    trainImputed = safeCopyPlpData(trainData),
    testImputed = if (is.null(testData)) NULL else safeCopyPlpData(testData)
  )
}

#' @title Run a complete-case analysis on PLP data
#' @description
#' This function restricts the training and test data to only keeping the rows that
#' have observed values for all requested target covariates and returns the result in standard imputation-result format
#' @param trainData The training PLP data object
#' @param testData The test PLP data object
#' @param targetCovariateIds The target covariate identifiers that must be observed for a row to be retained
#' @return A list containing the complete-case training and test data
#'
#' @export
completeCasePlp <- function(trainData, testData, targetCovariateIds) {
  # Standardize target variables
  targetCovariateIds <- sort(unique(targetCovariateIds))

  # Keep only the rows with observed valeus for all requested target variables
  getObservedRows <- function(plpData, covariateIds) {
    plpData$covariateData$covariates %>%
      dplyr::filter(.data$covariateId %in% covariateIds) %>%
      dplyr::distinct(.data$rowId, .data$covariateId) %>%
      collectIfNeeded() %>%
      dplyr::count(.data$rowId, name = "nObservedTargets") %>%
      dplyr::filter(.data$nObservedTargets == length(covariateIds)) %>%
      dplyr::pull(.data$rowId)
  }

  trainRows <- getObservedRows(trainData, targetCovariateIds)
  testRows <- getObservedRows(testData, targetCovariateIds)

  # Return the filtered datasets in the standard imputed data structure
  list(
    methodName = "completeCase",
    trainImputed = subsetPlpDataRows(trainData, trainRows),
    testImputed = subsetPlpDataRows(testData, testRows)
  )
}
