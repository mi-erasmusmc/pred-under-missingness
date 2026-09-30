#' @title Run one or more PLP imputation methods
#' @description
#' This function applies the requested imputation methods to the training data and,
#' when test data are supplied, applies the fitted imputers to the test data using
#' the corresponding fitted settings.
#' @param trainData A PLP data object used to fit the imputation methods
#' @param testData An optional PLP data object to impute using the fitted imputation
#' settings
#' @param imputationMethods A vector naming the imputation methods
#' @param imputationOverview An optional named list of imputation settings. If `NULL`,
#' the default settings from `createPLPImputationMethods()` are used.
#' @return A named list containing the imputed training data, optional imputed test data,
#' and the original and fitted settings for each method.
#' @export
runPlpImputers <- function(trainData,
                           testData = NULL,
                           imputationMethods,
                           imputationOverview = NULL) {
  # Validate the requested imputation methods
  if (missing(imputationMethods) || length(imputationMethods) == 0) {
    stop("Provide at least one imputation name.")
  }

  if (is.null(imputationOverview)) {
    imputationOverview <- createPLPImputationMethods()
  }

  passthroughMethods <- c("noImputation", "none")
  unknownMethods <- setdiff(imputationMethods, c(names(imputationOverview), passthroughMethods))
  if (length(unknownMethods) > 0) {
    stop("Unknown imputation method(s): ", paste(unknownMethods, collapse = ", "))
  }

  # Run each requested imputation method in turn
  results <- lapply(imputationMethods, function(methodName) {
    # Skip imputation if no imputation is requested
    if (methodName %in% passthroughMethods) {
      message("Skipping imputation for method '", methodName, "'.")
      return(noImputationPlp(
        trainData = trainData,
        testData = testData,
        methodName = methodName
      ))
    }

    imputerSettings <- imputationOverview[[methodName]]
    funName <- attr(imputerSettings, "fun")
    imputerFun <- getPLPImputer(funName)

    startTime <- Sys.time()

    # Fit the imputer on the train data
    message(
      "Running imputation method '", methodName,
      "'using function '", funName, "' on the training set."
    )
    trainImputed <- runImputerWithOptionalSeed(
      imputerFun = imputerFun,
      trainData = trainData,
      featureEngineeringSettings = imputerSettings,
      done = FALSE
    )

    message(
      "Finished training imputation for '", methodName,
      "' in ", formatImputationElapsed(startTime), "."
    )

    # Extract fitted imputer settings
    fittedSettings <- getFittedPLPImputerSettings(
      imputedTrainData = trainImputed,
      imputerSettings = imputerSettings
    )

    testImputed <- NULL
    if (!is.null(testData)) {
      testStartTime <- Sys.time()
      # Apply the fitted imputer to the test set
      message("Applying imputation method '", methodName, "' to the test set.")
      testImputed <- runImputerWithOptionalSeed(
        imputerFun = imputerFun,
        trainData = testData,
        featureEngineeringSettings = fittedSettings,
        done = TRUE
      )

      message(
        "Finished test-set imputation for '", methodName,
        "' in ", formatImputationElapsed(testStartTime), "."
      )
    }

    # Return imputation output for this workflow
    message(
      "Completed imputation workflow for '", methodName,
      "' in ", formatImputationElapsed(startTime), "."
    )

    list(
      methodName = methodName,
      functionName = funName,
      originalSettings = imputerSettings,
      fittedSettings = fittedSettings,
      trainImputed = trainImputed,
      testImputed = testImputed
    )
  })

  # Name results by imputation method
  names(results) <- imputationMethods
  results
}

#' @title Run an imputer with an optional PMM seed
#' @description
#' This function runs an imputation function and appllies a seed only when the imputer settings correspond
#' to PMM
#' @param imputerFun The imputation function to execute
#' @param trainData The PLP data object to impute
#' @param featureEngineeringSettings The imputer settings passed to the imputation function
#' @param done Logical parameter to indicate whether the imputer is being applied to
#' new data using the fitted settings
#' @return The results returned by imputerFun
#'
runImputerWithOptionalSeed <- function(imputerFun, trainData, featureEngineeringSettings, done = FALSE) {
  # Derive PMM seed
  seed <- getPmmExecutionSeed(featureEngineeringSettings, done = done)

  # Run imputer directly if no specific PMM seed is set
  if (is.null(seed)) {
    return(imputerFun(
      trainData = trainData,
      featureEngineeringSettings = featureEngineeringSettings,
      done = done
    ))
  }

  # Run imputer with requested seed
  withr::with_seed(
    seed = seed,
    code = imputerFun(
      trainData = trainData,
      featureEngineeringSettings = featureEngineeringSettings,
      done = done
    )
  )
}

#' @title Apply scenario-specific seeds to an imputation overview
#' @description This function creates a local copy of an imputation overview and updates the sparse xgboost imputer
#' when that method is present. It is used to allow externally defined imputation methods to use the
#' generated scenario-specific seed in the simulation workflow.
#' @param imputationOverview A named list of imputation settings
#' @param iterativeXgboostSeed An optional scenario-specific seed for the sparse xgboost imputer
#' @return A copy of `imputationOverview` with the sparse xgboost seed updated when applicable
applyScenarioSeedsToImputationOverview <- function(imputationOverview, iterativeXGBoostSeed = NULL) {
  # Work on a local copy of the imputation overview
  updatedOverview <- imputationOverview

  for (methodName in names(updatedOverview)) {
    settings <- updatedOverview[[methodName]]

    # Skip entries that are not feature-engineering settings objects
    if (!inherits(settings, "featureEngineeringSettings")) {
      next
    }

    funName <- attr(settings, "fun")

    # Overwrite the iterative xgboost seed so that externally created imputers
    # can still use the scenario-specific seed in the simulation
    if (identical(funName, "implementXgboostImputer") &&
        !is.null(iterativeXGBoostSeed)) {
      settings$seed <- iterativeXGBoostSeed
    }

    updatedOverview[[methodName]] <- settings
  }

  updatedOverview
}

#' @title Format an imputation runtime
#' @description
#' This function formats the elapsed time between two timestamp as a readable string in seconds
#' @param startTime The start time
#' @param endTime The end time. Defaults to the current time
#' @return A character string describing the elapsed time in seconds
#'
formatImputationElapsed <- function(startTime, endTime = Sys.time()) {
  # Covnert elapsed runtime to a string representing time in seconds
  sprintf("%.2f seconds", as.numeric(difftime(endTime, startTime, units = "secs")))
}
