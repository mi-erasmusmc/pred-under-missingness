#' @title Resolve a PLP imputation function
#' @description
#' This function looks up an imputation function by name and returns a wrapper
#' that copies the PLP data object before execution. It supports both functions provided
#' by the PatientLevelPrediction package and local imputation functions defined in this package
#' @param funName The name of the imputation function to resolve
#' @return A function with arguments `trainData`, `featureEngineeringSettings` and `done`
#'
#' @export
getPLPImputer <- function(funName) {
  # Get imputation method from the PLP namespace functions if provided
  if (exists(funName, envir = asNamespace("PatientLevelPrediction"), inherits = FALSE)) {
    baseFun <- utils::getFromNamespace(funName, "PatientLevelPrediction")

    # Wrap the PLP imputer so it always receives a safe copy of the train input data
    return(function(trainData, featureEngineeringSettings, done = FALSE) {
      baseFun(
        trainData = safeCopyPlpData(trainData),
        featureEngineeringSettings = featureEngineeringSettings,
        done = done
      )
    })
  }

  # If it doesn't exist in PLP, look in the package
  if (!exists(funName, mode = "function", inherits = TRUE)) {
    stop("Unknown PLP imputation function: ", funName)
  }

  localFun <- get(funName, mode = "function", inherits = TRUE)

  # Wrap the custom imputer so it follows the same safe copy interface
  function(trainData, featureEngineeringSettings, done = FALSE) {
    localFun(
      trainData = safeCopyPlpData(trainData),
      featureEngineeringSettings = featureEngineeringSettings,
      done = done
    )
  }
}

#' @title Map an imputation function to its metadata key
#' @description
#' This function returns the feature-engineering metadata key to store
#' fitted settings for a given PLP imputation function
#' @param funName The name of the imputation function
#' @return A string containing the metadata key for the imputer
#'
getPLPImputationMetaKey <- function(funName) {
  # Translate the imputation function name into the corresponding metadata key
  switch(
    funName,
    simpleImpute = "simpleImputer",
    iterativeImpute = "iterativeImputer",
    sklearnIterativeImpute = "sklearnIterativeImputer",
    implementMissForestImputer = "missForestImputer",
    implementXgboostImputer = "xgboostImputer",
    stop("Unknown PLP imputation function: ", funName)
  )
}

#' @title Extract fitted PLP imputer settings
#' @description
#' This function retrieves the fitted feature-engineering settings stored in the
#' imputed training data and returns them in the format needed for test-set imputation.
#' For PMM imputers, the orginal train and test seeds are restored after extraction
#' @param imputedTrainData The imputed training data returned by an imputer
#' @param imputerSettings The original imputer settings used for training-set imputation
#' @return A fitted `featureEngineeringSettings` object that can be reused to impute new data
#' @export
getFittedPLPImputerSettings <- function(imputedTrainData, imputerSettings) {
  # Identify the metadata entry where the fitted imputer stored its settings
  funName <- attr(imputerSettings, "fun")
  metaKey <- getPLPImputationMetaKey(funName)

  # Extract the fitted settings recorded during training imputation
  fittedSettings <- attr(imputedTrainData$covariateData, "metaData")$
    featureEngineering[[metaKey]]$settings$featureEngineeringSettings

  # Restore the original PMM execution seeds so train/test application remains reproducible
  if (isPmmImputerSettings(imputerSettings)) {
    fittedSettings$methodSettings$pmm$trainSeed <- imputerSettings$methodSettings$pmm$trainSeed
    fittedSettings$methodSettings$pmm$testSeed <- imputerSettings$methodSettings$pmm$testSeed
  }

  fittedSettings
}
