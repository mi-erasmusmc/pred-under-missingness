#' @title Set the PLP prediction models
#' @description This function creates named PLP prediction models required for the
#' simulation study, using default settings.
#' @param lassoSeed A random seed for the lasso prediction model
#' @param xgboostSeed A random seed for the XGBoost prediction model
#' @param transformerSeed A random seed for the transformer model
#' @param transformerDevice A parameter to specify the device to run the transformer model on
#' @return A named list containing three PLP model specifications
#' @export
createPLPPredictionModels <- function(lassoSeed = 12L,
                                      xgboostSeed = 12L,
                                      transformerSeed = 12L,
                                      transformerDevice = NULL) {

  # If no device is supplied, run on the default device (cpu)
  transformerEstimator <- if (is.null(transformerDevice)) {
    DeepPatientLevelPrediction::setEstimator(seed = transformerSeed)
  } else {
    DeepPatientLevelPrediction::setEstimator(seed = transformerSeed, device = transformerDevice)
  }

  # Return the model specifications as a named list
  list(
    lasso = PatientLevelPrediction::setLassoLogisticRegression(seed = lassoSeed),
    xgboost = PatientLevelPrediction::setGradientBoostingMachine(
      seed = xgboostSeed
    ),
    transformer = DeepPatientLevelPrediction::setDefaultTransformer(
      estimatorSettings = transformerEstimator
    )
  )
}

#' @title Extract PLP evaluation metrics
#' @description
#' Converts the evaluation statistics to a data frame and identifies the corresponding
#' imputation method and prediction model
#' @param evaluation A PLP evaluation object
#' @param imputation The imputation method used in the scenario
#' @param predictionModel The PLP model used in the scenario
#' @return A data frame containing the PLP evaluation statistics together with
#' `imputation` and `predictionModel` columns
getPlpMetrics <- function(evaluation, imputation, predictionModel) {
  # Standardize to correct format
  # Identify the imputation/prediction model combination that produced each metric
  evaluation$evaluationStatistics %>%
    dplyr::mutate(
      evaluation = as.character(.data$evaluation),
      metric = as.character(.data$metric),
      value = as.numeric(.data$value),
      imputation = imputation,
      predictionModel = predictionModel
    )
}

#' @title Fit and evaluate PLP models across imputed data sets
#' @description Fits multiple PLP models for each imputed data set, generates predictions for
#' training and test set, evaluates model performance, and collects both detailed model results and evaluation metrics
#' @param imputationResults Named list of imputation results. Each element contains a `trainImputed` object and
#' optionally a `testImputed` object
#' @param modelOverview A named list of PLP model specifications.
#' If `NULL` models are created with `createPLPPredictionModels()`
#' @param preprocessSettings Preprocessing setings. By default, covariates occurring in less than 0.1 percent of the
#' population are removed, normalization is enabled and redundant covariates are removed
#' @param hyperparameterSettings Hyperparameter search settings
#' @param analysisPrefix Prefix used when constructing `analysisId` supplied to [fitPlp()]
#' @param analysisPath Directory in which PLP analysis object are written
#' @return A list with model results for each imputation/prediction combination, recording the fitted model,
#' test predictions, evaluation object, success status, and any error message
#' @export
runPlpModels <- function(imputationResults,
                         modelOverview = NULL,
                         preprocessSettings = PatientLevelPrediction::createPreprocessSettings(
                           minFraction = 0.001,
                           normalize = TRUE,
                           removeRedundancy = TRUE
                         ),
                         hyperparameterSettings = PatientLevelPrediction::createHyperparameterSettings(),
                         analysisPrefix = "sim",
                         analysisPath = tempdir()) {
  # Construct the defulat set of prediction models
  if (is.null(modelOverview)) {
    modelOverview <- createPLPPredictionModels()
  }

  # Store the detailed results and evaluation metrics separately
  results <- list()
  metrics <- list()

  # Use a single index because each imputation can be evaluated with multiple
  # prediction models
  index <- 1

  dir.create(analysisPath, recursive = TRUE, showWarnings = FALSE)

  # Iterate over each imputed data set
  for (imputation in names(imputationResults)) {
    imputationResult <- imputationResults[[imputation]]

    # Fit each prediction model to the current imputed dataset
    for (predModel in names(modelOverview)) {
      # Catch errors such that it does not terminate the simulation study
      fitResult <- tryCatch(
        {
          # Work with copies of the PLP data, such that the original object is not modified
          trainData <- safeCopyPlpData(imputationResult$trainImputed)
          testData <- if (is.null(imputationResult$testImputed)) {
            NULL
          } else {
            safeCopyPlpData(imputationResult$testImputed)
          }

      # Preprocess the training data before fitting the prediction model
      trainData$covariateData <- PatientLevelPrediction::preprocessData(
        covariateData = trainData$covariateData,
        preprocessSettings = preprocessSettings
      )

      # Fit the PLP model
      # The analysisId contains the imputation and prediction model name
      model <- PatientLevelPrediction::fitPlp(trainData = trainData,
                                              modelSettings = modelOverview[[predModel]],
                                              hyperparameterSettings = hyperparameterSettings,
                                              analysisId = paste0(analysisPrefix, "_", imputation, "_", predModel),
                                              analysisPath = analysisPath
                                              )

      # Generate predictions for the test set
      predictionTest <- PatientLevelPrediction::predictPlp(
        plpModel = model,
        plpData = testData,
        population = testData$labels
        )

      # Store the predictions for the train set within the fitted model object
      # Convert them to a data frame before combining with test predictions
      predictionTrain <- as.data.frame(model$prediction)
      predictionTrain$evaluationType <- "Train"

      predictionTestEval <- as.data.frame(predictionTest)
      predictionTestEval$evaluationType <- "Test"

      combinedPrediction <- dplyr::bind_rows(
        predictionTrain,
        predictionTestEval
        )

      # Evaluate both train and test predictions
      evaluation <- PatientLevelPrediction::evaluatePlp(
        prediction = combinedPrediction,
        typeColumn = "evaluationType"
        )

      # Return a standardized object when a run has finished successfully
      list(
        success = TRUE,
        model = model,
        predictionTest = predictionTest,
        evaluation = evaluation,
        errorMessage = NA_character_
        )
    },
    error = function(e) {
      errorMessage <- conditionMessage(e)

      # Preserve the error text so it can be inspected
      warning(
            paste0(
              "PLP model failed for imputation '", imputation,
              "' and prediction model '", predModel,
              "': ", errorMessage
            ),
            call. = FALSE
          )
      # Use the same object structure as for a successful run
          list(
            success = FALSE,
            model = NULL,
            predictionTest = NULL,
            evaluation = NULL,
            errorMessage = errorMessage
          )
        }
      )

      # Store detailed information for this scenario
      results[[index]] <- list(
        imputation = imputation,
        predictionModel = predModel,
        plpModel = fitResult$model,
        predictionTest = fitResult$predictionTest,
        evaluation = fitResult$evaluation,
        success = fitResult$success,
        errorMessage = fitResult$errorMessage
      )

      # Exctract evaluation metrics when fitting and evaluation was successful
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

  # Return the model results and a combined metrics table
  list(
    modelResults = results,
    metrics = dplyr::bind_rows(metrics)
  )
}

#' @title Summarize PLP model run status
#' @description
#' Examines the output of `runPlpModels()` and determines whether the requested model runs completed
#' successfully, partially failed (e.g. for one of the imputation methods/mechanisms/ratios) or all failed
#' @param modelResults Result object returned by `runPlpModels()`
#' @return A list containing the completion status and the reported error messages
summariseModelRunStatus <- function(modelResults) {
  # Retain only model/imputation combinations whose execution failed
  failedRuns <- Filter(function(x) !isTRUE(x$success), modelResults$modelResults)

  # No failures measns the requested model runs completed successfully
  if (length(failedRuns) == 0) {
    return(list(
      status = "completed",
      modelErrors = NA_character_
    ))
  }

  # Format each error with its model and imputation identifiers to trace back the relevant model run
  modelErrors <- vapply(
    failedRuns,
    function(x) {
      paste0(
        x$predictionModel,
        "[",
        x$imputation,
        "]: ",
        x$errorMessage
      )
    },
    character(1)
  )

  # Distinguish between complete fails and cases where at least one model still completed successfully
  status <- if (length(failedRuns) == length(modelResults$modelResults)) {
    "failed"
  } else {
    "partial_failure"
  }

  # Collapse individual model errors into a single field
  list(status = status, modelErrors = paste(modelErrors, collapse = "|"))
}
