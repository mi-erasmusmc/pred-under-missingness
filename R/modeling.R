createPLPPredictionModels <- function(lassoSeed = 12L,
                                      xgboostSeed = 12L) {
  list(
    lasso = PatientLevelPrediction::setLassoLogisticRegression(seed = lassoSeed),
    xgboost = PatientLevelPrediction::setGradientBoostingMachine(
      seed = xgboostSeed
    )
  )
}

getPlpMetrics <- function(evaluation, imputation, predictionModel) {
  evaluation$evaluationStatistics %>%
    dplyr::mutate(
      evaluation = as.character(.data$evaluation),
      metric = as.character(.data$metric),
      value = as.numeric(.data$value),
      imputation = imputation,
      predictionModel = predictionModel
    )
}

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
  if (is.null(modelOverview)) {
    modelOverview <- createPLPPredictionModels()
  }

  results <- list()
  metrics <- list()
  index <- 1

  dir.create(analysisPath, recursive = TRUE, showWarnings = FALSE)

  for (imputation in names(imputationResults)) {
    imputationResult <- imputationResults[[imputation]]

    for (predModel in names(modelOverview)) {
      trainData <- safeCopyPlpData(imputationResult$trainImputed)
      testData <- imputationResult$testImputed

      trainData$covariateData <- PatientLevelPrediction::preprocessData(
        covariateData = trainData$covariateData,
        preprocessSettings = preprocessSettings
      )

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

  list(
    modelResults = results,
    metrics = dplyr::bind_rows(metrics)
  )
}
