resolveXgboostTuningMethodSettings <- function(methodName, imputationOverviewCurrent) {
  if (!methodName %in% names(imputationOverviewCurrent)) {
    return(NULL)
  }
  
  methodSettings <- imputationOverviewCurrent[[methodName]]
  if (!inherits(methodSettings, "featureEngineeringSettings")) {
    return(NULL)
  }
  
  if (!identical(attr(methodSettings, "fun"), "implementXgboostImputer")) {
    return(NULL)
  }
  
  methodSettings
}

maybeTuneXgboostInSimulation <- function(trainMissingData,
                                         imputationOverviewCurrent,
                                         imputationMethods,
                                         xgboostTuningGrid = NULL,
                                         xgboostValidationFraction = 0.20,
                                         xgboostMinHoldoutPerColumn = 1,
                                         xgboostTuningMetric = "rmse",
                                         xgboostTuningWorkers = 1L,
                                         xgboostSourceDir = NULL,
                                         xgboostVerbose = FALSE) {
  sparseMethods <- imputationMethods[grepl("^xgboost_", imputationMethods)]
  
  if (length(sparseMethods) == 0 || is.null(xgboostTuningGrid)) {
    return(list(
      imputationOverview = imputationOverviewCurrent, 
      tuning = list()
    ))
  }
  
  updatedOverview <- imputationOverviewCurrent
  tuningResults <- vector("list", length(sparseMethods))
  names(tuningResults) <- sparseMethods
  
  for(methodName in sparseMethods) {
    methodSettings <- resolveXgboostTuningMethodSettings(
      methodName = methodName,
      imputationOverviewCurrent = updatedOverview
    )
    
    if (is.null(methodSettings)) {
      next
    }
    
    if (isTRUE(xgboostVerbose)) {
      message(
        "Sparse xgboost tuning inside simulation for method '",
        methodName,
        "'."
      )
    }
    
    tuning <- 
      suppressWarnings(
        tuneXgboostImputer(
          trainData = trainMissingData,
          grid = xgboostTuningGrid,
          validationFraction = xgboostValidationFraction,
          minHoldoutPerColumn = xgboostMinHoldoutPerColumn,
          metric = xgboostTuningMetric,
          numWorkers= xgboostTuningWorkers,
          sourceDir = xgboostSourceDir,
          missingThreshold = methodSettings$missingThreshold,
          nrounds = methodSettings$nrounds,
          maxDepth = methodSettings$maxDepth,
          eta = methodSettings$eta,
          subsample = methodSettings$subsample,
          colsampleBytree = methodSettings$colsampleBytree,
          minChildWeight = methodSettings$minChildWeight,
          numThreads = methodSettings$numThreads,
          includeSparseZeroPredictors = methodSettings$includeSparseZeroPredictors,
          addMissingIndicator = methodSettings$addMissingIndicator,
          seed = methodSettings$seed,
          verbose = xgboostVerbose
          
        )
        
      )
    
    bestSettings <- tuning$bestSettings
    updatedOverview[[methodName]] <- bestSettings
    tuning$selectedMethodName <- methodName
    tuningResults[[methodName]] <- tuning
  }
  
  list(imputationOverview = updatedOverview, tuning = tuningResults)
}