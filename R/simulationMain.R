#' @title Run a missingness/imputation/prediction simulation study
#' @description
#'  This function runs a full simulation workflow on PLP data. For each 
#'  simulation run and missingness scenario, it creates the incomplete training and test data,
#'  applies the imputation methods, fits the prediction models and stores the resulting 
#'  evaluation outputs
#'  @param data A complete PLP data object
#'  @param population The study population used to create train/test splits
#'  @param targetCovariateId The target covariate identifier(s) for missingness simulation
#'  @param causeCovariateIds Optional covariate identifiers driving missingness
#'  @param completeInputCovariateIds Optional set of covariates that must be fully observed
#'  in the input data
#'  @param completeCase Logical parameter indicating whether complete-case analysis should be run
#'  @param mechanisms The missingness mechanisms to simulate
#'  @param missingnessRatios The missingness ratios to simulate
#'  @param patterns Optional explicit missingness patterns
#'  @param freq Optional pattern frequencies
#'  @param targetCovariateIdPerMech Optional mechanism-specific target covariate identifiers
#'  @param causeCovariateIdsPerMech Optional mechanism-specific cause covariate identifiers
#'  @param typePerMech Optional mechanism-specific score transformation types
#'  @param imputationOverview Optional imputation method settigng
#'  @param imputationMethods Optional list of imputation methods to run
#'  @param xgboostTuningGrid Candidate hyperparameters to tune for iterative XGBoost imputation
#' @param xgboostValidationFraction Fraction of observed target cells to use as validation set during tuning
#' @param minHoldoutPerColumn Minimum number of cells to hold out per imputatble covariate during tuning
#' @param xgboostTuningMetric Metric used for tuning
#' @param xgboostTuningWorkers The number of parallel workers to use for iterative XGBoost imputation tuning
#' @param xgboostSourceDir Source of R-package files
#' @param xgboostVerbose Logical parameter indicating whether to print extra tuning progress messages
#'  @param predictionModels The prediction models to run
#'  @param hyperparameterSettings The hyperparameter settings used for the prediction models
#'  @param runs The number of simulation runs
#'  @param outputFolder Folder directory to store evaluation objects
#' @param analysisPath Path to write model fit results to
#' @param appendIntermediateCsv Logical parameter to control whether scenario results are written as intermediate raw csv files
#' @param savePlpModelResults Logical parameter denoting whether PLP models should be saved on disk
#' @param plpResultsFolder Folder directory to store PLP model results
#' @param xgboostThreads Optional parameter to set number of threads for PLP xgboost prediction
#' @param transformerDevice Optional parameter to set the device to run the transformer model on
#' @param startSimulation The simulation run to resume the simulation 
#' @param startScenarioId The scenario within a simulation run to resume the simulation
#' @param seed The master seed used to derive run-specific seeds
#'  
#' @export
runMissingnessSimulation <- function(data,
                                     population,
                                     targetCovariateId,
                                     causeCovariateIds = NULL,
                                     completeInputCovariateIds = NULL,
                                     completeCase = FALSE,
                                     mechanisms = c("MCAR", "MAR", "MNAR"),
                                     missingnessRatios = seq(0, 0.8, by = 0.2),
                                     patterns = NULL,
                                     freq = NULL,
                                     targetCovariateIdPerMech = NULL,
                                     causeCovariateIdsPerMech = NULL,
                                     typePerMech = list(MAR = "RIGHT", MNAR = "RIGHT"),
                                     imputationOverview = NULL,
                                     imputationMethods = c(
                                       "simpleMean_noIndicator",
                                       "simpleMean_withIndicator",
                                       "simpleMedian_noIndicator",
                                       "simpleMedian_withIndicator",
                                       "iterativePMM_noIndicator",
                                       "iterativePMM_withIndicator",
                                       "sklearnIterative_noIndicator",
                                       "sklearnIterative_withIndicator",
                                       "xgboost_noIndicator",
                                       "xgboost_withIndicator"
                                     ),
                                     xgboostTuningGrid = NULL,
                                     xgboostValidationFraction = 0.20,
                                     xgboostMinHoldoutPerColumn = 1,
                                     xgboostTuningMetric = "rmse",
                                     xgboostTuningWorkers = 1L,
                                     xgboostSourceDir = NULL,
                                     xgboostVerbose = FALSE,
                                     predictionModels = c("lasso", "xgboost", "transformer"),
                                     hyperparameterSettings = PatientLevelPrediction::createHyperparameterSettings(),
                                     runs = 10,
                                     outputFolder = tempdir(),
                                     analysisPath = NULL,
                                     appendIntermediateCsv = TRUE,
                                     savePlpModelResults = FALSE,
                                     plpResultsFolder = NULL,
                                     xgboostThreads = NULL,
                                     transformerDevice = "cuda:0",
                                     startSimulation = 1L,
                                     startScenarioId = 1L,
                                     seed = 123L) {
  results <- list()
  count <- 1
  
  # Initialize output folders
  if (is.null(analysisPath)) {
    analysisPath <- file.path(outputFolder, "plpModels")
  }
  
  if (is.null(plpResultsFolder)) {
    plpResultsFolder <- file.path(outputFolder, "plpResults")
  }
  
  # Resolve the required complete data inputs
  if (is.null(completeInputCovariateIds)) {
    completeInputCovariateIds <- defaultSimulationInputs(
      targetCovariateId = targetCovariateId,
      causeCovariateIds = causeCovariateIds,
      targetCovariateIdPerMech = targetCovariateIdPerMech,
      causeCovariateIdsPerMech = causeCovariateIdsPerMech,
      patterns = patterns
    )
  }
  
  if (startSimulation < 1L) {
    stop("startSimulation must be at least 1.")
  }
  
  if (startSimulation > runs) {
    stop("startSimulation cannot be larger than runs.")
  }
  
  if (startScenarioId < 1L) {
    stop("startScenarioId must be at least 1.")
  }
  
  # Start the simulation
  for (simulation in seq_len(runs)) {
    if (simulation < startSimulation) {
      next
    }
    
    # Build the train/test split for each simulation run
    splitSeed <- makeSeed(seed, simulationId = simulation, seedLabel = "split")
    splitSettings <- PatientLevelPrediction::createDefaultSplitSetting(
      testFraction = 0.25,
      trainFraction = 0.75,
      nfold = 3,
      splitSeed = splitSeed,
      type = "stratified"
    )
    
    splitPlpData <- PatientLevelPrediction::splitData(
      plpData = data,
      population = population,
      splitSettings = splitSettings
    )
    
    trainData <- splitPlpData$Train
    testData <- splitPlpData$Test
    
    # Check whether the train and test data are complete
    checkPlpDataComplete(trainData, completeInputCovariateIds)
    checkPlpDataComplete(testData, completeInputCovariateIds)
    
    scenarioId <- 0L
    
    # If no scenarioId and startSimulation are given, start from scenario 1
    firstScenarioThisSimulation <- if (simulation == startSimulation) startScenarioId else 1L
    
    # Iterate over the requested missingness scenarios
    for (mech in mechanisms) {
      for (ratio in missingnessRatios) {
        mech <- toupper(mech)
        scenarioId <- scenarioId + 1L
        
        if (scenarioId < firstScenarioThisSimulation) {
          next
        }
        
        type <- resolvePerMechanismValue(typePerMech, mech, defaultValue = "RIGHT")
        scenarioTargetCovariateId <- resolvePerMechanismValue(
          targetCovariateIdPerMech,
          mech,
          defaultValue = targetCovariateId
        )
        
        scenarioCauseCovariateIds <- resolvePerMechanismValue(
          causeCovariateIdsPerMech,
          mech,
          defaultValue = causeCovariateIds
        )
        
        if (!is.null(patterns)) {
          patternsSettings <- patterns
        } else {
          patternsSettings <- list(
            buildMechanismPattern(
              mech = mech,
              targetCovariateId = scenarioTargetCovariateId,
              causeCovariateIds = scenarioCauseCovariateIds,
              type = type
            )
          )
        }
        
        scenarioTargetCovariateIds <- extractPatternIds(patternsSettings, "targetCovariateIds")
        targetVariables <- collapseIds(scenarioTargetCovariateIds)
        causeVariables <- collapseIds(extractPatternIds(patternsSettings, "causeCovariateIds"))
        
        # Create seeds from the master seed
        trainMissingSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "trainMissing")
        testMissingSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "testMissing")
        sklearnSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "sklearnImputer")
        iterativeXGBoostSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "iterativeXGBoostImputer")
        lassoSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "lasso")
        xgboostSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "xgboost")
        transformerSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "transformer")
        pmmTrainSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "pmmTrain")
        pmmTestSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "pmmTest")
        
        # Initialize PLP prediction models
        modelOverview <- createPLPPredictionModels(
          lassoSeed = lassoSeed,
          xgboostSeed = xgboostSeed,
          transformerSeed = transformerSeed,
          transformerDevice = transformerDevice
        )[predictionModels]
        
        # Initialize requested imputation methods
        if (is.null(imputationOverview)) {
          imputationOverviewCurrent <- createPLPImputationMethods(
            includeMissingIndicator = c(FALSE, TRUE),
            sklearnRandomState = sklearnSeed,
            iterativeXGBoostSeed = iterativeXGBoostSeed,
            pmmTrainSeed = pmmTrainSeed,
            pmmTestSeed = pmmTestSeed
          )
        } else {
          imputationOverviewCurrent <- applyScenarioSeedsToImputationOverview(
            imputationOverview = imputationOverview,
            iterativeXGBoostSeed = iterativeXGBoostSeed
          )
        }
        
        # Simulate missingness in the train and test data independently
        simulateTrain <- simulateMultivariateMissingness(
          data = trainData,
          patterns = patternsSettings,
          ratio = ratio,
          freq = freq,
          seed = trainMissingSeed
        )
        
        targetId <- targetCovariateId
        
        trainDroppedTarget <- simulateTrain$droppedPairs %>%
          dplyr::filter(.data$covariateId == targetId) %>%
          dplyr::summarise(nDropped = dplyr::n()) %>%
          dplyr::pull(.data$nDropped)
        
        trainMissingReport <- reportMissingness(
          simulateTrain$plpData,
          covariateIds = targetId
        )
        
        message(
          "Train target check: covariateId = ", targetId,
          ", droppedPairs = ", trainDroppedTarget,
          ", nObserved = ", trainMissingReport$nObserved,
          ", nMissing = ", trainMissingReport$nMissing,
          ", missingFraction = ", signif(trainMissingReport$missingFraction, 4)
        )
        
        
        simulateTest <- simulateMultivariateMissingness(
          data = testData,
          patterns = patternsSettings,
          ratio = ratio,
          freq = freq,
          seed = testMissingSeed
        )
        
        testDroppedTarget <- simulateTest$droppedPairs %>%
          dplyr::filter(.data$covariateId == targetId) %>%
          dplyr::summarise(nDropped = dplyr::n()) %>%
          dplyr::pull(.data$nDropped)
        
        testMissingReport <- reportMissingness(
          simulateTest$plpData,
          covariateIds = targetId
        )
        
        message(
          "Test target check: covariateId = ", targetId,
          ", droppedPairs = ", testDroppedTarget,
          ", nObserved = ", testMissingReport$nObserved,
          ", nMissing = ", testMissingReport$nMissing,
          ", missingFraction = ", signif(testMissingReport$missingFraction, 4)
        )
        
        trainMissingData <- simulateTrain$plpData
        testMissingData <- simulateTest$plpData
        
        # Optional: tune iterative xgboost imputation
        xgboostTuning <- maybeTuneXgboostInSimulation(
          trainMissingData = trainMissingData,
          imputationOverviewCurrent = imputationOverviewCurrent,
          imputationMethods = imputationMethods,
          xgboostTuningGrid = xgboostTuningGrid,
          xgboostValidationFraction = xgboostValidationFraction,
          xgboostMinHoldoutPerColumn = xgboostMinHoldoutPerColumn,
          xgboostTuningMetric = xgboostTuningMetric,
          xgboostTuningWorkers = xgboostTuningWorkers,
          xgboostSourceDir = xgboostSourceDir,
          xgboostVerbose = xgboostVerbose
        )
        
        imputationOverviewCurrent <- xgboostTuning$imputationOverview
        
        # Impute the train and test data using the requested methods
        imputationResults <- runPlpImputers(
          trainData = trainMissingData,
          testData = testMissingData,
          imputationMethods = imputationMethods,
          imputationOverview = imputationOverviewCurrent
        )
        
        if (isTRUE(completeCase)) {
          completeCaseResult <- completeCasePlp(
            trainData = trainMissingData,
            testData = testMissingData,
            targetCovariateIds = targetCovariateId
          )
          
          imputationResults <- c(list(completeCase = completeCaseResult), imputationResults)
        }
        
        # Run PLP prediction models
        modelResults <- runPlpModels(
          imputationResults = imputationResults,
          modelOverview = modelOverview,
          hyperparameterSettings = hyperparameterSettings,
          analysisPrefix = paste0("sim", simulation, "_", mech, "_ratio", ratio),
          analysisPath = analysisPath
        )
        
        # Collect evaluation tables for the scenario
        evaluationTables <- collectEvaluationTables(
          modelResults = modelResults,
          simulation = simulation,
          mechanism = mech,
          ratio = ratio,
          type = type,
          targetVariables = targetVariables,
          causeVariables = causeVariables
        )
        
        # Track the scenario progress
        modelRunStatus <- summariseModelRunStatus(modelResults)
        
        scenarioProgress <- tibble::tibble(
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
          iterativeXGBoostSeed = iterativeXGBoostSeed,
          pmmTrainSeed = pmmTrainSeed,
          pmmTestSeed = pmmTestSeed,
          lassoSeed = lassoSeed,
          xgboostSeed = xgboostSeed,
          transformerSeed = transformerSeed,
          status = modelRunStatus$status,
          modelErrors = modelRunStatus$modelErrors
        )
        
        # Append new rows of new scenarios
        if (isTRUE(appendIntermediateCsv)) {
          appendScenarioResults(
            evaluationTables = evaluationTables,
            scenarioStatus = modelRunStatus$status,
            modelErros = modelRunStatus$modelErrors,
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
        
        # Create a list of the results of the scenario
        results[[count]] <- list(
          simulation = simulation,
          mechanism = mech,
          ratio = ratio,
          type = type,
          targetVariables = targetVariables,
          causeVariables = causeVariables,
          seeds = list(
            seed = seed,
            split = splitSeed,
            trainMissing = trainMissingSeed,
            testMissing = testMissingSeed,
            sklearnImputer = sklearnSeed,
            iterativeXGBoostImputer = iterativeXGBoostSeed,
            pmmTrain = pmmTrainSeed,
            pmmTest = pmmTestSeed,
            lasso = lassoSeed,
            xgboost = xgboostSeed,
            transformer = transformerSeed
          ),
          trainMissingSummary = simulateTrain$patternSummary,
          testMissingSummary = simulateTest$patternSummary,
          xgboostTuning = xgboostTuning$tuning,
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
        
        # Save plp model results if requested
        if (isTRUE(savePlpModelResults)) {
          dir.create(plpResultsFolder, recursive = TRUE, showWarnings = FALSE)
          
          modelResultsFile <- file.path(
            plpResultsFolder,
            paste0(
              "simulation_", simulation,
              "_scenario", scenarioId,
              "_", mech,
              "_ratio_", ratio,
              ".rds"
            )
          )
          
          saveRDS(results[[count]]$modelResults, file = modelResultsFile)
          results[[count]]$modelResultsFile <- modelResultsFile
        }
        
        count <- count + 1
      }
    }
  }
  
  # Aggregate evaluation tables after all simulations
  if (isTRUE(appendIntermediateCsv)) {
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
  saveSummarisedResults(aggregatedTables = aggregatedTables, folder = outputFolder)
  
  list(
    results = results,
    allEvaluationTables = allEvaluationTables,
    aggregatedTables = aggregatedTables
  )
}

#' @title Run a multivariate missingness/imputation/prediction simulation study
#' @description
#'  This function runs a full simulation workflow on PLP data, 
#'  using fixed missingness patterns and target pattern proportions. For each 
#'  simulation run, it splits the PLP data in to a train and test set, applies the pattern missingness design, 
#'  runs the selected imputation methods and prediction models and stores the resulting evaluation outputs.
#'  @param data A complete PLP data object
#'  @param population The study population used to create train/test splits
#'  @param patterns A list of pattern specifications for multivariate missingness
#'  @param freq The desired final proportions for the missingness patterns
#'  @param completeInputCovariateIds Optional set of covariates that must be fully observed
#'  in the input data
#'  @param completeCase Logical parameter indicating whether complete-case analysis should be run
#'  @param imputationOverview Optional imputation method settigng
#'  @param imputationMethods Optional list of imputation methods to run
#'  @param xgboostTuningGrid Candidate hyperparameters to tune for iterative XGBoost imputation
#' @param xgboostValidationFraction Fraction of observed target cells to use as validation set during tuning
#' @param minHoldoutPerColumn Minimum number of cells to hold out per imputatble covariate during tuning
#' @param xgboostTuningMetric Metric used for tuning
#' @param xgboostTuningWorkers The number of parallel workers to use for iterative XGBoost imputation tuning
#' @param xgboostSourceDir Source of R-package files
#' @param xgboostVerbose Logical parameter indicating whether to print extra tuning progress messages
#'  @param predictionModels The prediction models to run
#'  @param hyperparameterSettings The hyperparameter settings used for the prediction models
#'  @param runs The number of simulation runs
#'  @param outputFolder Folder directory to store evaluation objects
#' @param analysisPath Path to write model fit results to
#' @param appendIntermediateCsv Logical parameter to control whether scenario results are written as intermediate raw csv files
#' @param savePlpModelResults Logical parameter denoting whether PLP models should be saved on disk
#' @param plpResultsFolder Folder directory to store PLP model results
#' @param xgboostThreads Optional parameter to set number of threads for PLP xgboost prediction
#' @param transformerDevice Optional parameter to set the device to run the transformer model on
#' @param startSimulation The simulation run to resume the simulation 
#' @param seed The master seed used to derive run-specific seeds
#' @param mechanismLabel The label used to identify the pattern-matched missingness mechanism in the output tables
#' @param typeLabel The label used to identify the pattern-matched missingness type in the output label
#' @export
runDynamicMissingnessSimulation <- function(data,
                                            population,
                                            patterns,
                                            freq,
                                            completeInputCovariateIds = NULL,
                                            completeCase = FALSE,
                                            imputationOverview = NULL,
                                            imputationMethods = c(
                                              "simpleMean_noIndicator",
                                              "simpleMean_withIndicator",
                                              "simpleMedian_noIndicator",
                                              "simpleMedian_withIndicator",
                                              "iterativePMM_noIndicator",
                                              "iterativePMM_withIndicator",
                                              "sklearnIterative_noIndicator",
                                              "sklearnIterative_withIndicator",
                                              "xgboost_noIndicator",
                                              "xgboost_withIndicator"
                                            ),
                                            xgboostTuningGrid = NULL,
                                            xgboostValidationFraction = 0.20,
                                            xgboostMinHoldoutPerColumn = 1,
                                            xgboostTuningMetric = "rmse",
                                            xgboostTuningWorkers = 1L,
                                            xgboostSourceDir = NULL,
                                            xgboostVerbose = FALSE,
                                            predictionModels = c("lasso", "xgboost", "transformer"),
                                            hyperparameterSettings = PatientLevelPrediction::createHyperparameterSettings(),
                                            runs = 10,
                                            outputFolder = tempdir(),
                                            analysisPath = NULL,
                                            appendIntermediateCsv = TRUE,
                                            savePlpModelResults = FALSE,
                                            plpResultsFolder = NULL,
                                            xgboostThreads = NULL,
                                            transformerDevice = "cuda:0",
                                            startSimulation = 1L,
                                            seed = 123L,
                                            mechanismLabel = "DYNAMIC_MAR",
                                            typeLabel = "DYNAMIC_PATTERN") {
  # Initialize output folders and simulation identifiers
  results <- list()
  count <- 1
  scenarioId <- 1
  ratio <- 1
  
  if (is.null(analysisPath)) {
    analysisPath <- file.path(outputFolder, "plpModels")
  }
  
  if (is.null(plpResultsFolder)) {
    plpResultsFolder <- file.path(outputFolder, "plpResults")
  }
  
  # Resolve require complete data inputs
  if (is.null(completeInputCovariateIds)) {
    completeInputCovariateIds <- defaultSimulationInputs(
      targetCovariateId = NULL,
      causeCovariateIds = NULL,
      patterns = patterns
    )
  }
  
  if (startSimulation < 1L) {
    stop("startSimulation must be at least 1.")
  }
  
  if (startSimulation > runs) {
    stop("startSimulation cannot be larger than runs.")
  }
  
  # Derive output labels for the target and cause covariates
  scenarioTargetCovariateIds <- extractPatternIds(patterns, "targetCovariateIds")
  targetVariables <- collapseIds(scenarioTargetCovariateIds)
  causeVariables <- collapseIds(extractPatternIds(patterns, "causeCovariateds"))
  
  # Run the simulation
  for (simulation in seq_len(runs)) {
    if (simulation < startSimulation) {
      next
    }
    
    # Build the train/test split for each simulation run
    splitSeed <- makeSeed(seed, simulationId = simulation, seedLabel = "split")
    splitSettings <- PatientLevelPrediction::createDefaultSplitSetting(
      testFraction = 0.25,
      trainFraction = 0.75,
      nfold = 3,
      splitSeed = splitSeed,
      type = "stratified"
    )
    
    splitPlpData <- PatientLevelPrediction::splitData(
      plpData = data,
      population = population,
      splitSettings = splitSettings
    )
    
    trainData <- splitPlpData$Train
    testData <- splitPlpData$Test
    
    # Check whether the train and test data are complete
    checkPlpDataComplete(trainData, completeInputCovariateIds)
    checkPlpDataComplete(testData, completeInputCovariateIds)
    
    # Create seeds from the master seed
    trainMissingSeed <- makeSeed(
      seed,
      simulationId = simulation,
      scenarioId = scenarioId,
      seedLabel = "trainMissing"
    )
    testMissingSeed <- makeSeed(
      seed,
      simulationId = simulation,
      scenarioId = scenarioId,
      seedLabel = "testMissing"
    )
    sklearnSeed <- makeSeed(
      seed,
      simulationId = simulation,
      scenarioId = scenarioId,
      seedLabel = "sklearnImputer"
    )
    iterativeXGBoostSeed <- makeSeed(
      seed,
      simulationId = simulation,
      scenarioId = scenarioId,
      seedLabel = "iterativeXGBoostImputer"
    )
    lassoSeed <- makeSeed(
      seed,
      simulationId = simulation,
      scenarioId = scenarioId,
      seedLabel = "lasso"
    )
    xgboostSeed <- makeSeed(
      seed,
      simulationId = simulation,
      scenarioId = scenarioId,
      seedLabel = "xgboost"
    )
    transformerSeed <- makeSeed(
      seed,
      simulationId = simulation,
      scenarioId = scenarioId,
      seedLabel = "transformer"
    )
    pmmTrainSeed <- makeSeed(
      seed,
      simulationId = simulation,
      scenarioId = scenarioId,
      seedLabel = "pmmTrain"
    )
    pmmTestSeed <- makeSeed(
      seed,
      simulationId = simulation,
      scenarioId = scenarioId,
      seedLabel = "pmmTest"
    )
    
    # Initialize PLP prediction models
    modelOverview <- createPLPPredictionModels(
      lassoSeed = lassoSeed,
      xgboostSeed = xgboostSeed,
      transformerSeed = transformerSeed,
      transformerDevice = transformerDevice
    )[predictionModels]
    
    # Initialize requested imputation methods
    if (is.null(imputationOverview)) {
      imputationOverviewCurrent <- createPLPImputationMethods(
        includeMissingIndicator = c(FALSE, TRUE),
        sklearnRandomState = sklearnSeed,
        iterativeXGBoostSeed = iterativeXGBoostSeed,
        pmmTrainSeed = pmmTrainSeed,
        pmmTestSeed = pmmTestSeed
      )
    } else {
      imputationOverviewCurrent <- applyScenarioSeedsToImputationOverview(
        imputationOverview = imputationOverview,
        iterativeXGBoostSeed = iterativeXGBoostSeed
      )
    }
    
    # Simulate missingness in the train and test data independently
    simulateTrain <- simulateMultivariatePatternMissingness(
      data = trainData,
      patterns = patterns,
      freq = freq,
      seed = trainMissingSeed
    )
    
    simulateTest <- simulateMultivariatePatternMissingness(
      data = testData,
      patterns = patterns,
      freq = freq,
      seed = testMissingSeed
    )
    trainMissingData <- simulateTrain$plpData
    testMissingData <- simulateTest$plpData
    # Optional: tune iterative xgboost imputation
    xgboostTuning <- maybeTuneXgboostInSimulation(
      trainMissingData = trainMissingData,
      imputationOverviewCurrent = imputationOverviewCurrent,
      imputationMethods = imputationMethods,
      xgboostTuningGrid = xgboostTuningGrid,
      xgboostValidationFraction = xgboostValidationFraction,
      xgboostMinHoldoutPerColumn = xgboostMinHoldoutPerColumn,
      xgboostTuningMetric = xgboostTuningMetric,
      xgboostTuningWorkers = xgboostTuningWorkers,
      xgboostSourceDir = xgboostSourceDir,
      xgboostVerbose = xgboostVerbose
    )
    
    imputationOverviewCurrent <- xgboostTuning$imputationOverview
    
    imputationResults <- runPlpImputers(
      trainData = trainMissingData,
      testData = testMissingData,
      imputationMethods = imputationMethods,
      imputationOverview = imputationOverviewCurrent
    )
    
    if (isTRUE(completeCase)) {
      completeCaseResult <- completeCasePlp(
        trainData = trainMissingData,
        testData = testMissingData,
        targetCovariateIds = scenarioTargetCovariateIds
      )
      
      imputationResults <- c(list(completeCase = completeCaseResult),
                             imputationResults)
    }
    
    # Run PLP prediction models
    modelResults <- runPlpModels(
      imputationResults = imputationResults,
      modelOverview = modelOverview,
      hyperparameterSettings = hyperparameterSettings,
      analysisPrefix = paste0("sim", simulation, "_dynamicPattern"),
      analysisPath = analysisPath
    )
    
    # Collect evaluation tables for the scenario
    evaluationTables <- collectEvaluationTables(
      modelResults = modelResults,
      simulation = simulation,
      mechanism = mechanismLabel,
      ratio = ratio,
      type = typeLabel,
      targetVariables = targetVariables,
      causeVariables = causeVariables
    )
    
    # Track scenario progress
    modelRunStatus <- summariseModelRunStatus(modelResults)
    
    scenarioProgress <- tibble::tibble(
      simulation = simulation,
      scenarioId = scenarioId,
      resultIndex = count,
      mechanism = mechanismLabel,
      ratio = ratio,
      type = typeLabel,
      targetVariables = targetVariables,
      causeVariables = causeVariables,
      splitSeed = splitSeed,
      trainMissingSeed = trainMissingSeed,
      testMissingSeed = testMissingSeed,
      sklearnSeed = sklearnSeed,
      iterativeXGBoostSeed = iterativeXGBoostSeed,
      pmmTrainSeed = pmmTrainSeed,
      pmmTestSeed = pmmTestSeed,
      lassoSeed = lassoSeed,
      xgboostSeed = xgboostSeed,
      transformerSeed = transformerSeed,
      status = modelRunStatus$status,
      modelErrors = modelRunStatus$modelErrors
    )
    
    # Append new rows of new scenarios
    if (isTRUE(appendIntermediateCsv)) {
      appendScenarioResults(
        evaluationTables = evaluationTables,
        scenarioStatus = modelRunStatus$status,
        modelErros = modelRunStatus$modelErrors,
        metrics = modelResults$metrics %>%
          dplyr::mutate(
            simulation = simulation,
            scenarioId = scenarioId,
            mechanism = mechanismLabel,
            ratio = ratio,
            type = typeLabel,
            targetVariables = targetVariables,
            causeVariables = causeVariables
          ),
        progress = scenarioProgress,
        folder = outputFolder
      )
    }
    
    # Create a list of the results of the scenario
    results[[count]] <- list(
      simulation = simulation,
      mechanism = mechanismLabel,
      ratio = ratio,
      type = typeLabel,
      targetVariables = targetVariables,
      causeVariables = causeVariables,
      freq = freq,
      seeds = list(
        seed = seed,
        split = splitSeed,
        trainMissing = trainMissingSeed,
        testMissing = testMissingSeed,
        sklearnImputer = sklearnSeed,
        iterativeXGBoostImputer = iterativeXGBoostSeed,
        pmmTrain = pmmTrainSeed,
        pmmTest = pmmTestSeed,
        lasso = lassoSeed,
        xgboost = xgboostSeed,
        transformerSeed = transformerSeed
      ),
      trainMissingSummary = simulateTrain$patternSummary,
      testMissingSummary = simulateTest$patternSummary,
      xgboostTuning = xgboostTuning$tuning,
      imputationResults = imputationResults,
      modelResults = modelResults,
      evaluationTables = evaluationTables,
      scenarioStatus = modelRunStatus$status,
      modelErrors = modelRunStatus$modelErrors,
      metrics = modelResults$metrics %>%
        dplyr::mutate(
          simulation = simulation,
          scenarioId = scenarioId,
          mechanism = mechanismLabel,
          ratio = ratio,
          type = typeLabel,
          targetVariables = targetVariables,
          causeVariables = causeVariables
        )
    )
    
    # Save plp model results if requested
    if (isTRUE(savePlpModelResults)) {
      dir.create(plpResultsFolder,
                 recursive = TRUE,
                 showWarnings = FALSE)
      
      modelResultsFile <- file.path(
        plpResultsFolder,
        paste0(
          "simulation_",
          simulation,
          "_scenario_",
          scenarioId,
          "_dynamicPattern.rds"
        )
      )
      
      saveRDS(results[[count]]$modelResults, file = modelResultsFile)
      results[[count]]$modelResultsFile <- modelResultsFile
    }
    
    count <- count + 1
    
  }
  
  # Aggregate evaluation tables after all simulations
  if (isTRUE(appendIntermediateCsv)) {
    allEvaluationTables <- loadRawEvaluationTables(outputFolder)
  } else {
    allEvaluationTables <- list(
      evaluationStatistics = dplyr::bind_rows(
        lapply(results, function(x)
          x$evaluationTables$evaluationStatistics)
      ),
      calibrationSummary = dplyr::bind_rows(
        lapply(results, function(x)
          x$evaluationTables$calibrationSummary)
      ),
      thresholdSummary = dplyr::bind_rows(
        lapply(results, function(x)
          x$evaluationTables$thresholdSummary)
      ),
      demographicSummary = dplyr::bind_rows(
        lapply(results, function(x)
          x$evaluationTables$demographicSummary)
      ),
      predictionDistribution = dplyr::bind_rows(
        lapply(results, function(x)
          x$evaluationTables$predictionDistribution)
      )
    )
  }
  
  aggregatedTables <- lapply(allEvaluationTables, summariseOverSims)
  saveSummarisedResults(aggregatedTables = aggregatedTables, folder = outputFolder)
  
  list(
    results = results,
    allEvaluationTables = allEvaluationTables,
    aggregatedTables = aggregatedTables
  )
}
