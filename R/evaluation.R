extractEvaluationObjects <- function(evaluation,
                                     simulation,
                                     mechanism,
                                     ratio,
                                     type,
                                     targetVariables,
                                     causeVariables,
                                     imputation,
                                     predictionModel) {
  if (is.null(evaluation)) {
    return(list(
      evaluationStatistics = NULL,
      calibrationSummary = NULL,
      thresholdSummary = NULL,
      demographicSummary = NULL,
      predictionDistribution = NULL
    ))
  }

  addTable <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }

    as.data.frame(x) %>%
      dplyr::mutate(
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

collectEvaluationTables <- function(modelResults,
                                    simulation,
                                    mechanism,
                                    ratio,
                                    type,
                                    targetVariables,
                                    causeVariables) {
  allTables <- list(
    evaluationStatistics = list(),
    calibrationSummary = list(),
    thresholdSummary = list(),
    demographicSummary = list(),
    predictionDistribution = list()
  )

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

  lapply(allTables, dplyr::bind_rows)
}

appendCsvTable <- function(table, filePath) {
  if (is.null(table) || nrow(table) == 0) {
    return(invisible(NULL))
  }

  dir.create(dirname(filePath), recursive = TRUE, showWarnings = FALSE)

  readr::write_csv(
    table,
    file = filePath,
    append = file.exists(filePath),
    col_names = !file.exists(filePath)
  )

  invisible(filePath)
}

evaluationStatsHelper <- function(table) {
  if (is.null(table) || nrow(table) == 0) {
    return(table)
  }

  table %>%
    dplyr::mutate(
      evaluation = as.character(.data$evaluation),
      metric = as.character(.data$metric),
      value = readr::parse_double(as.character(.data$value))
    )
}

appendScenarioResults <- function(evaluationTables,
                                  metrics = NULL,
                                  progress = NULL,
                                  folder) {
  rawFolder <- file.path(folder, "raw")

  for (tableName in names(evaluationTables)) {
    table <- evaluationTables[[tableName]]

    if (tableName == "evaluationStatistics" && !is.null(table) && nrow(table) > 0) {
      table <- evaluationStatsHelper(table)
    }

    appendCsvTable(table, file.path(rawFolder, paste0(tableName, ".csv")))
  }

  appendCsvTable(metrics, file.path(rawFolder, "metrics.csv"))
  appendCsvTable(progress, file.path(rawFolder, "scenarioProgress.csv"))
}

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

summariseOverSims <- function(table) {
  if (is.null(table) || nrow(table) == 0) {
    return(table)
  }

  if (all(c("evaluation", "metric", "value") %in% names(table))) {
    table <- evaluationStatsHelper(table)
  }

  groupVars <- intersect(
    c(
      "simulation",
      "mechanism",
      "ratio",
      "type",
      "imputation",
      "targetVariables",
      "causeVariables",
      "predictionModel",
      "evaluation",
      "evaluationType",
      "metric",
      "category",
      "threshold"
    ),
    names(table)
  )

  numericV <- table %>%
    dplyr::select(dplyr::where(is.numeric)) %>%
    names()

  numericV <- setdiff(numericV, c("simulation", groupVars))

  table %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(setdiff(groupVars, "simulation")))) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(numericV),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ stats::sd(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      nSimulation = dplyr::n_distinct(.data$simulation),
      .groups = "drop"
    )
}

saveSummarisedResults <- function(aggregatedTables, folder) {
  mainFolder <- file.path(folder, "master")
  specificFolder <- file.path(folder, "specifications")

  dir.create(mainFolder, recursive = TRUE, showWarnings = FALSE)
  dir.create(specificFolder, recursive = TRUE, showWarnings = FALSE)

  for (tableName in names(aggregatedTables)) {
    table <- aggregatedTables[[tableName]]

    if (is.null(table) || nrow(table) == 0) {
      next
    }

    readr::write_csv(table, file.path(mainFolder, paste0(tableName, ".csv")))

    if (!all(c("mechanism", "ratio", "type", "imputation", "predictionModel") %in% names(table))) {
      next
    }

    configurations <- table %>%
      dplyr::distinct(
        .data$mechanism,
        .data$ratio,
        .data$type,
        .data$imputation,
        .data$predictionModel
      )

    for (i in seq_len(nrow(configurations))) {
      config <- configurations[i, ]

      configTab <- table %>%
        dplyr::filter(
          .data$mechanism == config$mechanism,
          .data$ratio == config$ratio,
          .data$type == config$type,
          .data$imputation == config$imputation,
          .data$predictionModel == config$predictionModel
        )

      outputDirectory <- file.path(
        specificFolder,
        config$mechanism,
        paste0("ratio_", config$ratio),
        config$imputation,
        config$predictionModel
      )

      dir.create(outputDirectory, recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(configTab, file.path(outputDirectory, paste0(tableName, ".csv")))
    }
  }
}

collapseIds <- function(x) {
  if (is.null(x)) {
    return(NA_character_)
  }

  paste(sort(unique(x)), collapse = "_")
}

extractPatternIds <- function(patterns, field) {
  ids <- unlist(lapply(patterns, function(pattern) pattern[[field]]), use.names = FALSE)

  if (length(ids) == 0) {
    return(NULL)
  }

  ids
}

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
                                       "sklearnIterative_withIndicator"
                                     ),
                                     predictionModels = c("lasso", "xgboost"),
                                     hyperparameterSettings = PatientLevelPrediction::createHyperparameterSettings(),
                                     runs = 10,
                                     outputFolder = tempdir(),
                                     appendIntermediateCsv = TRUE,
                                     startSimulation = 1L,
                                     startScenarioId = 1L,
                                     seed = 123L) {
  results <- list()
  count <- 1

  if (is.null(completeInputCovariateIds)) {
    completeInputCovariateIds <- defaultRequiredCovariates(targetCovariateId, causeCovariateIds)
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

  for (simulation in seq_len(runs)) {
    if (simulation < startSimulation) {
      next
    }

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
            switch(
              mech,
              MCAR = list(targetCovariateIds = targetCovariateId, mechanism = "MCAR"),
              MAR = list(
                targetCovariateIds = targetCovariateId,
                causeCovariateIds = causeCovariateIds,
                mechanism = "MAR",
                type = type
              ),
              MNAR = list(
                targetCovariateIds = targetCovariateId,
                causeCovariateIds = targetCovariateId,
                mechanism = "MNAR",
                type = type
              )
            )
          )
        }

        targetVariables <- collapseIds(extractPatternIds(patternsSettings, "targetCovariateIds"))
        causeVariables <- collapseIds(extractPatternIds(patternsSettings, "causeCovariateIds"))

        trainMissingSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "trainMissing")
        testMissingSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "testMissing")
        sklearnSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "sklearnImputer")
        lassoSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "lasso")
        xgboostSeed <- makeSeed(seed, simulationId = simulation, scenarioId = scenarioId, seedLabel = "xgboost")

        modelOverview <- createPLPPredictionModels(
          lassoSeed = lassoSeed,
          xgboostSeed = xgboostSeed
        )[predictionModels]

        if (is.null(imputationOverview)) {
          imputationOverviewCurrent <- createPLPImputationMethods(
            includeMissingIndicator = c(FALSE, TRUE),
            sklearnRandomState = sklearnSeed
          )
        } else {
          imputationOverviewCurrent <- imputationOverview
        }

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

        if (isTRUE(completeCase)) {
          completeCaseResult <- completeCasePlp(
            trainData = trainMissingData,
            testData = testMissingData,
            targetCovariateIds = targetCovariateId
          )

          imputationResults <- c(list(completeCase = completeCaseResult), imputationResults)
        }

        modelResults <- runPlpModels(
          imputationResults = imputationResults,
          modelOverview = modelOverview,
          hyperparameterSettings = hyperparameterSettings,
          analysisPrefix = paste0("sim", simulation, "_", mech, "_ratio", ratio),
          analysisPath = file.path(tempdir(), "missingnessSimulation")
        )

        evaluationTables <- collectEvaluationTables(
          modelResults = modelResults,
          simulation = simulation,
          mechanism = mech,
          ratio = ratio,
          type = type,
          targetVariables = targetVariables,
          causeVariables = causeVariables
        )

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
          lassoSeed = lassoSeed,
          xgboostSeed = xgboostSeed,
          status = "completed"
        )

        if (isTRUE(appendIntermediateCsv)) {
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

        count <- count + 1
      }
    }
  }

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
