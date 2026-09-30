validateXgboostMetric <- function(metric) {
  metric <- tolower(as.character(metric))
  if(length(metric) != 1 || !metric %in% c("rmse","mae")) {
    stop("metric must be 'rmse' or 'mae'.")
  }
  metric
}
 builddXgboostTuningHoldout <-function(prepared, validationFraction = 0.20, minHoldOutPerColumn = 1L, seed = 42L) {
   validationFraction <- validateMissForestNumeric(
     validationFraction,
     name = "validationFraction",
     lower = 0,
     upper = 1,
     lowerInclusive = FALSE,
     upperInclusive = FALSE
   )

   minHoldOutPerColumn <- validateMissForestWholeNumber(
     minHoldOutPerColumn,
     name = "minHoldOutPerColumns"
   )

   observedTargets <- prepared$imputableObserved %>%
     dplyr::arrange(.data$covariateId, .data$rowId)

   if (nrow(observedTargets) == 0) {
     return(tibble::tibble(
       rowId = integer(),
       covariateId = numeric(),
       trueValue = numeric()
     ))
   }

   splitObservedTargets <- split(observedTargets, observedTargets$covariateId)
   originalSeed <- if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
     get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
   } else {
     NULL
   }

   on.exit({
     if(!is.null(originalSeed)) {
       .Random.seed <<- originalSeed
     }
   }, add = TRUE)
   set.seed(seed)

   holdoutList <- lapply(splitObservedTargets, function(targetRows) {
     nObserved <- nrow(targetRows)
     if (nObserved <= 1) {
       return(NULL)
     }

     nHoldout <- max(minHoldOutPerColumn, floor(nObserved * validationFraction))
     nHoldout <- min(nHoldout, nObserved - 1L)

     if (nHoldout <= 0) {
       return(NULL)
     }

     sampledIndex <- sample.int(nObserved, size = nHoldout, replace = FALSE)
     targetRows[sampledIndex, , drop = FALSE] %>%
       dplyr::transmute(
         rowId = .data$rowId,
         covariateId = .data$covariateId,
         trueValue = .data$covariateValue
       )
   })

   dplyr::bind_rows(holdoutList) %>%
     dplyr::arrange(.data$covariateId, .data$rowId)
 }

maskXgboostTuningCells <- function(trainData, holdoutCells) {
  if (nrow(holdoutCells) == 0) {
    return(safeCopyPlpData(trainData))
  }

  covariates <- collectIfNeeded(trainData$covariateData$covariates) %>%
    dplyr::anti_join(
      holdoutCells %>% dplyr::select(.data$rowId, .data$covariateId),
      by = c("rowId", "covariateId")
    )

  plpDataHelper(
    labels = trainData$labels,
    folds = trainData$folds,
    covariates = covariates,
    covariateRef = collectIfNeeded(trainData$covariateData$covariateRef),
    analysisRef = collectIfNeeded(trainData$covariateData$analysisRef),
    templatePLPData = trainData
  )
}

materializeXgboostTuningData <- function(trainData) {
  out <- list(
    labels = as.data.frame(trainData$labels),
    covariateData = list(
      covariates = collectIfNeeded(trainData$covariateData$covariates),
      covariateRef = collectIfNeeded(trainData$covariateData$covariateRef),
      analysisRef = collectIfNeeded(trainData$covariateData$analysisRef)
    )
  )

  if (!is.null(trainData$folds)) {
    out$folds <- as.data.frame(trainData$folds)
  }

  class(out) <- "plpData"
  attr(out, "metaData") <- attr(trainData, "metaData")
  out
}

extractXgboostImputedCells <- function(imputedData, holdoutCells) {
  if (nrow(holdoutCells) == 0) {
    return(holdoutCells %>% dplyr::mutate(predictedValue = numeric()))
  }

  imputedCovariates <- collectIfNeeded(imputedData$covariateData$covariates) %>%
    dplyr::select(.data$rowId, .data$covariateId, .data$covariateValue) %>%
    dplyr::distinct(.data$rowId, .data$covariateId, .keep_all = TRUE)

  holdoutCells %>%
    dplyr::left_join(
      imputedCovariates,
      by = c("rowId", "covariateId")
    ) %>%
    dplyr::rename(predictedValue = .data$covariateValue)
}

extractXgboostImputedStateCells <- function(imputedState, holdoutCells) {
  if (nrow(holdoutCells) == 0) {
    return(holdoutCells %>% dplyr::mutate(predictedValue = numeric()))
  }

  imputedCells <- imputedState$currentImputedValues %>%
    dplyr::select(.data$rowId, .data$covariateId, .data$covariateValue) %>%
    dplyr::distinct(.data$rowId, .data$covariateId, .keep_all = TRUE)

  holdoutCells %>% dplyr::left_join(
    imputedCells,
    by = c("rowId", "covariateId")
  ) %>% dplyr::rename(predictedValue = .data$covariateValue)
}

scoreXgboostTuningPredictions <- function(predictions, metric = "rmse") {
  metric <- validateXgboostMetric(metric)

  missingPredictionCount <- sum(is.na(predictions$predictedValue))
  completePredictions <- predictions %>%
    dplyr::filter(!is.na(.data$predictedValue))

  if (nrow(completePredictions) == 0) {
    return(list(
      score = Inf,
      missingPredictionCount = missingPredictionCount,
      observedpredictionCount = 0L
    ))
  }

  errors <- completePredictions$predictedValue - completePredictions$trueValue
  score <- if (metric == "rmse") {
    sqrt(mean(errors^2))
  } else {
    mean(abs(errors))
  }

  if (missingPredictionCount > 0) {
    score <- Inf
  }

  list(
    score = as.numeric(score),
    missingPredictionCount = missingPredictionCount,
    observedPredictionCount = nrow(completePredictions)
  )
}

normalizeXgboostTuningCandidates <- function(candidates = NULL, grid = NULL, defaultArgs = list()) {
  if (!is.null(candidates) && !is.null(grid)) {
    stop("Provide either candidates or grid, not both.")
  }

  if (is.null(candidates) && is.null(grid)) {
    stop("Provide candidate settings via candidates or grid.")
  }

  if (!is.null(candidates)) {
    if (inherits(candidates, "featureEngineeringSettings")) {
      candidates <- list(candidate_1 = candidates)
    }

    if (!is.list(candidates) || length(candidates) == 0) {
      stop("candidates must be a non-empty list of featureEngineeringSettings objects.")
    }

    if (is.null(names(candidates)) || any(names(candidates) == "")) {
      names(candidates) <- paste0("candidate_", seq_along(candidates))
    }

    invalidCandidates <- !vapply(
      candidates,
      function(x) inherits(x, "featureEngineeringSettings") &&
        identical(attr(x, "fun"), "implementXgboostImputer"),
      logical(1)
    )

    if (any(invalidCandidates)) {
      stop("All candidates must be xgboost featureEngineeringSettings objects.")
    }

    return(list(
      candidates =candidates,
      grid = NULL
    ))
  }

  grid <- as.data.frame(grid, stringsAsFactors = FALSE)
  if (nrow(grid)==0) {
    stop("grid must contain at least one row.")
  }

  candidateNames <- paste0("candidate_", seq_len(nrow(grid)))
  candidates <- lapply(seq_len(nrow(grid)), function(i) {
    candidateArgs <- utils::modifyList(defaultArgs, as.list(grid[i, , drop=FALSE]))
    suppressWarnings(do.call(createXgboostImputer, candidateArgs))
  })
  names(candidates) <- candidateNames

  list(
    candidates = candidates,
    grid = tibble::as_tibble(grid) %>%
      dplyr::mutate(candidate = candidateNames, .before = 1)
  )
}

evaluateXgboostTuningCandidate <- function(candidateName, candidateSettings, maskedTrainData, holdoutCells, metric = "rmse") {
  prepared <- prepareMissForestPlpData(maskedTrainData)

  fittedImputer <- fitXgboostImputer(
    prepared = prepared,
    settings = candidateSettings
  )

  imputedState <- applyXgboostImputer(
    prepared = prepared,
    fittedImputer = fittedImputer,
    settings = candidateSettings
  )

  predictions <- extractXgboostImputedStateCells(
    imputedState = imputedState,
    holdoutCells = holdoutCells
  )
  scoreInfo <- scoreXgboostTuningPredictions(
    predictions = predictions,
    metric = metric
  )

  tibble::tibble(
    candidate = candidateName,
    score = scoreInfo$score,
    metric = metric,
    holdoutCells = nrow(holdoutCells),
    predictedCells = scoreInfo$observedPredictionCount,
    missingPredictions = scoreInfo$missingPredictionCount
  )
}

runXgboostTuningCandidates <- function(candidateNames, candidates, maskedTrainData, holdoutCells, metric = "rmse", numWorkers = 1L, sourceDir = NULL, verbose = FALSE) {
  numWorkers <- validateMissForestWholeNumber(
    numWorkers,
    name = "numWorkers"
  )

  if (numWorkers == 1L || length(candidateNames) <= 1L) {
    return(lapply(seq_along(candidateNames), function(i) {
      candidateName <- candidateName[[i]]
      if (isTRUE(verbose)) {
        message(
          "Xgboost tuning: evaluate candidate ",
          i,
          "/",
          length(candidateName),
          " (",
          candidateName,
          ")."
        )
      }

      evaluateXgboostTuningCandidate(
        candidateName = candidateName,
        candidateSettings = candidate[[candidateName]],
        maskedTrainData = maskedTrainData,
        holdoutCells = holdoutCells,
        metric = metric
      )
    }))
  }

  projectRoot <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  while (!file.exists(file.path(projectRoot, "DESCRIPTION")) &&
         dirname(projectRoot) != projectRoot) {
    projectRoot <- dirname(projectRoot)
  }
  if (!file.exists(file.path(projectRoot, "DESCRIPTION"))) {
    stop(
      "Parallel tuning could not locate the package root containing DESCRIPTION from the current working directory."
    )
  }

  workerFiles <- list.files(
    file.path(projectRoot, "R"),
    pattern = "[.]R$",
    full.names = TRUE
  )

  missingFiles <- workerFiles[!file.exists(workerFiles)]
  if(length(missingFiles) > 0) {
    stop(
      "Parallel tuning could not find required source files from the current working directory: ",
      paste(missingFiles, collapse = ", ")
    )
  }

  if (isTRUE(verbose) && any(vapply(candidates, function(x) !is.null(x$numThreads) && x$numThreads > 1L, logical(1)))) {
    message(
      "Xgboost tuning: parallel candidate evaluation works best when each candidate uses numThreads = 1 to avoid oversubscription."
    )
  }

  clusterSize <- min(numWorkers, length(candidateNames))
  cl <- parallel::makeCluster(clusterSize)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterEvalQ(cl, {
    library(magrittr)
    invisible(loadNamespace("Andromeda"))
    invisible(loadNamespace("xgboost"))
    invisible(loadNamespace("Matrix"))
    NULL
  })
  parallel::clusterExport(
    cl,
    varlist = c("maskedTrainData", "holdoutCells", "candidates", "metric"),
    envir = environment()
  )

  parallel::clusterCall(cl, function(paths) {
    invisible(lapply(paths, source))
    NULL
  }, paths = workerFiles)

  if (isTRUE(verbose)) {
    message(
      "Xgboost tuning: evaluating ",
      length(candidateNames),
      " candidates in parallel across ",
      clusterSize,
      " worker(s)."
    )
  }

  parallel::parLapply(cl, candidateNames, function(candidateName) {
    evaluateXgboostTuningCandidate(
      candidateName = candidateName,
      candidateSettings = candidates[[candidateName]],
      maskedTrainData = maskedTrainData,
      holdoutCells = holdoutCells,
      metric = metric
    )
  })
}

#' @title Tune iterative xgboost imputer settings
#' @description
#' Evaluate one or more iterative xgboost imputation candidates on held-out
#' observed cells and return tuning diagnostics together with the selected
#' candidate.
#' @export
tuneXgboostImputer <- function(trainData, candidates = NULL, grid = NULL, validationFraction = 0.20, minHoldoutPerColumn = 1L, metric = "rmse", numWorkers = 1L, sourceDir = NULL,
                               missingThreshold = 0.95, nrounds = 50L, maxDepth = 6L, eta = 0.10, subsample = 0.80, colsampleBytree = 0.80, minChildWeight = 1,
                               numThreads = NULL, includeSparseZeroPredictors = TRUE, addMissingIndicator = FALSE, seed = 42L, verbose = FALSE, returnMaskedData, device = NULL, treeMethod = "hist", gamma = 0, samplingMethod = "uniform",
                               colsampleBylevel = 1, colsampleBynode = 1, maxDeltaStep = 0, lambda = 1, alpha = 0, maxLeaves = 0L,
                               maxBin = 256L, numParallelTree = 1L) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required for xgboost tuning.")
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for xgboost tuning.")
  }

  metric <- validateXgboostMetric(metric)
  returnMaskedData <- validateMissForestLogical(returnMaskedData, name = "returnMaskedData")
  verbose <- validateMissForestLogical(verbose, name = "verbose")

  normalized <- normalizeXgboostTuningCandidates(
    candidates = candidates,
    grid = grid,
    defaultArgs = list(
      missingThreshold = missingThreshold,
      nrounds = nrounds,
      device = device,
      treeMethod = treeMethod,
      maxDepth = maxDepth,
      eta =eta,
      gamma = gamma,
      subsample = subsample,
      samplingMethod = samplingMethod,
      colsampleBytree = colsampleBytree,
      colsampleBylevel = colsampleBylevel,
      colsampleBynode = colsampleBynode,
      minChildWeight = minChildWeight,
      maxDeltaStep = maxDeltaStep,
      lambda = lambda,
      alpha = alpha,
      maxLeaves = maxLeaves,
      maxBin = maxBin,
      numParallelTree = numParallelTree,
      numThreads = numThreads,
      includeSparseZeroPredictors = includeSparseZeroPredictors,
      addMissingIndicator = addMissingIndicator,
      seed = seed,
      verbose = FALSE
    )
  )

  prepared <- prepareMissForestPlpData(trainData)
  holdoutCells <- builddXgboostTuningHoldout(
    prepared = prepared,
    validationFraction = validationFraction,
    minHoldOutPerColumn = minHoldoutPerColumn,
    seed = seed
  )

  if (nrow(holdoutCells) == 0) {
    stop("No holdout cells could be created for tuning. Check target covariates and validation settings.")
  }

  maskedTrainData <- maskXgboostTuningCells(trainData, holdoutCells)
  maskedTrainData <- materializeXgboostTuningData(maskedTrainData)
  candidateNames <- names(normalized$candidates)
  # tuningRows <- vector("list", length(candidateNames))
  #
  # for(i in seq_along(candidateNames)) {
  #   candidateName <- candidateNames[[i]]
  #   candidateSettings <- normalized$candidates[[candidateName]]
  #
  #   if (isTRUE(verbose)) {
  #     message(
  #       "Sparse xgboost tuning: evaluating candidate ",
  #       i,
  #       "/",
  #       length(candidateNames),
  #       " (",
  #       candidateName,
  #       ")."
  #     )
  #   }
  #
  #   imputedTrainData <- implementXgboostImputer(
  #     trainData = maskedTrainData,
  #     featureEngineeringSettings = candidateSettings,
  #     done = FALSE
  #   )
  #
  #   predictions <- extractXgboostImputedCells(
  #     imputedData = imputedTrainData,
  #     holdoutCells = holdoutCells
  #   )
  #
  #   scoreInfo <- scoreXgboostTuningPredictions(
  #     predictions = predictions,
  #     metric = metric
  #   )
  #
  #   tuningRows[[i]] <- tibble::tibble(
  #     candidate = candidateName,
  #     score = scoreInfo$score,
  #     metric = metric,
  #     holdoutCells = nrow(holdoutCells),
  #     predictedCells = scoreInfo$observedPredictionCount,
  #     missingPredictions = scoreInfo$missingPredictionCount
  #   )
  # }

  tuningRows <- runXgboostTuningCandidates(
    candidateNames = candidateNames,
    candidates = normalized$candidates,
    maskedTrainData = maskedTrainData,
    holdoutCells = holdoutCells,
    metric = metric,
    numWorkers =numWorkers,
    sourceDir = sourceDir,
    verbose = verbose
  )

  tuningResults <- dplyr::bind_rows(tuningRows) %>%
    dplyr::arrange(.data$score, .data$candidate)

  if (!is.null(normalized$grid)) {
    tuningResults <- tuningResults %>%
      dplyr::left_join(normalized$grid, by = "candidate")
  }

  bestCandidate <- tuningResults$candidate[[1]]
  bestSettings <- normalized$candidates[[bestCandidate]]

  output <- list(
    bestCandidate = bestCandidate,
    bestScore = tuningResults$score[[1]],
    bestSettings = bestSettings,
    bestImputationOverview = stats::setNames(list(bestSettings), bestCandidate),
    tuningResults = tuningResults,
    holdoutCells = holdoutCells
  )

  if (isTRUE(returnMaskedData)) {
    output$maskedTrainData <- maskedTrainData
  }

  class(output) <- "xgboostImputerTuning"
  output
}
