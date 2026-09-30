# missForest Imputer
# Custom-made missForest-style imputation method outside PatientLevelPrediction package


buildMissForestCompleteGrid <- function(rowIds, covariateIds, observedCovariates) {
  if (length(rowIds) == 0 || length(covariateIds) == 0) {
    return(tibble::tibble(
      rowId = integer(),
      covariateId = integer(),
      covariateValue = numeric()
    ))
  }

  completeIds <- expand.grid(
    rowId = rowIds,
    covariateId = covariateIds,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  completeGrid <- tibble::as_tibble(completeIds) %>%
    dplyr::mutate(
      rowIdKey = normalizeMissForestIdKey(.data$rowId),
      covariateIdKey = normalizeMissForestIdKey(.data$covariateId)
    )

  rowIdKeys <- unique(completeGrid$rowIdKey)
  covariateIdKeys <- unique(completeGrid$covariateIdKey)

  observedSubset <- observedCovariates %>%
    dplyr::mutate(
      rowIdKey = normalizeMissForestIdKey(.data$rowId),
      covariateIdKey = normalizeMissForestIdKey(.data$covariateId)
    ) %>%
    dplyr::filter(.data$rowIdKey %in% rowIdKeys, .data$covariateIdKey %in% covariateIdKeys) %>%
    dplyr::distinct(.data$rowIdKey, .data$covariateIdKey, .keep_all = TRUE) %>%
    dplyr::select("rowIdKey", "covariateIdKey", "covariateValue")

  completeGrid %>%
    dplyr::left_join(observedSubset, by = c("rowIdKey", "covariateIdKey")) %>%
    dplyr::select("rowId", "covariateId", "covariateValue") %>%
    dplyr::arrange(.data$covariateId, .data$rowId)
}
#' @title Create missForest imputer settings
#' @description
#' Create a `featureEngineeringSettings` object for the package's experimental
#' missForest-style imputer.
#' @export
createMissForestImputer <- function(missingThreshold = 0.3,
                                    maxiter = 5L,
                                    tol = 1e-3,
                                    ntree = 100L,
                                    mtry = NULL,
                                    minNodeSizeNumeric = 5L,
                                    minNodeSizeBinary = 1L,
                                    sampleFraction = NULL,
                                    replace = TRUE,
                                    numThreads = NULL,
                                    decreasing = FALSE,
                                    addMissingIndicator = FALSE,
                                    seed = 42L,
                                    verbose = FALSE) {
  warning("MissForest imputation is experimental in this package.", call. = FALSE)

  featureEngineeringSettings <- list(
    missingThreshold = validateNumeric(
      missingThreshold,
      name = "missingThreshold",
      lower = 0,
      upper = 1,
      lowerInclusive = FALSE,
      upperInclusive = TRUE
    ),
    maxiter = validateInteger(maxiter, name = "maxiter"),
    tol = validateNumeric(tol, name = "tol", lower = 0),
    ntree = validateInteger(ntree, name = "ntree"),
    mtry = validateInteger(mtry, name = "mtry", allowNull = TRUE),
    minNodeSizeNumeric = validateInteger(minNodeSizeNumeric, name = "minNodeSizeNumeric"),
    minNodeSizeBinary = validateInteger(minNodeSizeBinary, name = "minNodeSizeBinary"),
    sampleFraction = validateNumeric(
      sampleFraction,
      name = "sampleFraction",
      lower = 0,
      upper = 1,
      lowerInclusive = FALSE,
      upperInclusive = TRUE,
      allowNull = TRUE
    ),
    replace = validateLogical(replace, name = "replace"),
    numThreads = validateInteger(numThreads, name = "numThreads", allowNull = TRUE),
    decreasing = validateLogical(decreasing, name = "decreasing"),
    addMissingIndicator = validateLogical(addMissingIndicator, name = "addMissingIndicator"),
    seed = validateInteger(seed, name = "seed", minValue = 0L),
    verbose = validateLogical(verbose, name = "verbose")
  )

  attr(featureEngineeringSettings, "fun") <- "implementMissForestImputer"
  class(featureEngineeringSettings) <- "featureEngineeringSettings"
  featureEngineeringSettings
}

implementMissForestImputer <- function(trainData, featureEngineeringSettings, done = FALSE) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is required for MissForest imputation.")
  }

  if(isTRUE(done) && is.null(featureEngineeringSettings$model)) {
    stop("Fitted MissForest settings must include a trained model when done = TRUE.")
  }

  if (isTRUE(featureEngineeringSettings$verbose)) {
    if (isTRUE(done)) {
      message("MissForest: applying fitted imputer to new data.")
    } else {
      message("MissForest: fitting imputer on training data.")
    }
  }

  prepared <- prepareTrainPlpData(trainData)

  if (isTRUE(featureEngineeringSettings$verbose)) {
    #totalMissingCells <- sum(is.na(prepared$imputableData))
    message(
      "MissForest: prepared ",
      #nrow(prepared$imputableData),
      length(prepared$rowIds),
      " training rows ",
      #ncol(prepared$imputableData),
      length(prepared$imputableCovariateIds),
      " imputable covariates, and",
      length(prepared$predictorCovariateIds),
      " predictor covariates. Target missing cells = ",
      #totalMissingCells,
      prepared$totalMissingCells,
      "."
    )

    if (length(prepared$imputableCovariateIds) == 0) {
      message(
        "MissForest: no imputable target covariates were identified after applying isBinary/missingMeansZero filters."
      )
    }
  }
  fittedImputer <- if (isTRUE(done)) {
    featureEngineeringSettings$model
  } else {
    fitMissForestImputer(prepared = prepared, settings = featureEngineeringSettings)
  }

  if (!isTRUE(done)) {
    fittedImputer <- addMetadata(
      prepared = prepared,
      fittedImputer = fittedImputer
    )
  }

  imputedState <- applyMissForestImputer(
    prepared = prepared,
    fittedImputer = fittedImputer,
    settings = featureEngineeringSettings
  )

  rebuilt <- rebuildTrainPlpData(
    trainData = trainData,
    prepared = prepared,
    #imputedDataFrame = imputedDataFrame,
    imputedState = imputedState,
    fittedImputer = fittedImputer
  )

  if (!isTRUE(done)) {
    fittedSettings <- featureEngineeringSettings
    fittedSettings$model <- fittedImputer

    metaData <- attr(rebuilt$covariateData, "metaData")
    if(is.null(metaData)) {
      metaData <- list()
    }

    if(is.null(metaData$featureEngineering)) {
      metaData$featureEngineering <- list()
    }

    metaData$featureEngineering[["missForestImputer"]] <- list(
      funct = "implementMissForestImputer",
      settings = list(featureEngineeringSettings = fittedSettings)
    )
    attr(rebuilt$covariateData, "metaData") <- metaData
  }
  rebuilt
}

prepareMissForestTrainingState <- function(prepared, settings) {
  rowIds <- sort(unique(prepared$rowIds))
  imputableIds <- sort(unique(prepared$imputableCovariateIds))

  missingInfo <- prepared$missingInfo
  keepCols <- missingInfo %>%
    dplyr::filter(.data$missing <= settings$missingThreshold, .data$missing < 1) %>%
    dplyr::pull(.data$covariateId) %>%
    unique() %>%
    sort()

  removedCols <- setdiff(imputableIds, keepCols)
  predictorCovariateIds <- setdiff(prepared$predictorCovariateIds, removedCols)

  observedCovariates <- prepared$covariates %>%
    dplyr::filter(.data$covariateId %in% predictorCovariateIds, !is.na(.data$covariateValue)) %>%
    dplyr::arrange(.data$covariateId, .data$rowId)

  targetMatrix <- subsetPreparedTargetMatrix(
    prepared = prepared,
    covariateIds = keepCols
  )
  missingMask <- is.na(targetMatrix)

  predictorDefaults <- prepared$covariateInfo %>%
    dplyr::filter(.data$covariateId %in% predictorCovariateIds) %>%
    dplyr::distinct(.data$covariateId, .keep_all = TRUE) %>%
    dplyr::transmute(
      covariateId = .data$covariateId,
      defaultValue = dplyr::if_else(.data$missingMeansZero == "Y",0,NA_real_)
    )

  predictorDefaults <- stats::setNames(
    predictorDefaults$defaultValue,
    predictorDefaults$covariateId
  )
  missingIndexList <- vector("list", length(keepCols))
  names(missingIndexList) <- as.character(keepCols)
  binaryCols <- logical(length(keepCols))
  names(binaryCols) <- as.character(keepCols)
  initialFill <- vector("list", length(keepCols))
  names(initialFill) <- as.character(keepCols)
  missingCounts <- integer(length(keepCols))
  names(missingCounts) <- as.character(keepCols)

  for (i in seq_along(keepCols)) {
    covariateId <- keepCols[i]
    missingRows <- rowIds[missingMask[, i]]

    missingCounts[i] <- length(missingRows)
    if (length(missingRows) > 0) {
      missingIndexList[[i]] <- data.frame(
        rowId = missingRows,
        covariateId = covariateId
      )
    } else {
      missingIndexList[[i]] <- data.frame(
        rowId = integer(),
        covariateId = integer()
      )
    }

    covariateValues <- targetMatrix[, i]
    covariateValues <- covariateValues[!is.na(covariateValues)]
    binaryCols[i] <- isBinaryCovariate(covariateValues)
    if (binaryCols[i]) {
      initialFill[[i]] <- computeBinaryMode(covariateValues)
    } else {
      fillValue <- mean(covariateValues, na.rm=TRUE)
      if (!is.finite(fillValue)) {
        fillValue <- 0
      }
      initialFill[[i]] <- fillValue
    }
  }

  list(
    rowIds = rowIds,
    originalColCount = length(imputableIds),
    keepCols = keepCols,
    removedCols = removedCols,
    predictorCovariateIds = predictorCovariateIds,
    predictorDefaults = predictorDefaults,
    observedCovariates = observedCovariates,
    missingIndex = if (length(missingIndexList)==0) {
      tibble::tibble(rowId = integer(), covariateId = integer())
    } else {
      dplyr::bind_rows(missingIndexList) %>%
        dplyr::arrange(.data$covariateId, .data$rowId)
    },
    missingCounts = missingCounts,
    binaryCols = binaryCols,
    initialFill = initialFill
  )
}

materializeMissForestPredictors <- function(rowIds, predictorIds, predictorDefaults, observedCovariates, imputedCovariates) {
  if (length(predictorIds) == 0) {
    return(data.frame(row.names = seq_along(rowIds)))
  }

  observedSubset <- observedCovariates %>%
    dplyr::filter(.data$rowId %in% rowIds, .data$covariateId %in% predictorIds)

  imputedSubset <- imputedCovariates %>%
    dplyr::filter(.data$rowId %in% rowIds, .data$covariateId %in% predictorIds)

  values <- dplyr::bind_rows(observedSubset, imputedSubset)
  defaultValues <- predictorDefaults[as.character(predictorIds)]
  defaultValues[is.na(defaultValues)] <- NA_real_
  predictorMatrix <- matrix(
    #NA_real_,
    rep(defaultValues, each = length(rowIds)),
    nrow = length(rowIds),
    ncol = length(predictorIds),
    dimnames = list(as.character(rowIds), as.character(predictorIds))
  )

  if (nrow(values) > 0) {
    rowIndex <- match(values$rowId, rowIds)
    colIndex <- match(values$covariateId, predictorIds)
    predictorMatrix[cbind(rowIndex, colIndex)] <- values$covariateValue
  }

  if (anyNA(predictorMatrix)) {
    stop("MissForst predictor matrix could not be fully materialized for the requested rows.")
  }

  as.data.frame(predictorMatrix, check.names = FALSE, stringsAsFactors = FALSE)
}

fitSingleImputationModelLong <- function(state, targetCovariateId, settings, fitModelOnly = FALSE) {
  targetObserved <- state$observedCovariates %>%
    dplyr::filter(.data$covariateId == .env$targetCovariateId) %>%
    dplyr::arrange(.data$rowId)

  missingRows <- state$missingIndex %>%
    dplyr::filter(.data$covariateId == .env$targetCovariateId) %>%
    dplyr::pull(.data$rowId)

  observedRows <- targetObserved$rowId
  if (length(missingRows) == 0 || length(observedRows) == 0) {
    return(list(model = NULL, predictions = NULL))
  }

  predictorIds <- setdiff(state$predictorCovariateIds, targetCovariateId)
  if(length(predictorIds) == 0L) {
    return(list(model = NULL, predictions = NULL))
  }

  trainPredictors <- materializeMissForestPredictors(
    rowIds = observedRows,
    predictorIds = predictorIds,
    predictorDefaults = state$predictorDefaults,
    observedCovariates = state$observedCovariates,
    imputedCovariates = state$currentImputedValues
  )

  response <- targetObserved$covariateValue[match(observedRows, targetObserved$rowId)]
  isBinaryTarget <- isTRUE(state$binaryCols[[as.character(targetCovariateId)]])
  targetName <- as.character(targetCovariateId)

  if (isBinaryTarget) {
    response <- factor(as.character(as.integer(response)), levels = c("0","1"))
    if(nlevels(droplevels(response)) <= 1) {
      constantValue <- as.character(levels(droplevels(response)))[1]
      predictions <- if (fitModelOnly) NULL else rep(as.numeric(constantValue), length(missingRows))
      return(list(model = list(kind = "constant_binary", value = constantValue),
             predictions = predictions))
    }
  } else if (length(unique(response)) <= 1) {
    constantValue <- as.numeric(response[1])
    predictions <- if(fitModelOnly) NULL else rep(constantValue, length(missingRows))
    return(list(
      model = list(kind = "constant_numeric", value = constantValue),
      predictions = predictions
    ))
  }

  trainingDataFrame <- trainPredictors
  trainingDataFrame[[targetName]] <- response

  fittedModel <- do.call(
    ranger::ranger,
    buildRangerArguments(
      trainingDataFrame = trainingDataFrame,
      targetName = targetName,
      settings = settings,
      isBinaryTarget = isBinaryTarget
    )
  )

  predictions <- NULL
  if (!fitModelOnly) {
    predictionDataFrame <- materializeMissForestPredictors(
      rowIds = missingRows,
      predictorIds = predictorIds,
      predictorDefaults = state$predictorDefaults,
      observedCovariates = state$observedCovariates,
      imputedCovariates = state$currentImputedValues
    )
    predictions <- predict(fittedModel, data= predictionDataFrame)$predictions
    if (isBinaryTarget) {
      predictions <- as.numeric(as.character(predictions))
    } else {
      predictions <- as.numeric(predictions)
    }
  }

  list(
    model = list(
      kind = if (isBinaryTarget) "ranger_binary" else "range_numeric",
      fit = fittedModel
    ),
    predictions = predictions
  )
}

buildRangerArguments <- function(trainingDataFrame, targetName, settings, isBinaryTarget) {
  predictorCount <- ncol(trainingDataFrame) - 1L
  localMtry <- if (is.null(settings$mtry)) {
    max(1L, floor(sqrt(max(1L, predictorCount))))
  } else {
    min(settings$mtry, max(1L, predictorCount))
  }

  rangerArguments <- list(
    dependent.variable.name = targetName,
    data = trainingDataFrame,
    num.trees = settings$ntree,
    mtry = localMtry,
    replace = settings$replace,
    write.forest = TRUE,
    seed = settings$seed,
    verbose = FALSE
  )

  if (!is.null(settings$sampleFraction)) {
    rangerArguments$sample.fraction <- settings$sampleFraction
  }

  if (!is.null(settings$numThreads)) {
    rangerArguments$num.threads <- settings$numThreads
  }

  if (isBinaryTarget) {
    rangerArguments$min.node.size <- settings$minNodeSizeBinary
    rangerArguments$probability <- FALSE
    rangerArguments$respect.unordered.factors <- "order"
  } else {
    rangerArguments$min.node.size <- settings$minNodeSizeNumeric
  }

  rangerArguments
}

predictSingleImputationModel <- function(modelObject, newDataFrame) {
  if (is.null(modelObject)) {
    return(NULL)
  }

  if (identical(modelObject$kind, "constant_numeric")) {
    return(rep(modelObject$value, nrow(newDataFrame)))
  }

  if (identical(modelObject$kind, "constant_binary")) {
    #return(factor(rep(modelObject$value, nrow(newDataFrame)), levels = c("0","1")))
    return(rep(as.numeric(modelObject$value), nrow(newDataFrame)))
  }

  predictions <- predict(modelObject$fit, data = newDataFrame)$predictions

  if (identical(modelObject$kind, "ranger_binary")) {
    #return(factor(as.character(predictions), levels = c("0", "1")))
    return(as.numeric(as.character(predictions)))
  }

  as.numeric(predictions)
}

fitMissForestImputer <- function(prepared, settings) {
  preparation <- prepareMissForestTrainingState(prepared, settings)

  if (isTRUE(settings$verbose)) {
    message(
      "MissForest: retained ",
      length(preparation$keepCols),
      " of ",
      preparation$originalColCount,
      " imputable columns after applying missingThreshold = ",
      settings$missingThreshold,
      "."
    )

    if (length(preparation$removedCols) > 0) {
      message(
        "MissForest: removed ",
        length(preparation$removedCols),
        " columns because they were fully missing or exceeded the missingness threshold."
      )
    }
  }

  if (length(preparation$keepCols) == 0) {
    if (!isTRUE(settings$verbose)) {
      message("MissForest: no columns remained for iterative fitting.")
    }
    return(list(
      keepCols = character(),
      #removedCols = names(x),
      removedCols = preparation$removedCols,
      binaryCols = logical(),
      initialFill = list(),
      visitOrder = integer(),
      models = list(),
      addMissingIndicator = isTRUE(settings$addMissingIndicator),
      indicatorCovariateIds = integer(),
      indicatorAnalysisIds = integer(),
      indicatorColumnNames = character()
    ))
  }

  currentImputedValues <- initializeLongImputations(
    missingIndex = preparation$missingIndex,
    initialFill = preparation$initialFill
  )

  visitOrder <- preparation$keepCols[order(preparation$missingCounts)]
  if (settings$decreasing) {
    visitOrder <- rev(visitOrder)
  }

  visitOrder <- visitOrder[preparation$missingCounts[as.character(visitOrder)] > 0]

  if (isTRUE(settings$verbose)) {
    visitRatios <- preparation$missingCounts[as.character(visitOrder)] / length(preparation$rowIds)
    message(
      "MissForest: ",
      length(visitOrder),
      " retained columns contain at least one missing value and will enter the iterative loop. Missing ratio(s) = ",
      if (length(visitRatios) == 0) {
        "none"
      } else {
        paste(formatMissForestRatio(visitRatios), collapse = ", ")
      },
      "."
    )
  }
  previousDelta <- Inf

  if (length(visitOrder) > 0) {
    for (iteration in seq_len(settings$maxiter)) {
      if (isTRUE(settings$verbose)) {
        message(sprintf("MissForest training iteration %s", iteration))
      }
      #previousImputedData <- xImputed

      previousImputedValues <- currentImputedValues
      state <- list(
        rowIds = preparation$rowIds,
        keepCols = preparation$keepCols,
        predictorCovariateIds = preparation$predictorCovariateIds,
        predictorDefaults = preparation$predictorDefaults,
        observedCovariates = preparation$observedCovariates,
        missingIndex = preparation$missingIndex,
        binaryCols = preparation$binaryCols,
        currentImputedValues = currentImputedValues
      )

      for (targetIndex in visitOrder) {
        targetCovariateId <- visitOrder[[targetIndex]]
        if (isTRUE(settings$verbose)) {
          message(
            "MissForest: training target ",
            targetIndex,
            "/",
            length(visitOrder),
            " (covariateId = ",
            targetCovariateId,
            ", missing rows = ",
            preparation$missingCounts[[as.character(targetCovariateId)]],
            ", missing ratio = ",
            formatMissForestRatio(
              preparation$missingCounts[[as.character(targetCovariateId)]] / length(preparation$rowIds)
            ),
            ")."
          )
        }
        fitResult <- fitSingleImputationModelLong(
          state = state,
          targetCovariateId = targetCovariateId,
          settings = settings,
          fitModelOnly = FALSE
        )

        if (!is.null(fitResult$predictions)) {
          currentImputedValues <- updateImputedValues(
            currentImputedValues = currentImputedValues,
            targetCovariateId = targetCovariateId,
            predictions = fitResult$predictions
          )
          state$currentImputedValues <- currentImputedValues
        }
      }

      currentDelta <- computeDeltaLong(
        previousImputedValues = previousImputedValues,
        currentImputedValues = currentImputedValues,
        missingIndex = preparation$missingIndex,
        binaryCols = preparation$binaryCols
      )

      if (isTRUE(settings$verbose)) {
        message(sprintf("MissForest delta = %.6f", currentDelta))
      }

      if (currentDelta <= settings$tol || currentDelta > previousDelta) {
        if (isTRUE(settings$verbose)) {
          if (currentDelta <= settings$tol) {
            message("MissForest: stopping because delta is below tolerance.")
          } else {
            message("MissForest: stopping because delta is increased.")
          }
        }
        break
      }
      previousDelta <- currentDelta
    }
  } else if (isTRUE(settings$verbose)) {
    message("MissForest: skipping iterative fitting because no retained column has missing values.")
  }

  fittedModels <- vector("list", length = length(preparation$keepCols))
  names(fittedModels) <- as.character(preparation$keepCols)
  finalState <- list(
    rowIds = preparation$rowIds,
    keepCols = preparation$keepCols,
    predictorCovariateIds = preparation$predictorCovariateIds,
    predictorDefaults = preparation$predictorDefaults,
    observedCovariates = preparation$observedCovariates,
    missingIndex = preparation$missingIndex,
    binaryCols = preparation$binaryCols,
    currentImputedValues = currentImputedValues
  )

  for(targetIndex in seq_along(preparation$keepCols)) {
    targetCovariateId <- preparation$keepCols[[targetIndex]]
    if (isTRUE(settings$verbose)) {
      message(
        "MissForest: fitting final model ",
        targetIndex,
        "/",
        length(preparation$keepCols),
        " for covariateId= ",
        targetCovariateId,
        "."
      )
    }

    fittedModels[[as.character(targetCovariateId)]] <- fitSingleImputationModelLong(
      state = finalState,
      targetCovariateId = targetCovariateId,
      settings = settings,
      fitModelOnly = TRUE
    )$model
  }

  list(
    keepCols = preparation$keepCols,
    removedCols = preparation$removedCols,
    predictorCovariateIds = preparation$predictorCovariateIds,
    predictorDefaults = preparation$predictorDefaults,
    binaryCols = preparation$binaryCols,
    initialFill = preparation$initialFill,
    visitOrder = visitOrder,
    models = fittedModels,
    addMissingIndicator = isTRUE(settings$addMissingIndicator),
    indicatorCovariateIds = integer(),
    indicatorAnalysisIds = integer(),
    indicatorColumnNames = character()
  )
}

applyMissForestImputer <- function(prepared, fittedImputer, settings) {
  if(length(fittedImputer$keepCols) == 0) {
    if (isTRUE(settings$verbose)) {
      message("MissForest: no fitted columns available for apply step.")
    }
    #return(data.frame(row.names = seq_len(nrow(x))))
    return(list(
      currentImputedValues = tibble::tibble(
        rowId = integer(),
        covariateId = integer(),
        covariateValue = numeric()
      )
    ))
  }

  observedCovariates <- prepared$covariates %>%
    dplyr::filter(.data$covariateId %in% fittedImputer$predictorCovariateIds, !is.na(.data$covariateValue)) %>%
    dplyr::arrange(.data$covariateId, .data$rowId)
  targetMatrix <- subsetPreparedTargetMatrix(
    prepared = prepared,
    covariateIds = fittedImputer$keepCols
  )
  missingMask <- is.na(targetMatrix)

  missingIndexList <- lapply(fittedImputer$keepCols, function(covariateId) {
    targetIndex <- match(covariateId, fittedImputer$keepCols)
    missingRows <- prepared$rowIds[missingMask[, targetIndex]]

    if (length(missingRows) == 0) {
      return(data.frame(rowId = integer(), covariateId = integer()))
    }

    data.frame(rowId = missingRows, covariateId = covariateId)
  })

  missingIndex <- dplyr::bind_rows(missingIndexList) %>%
    dplyr::arrange(.data$covariateId, .data$rowId)
  currentImputedValues <- initializeLongImputations(
    missingIndex = missingIndex,
    initialFill = fittedImputer$initialFill
  )

  if (isTRUE(settings$verbose)) {
    message(
      "MissForest: apply step received ",
      length(prepared$rowIds),
      " rows x ",
      length(fittedImputer$keepCols),
      " columns with ",
      nrow(missingIndex),
      " missing cells."
    )
  }

  if(length(fittedImputer$visitOrder) > 0) {
    state <- list(
      rowIds = prepared$rowIds,
      keepCols = fittedImputer$keepCols,
      predictorCovariateIds = fittedImputer$predictorCovariateIds,
      predictorDefaults = fittedImputer$predictorDefaults,
      observedCovariates = observedCovariates,
      missingIndex = missingIndex,
      binaryCols = fittedImputer$binaryCols,
      currentImputedValues = currentImputedValues
    )
    for (iteration in seq_len(settings$maxiter)) {
      if (isTRUE(settings$verbose)) {
        message(sprintf("MissForest apply iteration %s", iteration))
      }
      for (targetCovariateId in fittedImputer$visitOrder) {
        targetMissingRows <- missingIndex$covariateId == targetCovariateId
        if (!any(targetMissingRows)) {
          next
        }

        predictorIds <- setdiff(fittedImputer$predictorCovariateIds, targetCovariateId)
        if(length(predictorIds) == 0L) {
          next
        }

        predictionDataFrame <- materializeMissForestPredictors(
          rowIds = missingIndex$rowId[targetMissingRows],
          predictorIds = predictorIds,
          predictorDefaults = fittedImputer$predictorDefaults,
          observedCovariates = observedCovariates,
          imputedCovariates = currentImputedValues
        )

        predictions <- predictSingleImputationModel(
          modelObject = fittedImputer$models[[as.character(targetCovariateId)]],
          newDataFrame = predictionDataFrame
        )

        if (!is.null(predictions)) {
          currentImputedValues <- updateImputedValues(
            currentImputedValues = currentImputedValues,
            targetCovariateId = targetCovariateId,
            predictions = predictions
          )
          state$currentImputedValues <- currentImputedValues
        }
      }
    }
  } else if (isTRUE(settings$verbose)) {
    message("MissForest: apply step skipped iterative updates because the fitted visit order is empty")
  }

  list(currentImputedValues = currentImputedValues)

}






























