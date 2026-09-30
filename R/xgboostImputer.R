#' @title Validate an iterative xgboost choice argument
#' @description
#' This function checks that an iterative xgboost imputation argument is supplied
#' as a single character value and matches one of the allowed choices.
#' @param x The value to validate
#' @param name The argument name used in error messages
#' @param choices The allowed character values
#' @return The validated character value
validateXgboostChoice <- function(x, name, choices) {
  # Require exactly one character value
  if (length(x) != 1 || !is.character(x) || is.na(x)) {
    stop(name, " must be a single character value.")
  }

  # Restrict the argument to the supported iterative xgboost imputation options
  if (!x %in% choices) {
    stop(name, " must be one of: ", paste(choices, collapse = ", "), ".")
  }
  x
}

#' @title Validate a device setting
#' @description This helper checks that the iterative xgboost device argument is either
#' `NULL` or a single character value describing a supported compute device
#' @param device The device setting to validate
#' @return `NULL` or the validated device name string
validateXgboostDevice <- function(device) {
  # If no device is given, use the default
  if (is.null(device)) {
    return(NULL)
  }

  # Require exactly one non missing value
  if (length(device) != 1 || !is.character(device) || is.na(device)) {
    stop("device must be NULL or a single character value.")
  }

  device <- trimws(device)
}

#' @title Create iterative xgboost imputer settings
#' @description This function creates the settings for the iterative xgboost
#' imputation. It defines how target covariates are selected for imputation, how the
#' iterative procedure is run, and which hyperparameters are used to fit each conditional imputation model
#' @param missingThreshold The maximum allowed missingness  fraction for a target covariate to
#' remain eligible for imputation
#' @param nrounds The number of boosting rounds used for each xgboost model
#' @param maxDepth The maximum tree depth used by xgboost
#' @param eta The xgboost learning rate
#' @param subsample The fraction of training rows sampled for each boosting round
#' @param colsampleBytree The fraction of predictors sampled for each tree
#' @param minChildWeight The minimum child weight used by xgboost
#' @param maxiter The maximum number of imputation iterations
#' @param tol The stopping tolerance for the iterative procedure
#' @param numThreads The number of threads used by xgboost. If `NULL`, xgboost uses its default
#' @param decreasing Logical variable indicating whether variables should be
#' visited from highest to lowest missingness instead of the default increasing order
#' @param includeSparseZeroPredictors Logical variable indicating whether covariates with
#' `missingMeansZero = "Y"` should be included as sparse zero-default predictors
#' @param addMissingIndicator Logical variable indicating whether missingness indicators should
#' be added for imputed target covariates
#' @param seed The random seed passed to xgboost
#' @param verbose Logical variable indicating whether progress messages should be printed
#' @param device Optional xgboost device string such as `"cpu"` or `"cuda"`
#' @param treeMethod The xgboost tree construction method
#' @param gamma The xgboost minimum loss reduction parameter
#' @param samplingMethod The xgboost row-sampling method
#' @param colsampleBylevel The fraction of variables sampled at each tree level
#' @param colsampleByNode The fraction of variables sampled at each split
#' @param maxDeltaStep The xgboost maximum delta step
#' @param lambda The xgboost L2 regularization parameter
#' @param alpha The xgboost L1 regularization parameter
#' @param maxLeaves The maximum number of leaves per tree
#' @param maxBin The maximum number of histogram bins
#' @param numParallelTree The number of trees constructed in parallel per round
#' @return A `featureEngineeringSettings` object for `implementXgboostImputer()`
#' @export
createXgboostImputer <- function(missingThreshold = 0.95,
                                 nrounds = 20L,
                                 maxDepth = 6L,
                                 eta = 0.30,
                                 subsample = 0.70,
                                 colsampleBytree = 1,
                                 minChildWeight = 1,
                                 maxiter = 5L,
                                 tol = 1e-3,
                                 numThreads = 1L,
                                 decreasing = FALSE,
                                 includeSparseZeroPredictors = TRUE,
                                 addMissingIndicator = FALSE,
                                 seed = 42L,
                                 verbose = TRUE,
                                 device = "cpu",
                                 treeMethod = "hist",
                                 gamma = 0,
                                 samplingMethod = "uniform",
                                 colsampleBylevel = 1,
                                 colsampleBynode = 1,
                                 maxDeltaStep = 0,
                                 lambda = 1,
                                 alpha = 0,
                                 maxLeaves = 0L,
                                 maxBin = 256L,
                                 numParallelTree = 1L) {
  # Flag method as experimental
  warning("XGBoost imputation is experimental in this package", call. = FALSE)

  # Validate and store all imputer and xgboost hyperparameter settings
  featureEngineeringSettings <- list(
    missingThreshold = validateNumeric(
      missingThreshold,
      name = "missingThreshold",
      lower = 0,
      upper = 1,
      lowerInclusive = FALSE,
      upperInclusive = TRUE
    ),
    nrounds = validateInteger(nrounds, name = "nrounds"),
    device = validateXgboostDevice(
      device),
    treeMethod = validateXgboostChoice(
      treeMethod,
      name = "treeMethod",
      choices = c("auto","approx","exact","hist")
    ),
    maxDepth = validateInteger(maxDepth, name = "maxDepth"),
    eta = validateNumeric(
      eta,
      name = "eta",
      lower = 0,
      upper = 1,
      lowerInclusive = FALSE,
      upperInclusive = TRUE
    ),
    gamma = validateNumeric(
      gamma,
      name = "gamma",
      lower = 0,
      lowerInclusive = TRUE
    ),
    subsample = validateNumeric(
      subsample,
      name = "colsampleBytree",
      lower = 0,
      upper = 1,
      lowerInclusive = FALSE,
      upperInclusive = TRUE
    ),
    samplingMethod = validateXgboostChoice(
      samplingMethod,
      name = "samplingMethod",
      choices = c("uniform","gradient_based")
    ),
    colsampleBytree = validateNumeric(
      colsampleBytree,
      name = "colsampleBytree",
      lower = 0,
      upper = 1,
      lowerInclusive = FALSE,
      upperInclusive = TRUE
    ),
    colsampleBylevel = validateNumeric(
      colsampleBylevel,
      name = "colsampleBylevel",
      lower = 0,
      upper = 1,
      lowerInclusive = FALSE,
      upperInclusive = TRUE
    ),
    colsampleBynode = validateNumeric(
      colsampleBynode,
      name = "colsampleBynode",
      lower = 0,
      upper = 1,
      lowerInclusive = FALSE,
      upperInclusive = TRUE
    ),
    minChildWeight = validateNumeric(
      minChildWeight,
      name = "minChildWeight",
      lower = 0,
      lowerInclusive = TRUE
    ),
    maxDeltaStep = validateNumeric(
      maxDeltaStep,
      name = "maxDeltaStep",
      lower = 0,
      lowerInclusive = TRUE
    ),
    lambda = validateNumeric(
      lambda,
      name = "lambda",
      lower = 0,
      lowerInclusive = TRUE
    ),
    alpha = validateNumeric(
      alpha,
      name = "alpha",
      lower = 0,
      lowerInclusive = TRUE
    ),
    maxLeaves = validateInteger(
      maxLeaves,
      name= "maxLeaves",
      minValue = 0L
    ),
    maxBin = validateInteger(
      maxBin,
      name = "maxBin"
    ),
    numParallelTree = validateInteger(
      numParallelTree,
      name = "numParallelTree"
    ),
    maxiter = validateInteger(maxiter, name = "maxiter"),
    tol = validateNumeric(tol, name = "tol", lower = 0),
    numThreads = validateInteger(numThreads, name = "numThreads", allowNull = TRUE),
    decreasing = validateLogical(decreasing, name = "decreasing"),
    includeSparseZeroPredictors = validateLogical(includeSparseZeroPredictors, name = "includeSparseZeroPredictors"),
    addMissingIndicator = validateLogical(addMissingIndicator, name = "addMissingIndicator"),
    seed = validateInteger(seed, name = "seed", minValue = 0L),
    verbose = validateLogical(verbose, name = "verbose")
  )

  # Label the settings with the implementation function used by PLP
  attr(featureEngineeringSettings, "fun") <- "implementXgboostImputer"
  class(featureEngineeringSettings) <- "featureEngineeringSettings"
  featureEngineeringSettings
}

#' @title Build a sparse zero-default predictor matrix
#' @description This function contructs a sparse predictor matrix for covariates whose unobserved values
#' are interpreted as zer. Only nonzero observed values are materialized so the resulting matrix
#' stays sparse
#' @param rowIds The row ids to include
#' @param predictorIds The predictor covariate ids to include
#' @param observedCovariates A long-format table of observed covariate values
#' @param rowIdIndex Optional lookup index for `rowIds`
#' @param predictorIdIndex Optional lookup index for `predictorIds`
#' @return A sparse matrix with rows indexed by `rowIds` and columns by `predictorIds`
buildSparseZeroPredictorMatrix <- function(rowIds, predictorIds, observedCovariates, rowIdIndex = NULL, predictorIdIndex = NULL) {
  # Iterative xgboost relies on `Matrix` package to represent zero-default predictors efficiently
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for xgboost imputation.")
  }

  # Start from an all zero sparse matrix so unobserved values are treated as zero
  emptyMatrix <- Matrix::sparseMatrix(
    i = integer(),
    j = integer(),
    x = numeric(),
    dims = c(length(rowIds), length(predictorIds)),
    dimnames = list(as.character(rowIds), as.character(predictorIds))
  )

  if (length(predictorIds) == 0) {
    return(emptyMatrix)
  }

  # Keep only observed nonzero predictor values because zeros are already
  # represented implicitly by the sparse matrix structure
  predictorValues <- observedCovariates %>%
    dplyr::filter(
      .data$covariateId %in% predictorIds,
      !is.na(.data$covariateValue),
      .data$covariateValue != 0
    )

  if (nrow(predictorValues) == 0) {
    return(emptyMatrix)
  }

  # Build lookup indices once so long format row/covariate pairs can be mappedi nto matrix row and column positions
  if (is.null(rowIdIndex)) {
    rowIdIndex <- buildLookupIndex(rowIds)
  }

  if (is.null(predictorIdIndex)) {
    predictorIdIndex <- buildLookupIndex(predictorIds)
  }

  rowIndex <- unname(rowIdIndex[normalizeIdKey(predictorValues$rowId)])
  colIndex <- unname(predictorIdIndex[normalizeIdKey(predictorValues$covariateId)])
  valid <- !is.na(rowIndex) & !is.na(colIndex)

  if (!any(valid)) {
    return(emptyMatrix)
  }

  # Materialize only the valid observed nonzero entries into the sparse matrix
  Matrix::sparseMatrix(
    i = rowIndex[valid],
    j = colIndex[valid],
    x = predictorValues$covariateValue[valid],
    dims = c(length(rowIds), length(predictorIds)),
    dimnames = list(as.character(rowIds), as.character(predictorIds))
  )

}

#' @title Prepare the iterative xgboost imputation training state
#' @description This helper prepares the matrices, indices and storing objects
#' tneeded to fit the iterative xgboost needed to fit the iterative xgboost imputer on the training data
#' @param prepared A prepared PLP imputation object
#' @param settings The iterative xgboost imputer settings
#' @return A list containing the prepared training state for iterative xgboost imputation
prepareXgboostTrainingState <- function(prepared, settings) {
  rowIds <- sort(unique(prepared$rowIds))
  imputableIds <- sort(unique(prepared$imputableCovariateIds))
  missingInfo <- prepared$missingInfo

  # Retain only target covariates whose missingness ratio is low enough to be imputed
  keepCols <- missingInfo %>%
    dplyr::filter(.data$missing <= settings$missingThreshold, .data$missing < 1) %>%
    dplyr::pull(.data$covariateId) %>%
    unique() %>%
    sort()

  removedCols <- setdiff(imputableIds, keepCols)

  zeroDefaultPredictorIds <- numeric()

  if (isTRUE(settings$includeSparseZeroPredictors)) {
    # Use missingMeansZero covariates as sparse predictors, excluding the target covariates
    # that are being imputed
    zeroDefaultPredictorIds <- prepared$covariateInfo %>%
      dplyr::filter(
        .data$covariateId %in% prepared$predictorCovariateIds,
        .data$missingMeansZero == "Y",
        !(.data$covariateId %in% keepCols)
      ) %>%
      dplyr::pull(.data$covariateId) %>%
      unique() %>%
      sort()
  }

  # Build the dense target matrix for the covariates that will be imputed
  targetMatrix <-subsetPreparedTargetMatrix(
    prepared = prepared,
    covariateIds = keepCols
  )

  missingMask <- is.na(targetMatrix)
  filledTargetMatrix <- targetMatrix

  # Initialize per-target storage used iterative imputation
  keepColIndex <- buildLookupIndex(keepCols)
  binaryCols <- logical(length(keepCols))
  names(binaryCols) <- as.character(keepCols)
  initialFill <- vector("list", length(keepCols))
  names(initialFill) <- as.character(keepCols)
  missingCounts <- integer(length(keepCols))
  names(missingCounts) <- as.character(keepCols)

  missingIndexList <- vector("list", length(keepCols))
  missingRowsByCol <- vector("list", length(keepCols))
  observedRowsByCol <- vector("list", length(keepCols))
  names(missingRowsByCol) <- as.character(keepCols)
  names(observedRowsByCol) <- as.character(keepCols)

  for (i in seq_along(keepCols)) {
    observedValues <- targetMatrix[, i]
    observedValues <- observedValues[!is.na(observedValues)]

    # Decide whether the current target covariate is binary and choose initial fill
    binaryCols[i] <- isBinaryCovariate(observedValues)
    initialFill[[i]] <- if (binaryCols[i]) {
      # Initialize with mode for binary targets
      computeBinaryMode(observedValues)
    } else {
      # Initialize with mean for continuous targets
      fillValue <- mean(observedValues, na.rm=TRUE)
      if (!is.finite(fillValue)) {
        fillValue <- 0
      }
      fillValue
    }

    # Keep track of which rows are observed and missing for this target covariate
    missingRows <- rowIds[missingMask[, i]]
    observedRowsByCol[[i]] <- which(!missingMask[, i])
    missingCounts[i] <- length(missingRows)
    missingRowsByCol[[i]] <- which(missingMask[, i])

    # Fill missing target values with the starting values used in the first iteration
    filledTargetMatrix[missingMask[, i], i] <- initialFill[[i]]

    missingIndexList[[i]] <- if(length(missingRows) == 0) {
      data.frame(rowId = integer(), covariateId = numeric())
    } else {
      data.frame(rowId = missingRows, covariateId = keepCols[i])
    }
  }

  list(
    rowIds = rowIds,
    originalColCount = length(imputableIds),
    keepCols = keepCols,
    removedCols = removedCols,
    zeroDefaultPredictorIds = zeroDefaultPredictorIds,
    sparsePredictorMatrix = buildSparseZeroPredictorMatrix(
      rowIds = rowIds,
      predictorIds = zeroDefaultPredictorIds,
      observedCovariates = prepared$covariates,
      rowIdIndex = prepared$rowIdIndex
    ),
    targetMatrix = targetMatrix,
    currentFilledTargetMatrix = filledTargetMatrix,
    missingMask = missingMask,
    keepColIndex = keepColIndex,
    missingIndex = if(length(missingIndexList) == 0) {
      tibble::tibble(rowId = integer(), covariateId = numeric())
    } else {
      dplyr::bind_rows(missingIndexList) %>%
        dplyr::arrange(.data$covariateId, .data$rowId)
    },
    missingRowsByCol = missingRowsByCol,
    observedRowsByCol = observedRowsByCol,
    missingCounts = missingCounts,
    binaryCols = binaryCols,
    initialFill = initialFill
  )
}

#' @title Prepare the iterative xboost apply state
#' @description This helper prepares matrices, indices and storage objects needed
#' to apply a fitted iterative xgboost imputer to new data
#' @param prepared A prepared PLP imputation object
#' @param fittedImputer Object with fitted iterative xgboost imputation settings
#' @return A list containing the prepared apply state for iterative xgboost imputation
prepareXgboostApplyState <- function(prepared, fittedImputer) {
  rowIds <- sort(unique(prepared$rowIds))
  keepCols <- fittedImputer$keepCols

  # Rebuild the target matrix for the covariates retained during training
  targetMatrix <- subsetPreparedTargetMatrix(
    prepared = prepared,
    covariateIds = keepCols
  )

  missingMask <- is.na(targetMatrix)
  filledTargetMatrix <- targetMatrix

  # Recreate the target specific storage using the fitted training time column set
  keepColIndex <- buildLookupIndex(keepCols)
  missingIndexList <- vector("list", length(keepCols))
  missingRowsByCol <- vector("list", length(keepCols))
  names(missingRowsByCol) <- as.character(keepCols)
  missingCounts <- integer(length(keepCols))
  names(missingCounts) <- as.character(keepCols)

  for (i in seq_along(keepCols)) {
    # Reuse the training time initial fill for the current target covariate
    fillValue <- fittedImputer$initialFill[[as.character(keepCols[i])]]
    filledTargetMatrix[missingMask[, i], i] <- fillValue

    # Keep track of which rows are missing for this target in the new data
    missingRows <- rowIds[missingMask[, i]]
    missingCounts[i] <- length(missingRows)
    missingRowsByCol[[i]] <- which(missingMask[, i])
    missingIndexList[[i]] <- if(length(missingRows) == 0) {
      data.frame(rowId = integer(), covariateId = numeric())
    } else {
      data.frame(rowId = missingRows, covariateId = keepCols[i])
    }
  }

  list(
    rowIds = rowIds,
    keepCols = keepCols,
    sparsePredictorMatrix = buildSparseZeroPredictorMatrix(
      rowIds = rowIds,
      predictorIds = fittedImputer$zeroDefaultPredictorIds,
      observedCovariates = prepared$covariates,
      rowIdIndex = prepared$rowIdIndex
    ),
    targetMatrix = targetMatrix,
    currentFilledTargetMatrix = filledTargetMatrix,
    missingMask = missingMask,
    keepColIndex = keepColIndex,
    missingIndex = if (length(missingIndexList) == 0) {
      tibble::tibble(rowId = integer(), covariateId = numeric())
    } else {
      dplyr::bind_rows(missingIndexList) %>%
        dplyr::arrange(.data$covariateId, .data$rowId)
    },
    missingRowsByCol = missingRowsByCol,
    missingCounts = missingCounts,
    binaryCols = fittedImputer$binaryCols
  )
}

#' @title Build a design matrix for the iterative xgboost imputation
#' @description This function constructs the predictor matrix used to fit or apply one iterative
#' xgboost imputation model. It combines sparse zero-default predictors with the current filled values of the other
#' imputable dense target covariates
#' @param state An iterative xgboost imputation training or apply state object
#' @param rowIdIndex The row idices to include in the design matrix
#' @param targetIndex The column index of the target covariate currently being modeled
#' @return A sparse design matrix, or `NULL` when no predictors are available
buildXgboostDesignMatrix <- function(state, rowIndex, targetIndex) {
  # Iterative xgboost relies on `Matrix` to build the combined sparse design matrix
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for xgboost imputation.")
  }

  # Start with the sparse zero-default predictors that are always available regardless
  # of which dense target is currently being modeled
  sparseBlock <- state$sparsePredictorMatrix[rowIndex, , drop=FALSE]
  denseBlock <- NULL

  if (length(state$keepCols) > 1L) {
    # use the other dense target covariates as predictors, excluding the current target so the
    # model does not predict a variable from itself
    denseIndex <- seq_along(state$keepCols)
    denseIndex <- denseIndex[denseIndex != targetIndex]
    denseBlock <- Matrix::Matrix(
      state$currentFilledTargetMatrix[rowIndex, denseIndex, drop=FALSE],
      sparse = TRUE
    )
    colnames(denseBlock) <- paste0("dense_", state$keepCols[denseIndex])
  }

  # Return `NULL` when there are no sparse zero predictors and no dense target
  # predictors
  if (ncol(sparseBlock) == 0 && is.null(denseBlock)) {
    return(NULL)
  }

  # Return dense predictors if no sparse zero predictors exist
  if (ncol(sparseBlock) == 0) {
    return(denseBlock)
  }

  # Return sparse zero predictors if no dense predictors exist
  if (is.null(denseBlock) || ncol(denseBlock) == 0) {
    return(sparseBlock)
  }

  # Else combine both blocks into one matrix
  cbind(sparseBlock, denseBlock)
}

#' @title Build iterative xgboost imputation training argumnets
#' @description
#' This function constructs the xgboost parameter list used to fit one imputation model,
#' adapting the objective and evaluation metric to the target covariate type
#' @param settings The iterative xgboost imputer settings
#' @param isBinaryTarget Logical variable indicating whether the current target covariate is binary
#' @return A named list of xgboost training parameters
buildXgboostArguments <- function(settings, isBinaryTarget) {
  # Start form the shared xgboost hyperparametrs stored in the imputer settings
  params <- list(
    tree_method = settings$treeMethod,
    eta = settings$eta,
    gamma = settings$gamma,
    max_depth = settings$maxDepth,
    max_delta_step = settings$maxDeltaStep,
    subsample = settings$subsample,
    sampling_method = settings$samplingMethod,
    colsample_bytree = settings$colsampleBytree,
    colsample_bylevel = settings$colsampleBylevel,
    colsample_bynode = settings$colsampleBynode,
    min_child_weight = settings$minChildWeight,
    lambda = settings$lambda,
    alpha = settings$alpha,
    max_leaves= settings$maxLeaves,
    max_bin = settings$maxBin,
    num_parallel_tree = settings$numParallelTree,
    verbosity = 0,
    seed = settings$seed
  )

  # Add optional compute settings only when they were explicitly configured
  if (!is.null(settings$device)) {
    params$device <- settings$device
  }

  if (!is.null(settings$numThreads)) {
    params$nthread <- settings$numThreads
  }

  # Choose the xgboost objective and metric to match the target type
  if (isBinaryTarget) {
    params$objective <- "binary:logistic"
    params$eval_metric <- "logloss"
  } else {
    params$objective <- "reg:squarederror"
    params$eval_metric <- "rmse"
  }
  params
}

#' @title Fit one iterative xgboost imputation model
#' @description
#' This function fits the iterative xgboost imputation model for one covariate and, unless
#' `fitModelOnly = TRUE`, also predicts the missing value for that target
#' @param state An iterative xgboost training state object
#' @param targetCovariateId The covariate id currently being imputed
#' @param settings The iterative xboost imputer settings
#' @param fitModelOnly Logical variable indicating whether only the fitted model should
#' be returned without generating predictions for missing values
#' @return A list with elements `model` and `predictions`
fitSingleXgboostModel <- function(state, targetCovariateId, settings, fitModelOnly = FALSE) {
  # Map the target covariate id back to its column position in the prepared state
  targetIndex <- unname(state$keepColIndex[[normalizeIdKey(targetCovariateId)]])
  observedRows <- state$observedRowsByCol[[targetIndex]]
  missingRows <- state$missingRowsByCol[[targetIndex]]

  # If the target has no observed values, there is nothing to train on
  if (length(observedRows) == 0) {
    return(list(model = NULL, predictions = NULL))
  }

  response <- state$targetMatrix[observedRows, targetIndex]
  isBinaryTarget <- isTRUE(state$binaryCols[[as.character(targetCovariateId)]])

  if (isBinaryTarget) {
    # Convert binary targets to numeric 0/1 labels for xgboost classification
    response <- as.numeric(as.integer(response))
    # If the target is constant return a constant model instead of fitting xgboost
    if (length(unique(response)) <= 1) {
      constantValue <- as.numeric(response[1])
      predictions <- if (fitModelOnly || length(missingRows) == 0) NULL else rep(constantValue, length(missingRows))
      return(list(model = list(kind = "constant_binary", value = constantValue),
                  predictions = predictions
                  ))
    }
  } else if (length(unique(response)) <= 1) {
    # Do the same for non binary targets with no variation
    constantValue <- as.numeric(response[1])
    predictions <- if (fitModelOnly || length(missingRows) == 0) NULL else rep(constantValue, length(missingRows))
    return(list(
      model = list(kind = "constant_numeric", value = constantValue),
      predictions = predictions
    ))
  }

  # Build the predictor matrix from sparse zero-default covariates and the current
  # filled valus of the other dense target covariates
  trainingDesign <- buildXgboostDesignMatrix(
    state = state,
    rowIndex = observedRows,
    targetIndex = targetIndex
  )

  # If no predictors are available, fall back to a constant imputation rule
  if(is.null(trainingDesign) || ncol(trainingDesign) ==0) {
    constantValue <- if (isBinaryTarget) computeBinaryMode(response) else mean(response, na.rm = TRUE)
    predictions <- if(fitModelOnly || length(missingRows) == 0) NULL else rep(constantValue, length(missingRows))
    return(list(
      model = list(
        kind = if (isBinaryTarget) "constant_binary" else "constant_numeric",
        value = constantValue
      ),
      predictions = predictions
    ))
  }

  # Fit XGboost model for the current target covariate
  dtrain <- xgboost::xgb.DMatrix(data = trainingDesign, label = response)
  fit <- xgboost::xgb.train(
    params = buildXgboostArguments(settings = settings, isBinaryTarget = isBinaryTarget),
    data = dtrain,
    nrounds = settings$nrounds,
    verbose = 0
  )

  predictions <- NULL
  if (!fitModelOnly && length(missingRows) > 0) {
    # Predict only the rows where the target covariate is missing
    predictionDesign <- buildXgboostDesignMatrix(
      state = state,
      rowIndex = missingRows,
      targetIndex = targetIndex
    )
    predictions <- predict(fit, newdata = xgboost::xgb.DMatrix(data = predictionDesign))
    # Convert predictions back to the imputed target scale
    if (isBinaryTarget) {
      predictions <- as.numeric(predictions >= 0.5)
    } else {
      predictions <- as.numeric(predictions)
    }
  }

  list(
    model = list(
      kind = if (isBinaryTarget) "xgboost_binary" else "xgboost_numeric",
      fit = fit
    ),
    predictions = predictions
  )
}

#' @title Predict from one iterative xgboost model
#' @description This functin generates predictions for one fitted iterative xgboost imputation
#' model, handling both fitted xgboost models and constant fallback models
#' @param modelObject A fitted model object returned by `fitSingleXgboostModel()`
#' @param newData A predictor matrix for the rows to be imputed
#' @param rowCount The number of rows to predict. Used when `newData` is `NULL` or  when
#' constant predictions are returned
#' @return A numeric vector of predicted values, or `NULL` when no model is available
predictSingleXgboostModel <- function(modelObject, newdata, rowCount = NULL) {
  # Return NULL when no fitted model is available for the target covariate
  if (is.null(modelObject)) {
    return(NULL)
  }

  # Derive the number of rows when it was not supplied explicitly
  if (is.null(rowCount)) {
    rowCount <- if (is.null(newdata)) 0L else nrow(newdata)
  }

  # constant fallback models return the same value for every requested row
  if (identical(modelObject$kind, "constant_numeric")) {
    return(rep(as.numeric(modelObject$value), rowCount))
  }

  if (identical(modelObject$kind, "constant_binary")) {
    return(rep(as.numeric(modelObject$value), rowCount))
  }

  # Otherwise generate predictions from the fitted xgboost model
  predictions <- predict(modelObject$fit, newdata = xgboost::xgb.DMatrix(data= newdata))
  # COnvert binary model probabilities to 0/1 imputations
  if (identical(modelObject$kind, "xgboost_binary")) {
    return(as.numeric(predictions >= 0.5))
  }

  as.numeric(predictions)
}

#' @title Fit the iterative xgboost imputer
#' @description This function fits the full iterative xgboost imputer on prepared
#' training data. It selects the target covariates to impute, initializes the missing values,
#' performs the iterative update loop, and then fits the final per-covariate models during application
#' @param prepared A prepared PLP imputation object
#' @param settings The iterative xgboost imputer settings
#' @return A fitted iterative xgboost imputer object
fitXgboostImputer <- function(prepared, settings) {
  # prepare the target matrices, sparse predictors, and per column tracking
  state <- prepareXgboostTrainingState(prepared, settings)

  if (isTRUE(settings$verbose)) {
    message(
      "Iterative xgboost: retained ",
      length(state$keepCols),
      " of ",
      state$originalColCount,
      " target covariates after applying missingThreshold = ",
      settings$missingThreshold,
      "."
    )
  }

  # Return an empty fitted object when no target covariates remain eligible
  if (length(state$keepCols) == 0) {
    return(list(
      keepCols = numeric(),
      removedCols = state$removedCols,
      zeroDefaultPredictorIds = numeric(),
      binaryCols = logical(),
      initialFill = list(),
      visitOrder = numeric(),
      modelOrder = numeric(),
      model = list(),
      addMissingIndicator = isTRUE(settings$addMissingIndicator),
      indicatorCovariateIds = numeric(),
      indicatorAnalysisIds = numeric(),
      indicatorColumnNames = character()
    ))
  }

  # Visit columns in order of missingness, optionally reversing that order
  visitOrder <- state$keepCols[order(state$missingCounts)]
  if (isTRUE(settings$decreasing)) {
    visitOrder <- rev(visitOrder)
  }
  visitOrder <- visitOrder[state$missingCounts[as.character(visitOrder)] > 0]

  # Fit final model for all retained columns, including columns with no missing
  # values, after the iterative updates have stabilized
  modelOrder <- c(visitOrder, setdiff(state$keepCols, visitOrder))

  # Initialize the long format missing cells with their per-column starting values
  currentImputedValues <- initializeLongImputations(
    missingIndex = state$missingIndex,
    initialFill = state$initialFill
  )

  if(isTRUE(settings$verbose)) {
    visitRatios <- state$missingCounts[as.character(visitOrder)] / length(state$rowIds)
    message(
      "Xgboost: ",
      length(visitOrder),
      " retained columns contain at least one missing value and will enter the iterative loop. Missing ratio(s) = ",
      if (length(visitRatios) == 0) {
        "none"
      } else {
        paste(formatRatio(visitRatios), collapse = ", ")
      },
      "."
    )
  }

  previousDelta <- Inf

  if (length(visitOrder) > 0) {
    for (iteration in seq_len(settings$maxiter)) {
      if (isTRUE(settings$verbose)) {
        message(sprintf("Xgboost training iteration %s/%s", iteration, settings$maxiter))
      }

      # Compare each iteration against the previous imputed values to decide whether the iterative
      # updates have converged
      previousImputedValues <- currentImputedValues

      for (targetIndex in seq_along(visitOrder)) {
        targetCovariateId <- visitOrder[[targetIndex]]
        targetMissingCount <- state$missingCounts[[as.character(targetCovariateId)]]

        if (isTRUE(settings$verbose)) {
          message(
            "Xgboost: training target ",
            targetIndex,
            "/",
            length(visitOrder),
            " (covariateId = ",
            targetCovariateId,
            ", missing rows = ",
            targetMissingCount,
            ", missing ratio = ",
            formatRatio(targetMissingCount / length(state$rowIds)),
            ")."
          )
        }

        # Fit the model for the current target and update only its missing rows
        # in the current filled target matrix
        fitResult <- fitSingleXgboostModel(
          state = state,
          targetCovariateId = targetCovariateId,
          settings = settings,
          fitModelOnly = FALSE
        )

        if (!is.null(fitResult$predictions) && targetMissingCount > 0) {
          targetIndex <- state$keepColIndex[[normalizeIdKey(targetCovariateId)]]
          rowMask <- state$missingRowsByCol[[targetIndex]]
          state$currentFilledTargetMatrix[rowMask, targetIndex] <- fitResult$predictions
          currentImputedValues <- updateImputedValues(
            currentImputedValues = currentImputedValues,
            targetCovariateId = targetCovariateId,
            predictions = fitResult$predictions
          )
        }
      }

      # Stop when the iterative updates converge or begin to diverge
      currentDelta <- computeDeltaLong(
        previousImputedValues = previousImputedValues,
        currentImputedValues = currentImputedValues,
        missingIndex = state$missingIndex,
        binaryCols = state$binaryCols
      )

      if (isTRUE(settings$verbose)) {
        message(sprintf("Xgboost delta = %.6f", currentDelta))
      }

      if (currentDelta <= settings$tol || currentDelta > previousDelta) {
        if (isTRUE(settings$verbose)) {
          if (currentDelta <= settings$tol) {
            message("Xgboost: stopping because delta is below tolerance.")
          } else {
            message("Xgboost: stopping because delta increased.")
          }
        }
        break
      }
      previousDelta <- currentDelta
    }
  } else if (isTRUE(settings$verbose)) {
    message("Xgboost: skipping iterative fitting because no retained column has missing values.")
  }

  fittedModels <- vector("list", length = length(state$keepCols))
  names(fittedModels) <- as.character(state$keepCols)

  for(targetIndex in seq_along(modelOrder)) {
    targetCovariateId <- modelOrder[[targetIndex]]

    if (isTRUE(settings$verbose)) {
      message(
        "Iterative xgboost: fitting final model ",
        targetIndex,
        "/",
        length(modelOrder),
        " for covariateId = ",
        targetCovariateId,
        "."
      )
    }

    # Refit the final model for each retained target using the stabilized filled target matrix from the iterative procedure
    fittedModels[[as.character(targetCovariateId)]] <- fitSingleXgboostModel(
      state = state,
      targetCovariateId = targetCovariateId,
      settings = settings,
      fitModelOnly = TRUE
    )$model

  }

  list(
    keepCols = state$keepCols,
    removedCols = state$removedCols,
    zeroDefaultPredictorIds = state$zeroDefaultPredictorIds,
    binaryCols = state$binaryCols,
    initialFill = state$initialFill,
    visitOrder = visitOrder,
    modelOrder = modelOrder,
    models = fittedModels,
    addMissingIndicator = isTRUE(settings$addMissingIndicator),
    indicatorCovariateIds = numeric(),
    indicatorAnalysisIds = numeric(),
    indicatorColumnNames = character()
  )
}

#' @title Apply the iterative xgboost imputer
#' @description This function applies a fitted iterative xgboost imputer to
#' prepared data. It rebuilds the predictor state for the new data, initializes missing
#' target values with the stored fill values, and iteratively updates them using the fitted per-covariate
#' models
#' @param prepared A prepared PLp imputation object
#' @param fittedImputer A fitted iterative xgboost imputer object
#' @param settings The iterative xgboost imputer settings
#' @return A list containing the long format imputed values for the missing cells
applyXgboostImputer <- function(prepared, fittedImputer, settings) {
  if (length(fittedImputer$keepCols) == 0) {
    return(list(
      currentImputedValues = tibble::tibble(
        rowId = integer(),
        covariateId = numeric(),
        covariateValue = numeric()
      )
    ))
  }

  # Prepare initialized imputation values and apply state
  state <- prepareXgboostApplyState(prepared, fittedImputer)
  currentImputedValues <- initializeLongImputations(missingIndex = state$missingIndex,
                                                              initialFill = fittedImputer$initialFill)

  if (isTRUE(settings$verbose)) {
    message(
      "Iterative xgboost: apply step received ",
      length(state$rowIds),
      " rows x ",
      length(state$keepCols),
      " columns with ",
      nrow(state$missingIndex),
      " missing cells."
    )
  }

  # Apply iterative procedure
  if (length(fittedImputer$visitOrder) > 0) {
    for (iteration in seq_len(settings$maxiter)) {
      if (isTRUE(settings$verbose)) {
        message(sprintf("Iterative xgboost apply iteration %s", iteration))
      }

      for (targetIndex in seq_along(fittedImputer$visitOrder)) {
        targetCovariateId <- fittedImputer$visitOrder[[targetIndex]]
        targetColIndex <- state$keepColIndex[[normalizeIdKey(targetCovariateId)]]
        targetRows <- state$missingRowsByCol[[targetColIndex]]

        if (length(targetRows) == 0) {
          next
        }

        if (isTRUE(settings$verbose)) {
          message(
            "Iterative xgboost: applying target ",
            targetIndex,
            "/",
            targetCovariateId,
            ", missing rows = ",
            length(targetRows),
            ", missing ratio = ",
            formatRatio(length(targetRows) / length(state$rowIds)),
            ")."
          )
        }

        # Create design matrix with target covariates
        predictionDesign <- buildXgboostDesignMatrix(
          state = state,
          rowIndex = targetRows,
          targetIndex = targetIndex
        )
        # Impute values of a single target variable
        predictions <- predictSingleXgboostModel(
          modelObject = fittedImputer$models[[as.character(targetCovariateId)]],
          newdata = predictionDesign,
          rowCount = length(targetRows)
        )

        # Update previous imputed values
        if (!is.null(predictions)) {
          state$currentFilledTargetMatrix[targetRows, targetColIndex] <- predictions
          currentImputedValues <- updateImputedValues(
            currentImputedValues = currentImputedValues,
            targetCovariateId = targetCovariateId,
            predictions = predictions
          )
        }
      }
    }

  } else if (isTRUE(settings$verbose)) {
    message(
      "Iterative xgboost: apply step skipped iterative updates becasue the fitted visit order is empty."
    )

  }

  list(currentImputedValues = currentImputedValues)
}

#' @title Implement iterative xgboost imputation to work with PLP
#' @description This function imputes the PLP data object using iterative imputation
#' @param trainData A PLP training data object
#' @param featureEngineeringSettings Iterative XGBoost imputation sttings
implementXgboostImputer <- function(trainData, featureEngineeringSettings, done = FALSE) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required for xgboost imputation.")
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for xgboost imputation.")
  }

  if (isTRUE(done) && is.null(featureEngineeringSettings$model)) {
    stop("Fitted iterative xgboost settings must include a trained model when done = TRUE.")
  }

  if (isTRUE(featureEngineeringSettings$verbose)) {
    message(
      if (isTRUE(done)) {
        "Iterative xgboost: applying fitted imputer to new data."
      } else {
        "Iterative xgboost: fitted imputer on training data."
      }
    )
  }

  # Prepare train data
  prepared <- prepareTrainPlpData(trainData)

  if (isTRUE(featureEngineeringSettings$verbose)) {
    message(
      "Iterative xgboost: prepared ",
      length(prepared$rowIds),
      " training rows, ",
      length(prepared$imputableCovariateIds),
      " imputation target covariates, and ",
      length(prepared$predictorCovariateIds),
      " predictor covariates. Target missing cells = ",
      prepared$totalMissingCells,
      "."
    )
  }

  fittedImputer <- if (isTRUE(done)) {
    featureEngineeringSettings$model
  } else {
    fitXgboostImputer(prepared = prepared, settings = featureEngineeringSettings)
  }

  if (!isTRUE(done)) {
    fittedImputer <- addMetadata(
      prepared = prepared,
      fittedImputer = fittedImputer
    )
  }

  # Apply fitted imputation model
  imputedState <- applyXgboostImputer(
    prepared = prepared,
    fittedImputer = fittedImputer,
    settings = featureEngineeringSettings
  )

  # Reconstruct the PLP training data object with imputed values
  rebuilt <- rebuildTrainPlpData(
    trainData = trainData,
    prepared = prepared,
    imputedState = imputedState,
    fittedImputer = fittedImputer
  )

  if (!isTRUE(done)) {
    fittedSettings <-featureEngineeringSettings
    fittedSettings$model <- fittedImputer

    metaData <- attr(rebuilt$covariateData, "metaData")
    if (is.null(metaData)) {
      metaData <- list()
    }
    if(is.null(metaData$featureEngineering)) {
      metaData$featureEngineering <- list()
    }

    # Label the imputation as PLP settings
    metaData$featureEngineering[["xgboostImputer"]] <- list(
      funct = "implementXgboostImputer",
      settings = list(featureEngineeringSettings = fittedSettings)
    )
    attr(rebuilt$covariateData, "metaData") <- metaData
  }
  rebuilt
}






















