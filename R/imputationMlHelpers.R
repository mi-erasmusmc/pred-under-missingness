#' @title Validate an integer
#' @description
#'  This function validates that the input is a single finite integer within requested bounds
#'  @param x The value to validate
#'  @param name The argument name used in error messages
#'  @param minValue The minimum allowed integer value
#'  @param allowNull Logical parameter indicating whether `NULL` is allowed
#'  @return The validated integer value
#'
validateInteger <- function(x, name, minValue = 1L, allowNull = FALSE) {
  # Allow only NULL when requested
  if (is.null(x)) {
    if (isTRUE(allowNull)) {
      return(NULL)
    }
    stop(name, " must not be NULL.")
  }

  # Validate that the input is a single finite integer within bounds
  if (length(x) != 1 || !is.numeric(x) || !is.finite(x) || x != as.integer(x) || x < minValue) {
    stop(name, " must be a single integer >= ", minValue, ".")
  }

  as.integer(x)
}

#' @title Validate a numeric value
#' @description
#' This function validates that the input is a single finite numeric value within the requested bounds
#' @param x The value to validate
#' @param name The argument name used in the error messages
#' @param lower The lower bound
#' @param upper The upper bound
#' @param lowerInclusive Logical parameter indicating whether the lower bound is inclusive
#' @param upperInclusive Logical parameter indicating whether the upper bound is inclusive
#' @param allowNull Logical parameter indicating whether `NULL` is allowed
#' @return The validated numeric value
#'
validateNumeric <- function(x, name, lower = -Inf, upper = Inf, lowerInclusive = TRUE, upperInclusive = TRUE, allowNull = FALSE) {
  # Allow NULL only when requested
  if (is.null(x)) {
    if (isTRUE(allowNull)) {
      return(NULL)
    }
    stop(name, " must not be NULL.")
  }

  # Validate that the input is a single finite numeric value
  if (length(x) != 1 || !is.numeric(x) || !is.finite(x)) {
    stop(name, " must be a single finite numeric value.")
  }

  # Check that the value is within the required bounds
  lowerOk <- if (isTRUE(lowerInclusive)) x >= lower else x > lower
  upperOk <- if (isTRUE(upperInclusive)) x <= upper else x < upper

  if (!lowerOk || !upperOk) {
    stop(name, " is outside the allowed range.")
  }

  as.numeric(x)
}

#' @title Validate a logical value
#' @description
#' This function validates that the input is a single logical value
#' @param x The value to validate
#' @param name The argument name used in the error messages
#' @return The validated logical value
#'
validateLogical<- function(x, name) {
  # Validate that the input is a logical value
  if (length(x) != 1 || !is.logical(x) || is.na(x)) {
    stop(name, " must be TRUE or FALSE.")
  }
  x
}

#' @title Normalize identifier keys
#' @description
#' This function converts identifiers to a consistent character key for matching
#' and look up
#' @param x A vector of identifiers
#' @return A vector of normalized identifier keys
#'
normalizeIdKey <- function(x) {
  # Convert identifiers to trimmed character keys
  trimws(as.character(x))
}

#' @title Normalize flag values
#' @description
#' This function standardizes flag values to "Y" and "N" values
#' @param x A vector containing flag values
#' @param default The default value used when input values are missing or unrecognized
#' @return A character vector containing normalized flag values

normalizeBinaryFlag <- function(x, default = "N") {
  # Return empty vector when no input is supplied
  if (is.null(x)) {
    return(character(x))
  }

  # Standardize logical and numeric flags directly
  normalized <- rep(default, length(x))
  missingMask <- is.na(x)

  if (is.logical(x)) {
    normalized[!missingMask & x] <- "Y"
    normalized[!missingMask & !x] <- "N"
    return(normalized)
  }

  if (is.numeric(x)) {
    normalized[!missingMask & x != 0] <- "Y"
    normalized[!missingMask & x == 0] <- "N"
    return(normalized)
  }

  # Normalize character vectors
  xChar <- toupper(trimws(as.character(x)))
  normalized[!missingMask & xChar %in% c("Y", "YES", "TRUE", "T", "1")] <- "Y"
  normalized[!missingMask & xChar %in% c("N", "NO", "FALSE", "F", "0")] <- "N"
  normalized
}

#' @title Format a ratio for messages
#' @description
#' This function formats numeric ratios as short strings for progress messages and logs
#' @param x A numeric ratio
#' @return A formatted character string
#'

formatRatio <- function(x) {
  # Format a ratio using a consistent number of decimal places
  sprintf("%.3f", as.numeric(x))
}

#' @title Build a lookup index
#' @description
#' This function creates a named position index for a vector of identifiers so values
#' can be matched efficiently by key.
#' @param ids A vector of identifiers
#' @return An integer index named by normalized identifier keys
#'
buildLookupIndex <- function(ids) {
  index <- seq_along(ids)
  names(index) <- normalizeIdKey(ids)
  index
}

#' @title Build a target matrix from long covariate data
#' @description This helper materializes a row-by-covariate matrix from PLP covariates stored in long format
#' @param rowIds The row identifiers to include
#' @param covariateIds The covariate identifiers to include
#' @param observedCovariates A long-format table of observed covariate values
#' @param rowIdIndex Optional lookup index for `rowIds`
#' @param covariateIdIndex Optional lookup index for `covariateIds`
#' @return A matrix with rows indexed by `rowIds` and columns indexed by `covariateIds`
buildTargetMatrix <- function(rowIds, covariateIds, observedCovariates, rowIdIndex = NULL, covariateIdIndex = NULL) {
  # Return an empty matrix when there are no rows or no target covariates
  if (length(rowIds) == 0 || length(covariateIds) == 0) {
    return(matrix(
      numeric(),
      nrow = length(rowIds),
      ncol = length(covariateIds),
      dimnames = list(as.character(rowIds), as.character(covariateIds))
    ))
  }

  # Build lookup indices when they are not supplied
  if(is.null(rowIdIndex)) {
    rowIdIndex <- buildLookupIndex(rowIds)
  }

  if (is.null(covariateIdIndex)) {
    covariateIdIndex <- buildLookupIndex(covariateIds)
  }

  # Fill the matrix using the observed long-format covariate values
  targetMatrix <- matrix(
    NA_real_,
    nrow = length(rowIds),
    ncol = length(covariateIds),
    dimnames = list(as.character(rowIds), as.character(covariateIds))
  )

  rowIndex <- unname(rowIdIndex[normalizeIdKey(observedCovariates$rowId)])
  colIndex <- unname(covariateIdIndex[normalizeIdKey(observedCovariates$covariateId)])
  valid <- !is.na(rowIndex) & !is.na(colIndex)

  if (any(valid)) {
    targetMatrix[cbind(rowIndex[valid], colIndex[valid])] <- observedCovariates$covariateValue[valid]
  }

  targetMatrix
}

#' @title Subset a prepared target matrix
#' @description
#' This function extracts a covariate subset from a previously prepared target matrix while preserving the
#' requested row and column layout
#' @param prepared A prepared imputation data object
#' @param covariateIds The covariate identifiers to extract
#' @return A numeric matrix containing the requested covariate columns
#'
subsetPreparedTargetMatrix <- function(prepared, covariateIds) {
  # Return an empty matrix when no covariate columns are requested
  if (length(covariateIds) == 0) {
    return(matrix(
      numeric(),
      nrow = length(prepared$rowIds),
      ncol = 0,
      dimnames = list(as.character(prepared$rowIds), character())
    ))
  }

  # Map the requested covariates back to the prepared target matrix
  colIndex <- unname(
    prepared$imputableCovariateIdIndex[normalizeIdKey(covariateIds)]
  )
  valid <- !is.na(colIndex)

  # Copy the available columns into the requested output layer
  targetMatrix <- matrix(
    NA_real_,
    nrow = length(prepared$rowIds),
    ncol = length(covariateIds),
    dimnames = list(as.character(prepared$rowIds), as.character(covariateIds))
  )

  if (any(valid)) {
    targetMatrix[, valid] <- prepared$imputableTargetMatrix[, colIndex[valid], drop = FALSE]
  }
  targetMatrix
}

#' @title Add imputation metadata
#' @description This function augments a fitted imputer object with metadata
#' needed to rebuild missing-indicator covariates in the PLP output
#' @param prepared A prepared imputation data object
#' @param fittedImputer A fitted imputer object
#' @return The fitted imputer with missing indicator metadata added when needed
addMetadata <- function(prepared, fittedImputer) {
  # Return unchanged fitted imputer when no indicators are required
  if (!isTRUE(fittedImputer$addMissingIndicator) || length(fittedImputer$keepCols) == 0) {
    return(fittedImputer)
  }

  # Derive new covariate and analysis identifiers for missingness indicators.
  maxCovariateId <- max(prepared$covariateRef$covariateId, na.rm = TRUE)
  maxAnalysisId <- max(prepared$analysisRef$analysisId, na.rm = TRUE)
  keptCount <- length(fittedImputer$keepCols)

  fittedImputer$indicatorCovariateIds <- seq(maxCovariateId + 1, length.out = keptCount)
  fittedImputer$indicatorAnalysisIds <- seq(maxAnalysisId + 1L, length.out = keptCount)
  fittedImputer$indicatorColumnNames <- paste0(fittedImputer$keepCols, "_missing")
  fittedImputer
}

#' @title Prepare PLP data for training imputation
#' @description
#' This function materializes and annotates the PLP covariate data needed by the
#' iterative ML imputers
#' @param trainData A PLP dataobject
#' @return A prepared object containing the covariate matrices, metadata, and missingness
#' summaries used during imputer fitting and application
#'
prepareTrainPlpData <- function(trainData) {
  # Materialize the covariates and reference tables from the PLP data object
  covariates <- collectIfNeeded(trainData$covariateData$covariates)
  covariateRef <- collectIfNeeded(trainData$covariateData$covariateRef)
  analysisRef <- collectIfNeeded(trainData$covariateData$analysisRef)
  rowIds <- sort(unique(trainData$labels$rowId))
  rowIdIndex <- buildLookupIndex(rowIds)
  totalRows <- length(rowIds)

  # Validate that rowId/covariateId pairs are unique
  duplicatePairs <- covariates %>%
    dplyr::count(.data$rowId, .data$covariateId, name = "n") %>%
    dplyr::filter(.data$n > 1)

  if (nrow(duplicatePairs) > 0) {
    stop("Machine learning based imputation requires unique rowId/covariateId pairs.")
  }

  # Derive covariate meta data such as binary status and missingness behaviour
  covariateInfo <- covariateRef %>%
    dplyr::left_join(analysisRef, by = "analysisId")

  observedCounts <- covariates %>%
    dplyr::distinct(.data$rowId, .data$covariateId) %>%
    dplyr::count(.data$covariateId, name = "nObserved")

  observedValueSummary <- covariates %>%
    dplyr::group_by(.data$covariateId) %>%
    dplyr::summarise(
      inferredIsBinary = isBinaryCovariate(.data$covariateValue),
      .groups = "drop"
    )

  if (!"missingMeansZero" %in% names(covariateInfo)) {
    covariateInfo$missingMeansZero <- NA_character_
  }

  if (!"isBinary" %in% names(covariateInfo)) {
    covariateInfo$isBinary <- NA_character_
  }

  covariateInfo <- covariateInfo %>%
    dplyr::left_join(observedCounts, by = "covariateId") %>%
    dplyr::left_join(observedValueSummary, by = "covariateId") %>%
    dplyr::mutate(
      isBinary = normalizeBinaryFlag(.data$isBinary, default = NA_character_),
      inferredIsBinary = dplyr::coalesce(.data$inferredIsBinary, FALSE),
      isBinary = dplyr::coalesce(
        .data$isBinary,
        dplyr::if_else(.data$inferredIsBinary, "Y", "N")
      ),
      missingMeansZero = normalizeBinaryFlag(
        .data$missingMeansZero,
        default = NA_character_
      ),
      missingMeansZero = dplyr::coalesce(
        .data$missingMeansZero,
        dplyr::if_else(.data$isBinary == "Y", "Y", "N")
      ),
      nObserved = dplyr::coalesce(.data$nObserved, 0L),
      observedFraction = if (totalRows > 0) .data$nObserved / totalRows else 0
    )

  # Identify imputable target covariateIds and build the target matrix
  imputableCovariateIds <- covariateInfo %>%
    dplyr::filter(
      .data$isBinary == "N",
      .data$missingMeansZero == "N",
      .data$nObserved > 0L
    ) %>%
    dplyr::pull(.data$covariateId) %>%
    unique()

  imputableCovariateIds <- sort(imputableCovariateIds)
  imputableCovariateIdIndex <- buildLookupIndex(imputableCovariateIds)

  predictorCovariateIds <- covariates %>%
    dplyr::pull(.data$covariateId) %>%
    unique() %>%
    sort()

  imputableTargetMatrix <- buildTargetMatrix(
    rowIds = rowIds,
    covariateIds = imputableCovariateIds,
    observedCovariates = covariates,
    rowIdIndex = rowIdIndex,
    covariateIdIndex = imputableCovariateIdIndex
  )

  imputableObserved <- covariates %>%
    dplyr::filter(.data$covariateId %in% imputableCovariateIds, !is.na(.data$covariateValue)) %>%
    dplyr::distinct(.data$rowId, .data$covariateId, .keep_all = TRUE) %>%
    dplyr::arrange(.data$covariateId, .data$rowId)

  missingMask <- is.na(imputableTargetMatrix)
  totalMissingCells <- sum(missingMask)

  # Summarize the missingness structure
  missingInfo <- if(length(imputableCovariateIds) == 0) {
    tibble::tibble(
      covariateId = numeric(),
      counts = integer(),
      missing = numeric()
    )
  } else {
    tibble::tibble(
      covariateId = imputableCovariateIds,
      counts = as.integer(colSums(!missingMask)),
      missing = as.numeric(colMeans(missingMask))
    )
  }

  list(
    rowIds = rowIds,
    rowIdIndex = rowIdIndex,
    covariates = covariates,
    covariateRef = covariateRef,
    analysisRef = analysisRef,
    covariateInfo = covariateInfo,
    imputableCovariateIds = imputableCovariateIds,
    predictorCovariateIds = predictorCovariateIds,
    imputableCovariateIdIndex = imputableCovariateIdIndex,
    imputableTargetMatrix = imputableTargetMatrix,
    imputableObserved = imputableObserved,
    missingInfo = missingInfo,
    totalMissingCells = totalMissingCells
  )
}

#' @title Rebuild PLP data after imputation
#' @description
#' This function reconstructs a PLP data object from the prepared
#' data, imputed values, and fitted imputer mechanisms
#' @param trainData The original PLP data object
#' @param prepared A prepared imputation-data object
#' @param imputedState The imputed values produced by the imputer
#' @param fittedImputer A fitted imputer object
#' @return A rebuilt PLP data object containing observed, imputed and optional missing-
#' indicator covariates
#'
rebuildTrainPlpData <- function(trainData, prepared,
                                     imputedState,
                                     fittedImputer) {
  # Separate unchanged, observed, and newly imputed covariates
  keptImputableIds <- fittedImputer$keepCols
  removedImputableIds <- fittedImputer$removedCols

  unchangedCovariates <- prepared$covariates %>%
    dplyr::filter(!(.data$covariateId %in% prepared$imputableCovariateIds))

  observedImputableCovariates <- prepared$covariates %>%
    dplyr::filter(.data$covariateId %in% keptImputableIds, !is.na(.data$covariateValue))

  imputedCovariates <- imputedState$currentImputedValues %>%
    dplyr::transmute(
      rowId = .data$rowId,
      covariateId = .data$covariateId,
      covariateValue = .data$covariateValue
    )
  indicatorCovariates <- tibble::tibble(
    rowId = integer(),
    covariateId = integer(),
    covariateValue = numeric()
  )

  newCovariateRef <- prepared$covariateRef[0, , drop=FALSE]
  newAnalysisRef <- prepared$analysisRef[0, , drop=FALSE]

  # Build missing indicator covariates and reference rows when requested
  if (isTRUE(fittedImputer$addMissingIndicator) && length(fittedImputer$indicatorCovariateId) > 0) {

    indicatorMap <- stats::setNames(
      fittedImputer$indicatorCovariateIds,
      fittedImputer$keepCols
    )

    indicatorCovariates <- imputedState$currentImputedValues %>%
      dplyr::transmute(
        rowId = .data$rowId,
        covariateId = unname(indicatorMap[as.character(.data$covariateId)]),
        covariateValue = 1
      )

    newCovariateRef <- buildMissingIndicatorCovariateRef(
      covariateRef = prepared$covariateRef,
      sourceCovariateIds = keptImputableIds,
      indicatorCovariateIds = fittedImputer$indicatorCovariateIds,
      indicatorAnalysisIds = fittedImputer$indicatorAnalysisIds
    )

    newAnalysisRef <- buildMissingIndicatorAnalysisRef(
      analysisRef = prepared$analysisRef,
      sourceCovariateIds = keptImputableIds,
      sourceInfo = prepared$covariateInfo,
      indicatorAnalysisIds = fittedImputer$indicatorAnalysisIds
    )
  }

  # Reassemble the PLP covariate and reference tables in output form
  keptBaseIds <- sort(unique(c(setdiff(prepared$covariateRef$covariateId, removedImputableIds), keptImputableIds)))
  updatedCovariateRef <- prepared$covariateRef %>%
    dplyr::filter(.data$covariateId %in% keptBaseIds) %>%
    dplyr::bind_rows(newCovariateRef)

  updatedAnalysisRef <- prepared$analysisRef %>%
    dplyr::semi_join(updatedCovariateRef, by = "analysisId") %>%
    dplyr::bind_rows(newAnalysisRef) %>%
    dplyr::distinct(.data$analysisId, .keep_all = TRUE)

  updatedCovariates <- dplyr::bind_rows(
    unchangedCovariates,
    observedImputableCovariates,
    imputedCovariates,
    indicatorCovariates
  ) %>% dplyr::arrange(.data$rowId, .data$covariateId)

  plpDataHelper(
    labels = trainData$labels,
    folds = trainData$folds,
    covariates = updatedCovariates,
    covariateRef = updatedCovariateRef,
    analysisRef = updatedAnalysisRef,
    templatePLPData = trainData
  )
}

#' @title Initialize long format imputations
#' @description
#' This function creates the initial long format tale of imputed values for
#' missing cells using the supplied per covariate fill values
#' @param missingIndex A long format table identifying missing row/covariate pairs
#' @param initialFill A list of initial fill values by covariate
#' @return A long format table containing the initialized imputed values
#'
initializeLongImputations <- function(missingIndex, initialFill) {
  # Return an empty table when there are no missing cells to initialize
  if (nrow(missingIndex) == 0) {
    return(tibble::tibble(
      rowId = integer(),
      covariateId = integer(),
      covariateValue = numeric()
    ))
  }

  # Fill each missing cell with its covairate specific starting value
  covariateKeys <- as.character(missingIndex$covariateId)
  fillValues <- vapply(covariateKeys, function(key) initialFill[[key]], numeric(1))

  tibble::tibble(
    rowId = missingIndex$rowId,
    covariateId = missingIndex$covariateId,
    covariateValue = fillValues
  )
}

#' @title Update imputed values for one target covariate
#' @description
#' This function replaces the current imputed values for one covariate with a new set of predictions
#' @param currentImputedValues The current long format imputed value table
#' @param targetCovariateId The covariate identifier being updated
#' @param predictions The new predcited values for the missing cells
#' @return The updated imputed value table
updateImputedValues <- function(currentImputedValues, targetCovariateId, predictions) {
  # Return the current values when there is nothing to update
  targetRows <- currentImputedValues$covariateId == targetCovariateId
  if (!any(targetRows) || is.null(predictions)) {
    return(currentImputedValues)
  }

  # Replace the current values for the requested target covariate
  currentImputedValues$covariateValue[targetRows] <- as.numeric(predictions)
  currentImputedValues
}

#' @title Compute the iteration delta for long format imputation
#' @description This function measures the change between two iterations of imputation
#' and returns the stopping criterion used by the iterative ML imputers
#' @param previousImputedValues The imputed values from the previous iteration
#' @param currentImputedValues The imputed values from the current iteration
#' @param missingIndex A long format table identifying all missing row/covariate pairs
#' @param binaryCols A logical vector indicating which target covariates are binary
#' @return A numeric delta summarizing the change between iterations
computeDeltaLong <- function(previousImputedValues, currentImputedValues, missingIndex, binaryCols) {
  # Return zero when there are no missing cells to compare
  if (nrow(missingIndex) == 0) {
    return(0)
  }

  # Compute separate change measures for numeric and binary covariates
  binaryMask <- binaryCols[as.character(missingIndex$covariateId)]
  binaryMask[is.na(binaryMask)] <- FALSE
  numericMask <- !binaryMask

  numericDelta <- 0
  binaryDelta <- 0

  if (any(numericMask)) {
    previousValues <- previousImputedValues$covariateValue[numericMask]
    currentValues <- currentImputedValues$covariateValue[numericMask]
    denominator <- sum(currentValues^2)
    if(!is.finite(denominator) || denominator <= 0) {
      denominator <- .Machine$double.eps
    }
    numericDelta <- sum((currentValues - previousValues)^2) / denominator
  }

  if(any(binaryMask)) {
    previousValues <- previousImputedValues$covariateValue[binaryMask]
    currentValues <- currentImputedValues$covariateValue[binaryMask]
    binaryDelta <- mean(previousValues != currentValues)
  }

  # Use the largest as the stopping criterion
  max(numericDelta, binaryDelta)
}

#' @title Build missing indicator covariate reference rows
#' @description This function creates new `covariateRef` rows for missingness indicator
#' covaraites derived from imputatble source covariates
#' @param covariateRef The original covariate reference table
#' @param sourceCovariateIds The source covariate identifiers
#' @param indicatorCovariateIds The new indicator covariate identifiers
#' @param indicatorAnalysisIds The new indicator analysis covariates
#' @return A covariate reference table for the new indicator covariates
buildMissingIndicatorCovariateRef <- function(covariateRef, sourceCovariateIds, indicatorCovariateIds, indicatorAnalysisIds) {
  # Copy the source covariate metadata and replace the identifying fields
  indicatorRows <- covariateRef[rep(NA_integer_, length(sourceCovariateIds)), , drop=FALSE]
  sourceLookup <- covariateRef %>% dplyr::filter(.data$covariateId %in% sourceCovariateIds) %>%
    dplyr::distinct(.data$covariateId, .keep_all = TRUE)

  indicatorRows$covariateId <- indicatorCovariateIds

  if ("analysisId" %in% names(indicatorRows)) {
    indicatorRows$analysisId <- indicatorAnalysisIds
  }

  # Rename the indicator covariates so they are clearly marked as missingness flags
  if ("covariateName" %in% names(indicatorRows)) {
    sourceNames <- sourceLookup$covariateName[match(sourceCovariateIds, sourceLookup$covariateId)]
    indicatorRows$covariateName <- paste0(sourceNames, "_missing")
  }

  indicatorRows
}

#' @title Build missing indicator analysis reference rows
#' @description This function creates a new `analysisRef` rows for missingness indicator
#' covariates derived from imputable covariates
#' @param analysisRef The original analysis reference table
#' @param sourceCovariateIds The source covariate identifiers
#' @param sourceInfo A table containing metadata for the source covariates
#' @param indicatorAnalysisIds The new indicator analysis identifiers
#' @return An analysis reference table for the new indicator covariates
buildMissingIndicatorAnalysisRef <- function(analysisRef, sourceCovariateIds, sourceInfo, indicatorAnalysisIds) {
  # Create analysis reference rows that describe the derived missingness indicators
  indicatorRows <- analysisRef[rep(NA_integer_, length(sourceCovariateIds)), , drop=FALSE]
  sourceLookup <- sourceInfo %>%
    dplyr::filter(.data$covariateId %in% sourceCovariateIds) %>%
    dplyr::distinct(.data$covariateId, .keep_all = TRUE)

  indicatorRows$analysisId <- indicatorAnalysisIds

  # Mark the missing indicator analysis with the appropriate binary/missingness metadata
  if ("analysisName" %in% names(indicatorRows)) {
    sourceNames <- sourceLookup$covariateName[match(sourceCovariateIds, sourceLookup$covariateId)]
    indicatorRows$analysisName <- paste0("Missing indicator for ", sourceNames)
  }

  if ("isBinary" %in% names(indicatorRows)) {
    indicatorRows$isBinary <- "Y"
  }

  if ("missingMeansZero" %in% names(indicatorRows)) {
    indicatorRows$missingMeansZero <- "Y"
  }

  indicatorRows
}

#'@title Check whether a covariate is binary
#' @description This function determines whether the observed values for a covariate are binary
#' @param x A vector of observed covariate values
#' @return A logical value indicating whether the values are binary
isBinaryCovariate <- function(x) {
  # Use only unique observed values and check whether they are binary
  observedValues <- unique(x[!is.na(x)])
  length(observedValues) > 0 && all(observedValues %in% c(0,1))
}

#' @title Compute the binary mode
#' @description
#' This function returns the mode value of a binary vector
#' @param x A vector of observed binary covariate values
#' @return The modal binary value
#'
computeBinaryMode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(0)
  }

  counts <- table(x)
  as.numeric(names(counts)[which.max(counts)])
}
