#' @title Simulate univariate/multivariate missingenss
#' @description
#' This function simulates univariate or multivariate missingness in PLP data using a list
#' of pattern specifications. The overall amount of missingness is controlled by `ratio`, and
#' the missing rows are distributed across patterns according to `freq`
#' @param data A complete PLP data object
#' @param patterns A list of missingness pattern specifications
#' @param ratio The overall fraction of rows to assign to the missingness patterns
#' @param freq The relative frequences used to distribute missingness across the patterns
#' @param seed An optional random seed
#' @return A list containing the PLP data with simulated missingness, the row assignments,
#' the dropped row-covariate pairs and the pattern summary
#' @export
simulateMultivariateMissingness <- function(data,
                                            patterns,
                                            ratio,
                                            freq = NULL,
                                            seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Validate pattern specifications
  if (!is.list(patterns) || length(patterns) == 0) {
    stop("List of patterns cannot be empty.", call. = FALSE)
  }

  nrOfPatterns <- length(patterns)

  if (is.null(freq)) {
    freq <- rep(1 / nrOfPatterns, nrOfPatterns)
  }

  if (length(freq) != nrOfPatterns) {
    stop("Frequency list must have the same length as the number of patterns.", call. = FALSE)
  }

  if (any(freq < 0)) {
    stop("Frequencies must be nonnegative.", call. = FALSE)
  }

  if (sum(freq) == 0) {
    stop("Frequencies cannot sum to 0.", call. = FALSE)
  }

  # Normalize frequencies
  freq <- freq / sum(freq)

  # Select the rowIds and covariate data
  rowIds <- data$labels %>% dplyr::pull(.data$rowId)
  covariatesDf <- data$covariateData$covariates %>%
    collectIfNeeded() %>%
    dplyr::select(.data$rowId, .data$covariateId, .data$covariateValue)

  # Randomly sample which rows get assigned to which pattern
  assignRows <- tibble::tibble(
    rowId = rowIds,
    patternId = sample.int(
      n = nrOfPatterns,
      size = length(rowIds),
      replace = TRUE,
      prob = freq
    )
  )

  droppedPairs <- vector("list", nrOfPatterns)
  patternList <- vector("list", nrOfPatterns)

  # For each pattern, decide which rowIds will be made missing
  for (i in seq_along(patterns)) {
    pattern <- patterns[[i]]

    if (is.null(pattern$targetCovariateIds)) {
      stop("Pattern ", i, " must have targetCovariateIds defined.", call. = FALSE)
    }

    targetCovariateIds <- unique(pattern$targetCovariateIds)
    # Check whether the data is complete initially
    observedCovariateCheck(covariatesDf, targetCovariateIds)

    causeCovariateIds <- pattern$causeCovariateIds
    causeWeights <- pattern$causeWeights
    type <- if (is.null(pattern$type)) "RIGHT" else toupper(pattern$type)
    gamma2 <- if (is.null(pattern$gamma2)) 1 else pattern$gamma2
    standardizeScore <- if (is.null(pattern$standardizeScore)) TRUE else pattern$standardizeScore
    mechanism <- if (is.null(pattern$mechanism)) {
      if (is.null(causeCovariateIds)) "MCAR" else "CUSTOM"
    } else {
      toupper(pattern$mechanism)
    }

    # Select the rowIds assigned to the current pattern
    assignedRowIds <- assignRows %>%
      dplyr::filter(.data$patternId == i) %>%
      dplyr::pull(.data$rowId)

    if (length(assignedRowIds) == 0) {
      # Initialize a list of row-covariate pairs to drop
      droppedPairs[[i]] <- tibble::tibble(
        rowId = integer(),
        covariateId = integer(),
        patternId = integer()
      )

      patternList[[i]] <- tibble::tibble(
        patternId = i,
        mechanism = mechanism,
        assignedRows = 0L,
        droppedRows = 0L,
        droppedCells = 0L
      )
      next
    }

    # Select the target rowId-covariate pairs assigned to the current pattern
    targetPairs <- covariatesDf %>%
      dplyr::filter(
        .data$rowId %in% assignedRowIds,
        .data$covariateId %in% targetCovariateIds
      ) %>%
      dplyr::distinct(.data$rowId, .data$covariateId)

    # If MCAR: no cause covariates are given
    if (is.null(causeCovariateIds)) {
      # Rows assigned to the pattern get assigned the target ratio as probability
      probData <- tibble::tibble(rowId = assignedRowIds, prob = ratio)
    } else {
      # If MAR/MNAR and no weights are given: each covariate is equally likely to be made missing
      if (is.null(causeWeights)) {
        causeWeights <- rep(1, length(causeCovariateIds))
      }

      # Calculate individual missingness scores
      causeData <- calcScore(
        covariates = covariatesDf,
        rowIds = assignedRowIds,
        causeCovariateIds = causeCovariateIds,
        weights = causeWeights,
        standardize = standardizeScore
      ) %>%
        # Transform the score using the given type transformation
        dplyr::mutate(z = transformScore(.data$score, type))

      # Compute gamma 1
      gamma1 <- computeGamma1(
        z = causeData$z,
        targetProb = ratio,
        gamma2 = gamma2
      )

      # Use the score and gamma1 and the logit-inverse transformation to obtain
      # a missingness probability
      probData <- causeData %>%
        dplyr::mutate(prob = logitInverse(gamma1 + gamma2 * .data$z)) %>%
        dplyr::select(.data$rowId, .data$prob)
    }

    # Use the calculated probability to determine which rowIds will be
    # made missing
    dropRowIds <- probData %>%
      dplyr::mutate(drop = stats::rbinom(dplyr::n(), size = 1, prob = .data$prob) == 1) %>%
      dplyr::filter(.data$drop) %>%
      dplyr::pull(.data$rowId)

    # Select the coresponding row-covariate pairs based on the rows that will be made missing
    droppedPairsSelection <- targetPairs %>%
      dplyr::filter(.data$rowId %in% dropRowIds) %>%
      dplyr::mutate(patternId = i)

    # Summarise pattern assignment and dropped row-covariate pairs in the current pattern
    droppedPairs[[i]] <- droppedPairsSelection
    patternList[[i]] <- tibble::tibble(
      patternId = i,
      mechanism = mechanism,
      assignedRows = length(assignedRowIds),
      droppedRows = length(unique(dropRowIds)),
      droppedCells = nrow(droppedPairsSelection)
    )
  }

  # Summarise pattern assignment and dropped row-covariate pairs across all patterns
  droppedPairsSelection <- dplyr::bind_rows(droppedPairs)
  patternSummary <- dplyr::bind_rows(patternList)

  # Create the new covariate data after inducing missingness
  if (nrow(droppedPairsSelection) == 0) {
    newCovariates <- covariatesDf
  } else {
    newCovariates <- covariatesDf %>%
      dplyr::anti_join(
        droppedPairsSelection %>% dplyr::select("rowId", "covariateId"),
        by = c("rowId", "covariateId")
      )
  }

  # Rebuild the PLP data object
  newPlpData <- plpDataHelper(
    labels = data$labels,
    folds = data$folds,
    covariates = newCovariates,
    covariateRef = data$covariateData$covariateRef %>% collectIfNeeded(),
    analysisRef = data$covariateData$analysisRef %>% collectIfNeeded(),
    templatePLPData = data
  )

  # Return the incomplete PLP dat with pattern-level summary objects
  list(
    plpData = newPlpData,
    droppedPairs = droppedPairsSelection,
    assignments = assignRows,
    patternSummary = patternSummary,
    patterns = patterns,
    freq = freq,
    ratio = ratio
  )
}

#' @title Simulate multivariate missingness while matching final pattern fractions
#'
#' @description This function treats `freq` as the final row-pattern proportions.
#' Rows are assigned to the missingness patterns without replacement, using
#' MAR/MNAR- style weights derived from each pattern's score model.
#' This keeps the realized pattern shares close to the requested values while
#' still letting the specified cause covariates influence who becomes missing.
#'
#' @param data a PLP data object with complete target covariates
#' @param patterns A non-empty list of pattern specifications.
#' @param freq Desired final pattern fractions. Values are normalized to sum to 1.
#' @param seed Optional random seed
#'
#' @return A list with the simulated PLP data, row assignments, dropped pairs and
#' pattern summary.
#' @export
simulateMultivariatePatternMissingness <- function(data, patterns, freq, seed) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Validate pattern specifications
  if (!is.list(patterns) || length(patterns) == 0) {
    stop("List of patterns cannot be empty.", call. = FALSE)
  }

  nrOfPatterns <- length(patterns)
  if (missing(freq) || is.null(freq)) {
    stop("Provide freq for the desired final pattern shares.", .call = FALSE)
  }

  if (length(freq) != nrOfPatterns) {
    stop("Frequency list must have the same length as the number of patterns.", .call = FALSE)
  }

  if (any(!is.finite(freq)) || any(freq < 0)) {
    stop("Frequencies must be finite and nonnegative.", .call = FALSE)
  }

  if (sum(freq) == 0) {
    stop("Frequencies cannot sum to 0.", .call = FALSE)
  }

  # Normalize frequencies
  freq <- freq / sum(freq)

  # Extract the rows that can be assigned and the observed covariate values
  # used to compute pattern-specific selection weights
  rowIds <- data$labels %>% dplyr::pull(.data$rowId)
  covariatesDf <- data$covariateData$covariates %>%
    collectIfNeeded() %>%
    dplyr::select("rowId", "covariateId", "covariateValue")

  # Convert the requested pattern shares into integer row counts
  targetCounts <- allocatePatternCounts(length(rowIds), freq)
  names(targetCounts) <- paste0("pattern", seq_len(nrOfPatterns))

  # First assign patterns with more target covariates so the most restrictive patterns
  # claim rows before the simpler patterns do
  targetLengths <- vapply(
    patterns,
    function(pattern) length(unique(pattern$targetCovariateIds)),
    integer(1)
  )

  assignmentOrder <- c(order(targetLengths, decreasing = TRUE), integer())
  assignmentOrder <- unique(assignmentOrder)

  assignments <- vector("list", nrOfPatterns)
  droppedPairs <- vector("list", nrOfPatterns)
  patternList <- vector("list", nrOfPatterns)
  remainingRowIds <- rowIds

  # Loop over all patterns in the sorted order
  for (i in assignmentOrder) {
    pattern <- patterns[[i]]

    if (is.null(pattern$targetCovariateIds)) {
      stop("Pattern ", i, " must have targetCovariateIds defined.", call. = FALSE)
    }

    targetCovariateIds <- unique(pattern$targetCovariateIds)
    observedCovariateCheck(covariatesDf, targetCovariateIds)

    # Extract mechanism from the current pattern
    mechanism <- if (is.null(pattern$mechanism)) {
      if (is.null(pattern$causeCovariateIds)) "MCAR" else "CUSTOM"
    } else {
      toupper(pattern$mechanism)
    }

    # Select the rows for this pattern that are not yet used by another pattern
    # assignment
    assignedRowIds <- selectPatternRows(
      remainingRowIds = remainingRowIds,
      targetCount = targetCounts[[i]],
      covariatesDf = covariatesDf,
      pattern = pattern
    )
    print(length(assignedRowIds))
    assignments[[i]] <- tibble::tibble(
      rowId = assignedRowIds,
      patternId = i
    )

    # Remove from the remaining rows so they cannot be used for another pattern
    # assignment
    remainingRowIds <- setdiff(remainingRowIds, assignedRowIds)

    # Continue the process after all rows have been assigned to a pattern
    if (length(assignedRowIds) == 0 || length(targetCovariateIds) == 0) {
      # Initialize a list of row-covariate pairs to drop
      droppedPairs[[i]] <- tibble::tibble(
        rowId = integer(),
        covariateId = integer(),
        patternId = integer()
      )

      patternList[[i]] <- tibble::tibble(
        patternId = i,
        mechanism = mechanism,
        assignedRows = length(assignedRowIds),
        droppedRows = 0L,
        droppedCells = 0L
      )
      next
    }

    # Select the row-covariate target pairs that are made missing
    targetPairs <- covariatesDf %>%
      dplyr::filter(
        .data$rowId %in% assignedRowIds,
        .data$covariateId %in% targetCovariateIds
      ) %>%
      dplyr::distinct(.data$rowId, .data$covariateId) %>%
      dplyr::mutate(patternId = i)

    droppedPairs[[i]] <- targetPairs
    patternList[[i]] <- tibble::tibble(
      patternId = i,
      mechanism = mechanism,
      assignedRows = length(assignedRowIds),
      droppedRows = length(unique(assignedRowIds)),
      droppedCells = nrow(targetPairs)
    )
  }

  # Summarise the missing row-covariate pairs across all patterns
  assignRows <- dplyr::bind_rows(assignments) %>%
    dplyr::arrange(.data$rowId)
  droppedPairsSelection <- dplyr::bind_rows(droppedPairs)
  patternSummary <- dplyr::bind_rows(patternList) %>%
    dplyr::mutate(targetCount = unname(targetCounts[patternId]))

  newCovariates <- if (nrow(droppedPairsSelection) == 0) {
    covariatesDf
  } else {
    covariatesDf %>%
      dplyr::anti_join(
        droppedPairsSelection %>% dplyr::select(.data$rowId, .data$covariateId), by = c("rowId", "covariateId")
      )
  }

  # Rebuild the PLP data object
  newPlpData <- plpDataHelper(
    labels = data$labels,
    folds = data$folds,
    covariates = newCovariates,
    covariateRef = data$covariateData$covariateRef %>% collectIfNeeded(),
    analysisRef = data$covariateData$analysisRef %>% collectIfNeeded(),
    templatePLPData = data
  )

  # Summarise the pattern assignments and dropped row-covariate pairs
  list(
    plpData = newPlpData,
    droppedPairs = droppedPairsSelection,
    assignments = assignRows,
    patternSummary = patternSummary,
    patterns = patterns,
    freq = freq,
    targetCounts = targetCounts
  )
}
