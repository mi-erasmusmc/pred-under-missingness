logitInverse <- function(x) {
  1 / (1 + exp(-x))
}

computeGamma1 <- function(z, targetProb, gamma2 = 1) {
  stats::uniroot(
    function(gamma1) mean(logitInverse(gamma1 + gamma2 * z)) - targetProb,
    c(-50, 50),
    extendInt = "yes"
  )$root
}

makeSeed <- function(masterSeed, simulationId = 0, scenarioId = 0, seedLabel = "split") {
  seedOffsets <- c(
    split = 1,
    trainMissing = 2,
    testMissing = 3,
    sklearnImputer = 4,
    lasso = 5,
    xgboost = 6
  )

  if (!seedLabel %in% names(seedOffsets)) {
    stop("Unknown seed label: ", seedLabel)
  }

  as.numeric(masterSeed + simulationId * 10000 + scenarioId * 100 + seedOffsets[[seedLabel]])
}

calcScore <- function(covariates,
                      rowIds,
                      causeCovariateIds,
                      weights = NULL,
                      standardize = TRUE) {
  if (is.null(weights)) {
    weights <- rep(1, length(causeCovariateIds))
  }

  if (length(causeCovariateIds) != length(weights)) {
    stop(
      "Weights and number of variables causing missingness must have the same length."
    )
  }

  covariates <- collectIfNeeded(covariates)

  weightTable <- tibble::tibble(
    covariateId = causeCovariateIds,
    weight = weights
  )

  causeData <- covariates %>%
    dplyr::filter(
      .data$rowId %in% rowIds,
      .data$covariateId %in% causeCovariateIds
    ) %>%
    dplyr::inner_join(weightTable, by = "covariateId") %>%
    dplyr::mutate(weightedVal = .data$covariateValue * .data$weight) %>%
    dplyr::group_by(.data$rowId) %>%
    dplyr::summarise(
      score = sum(.data$weightedVal),
      .groups = "drop"
    )

  causeData <- tibble::tibble(rowId = rowIds) %>%
    dplyr::left_join(causeData, by = "rowId") %>%
    dplyr::mutate(score = dplyr::coalesce(.data$score, 0))

  if (isTRUE(standardize) && length(unique(causeData$score)) > 1) {
    causeData$score <- as.numeric(scale(causeData$score))
  }

  causeData
}

transformScore <- function(score, type = "RIGHT") {
  type <- toupper(type)
  meanScore <- mean(score, na.rm = TRUE)

  transformed <- switch(
    type,
    RIGHT = score - meanScore,
    LEFT = meanScore - score,
    MID = -abs(score - meanScore),
    TAIL = abs(score - meanScore),
    stop("Type must be either RIGHT, LEFT, MID or TAIL.")
  )

  if (length(unique(transformed)) > 1) {
    transformed <- as.numeric(scale(transformed))
  }

  transformed
}

observedCovariateCheck <- function(covariatesDf, targetCovariateIds) {
  observedRows <- covariatesDf %>%
    dplyr::filter(.data$covariateId %in% targetCovariateIds) %>%
    dplyr::distinct(.data$rowId, .data$covariateId)

  nrRows <- length(unique(covariatesDf$rowId))

  observedCheck <- observedRows %>%
    dplyr::group_by(.data$covariateId) %>%
    dplyr::summarise(
      nrObserved = dplyr::n_distinct(.data$rowId),
      .groups = "drop"
    ) %>%
    dplyr::mutate(nrMissing = nrRows - .data$nrObserved)

  if (any(observedCheck$nrMissing > 0)) {
    stop("At least one target covariate contains missing values before simulation.")
  }

  TRUE
}

checkPlpDataComplete <- function(plpData,
                                 requiredCovariateIds,
                                 treatMissingBinaryAsZero = TRUE) {
  requiredCovariateIds <- sort(unique(requiredCovariateIds))
  totalRows <- dplyr::n_distinct(plpData$labels$rowId)

  covInfo <- plpData$covariateData$covariateRef %>%
    dplyr::filter(.data$covariateId %in% requiredCovariateIds) %>%
    dplyr::inner_join(plpData$covariateData$analysisRef, by = "analysisId") %>%
    collectIfNeeded()

  if (!"missingMeansZero" %in% names(covInfo)) {
    covInfo$missingMeansZero <- NA_character_
  }

  missingIds <- setdiff(requiredCovariateIds, covInfo$covariateId)
  if (length(missingIds) > 0) {
    stop(
      "Required covariates not found in covariateRef/analysisRef: ",
      paste(missingIds, collapse = ", ")
    )
  }

  observed <- plpData$covariateData$covariates %>%
    dplyr::filter(.data$covariateId %in% requiredCovariateIds) %>%
    dplyr::distinct(.data$rowId, .data$covariateId) %>%
    collectIfNeeded() %>%
    dplyr::count(.data$covariateId, name = "nObserved")

  covInfo <- covInfo %>%
    dplyr::mutate(
      missingMeansZero = dplyr::coalesce(
        .data$missingMeansZero,
        dplyr::if_else(.data$isBinary == "Y" & treatMissingBinaryAsZero, "Y", "N")
      )
    )

  check <- tibble::tibble(covariateId = requiredCovariateIds) %>%
    dplyr::left_join(
      covInfo %>% dplyr::select(.data$covariateId, .data$isBinary, .data$missingMeansZero),
      by = "covariateId"
    ) %>%
    dplyr::left_join(observed, by = "covariateId") %>%
    dplyr::mutate(
      nObserved = dplyr::coalesce(.data$nObserved, 0L),
      requireCompleteness = .data$missingMeansZero != "Y",
      nMissing = dplyr::if_else(.data$requireCompleteness, totalRows - .data$nObserved, 0L)
    )

  invalidBinary <- check %>%
    dplyr::filter(.data$isBinary == "Y") %>%
    dplyr::pull(.data$covariateId)

  if (length(invalidBinary) > 0) {
    badBinary <- plpData$covariateData$covariates %>%
      dplyr::filter(.data$covariateId %in% invalidBinary) %>%
      dplyr::filter(!is.na(.data$covariateValue), !(.data$covariateValue %in% c(0, 1))) %>%
      dplyr::distinct(.data$covariateId) %>%
      collectIfNeeded()

    if (nrow(badBinary) > 0) {
      stop(
        "These binary covariates contain values other than 0/1: ",
        paste(badBinary$covariateId, collapse = ", ")
      )
    }
  }

  if (any(check$nMissing > 0)) {
    stop(
      "Simulation input is not complete for the required dense covariates: ",
      paste(check$covariateId[check$nMissing > 0], collapse = ", ")
    )
  }

  invisible(check)
}

reportMissingness <- function(plpData, covariateIds = NULL) {
  covariates <- plpData$covariateData$covariates %>%
    collectIfNeeded() %>%
    dplyr::select(.data$rowId, .data$covariateId)

  if (!is.null(covariateIds)) {
    covariates <- covariates %>%
      dplyr::filter(.data$covariateId %in% covariateIds)
  } else {
    covariateIds <- covariates %>%
      dplyr::distinct(.data$covariateId) %>%
      dplyr::pull(.data$covariateId)
  }

  totalRows <- dplyr::n_distinct(plpData$labels$rowId)

  tibble::tibble(covariateId = sort(unique(covariateIds))) %>%
    dplyr::left_join(
      covariates %>%
        dplyr::distinct(.data$rowId, .data$covariateId) %>%
        dplyr::count(.data$covariateId, name = "nObserved"),
      by = "covariateId"
    ) %>%
    dplyr::mutate(
      nObserved = dplyr::coalesce(.data$nObserved, 0L),
      nMissing = totalRows - .data$nObserved,
      missingFraction = .data$nMissing / totalRows
    ) %>%
    dplyr::arrange(dplyr::desc(.data$nMissing), .data$covariateId)
}

reportMissingFraction <- function(beforePLPData, afterPLPData, covariateIds = NULL) {
  before <- beforePLPData$covariateData$covariates %>%
    collectIfNeeded() %>%
    dplyr::select(.data$rowId, .data$covariateId) %>%
    dplyr::distinct()

  after <- afterPLPData$covariateData$covariates %>%
    collectIfNeeded() %>%
    dplyr::select(.data$rowId, .data$covariateId) %>%
    dplyr::distinct()

  if (!is.null(covariateIds)) {
    before <- before %>% dplyr::filter(.data$covariateId %in% covariateIds)
    after <- after %>% dplyr::filter(.data$covariateId %in% covariateIds)
  } else {
    covariateIds <- sort(unique(before$covariateId))
  }

  tibble::tibble(covariateId = sort(unique(covariateIds))) %>%
    dplyr::left_join(
      before %>% dplyr::count(.data$covariateId, name = "nObservedBefore"),
      by = "covariateId"
    ) %>%
    dplyr::left_join(
      after %>% dplyr::count(.data$covariateId, name = "nObservedAfter"),
      by = "covariateId"
    ) %>%
    dplyr::mutate(
      nObservedBefore = dplyr::coalesce(.data$nObservedBefore, 0L),
      nObservedAfter = dplyr::coalesce(.data$nObservedAfter, 0L),
      nDropped = .data$nObservedBefore - .data$nObservedAfter,
      droppedFraction = dplyr::if_else(
        .data$nObservedBefore > 0,
        .data$nDropped / .data$nObservedBefore,
        NA_real_
      )
    ) %>%
    dplyr::arrange(dplyr::desc(.data$droppedFraction), .data$covariateId)
}

simulateMultivariateMissingness <- function(data,
                                            patterns,
                                            ratio,
                                            freq = NULL,
                                            seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

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

  freq <- freq / sum(freq)

  rowIds <- data$labels %>% dplyr::pull(.data$rowId)
  covariatesDf <- data$covariateData$covariates %>%
    collectIfNeeded() %>%
    dplyr::select(.data$rowId, .data$covariateId, .data$covariateValue)

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

  for (i in seq_along(patterns)) {
    pattern <- patterns[[i]]

    if (is.null(pattern$targetCovariateIds)) {
      stop("Pattern ", i, " must have targetCovariateIds defined.", call. = FALSE)
    }

    targetCovariateIds <- unique(pattern$targetCovariateIds)
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

    assignedRowIds <- assignRows %>%
      dplyr::filter(.data$patternId == i) %>%
      dplyr::pull(.data$rowId)

    if (length(assignedRowIds) == 0) {
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

    targetPairs <- covariatesDf %>%
      dplyr::filter(
        .data$rowId %in% assignedRowIds,
        .data$covariateId %in% targetCovariateIds
      ) %>%
      dplyr::distinct(.data$rowId, .data$covariateId)

    if (is.null(causeCovariateIds)) {
      probData <- tibble::tibble(rowId = assignedRowIds, prob = ratio)
    } else {
      if (is.null(causeWeights)) {
        causeWeights <- rep(1, length(causeCovariateIds))
      }

      causeData <- calcScore(
        covariates = covariatesDf,
        rowIds = assignedRowIds,
        causeCovariateIds = causeCovariateIds,
        weights = causeWeights,
        standardize = standardizeScore
      ) %>%
        dplyr::mutate(z = transformScore(.data$score, type))

      gamma1 <- computeGamma1(
        z = causeData$z,
        targetProb = ratio,
        gamma2 = gamma2
      )

      probData <- causeData %>%
        dplyr::mutate(prob = logitInverse(gamma1 + gamma2 * .data$z)) %>%
        dplyr::select(.data$rowId, .data$prob)
    }

    dropRowIds <- probData %>%
      dplyr::mutate(drop = stats::rbinom(dplyr::n(), size = 1, prob = .data$prob) == 1) %>%
      dplyr::filter(.data$drop) %>%
      dplyr::pull(.data$rowId)

    droppedPairsSelection <- targetPairs %>%
      dplyr::filter(.data$rowId %in% dropRowIds) %>%
      dplyr::mutate(patternId = i)

    droppedPairs[[i]] <- droppedPairsSelection
    patternList[[i]] <- tibble::tibble(
      patternId = i,
      mechanism = mechanism,
      assignedRows = length(assignedRowIds),
      droppedRows = length(unique(dropRowIds)),
      droppedCells = nrow(droppedPairsSelection)
    )
  }

  droppedPairsSelection <- dplyr::bind_rows(droppedPairs)
  patternSummary <- dplyr::bind_rows(patternList)

  if (nrow(droppedPairsSelection) == 0) {
    newCovariates <- covariatesDf
  } else {
    newCovariates <- covariatesDf %>%
      dplyr::anti_join(
        droppedPairsSelection %>% dplyr::select(.data$rowId, .data$covariateId),
        by = c("rowId", "covariateId")
      )
  }

  newPlpData <- plpDataHelper(
    labels = data$labels,
    folds = data$folds,
    covariates = newCovariates,
    covariateRef = data$covariateData$covariateRef %>% collectIfNeeded(),
    analysisRef = data$covariateData$analysisRef %>% collectIfNeeded(),
    templatePLPData = data
  )

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
