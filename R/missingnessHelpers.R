#' @title Compute the inverse logit
#' @description
#' This function maps numeric values to probabilities on a (0,1) scale using the
#' inverse logit transformation
#' @param x A numeric value or vector
#' @return A numeric vector of probabilities
#' @export
logitInverse <- function(x) {
  # Transform the input from log odds to probabilities
  1 / (1 + exp(-x))
}

#' @title Compute the intercept for a target missingness probability
#' @description
#' This function computes the intercept term needed to match a target missingness
#' probability for a given score vector and slope parameter
#' @param z A numeric score vector
#' @param targetProb The target missingness probability
#' @param gamma2 The slope parameter applied to the score
#' @return A numeric intercept value
#' @export
computeGamma1 <- function(z, targetProb, gamma2 = 1) {
  # Solve for the intercept that matches the requested average probability
  stats::uniroot(
    function(gamma1) mean(logitInverse(gamma1 + gamma2 * z)) - targetProb,
    c(-50, 50),
    extendInt = "yes"
  )$root
}

#' @title Derive a reproducable seed
#' @description
#' This function constructs a seed for one simulation component from the master seed, simulation identifier,
#' scenario identifier and seed label
#' @param masterSeed The master seed for the full simulation study
#' @param simulationId The simulation run identifier
#' @param scenarioId The scenario identifier
#' @param seedLabel A label identifying the simulation component
#' @return A numeric seed value
#'
#' @export
makeSeed <- function(masterSeed, simulationId = 0, scenarioId = 0, seedLabel = "split") {
  # Map each simulation component to a fixed reference value
  seedOffsets <- c(
    split = 1,
    trainMissing = 2,
    testMissing = 3,
    sklearnImputer = 4,
    lasso = 5,
    xgboost = 6,
    transformer = 7,
    iterativeXGBoostImputer = 8,
    pmmTrain = 9,
    pmmTest = 10
  )

  if (!seedLabel %in% names(seedOffsets)) {
    stop("Unknown seed label: ", seedLabel)
  }

  # Combine the seed components deterministically
  as.numeric(masterSeed + simulationId * 10000 + scenarioId * 100 + seedOffsets[[seedLabel]])
}

#' @title Compute a missingness score
#' @description
#' This function computes a score for each individual used to drive missingness assignment from one or more
#' cause variables
#' @param covariates A long format covariate table
#' @param rowIds The row identifiers to score
#' @param causeCovariateIds The covariate identifiers used to build the score
#' @param weights Optional weights for the cause covariates
#' @param standardize A logical paramter indicating whether the score should be
#' standardized
#' @return A data frame or vector containing scores for each individual
#' @export
calcScore <- function(covariates,
                      rowIds,
                      causeCovariateIds,
                      weights = NULL,
                      standardize = TRUE) {
  # Use equal weights when no covariate specific weights are supplied
  if (is.null(weights)) {
    weights <- rep(1, length(causeCovariateIds))
  }

  # Each covariate must have exactly one weight
  if (length(causeCovariateIds) != length(weights)) {
    stop(
      "Weights and number of variables causing missingness must have the same length."
    )
  }

  covariates <- collectIfNeeded(covariates)

  # Store the cause covariate weights in a lookup table for joining
  weightTable <- tibble::tibble(
    covariateId = causeCovariateIds,
    weight = weights
  )

  # Select the cause variables driving missingness
  causeData <- covariates %>%
    dplyr::filter(
      .data$rowId %in% rowIds,
      .data$covariateId %in% causeCovariateIds
    ) %>%
    # Join with corresponding weights
    dplyr::inner_join(weightTable, by = "covariateId") %>%
    # Multiply each observed cause covariate value by its weight
    dplyr::mutate(weightedVal = .data$covariateValue * .data$weight) %>%
    dplyr::group_by(.data$rowId) %>%
    # Sum the weighted values to obtain one missingness score per row
    dplyr::summarise(
      score = sum(.data$weightedVal),
      .groups = "drop"
    )

  causeData <- tibble::tibble(rowId = rowIds) %>%
    dplyr::left_join(causeData, by = "rowId") %>%
    # Rows with no observed cause covariate get a score of zero
    dplyr::mutate(score = dplyr::coalesce(.data$score, 0))

  # Standardize the score if requested, so mechanisms depend on relative
  # row ranking
  if (isTRUE(standardize) && length(unique(causeData$score)) > 1) {
    causeData$score <- as.numeric(scale(causeData$score))
  }

  causeData
}

#' @title Transform a missingness score
#' @description
#' This function transforms a row level score so different missingness shapes can be
#' created, such as higher missingness in the right tail, left tail, middle or
#' extremes of the distribution
#' @param score A numeric score vector
#' @param type The score transformation type: `RIGHT`, `LEFT`, `MID`, `TAIL`
#' @return A transformed numeric score vector
#' @export
transformScore <- function(score, type = "RIGHT") {
  # Normalize the transformation type label and calculate the mean score
  type <- toupper(type)
  meanScore <- mean(score, na.rm = TRUE)

  transformed <- switch(
    type,
    # Higher original scores become more likely to be missing
    RIGHT = score - meanScore,
    # Lower original scores become more likely to be missing
    LEFT = meanScore - score,
    # Scores near the center become more likely to be missing
    MID = -abs(score - meanScore),
    # Scores in the tails become more likely to be missing
    TAIL = abs(score - meanScore),
    stop("Type must be either RIGHT, LEFT, MID or TAIL.")
  )

  if (length(unique(transformed)) > 1) {
    transformed <- as.numeric(scale(transformed))
  }
  transformed
}

#' @title Check that target covariates are observed
#' @description
#' This function checks whether the target covariates are complete before simulating
#' missingness
#' @param covariatesDf A long format covariate table
#' @param targetCovariateIds The covariate ids to check
#' @return `TRUE` if all target covariates are observed, otherwise it throws an error
#' @export
observedCovariateCheck <- function(covariatesDf, targetCovariateIds) {
  # Keep one row
  observedRows <- covariatesDf %>%
    dplyr::filter(.data$covariateId %in% targetCovariateIds) %>%
    dplyr::distinct(.data$rowId, .data$covariateId)

  # Identify the total number of rows each covariate should have
  nrRows <- length(unique(covariatesDf$rowId))

  observedCheck <- observedRows %>%
    dplyr::group_by(.data$covariateId) %>%
    dplyr::summarise(
      nrObserved = dplyr::n_distinct(.data$rowId),
      .groups = "drop"
    ) %>%
    # Any gap between the total row count and the observed row count means the
    # covariate is already missing before simulation
    dplyr::mutate(nrMissing = nrRows - .data$nrObserved)

  # Stop early when a covariate is not fully observed
  if (any(observedCheck$nrMissing > 0)) {
    stop("At least one target covariate contains missing values before simulation.")
  }

  TRUE
}

#' @title Check that the required PLP covariates are complete
#' @description This helper verifies that the required covariates in a PLP data object
#' are suitable as complete simulation inputs. Dense covariates must be fully observed,
#' whereas binary covariates can be treated as complete when missingness is interpreted as zero
#' @param plpData A PLP data object
#' @param requiredCovariateIds The covariate ids that must be checked
#' @param treatMissingBinaryAsZero Logical parameter indicating whether binary covariates with `missingMeansZero`
#' not explicitly set should be treated as complete by default
#' @return Invisibly returns a dataframe summarising the completeness check. Throws an error
#' when the required covariates are missing or invalid
#' @export
checkPlpDataComplete <- function(plpData,
                                 requiredCovariateIds,
                                 treatMissingBinaryAsZero = TRUE) {
  requiredCovariateIds <- sort(unique(requiredCovariateIds))
  # Determine the expected number of rows
  totalRows <- dplyr::n_distinct(plpData$labels$rowId)

  # Collect the covariate metadata needed to decide which covariates must be
  # fully observed and which can be treated as missing means zero
  covInfo <- plpData$covariateData$covariateRef %>%
    dplyr::filter(.data$covariateId %in% requiredCovariateIds) %>%
    dplyr::inner_join(plpData$covariateData$analysisRef, by = "analysisId") %>%
    collectIfNeeded()

  if (!"missingMeansZero" %in% names(covInfo)) {
    covInfo$missingMeansZero <- NA_character_
  }

  # Stop early if required covariates are missing from the PLP reference tables
  missingIds <- setdiff(requiredCovariateIds, covInfo$covariateId)
  if (length(missingIds) > 0) {
    stop(
      "Required covariates not found in covariateRef/analysisRef: ",
      paste(missingIds, collapse = ", ")
    )
  }

  # Count number of observed values each required covariate
  observed <- plpData$covariateData$covariates %>%
    dplyr::filter(.data$covariateId %in% requiredCovariateIds) %>%
    dplyr::distinct(.data$rowId, .data$covariateId) %>%
    collectIfNeeded() %>%
    dplyr::count(.data$covariateId, name = "nObserved")

  covInfo <- covInfo %>%
    dplyr::mutate(
      # When missing means zero is absent, optionally treat binary covariates as complete by
      # default because unobserved values are interpereted as zero
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
      # Only the continuous covariates require explicit observation in every row
      requireCompleteness = .data$missingMeansZero != "Y",
      nMissing = dplyr::if_else(.data$requireCompleteness, totalRows - .data$nObserved, 0L)
    )

  # For covariates market as binary, verify that all observed values are binary
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

  # Stop when any required continuous variable is incomplete before simulation
  if (any(check$nMissing > 0)) {
    stop(
      "Simulation input is not complete for the required dense covariates: ",
      paste(check$covariateId[check$nMissing > 0], collapse = ", ")
    )
  }

  invisible(check)
}



#' @title Report missingness in PLP data
#' @description
#' This function summarised the observed missingness for the selected covariates
#' in PLP data
#' @param plpData A PLP data object
#' @param covariateIds Optional covariate identifiers to summarise
#' @return A data frame describing the missingness per covariate
#' @export
reportMissingness <- function(plpData, covariateIds = NULL) {
  # Select covariate data
  covariates <- plpData$covariateData$covariates %>%
    collectIfNeeded() %>%
    dplyr::select(.data$rowId, .data$covariateId)

  # If covariate ids are given, filter on those ids
  if (!is.null(covariateIds)) {
    covariates <- covariates %>%
      dplyr::filter(.data$covariateId %in% covariateIds)
  } else {
    # Else filter on all unique covariate ids
    covariateIds <- covariates %>%
      dplyr::distinct(.data$covariateId) %>%
      dplyr::pull(.data$covariateId)
  }

  totalRows <- dplyr::n_distinct(plpData$labels$rowId)

  # Count the number of observed and missing values for each selected covariate
  tibble::tibble(covariateId = sort(unique(covariateIds))) %>%
    dplyr::left_join(
      covariates %>%
        dplyr::distinct(.data$rowId, .data$covariateId) %>%
        dplyr::count(.data$covariateId, name = "nObserved"),
      by = "covariateId"
    ) %>%
    # Return the missingness summary as a table
    dplyr::mutate(
      nObserved = dplyr::coalesce(.data$nObserved, 0L),
      nMissing = totalRows - .data$nObserved,
      missingFraction = .data$nMissing / totalRows
    ) %>%
    dplyr::arrange(dplyr::desc(.data$nMissing), .data$covariateId)
}

#' @title Report the change in missingenss
#' @description
#' This function compares the missingness before and after simulation
#' @param beforePLPData The original PLP data object
#' @param afterPLPData The PLP data object after missingness has been simulated
#' @param covariateIds Optional covariate ids to compare
#' @return A data frame describing the change in missingness per covariate
#' @export
reportMissingFraction <- function(beforePLPData, afterPLPData, covariateIds = NULL) {
  # Select the covariate data
  before <- beforePLPData$covariateData$covariates %>%
    collectIfNeeded() %>%
    dplyr::select(.data$rowId, .data$covariateId) %>%
    dplyr::distinct()

  after <- afterPLPData$covariateData$covariates %>%
    collectIfNeeded() %>%
    dplyr::select(.data$rowId, .data$covariateId) %>%
    dplyr::distinct()

  # Filter on given covariate ids if supplied
  if (!is.null(covariateIds)) {
    before <- before %>% dplyr::filter(.data$covariateId %in% covariateIds)
    after <- after %>% dplyr::filter(.data$covariateId %in% covariateIds)
  } else {
    # Else filter on all unique covariate ids
    covariateIds <- sort(unique(before$covariateId))
  }

  # Count number of observed values per covariate id
  tibble::tibble(covariateId = sort(unique(covariateIds))) %>%
    dplyr::left_join(
      before %>% dplyr::count(.data$covariateId, name = "nObservedBefore"),
      by = "covariateId"
    ) %>%
    dplyr::left_join(
      after %>% dplyr::count(.data$covariateId, name = "nObservedAfter"),
      by = "covariateId"
    ) %>%
    # Quantify the change in missingness
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
    # Return the summary as a table
    dplyr::arrange(dplyr::desc(.data$droppedFraction), .data$covariateId)
}

#' @title Compute softmax probabilities
#' @description
#' This function converts a numeric score vector into normalized probabilities
#' that sum to one
#' @param x A numeric vector
#' @return A numeric vector of softmax probabilities
softmax <- function(x) {
  # Stabilize the scores and normalize them to probabilities
  shifted <- x - max(x)
  expShifted <- exp(shifted)
  expShifted / sum(expShifted)
}
