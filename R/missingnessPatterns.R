#' @title Normalize target missing fractions
#' @description
#' This function standardizes the requested covariate level missing fractions into a names numeric vector
#' aligned with the supplied covariate identifiers
#' @param missingFractions A named vector of missingness fractions
#' @param covariateIds The covariate identifiers to align to
#' @return A numeric vector of missingness fractions
normalizeMissingFractionTargets <- function(missingFractions, covariateIds) {
  if (is.null(missingFractions)) {
    return(NULL)
  }

  if (!is.numeric(missingFractions)) {
    stop("missingFractions must be numeric.", call. = FALSE)
  }

  if (is.null(names(missingFractions))) {
    if (length(missingFractions) != length(covariateIds)) {
      stop(
        "Unnamed missingFractions must have the same length as covariateIds.",
        call. = FALSE
      )
    }

    names(missingFractions) <- as.character(covariateIds)
  } else {
    missingNames <- setdiff(as.character(covariateIds), names(missingFractions))
    if (length(missingNames) > 0) {
      stop(
        "missingFractions is missing values for covariateIds: ",
        paste(missingNames, collapse = ", "),
        call. = FALSE
      )
    }
    # Align the requested fractions to the target covariate ids
    missingFractions <- missingFractions[as.character(covariateIds)]
  }

  # Validate that all fractions are finite and lie between 0 and 1
  if (any(!is.finite(missingFractions)) || any(missingFractions < 0) || any(missingFractions > 1)) {
    stop("missingFractions must be finite probabilities between 0 and 1.", call. = FALSE)
  }

  missingFractions
}

#' @title Normalize pattern fractions
#' @description
#' This function standardizes the requested pattern-frequency vector into a numeric vector aligned with
#' the supplied number of patterns
#' @param patternFractions A named vector of pattern fractions
#' @param nrOfPatterns The number of missingness patterns
#' @return A numeric vector of pattern fractions
normalizePatternFractionTargets <- function(patternFractions, nrOfPatterns) {
  # Return null when no pattern fractions are supplied
  if (is.null(patternFractions)) {
    return(rep(NA_real_, nrOfPatterns))
  }

  if (!is.numeric(patternFractions)) {
    stop("patternFractions must be numeric.", .call = FALSE)
  }

  if (is.null(names(patternFractions))) {
    if (length(patternFractions) != nrOfPatterns) {
      stop("Unnamed patternFractions must have the same length as patterns.",
           call. = FALSE)
    }
    output <- patternFractions
  } else {
    output <- rep(NA_real_, nrOfPatterns)
    validNames <- pate0("pattern", seq_len(nrOfPatterns))
    unknownNames <- setdiff(names(patternFractions), validNames)

    if (length(unknownNames) > 0) {
      stop(
        "Names patternFractions entries must use names like pattern1, pattern2, ...",
        call. = FALSE
      )
    }

    # Align the requested fractions to the patterns
    output[match(names(patternFractions), validNames)] <- patternFractions

  }

  # Validate that all fractions are finite and nonnegative
  if (any(!is.na(output) & (!is.finite(output) | output < 0 | output > 1))) {
    stop("patternFractions must be finite probabilities between 0 and 1.")
  }
  output
}

#' @title Build a pattern missingness matrix
#' @description
#' This function converts a list of pattern specifications into a binary matrix
#' indicating which covariates are set missing by each pattern
#' @param patterns A list of missingness patterns
#' @param covariateIds The covariate ids to represent in the matrix
#' @return A binary matrix with covariates in rows and patterns in columns
buildPatternMissingnessMatrix <- function(patterns, covariateIds) {
  # Initialize a zero matrix
  patternMatrix <- matrix(
    0,
    nrow = length(covariateIds),
    ncol = length(patterns),
    dimnames = list(as.character(covariateIds), paste0("pattern", seq_along(patterns)))
  )

  for (i in seq_along(patterns)) {
    targetCovariateIds <- unique(patterns[[i]]$targetCovariateIds)

    if (is.null(targetCovariateIds)) {
      stop("Pattern ", i, " must have targetCovariateIds defined.", call. = FALSE)
    }

    # Mark each covariate pattern combination that becomes missing under the pattern.
    matchedIds <- intersect(as.character(targetCovariateIds), as.character(covariateIds))
    patternMatrix[matchedIds, i] <- 1
  }
  # Return binary pattern incidence matrix
  patternMatrix
}

#' @title Allocate integer row counts to patterns
#' @description
#' This function converts target pattern shares into integer row counts of which the total
#' matches the requested number of rows
#' @param nRows The total number of rows to allocate
#' @param freq The target pattern shares
#' @return An integer vector of allocated row counts
allocatePatternCounts <- function(nRows, freq) {
  # Scale the target shares to the requested number of rows
  rawCounts <- nRows * freq
  # Floor the number such that it never exceeds the number of rows
  counts <- floor(rawCounts)
  # Determine how many rows still need to be assigned
  remainder <- nRows - sum(counts)

  # Assign the remaining rows to patterns whose counts had the largest fractions
  # So the final allocation stays as close as possible to the requested shares
  if (remainder > 0) {
    addOrder <- order(rawCounts - counts, decreasing = TRUE)
    counts[addOrder[seq_len(remainder)]] <- counts[addOrder[seq_len(remainder)]] + 1L
  }

  # Return integer pattern counts that sum up to exactly nRows
  as.integer(counts)
}

#' @title Compute pattern selection weights
#' @description This function computes row level selection weights for one missingness pattern based on its mechanism
#' and cause covariates
#' @param covariatesDf A long format covariate table
#' @param rowIds The candidate row identifiers
#' @param pattern A missingness pattern specification
#' @return A data frame containing the row identifiers and their selection weights
computePatternSelectionWeights <- function(covariatesDf, rowIds, pattern) {
  # Return an empty table when there are no rows to assign
  if (length(rowIds) == 0) {
    return(tibble::tibble(
      rowId = integer(),
      weight = numeric()
    ))
  }

  # Extract the pattern settings that control how weights are built
  causeCovariateIds <- pattern$causeCovariateIds
  causeWeights <- pattern$causeWeights
  type <- if (is.null(pattern$type)) "RIGHT" else toupper(pattern$type)
  gamma2 <- if (is.null(pattern$gamma2)) 1 else pattern$gamma2
  standardizeScore <- if (is.null(pattern$standardizeScore)) TRUE else pattern$standardizeScore

  # If the pattern has no cause covariates, row assignment is effectively unweighted
  # and every candidate row gets the same selection weight
  if (is.null(causeCovariateIds)) {
    return(tibble::tibble(
      rowId = rowIds,
      weight = rep(1, length(rowIds))
    ))
  }
  # Use euqal covariate weights if no weights are supplied by the pattern
  if (is.null(causeWeights)) {
    causeWeights <- rep(1, length(causeCovariateIds))
  }

  # Calculate the row level score
  causeData <- calcScore(
    covariates = covariatesDf,
    rowIds = rowIds,
    causeCovariateIds = causeCovariateIds,
    weights = causeWeights,
    standardize = standardizeScore
  ) %>%
    dplyr::mutate(
      # Transform the score such that the requested pattern shape is used
      z = transformScore(.data$score, type),
      # Convert the transfermed score into selection probabilities
      weight = logitInverse(gamma2 * .data$z)
    ) %>%
    dplyr::select("rowId", "weight")

  # if negative or infinite weight --> set equal to zero
  causeData$weight[!is.finite(causeData$weight) | causeData$weight < 0] <- 0

  # If computed weights sum to zero, fall back to equal weights
  if (sum(causeData$weight) <= 0) {
    causeData$weight <- rep(1, nrow(causeData))
  }

  causeData
}

#' @title Select rows for one missingness pattern
#' @description This functin samples rows for one missingness pattern from the
#' rows that are still available for the assignment. Sampling is performed without
#' replacement using the pattern specific row weights
#' @param remainingRowIds The row identifiers still available for assignment
#' @param targetCount The number of rows to assign to the pattern
#' @param covariatesDf A long format covariate table
#' @param pattern A missingness pattern specification
#' @return A vector of selected row identifiers
selectPatternRows <- function(remainingRowIds, targetCount, covariatesDf, pattern) {
  # Return an empty selection when there are no rows to assign
  if (targetCount == 0L) {
    return(integer())
  }

  # Stop early if the target count exceeds the number of remaining rows
  if (targetCount > length(remainingRowIds)) {
    stop("Pattern target count exceeds the number of remaining rows.", call.=FALSE)
  }

  # Compute one selection weight per row still available for the assignment
  weights <- computePatternSelectionWeights(
    covariatesDf = covariatesDf,
    rowIds = remainingRowIds,
    pattern = pattern
  )

  # Stop if there are insufficient weights
  if (nrow(weights) != length(remainingRowIds)) {
    stop("Unable to compute selection weights for all remaining rows.", call. = FALSE)
  }
  # Name the weights by the row id
  prob <- stats::setNames(weights$weight, weights$rowId)

  # Sample rows without replacement so a row can only be assigned to one pattern
  sample(
    x = remainingRowIds,
    size = targetCount,
    replace = FALSE,
    prob = prob[as.character(remainingRowIds)]
  )
}

#' @title Solve row-pattern frequencies for a multivariate missingness design
#'
#' @description Given a list of row-level missingness patterns, this function solves for the `freq`
#' vector used by `simulateMultivariateMissingness()` so that the resulting variable-specific
#' missingness fractions match the required margin whenever that is feasible
#'
#' @param patterns A list of pattern specifications, each with `targetCovariateIds`.
#' @param covariateIds Optional vector of covariateIds to calibrate. Defaults to the union of all
#' `targetCovariateIds` in `patterns`.
#' @param missingnessFractions A vector with the desired variable-level missingness fractions.
#' @param patternFractions A vector with fixed pattern-level frequencies. Use `NA` for patterns
#' that should be solved. Named input should use names like `pattern1`, `pattern2`, etc.
#' @param tol Numerical tolerance used when checking feasibility.
#'
#' @return A list with the solved `freq` vector, the implied missingness fractions,
#' and the binary pattern matrix.
solvePatternFrequencies <- function(patterns, covariateIds = NULL, missingFractions = NULL, patternFractions = NULL, tol = 1e-8) {
  # Stop if no pattern is supplied at all
  if (!is.list(patterns) || length(patterns) == 0) {
    stop("Patterns must be a non-empty list.", call. = FALSE)
  }

  # If no covariate ids are supplied, solve over the union of all target covariates that appear in the
  # pattern specifications
  if (is.null(covariateIds)) {
    covariateIds <- unique(unlist(lapply(patterns, function(x) x$targetCovariateIds), use.names = FALSE))
  }

  covariateIds <- unique(covariateIds)

  # At least one covariate is needed to define the target missingness margins
  if (length(covariateIds) == 0) {
    stop("Provide at least one covariateId to calibrate.", call. = FALSE)
  }

  # Either missing fractions or pattern fractions should be supplied
  if (is.null(missingFractions) && is.null(patternFractions)) {
    stop("Provide missingFractions, patternFractions or both.", call. = FALSE)
  }

  # Represent each pattern as a binary matrix to see which covariates it makes missing
  patternMatrix <- buildPatternMissingnessMatrix(patterns, covariateIds)
  # Align the requested targets to the covariate and pattern order used
  missingFractions <- normalizeMissingFractionTargets(missingFractions, covariateIds)
  patternFractions <- normalizePatternFractionTargets(patternFractions, length(patterns))

  # Split the pattern frequencies into fixed entries and entries that still need to be solved
  knownIndex <- which(!is.na(patternFractions))
  unknownIndex <- which(is.na(patternFractions))
  knownTotal <- sum(patternFractions[knownIndex])

  # Known frequencies cannot sum to more than one
  if (knownTotal > 1 + tol) {
    stop("Known patternFractions cannot sum to more than 1.", call. = FALSE)
  }

  # Define the remaining mass that can still be distributed across free patterns
  remainingMass <- max(0, 1- knownTotal)

  if (!is.null(missingFractions)) {
    # Compute the missingness already forced by the fixed pattern frequencies
    fixedContribution <- if(length(knownIndex) == 0) {
      rep(0, length(covariateIds))
    } else {
      as.numeric(patternMatrix[, knownIndex, drop = FALSE] %*% patternFractions[knownIndex])
    }

    # Use the free patterns to derive the minimum and maximum missingness
    # that can still be achieved for each covariate
    unknownMatrix <- patternMatrix[, unknownIndex, drop = FALSE]
    minPossible <- if (length(unknownIndex) == 0) {
      fixedContribution
    } else {
      fixedContribution + remainingMass * apply(unknownMatrix, 1, min)
    }
    maxPossible <- if(length(unknownIndex) == 0) {
      fixedContribution
    } else {
      fixedContribution + remainingMass * apply(unknownMatrix, 1, max)
    }

    # Stop when missing fraction lies outside a feasible range
    infeasible <- missingFractions < (minPossible - tol) | missingFractions > (maxPossible + tol)
    if (any(infeasible)) {
      badIds <- names(missingFractions)[infeasible]
      details <- vapply(
        badIds,
        function(id) {
          index <- match(id, names(missingFractions))
          paste0(
            id,
            " requested =",
            formatC(missingFractions[[index]], digits = 4, format = "f"),
            ", feasible range =[",
            formatC(minPossible[[index]], digits = 4, format = "f"),
            ", ",
            formatC(maxPossible[[index]], digits = 4, format = "f"),
            "]"
          )
        },
        character(1)
      )

      stop(
        "No feasible pattern frequencies satisfy the requested missingness fractions.",
        paste(details, collapse = "; "),
        call. = FALSE
      )
    }
  }

  # If all pattern frequencies are already fixed, no optimization is needed
  if (length(unknownIndex) == 0) {
    freq <- patternFractions
  } else if (is.null(missingFractions)) {
    # Free pattern frequencies can only be solved if covariate targets are available
    stop(
      "missingFractions are required when some patternFractions are left unspecified.",
      call. = FALSE
    )
  } else {
    # Optimize the free pattern frequencies using softmax so they remain
    # nonnegative and automatically sum to the remaining probability mass
    objective <- function(par) {
      weights <- softmax(par)
      freq <- patternFractions
      freq[unknownIndex] <- remainingMass * weights
      residuals <- as.numeric(patternMatrix %*% freq) - missingFractions
      sum(residuals^2)
    }

    fit <- stats::optim(
      par = rep(0, length(unknownIndex)),
      fn = objective,
      method = "BFGS"
    )

    freq <- patternFractions
    freq[unknownIndex] <- remainingMass * softmax(fit$par)
  }

  # Translate the solved pattern frequencies back into covariate level
  # missingness fractions
  achievedMissingFractions <- as.numeric(patternMatrix %*% freq)
  names(achievedMissingFractions) <- rownames(patternMatrix)
  maxAbsError <- if (is.null(missingFractions)) {
    0
  } else {
    max(abs(achievedMissingFractions - missingFractions))
  }

  # Verify that the final solution matches the requested targets closely enough
  if (!is.null(missingFractions) && maxAbsError > tol) {
    stop(
      "Unable to match the requested missing fractions with tolerance. ",
      "Closest max absolute deviation = ",
      formatC(maxAbsError, digits = 4, format = "f"),
      call. = FALSE
    )
  }

  list(
    freq = stats::setNames(freq, colnames(patternMatrix)),
    achievedMissingFractions = achievedMissingFractions,
    requestedMissingFractions = requestedMissingFractions,
    patternMatrix = patternMatrix,
    maxAbsError = maxAbsError
  )
}
