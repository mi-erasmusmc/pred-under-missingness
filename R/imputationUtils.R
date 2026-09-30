#' @title Create a safe copy of PLP covariate data
#' @description This function copies a PLP `covariateData` object into a
#' new in-memory object so imputation and modeling steps can modify it without mutating the original input
#' @param covariateData A PLP covariate data object
#' @return A copied covariate data PLP object with the original class and metadata preserved
#' @export
# TODO Check if this needs additional changes such that parallelization of XGBoost works
safeCopyCovariateData <- function(covariateData) {
  # Create a new Andromeda object
  out <- Andromeda::andromeda(
    covariates = covariateData$covariates %>% collectIfNeeded(),
    covariateRef = covariateData$covariateRef %>% collectIfNeeded(),
    analysisRef = covariateData$analysisRef %>% collectIfNeeded()
  )

  # Preserve the optional time reference when it is present
  if (!is.null(covariateData$timeRef)) {
    out$timeRef <- covariateData$timeRef %>% collectIfNeeded()
  }

  # Restore the original class and metadata on the copied object
  class(out) <- class(covariateData)
  attr(out, "metaData") <- attr(covariateData, "metaData")
  out
}

#' @title Create a safe copy of PLP data
#' @description
#' This function copies a PLP data object into a new in-memory object
#' so workflows can modify it without changing the original input
#' @param x A PLP data object
#' @return A copied PLP data object with the original class and metadata preserved
#' @export
safeCopyPlpData <- function(x) {
  # Copy the labels and optional folds into plain in-memory data frames
  out <- list(labels = as.data.frame(x$labels))

  if (!is.null(x$folds)) {
    out$folds <- as.data.frame(x$folds)
  }

  # Copy the covariate data into a separate object
  out$covariateData <- safeCopyCovariateData(x$covariateData)

  # Restore the originaal class and metadata on the copied PLP object
  class(out) <- class(x)
  attr(out, "metaData") <- attr(x, "metaData")
  out
}

#' @title Subset a PLP data object by row identifiers
#' @description
#' This function returns a PLP data object containing only the
#' requested row identfiers while preserving the original associated labels,
#' folds and covariate reference tables
#' @param plpData A PLP data object
#' @param rowIds The row identifiers to retain
#' @return A PLP data object restricted to the requested rows
#' @export
subsetPlpDataRows <- function(plpData, rowIds) {
  # Standardize the requested row identifiers before subsetting
  rowIds <- sort(unique(rowIds))

  # Subset the labels and observed covariates to the requested rows
  labels <- as.data.frame(plpData$labels) %>%
    dplyr::filter(.data$rowId %in% rowIds)

  covariates <- plpData$covariateData$covariates %>%
    dplyr::filter(.data$rowId %in% rowIds) %>%
    collectIfNeeded()

  # Keep the covariaet and analysis reference tables unchanged
  covariateRef <- plpData$covariateData$covariateRef %>% collectIfNeeded()
  analysisRef <- plpData$covariateData$analysisRef %>% collectIfNeeded()

  # Subset the folds table when fold information is available
  folds <- NULL
  if (!is.null(plpData$folds)) {
    folds <- as.data.frame(plpData$folds) %>%
      dplyr::filter(.data$rowId %in% rowIds)
  }

  plpDataHelper(
    labels = labels,
    folds = folds,
    covariates = covariates,
    covariateRef = covariateRef,
    analysisRef = analysisRef,
    templatePLPData = plpData
  )
}
