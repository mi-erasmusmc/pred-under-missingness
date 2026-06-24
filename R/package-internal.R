collectIfNeeded <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }

  if (inherits(x, "data.frame")) {
    return(x)
  }

  tryCatch(
    dplyr::collect(x),
    error = function(e) as.data.frame(x)
  )
}

defaultRequiredCovariates <- function(targetCovariateId, causeCovariateIds = NULL) {
  sort(unique(stats::na.omit(c(targetCovariateId, causeCovariateIds))))
}
