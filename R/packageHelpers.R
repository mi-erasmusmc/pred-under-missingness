#' @title Materialize an object when possible
#' @description
#' This function calls `dplyr::collect()` on database objects when possible, otherwise
#' it returns the input unchanged
#' @param x An object that may require collection into memory
#' @return The collected object when collection is supported, otherwise the original input
#' @export
collectIfNeeded <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }

  # Leave ordinary data frames unchanged
  if (inherits(x, "data.frame")) {
    return(x)
  }

  # Try to collect the object, if it fails return the original data frame
  tryCatch(
    dplyr::collect(x),
    error = function(e) as.data.frame(x)
  )
}

#' @title Build the default required covariate set
#' @description This function combines target and cause covariate identifiers into one sorted unique set, removing missing values
#' @param targetCovariateId The target covariate id or ids
#' @param causeCovariateIds Optional cause covariate ids
#' @return A sorted vector of unique required covariate ids
#' @export
defaultRequiredCovariates <- function(targetCovariateId, causeCovariateIds = NULL) {
  # Combine target and cause covariates and drop missing values and duplicates
  sort(unique(stats::na.omit(c(targetCovariateId, causeCovariateIds))))
}
