#' @title Build and execute a cohort definition from JSON
#' @description This function reads a cohort definition from a JSON file, converts it to a CirceR
#' cohort expression, renders and translates the corresponding SQL for the target database, and executes it using the supplied connection
#' @param connection A `DatabaseConnector` connection object
#' @param jsonPath The path to the cohort-definition JSON file
#' @param cohortId The cohort definition id to assign in the target cohort table
#' @param cdmSchema The CDM database schema
#' @param vocabularySchema The vocabulary database schema
#' @param cohortSchema The schema where the cohort table will be written
#' @param cohortTable The name of the target cohort table
#' @param targetDialect The SQL dialect used by the target database
#' @return Invisibly returns the result of `DatabaseConnector::executeSql()`
#' @export
buildAndExecuteCohort <- function(connection,
                                  jsonPath,
                                  cohortId,
                                  cdmSchema,
                                  vocabularySchema,
                                  cohortSchema,
                                  cohortTable,
                                  targetDialect = "postgresql") {
  # Read the JSON cohort definition
  cohortJSON <- readr::read_file(jsonPath)
  cohortExpression <- CirceR::cohortExpressionFromJson(cohortJSON)

  # Build the cohort SQL
  cohortSQL <- CirceR::buildCohortQuery(
    expression = cohortExpression,
    options = CirceR::createGenerateOptions()
  )

  # Render the SQL by inserting the study-specific schemas, table name, and cohort id
  renderedSQL <- SqlRender::render(
    cohortSQL,
    cdm_database_schema = cdmSchema,
    vocabulary_database_schema = vocabularySchema,
    target_database_schema = cohortSchema,
    target_cohort_table = cohortTable,
    target_cohort_id = cohortId
  )

  # Translate the sql to dialect used by the connected database
  translatedSQL <- SqlRender::translate(
    renderedSQL,
    targetDialect = targetDialect
  )

  # Execute the translated cohort sql against the target database connection
  DatabaseConnector::executeSql(connection, translatedSQL)
}

#' @title Rebuild a PLP data object from its components
#' @description This helper constructs a new `plpData` object from supplied labels,
#' optional folds, and covariate tables while preserving the backing structure, class,
#' and metadata of a template PLP data object
#' @param labels A labels table fro the rebuilt PLP data object
#' @param folds An optional folds table
#' @param covariates A covariates table
#' @param covariateRef A covariate reference table
#' @param analysisRef An analysis reference table
#' @param templatePLPData A template `plpData` object whose structure and metadata should be preserved
#' @return A rebuilt `plpData` object
#' @export
plpDataHelper <- function(labels,
                          folds = NULL,
                          covariates,
                          covariateRef,
                          analysisRef,
                          templatePLPData) {
  # Copy the template covariateData container so the rebuilt object keeps the
  # same backing structure as the original PLP data
  covariateData <- Andromeda::copyAndromeda(templatePLPData$covariateData)

  # Replace the original covariate tables with the supplied versions
  covariateData$covariates <- covariates
  covariateData$covariateRef <- covariateRef
  covariateData$analysisRef <- analysisRef

  # Preserve the optional time reference from the template when it is present
  if (!is.null(templatePLPData$covariateData$timeRef)) {
    covariateData$timeRef <- templatePLPData$covariateData$timeRef
  }

  # Restore the original covariate data class and meta data
  class(covariateData) <- class(templatePLPData$covariateData)
  attr(covariateData, "metaData") <- attr(templatePLPData$covariateData, "metaData")

  # Build the outer PLP data object from the supplied labels and folds
  out <- list(labels = as.data.frame(labels))

  if (!is.null(folds)) {
    out$folds <- as.data.frame(folds)
  }

  out$covariateData <- covariateData

  # Restore the PLP class and top-level metadata from the template object
  class(out) <- "plpData"
  attr(out, "metaData") <- attr(templatePLPData, "metaData")
  out
}

#' @title Build a subset of the population given selected row ids
#' @description This helper filters PLP population tbale to the requested row ids while
#' preserving the original metadata
#' @param population A population table
#' @param selectedRowIds The row ids to retain
#' @return A filtered population table
#' @export
populationSubset <- function(population, selectedRowIds) {
  # Keep only the rows of the population corresponding to the row ids
  filteredPopulation <- as.data.frame(population) %>%
    dplyr::filter(.data$rowId %in% selectedRowIds)

  # Append original metadata to the subset
  attr(filteredPopulation, "metaData") <- attr(population, "metaData")
  filteredPopulation
}

#' @title Retain selected measurement covariates
#' @description This function filters the measurement covariates to a requested st of measurement concept
#' ids while leaving other types of covariates untouched
#' @param covariateRef A covariate reference table
#' @param covariateValues A covariates values table
#' @param analysisRef An analysis reference table
#' @param measurementConceptIds The measurement concept ids to retain
#' @return A list containing the filtered `covariateRef`, `covariateValues`, `analysisRef`, and
#' the retained `measurementRef`
#' @export
keepMeasurements <- function(covariateRef,
                             covariateValues,
                             analysisRef,
                             measurementConceptIds) {
  # If no measurment filtering is requested, return the inputs unchanged and provide
  # an empty measurement reference table
  if (length(measurementConceptIds) == 0) {
    return(list(
      covariateRef = covariateRef,
      covariateValues = covariateValues,
      analysisRef = analysisRef,
      measurementRef = covariateRef[0, , drop = FALSE]
    ))
  }

  # keep only the measurements whose concept ids are explicitly selected
  allowedMeasurementRef <- covariateRef %>%
    dplyr::filter(
      !is.na(.data$conceptId),
      .data$conceptId %in% measurementConceptIds
    ) %>%
    dplyr::distinct(.data$covariateId, .keep_all = TRUE)

  # Identify the full measurement analyses corresponding to the retained measurement covariates
  measurementAnalysisIds <- allowedMeasurementRef %>%
    dplyr::pull(.data$analysisId) %>%
    unique()

  measurementCovariateIds <- covariateRef %>%
    dplyr::filter(.data$analysisId %in% measurementAnalysisIds) %>%
    dplyr::pull(.data$covariateId) %>%
    unique()

  # Drop non selected covaraites from the measurement analyses, but keep all covariates that do not belong to the measurement
  # analyses being filtered
  retainedCovariateIds <- union(
    setdiff(covariateRef$covariateId, measurementCovariateIds),
    allowedMeasurementRef$covariateId
  )

  filteredCovariateRef <- covariateRef %>%
    dplyr::filter(.data$covariateId %in% retainedCovariateIds)

  filteredAnalysisRef <- analysisRef %>%
    dplyr::semi_join(filteredCovariateRef, by = "analysisId")

  filteredCovariateValues <- covariateValues %>%
    dplyr::filter(.data$covariateId %in% retainedCovariateIds)

  list(
    covariateRef = filteredCovariateRef,
    covariateValues = filteredCovariateValues,
    analysisRef = filteredAnalysisRef,
    measurementRef = allowedMeasurementRef
  )
}

#' @title Build a PLP data object for a subset population
#' @description This function constructs a `plpData` object for a selected set of row ids
#' by subsetting the popultaion, covariates and reference tables from a larger PLP data source
#' @param selectedRowids The row ids to include in the object
#' @param population A population table
#' @param covariateValues A covariate-values table
#' @param covariateRef A covariate reference table
#' @param analysisRef An analysis reference table
#' @param measurementConceptIds Optional measurement concept ids to retain
#' @param templatePLPData A template `plpData` object whose structure and metadata should
#' be preserved
#' @param includeFolds Logical variable indicating whether default folds table should added
#' to the output
#' @return A list containing the rebuilt `plpData` object and the filtered population table
#' @export
buildPopulationPLPData <- function(selectedRowIds,
                                   population,
                                   covariateValues,
                                   covariateRef,
                                   analysisRef,
                                   measurementConceptIds = integer(),
                                   templatePLPData,
                                   includeFolds = TRUE) {
  selectedRowIds <- sort(unique(selectedRowIds))

  # materialize the objects so the filtering works on in memory tables
  population <- as.data.frame(population)
  covariateValues <- collectIfNeeded(covariateValues)
  covariateRef <- collectIfNeeded(covariateRef)
  analysisRef <- collectIfNeeded(analysisRef)

  # Restrict the study population to the requested rows
  finalPopulation <- populationSubset(population, selectedRowIds)

  # Keep only the covariate values observed in the selected population
  analysisCovariates <- covariateValues %>%
    dplyr::filter(.data$rowId %in% selectedRowIds)

  # Restrict the reference tables to the covariates that are actually present after
  # row-level subsetting
  selectedCovariateIds <- analysisCovariates %>%
    dplyr::distinct(.data$covariateId) %>%
    dplyr::pull(.data$covariateId)

  selectedCovariateRef <- covariateRef %>%
    dplyr::filter(.data$covariateId %in% selectedCovariateIds)

  selectedAnalysisRef <- analysisRef %>%
    dplyr::semi_join(selectedCovariateRef, by = "analysisId")

  # If given, filter measurement covariates to the requested measurements
  if (length(measurementConceptIds) > 0) {
    filteredSet <- keepMeasurements(
      covariateRef = selectedCovariateRef,
      covariateValues = analysisCovariates,
      analysisRef = selectedAnalysisRef,
      measurementConceptIds = measurementConceptIds
    )

    analysisCovariates <- filteredSet$covariateValues
    selectedCovariateRef <- filteredSet$covariateRef
    selectedAnalysisRef <- filteredSet$analysisRef
  }

  # Add default folds table if folds should be included
  folds <- NULL
  if (isTRUE(includeFolds)) {
    folds <- data.frame(
      rowId = selectedRowIds,
      index = 1L
    )
  }

  # Rebuild the plp dat aobject and return it with the fitlered population
  list(
    plpData = plpDataHelper(
      labels = finalPopulation,
      folds = folds,
      covariates = analysisCovariates,
      covariateRef = selectedCovariateRef,
      analysisRef = selectedAnalysisRef,
      templatePLPData = templatePLPData
    ),
    population = finalPopulation
  )
}
