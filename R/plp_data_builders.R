buildAndExecuteCohort <- function(connection,
                                  jsonPath,
                                  cohortId,
                                  cdmSchema,
                                  vocabularySchema,
                                  cohortSchema,
                                  cohortTable,
                                  targetDialect = "postgresql") {
  cohortJSON <- readr::read_file(jsonPath)
  cohortExpression <- CirceR::cohortExpressionFromJson(cohortJSON)

  cohortSQL <- CirceR::buildCohortQuery(
    expression = cohortExpression,
    options = CirceR::createGenerateOptions()
  )

  renderedSQL <- SqlRender::render(
    cohortSQL,
    cdm_database_schema = cdmSchema,
    vocabulary_database_schema = vocabularySchema,
    target_database_schema = cohortSchema,
    target_cohort_table = cohortTable,
    target_cohort_id = cohortId
  )

  translatedSQL <- SqlRender::translate(
    renderedSQL,
    targetDialect = targetDialect
  )

  DatabaseConnector::executeSql(connection, translatedSQL)
}

plpDataHelper <- function(labels,
                          folds = NULL,
                          covariates,
                          covariateRef,
                          analysisRef,
                          templatePLPData) {
  covariateData <- Andromeda::copyAndromeda(templatePLPData$covariateData)

  covariateData$covariates <- covariates
  covariateData$covariateRef <- covariateRef
  covariateData$analysisRef <- analysisRef

  if (!is.null(templatePLPData$covariateData$timeRef)) {
    covariateData$timeRef <- templatePLPData$covariateData$timeRef
  }

  class(covariateData) <- class(templatePLPData$covariateData)
  attr(covariateData, "metaData") <- attr(templatePLPData$covariateData, "metaData")

  out <- list(labels = as.data.frame(labels))

  if (!is.null(folds)) {
    out$folds <- as.data.frame(folds)
  }

  out$covariateData <- covariateData

  class(out) <- "plpData"
  attr(out, "metaData") <- attr(templatePLPData, "metaData")
  out
}

populationSubset <- function(population, selectedRowIds) {
  filteredPopulation <- as.data.frame(population) %>%
    dplyr::filter(.data$rowId %in% selectedRowIds)

  attr(filteredPopulation, "metaData") <- attr(population, "metaData")
  filteredPopulation
}

keepMeasurements <- function(covariateRef,
                             covariateValues,
                             analysisRef,
                             measurementConceptIds) {
  if (length(measurementConceptIds) == 0) {
    return(list(
      covariateRef = covariateRef,
      covariateValues = covariateValues,
      analysisRef = analysisRef,
      measurementRef = covariateRef[0, , drop = FALSE]
    ))
  }

  allowedMeasurementRef <- covariateRef %>%
    dplyr::filter(
      !is.na(.data$conceptId),
      .data$conceptId %in% measurementConceptIds
    ) %>%
    dplyr::distinct(.data$covariateId, .keep_all = TRUE)

  measurementAnalysisIds <- allowedMeasurementRef %>%
    dplyr::pull(.data$analysisId) %>%
    unique()

  measurementCovariateIds <- covariateRef %>%
    dplyr::filter(.data$analysisId %in% measurementAnalysisIds) %>%
    dplyr::pull(.data$covariateId) %>%
    unique()

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

buildPopulationPLPData <- function(selectedRowIds,
                                   population,
                                   covariateValues,
                                   covariateRef,
                                   analysisRef,
                                   measurementConceptIds = integer(),
                                   templatePLPData,
                                   includeFolds = TRUE) {
  selectedRowIds <- sort(unique(selectedRowIds))

  population <- as.data.frame(population)
  covariateValues <- collectIfNeeded(covariateValues)
  covariateRef <- collectIfNeeded(covariateRef)
  analysisRef <- collectIfNeeded(analysisRef)

  finalPopulation <- populationSubset(population, selectedRowIds)

  analysisCovariates <- covariateValues %>%
    dplyr::filter(.data$rowId %in% selectedRowIds)

  selectedCovariateIds <- analysisCovariates %>%
    dplyr::distinct(.data$covariateId) %>%
    dplyr::pull(.data$covariateId)

  selectedCovariateRef <- covariateRef %>%
    dplyr::filter(.data$covariateId %in% selectedCovariateIds)

  selectedAnalysisRef <- analysisRef %>%
    dplyr::semi_join(selectedCovariateRef, by = "analysisId")

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

  folds <- NULL
  if (isTRUE(includeFolds)) {
    folds <- data.frame(
      rowId = selectedRowIds,
      index = 1L
    )
  }

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
