createPLPImputationMethods <- function(includeMissingIndicator = c(FALSE, TRUE),
                                       sklearnRandomState = 432) {
  methods <- list()

  for (addIndicator in includeMissingIndicator) {
    suffix <- if (addIndicator) "_withIndicator" else "_noIndicator"

    methods[[paste0("simpleMean", suffix)]] <-
      PatientLevelPrediction::createSimpleImputer(
        method = "mean",
        missingThreshold = 0.95,
        addMissingIndicator = addIndicator
      )

    methods[[paste0("simpleMedian", suffix)]] <-
      PatientLevelPrediction::createSimpleImputer(
        method = "median",
        missingThreshold = 0.95,
        addMissingIndicator = addIndicator
      )

    methods[[paste0("iterativePMM", suffix)]] <-
      PatientLevelPrediction::createIterativeImputer(
        missingThreshold = 0.95,
        method = "pmm",
        methodSettings = list(
          pmm = list(
            k = 5,
            iterations = 5,
            alpha = 1
          )
        ),
        addMissingIndicator = addIndicator
      )

    if ("createSklearnIterativeImputer" %in% getNamespaceExports("PatientLevelPrediction")) {
      methods[[paste0("sklearnIterative", suffix)]] <-
        PatientLevelPrediction::createSklearnIterativeImputer(
          missingThreshold = 0.95,
          methodSettings = list(
            maxIter = 10,
            initialStrategy = "mean",
            randomState = sklearnRandomState
          ),
          addMissingIndicator = addIndicator
        )
    }
  }

  methods
}

safeCopyCovariateData <- function(covariateData) {
  out <- Andromeda::andromeda(
    covariates = covariateData$covariates %>% collectIfNeeded(),
    covariateRef = covariateData$covariateRef %>% collectIfNeeded(),
    analysisRef = covariateData$analysisRef %>% collectIfNeeded()
  )

  if (!is.null(covariateData$timeRef)) {
    out$timeRef <- covariateData$timeRef %>% collectIfNeeded()
  }

  class(out) <- class(covariateData)
  attr(out, "metaData") <- attr(covariateData, "metaData")
  out
}

safeCopyPlpData <- function(x) {
  out <- list(labels = as.data.frame(x$labels))

  if (!is.null(x$folds)) {
    out$folds <- as.data.frame(x$folds)
  }

  out$covariateData <- safeCopyCovariateData(x$covariateData)

  class(out) <- class(x)
  attr(out, "metaData") <- attr(x, "metaData")
  out
}

getPLPImputer <- function(funName) {
  baseFun <- utils::getFromNamespace(funName, "PatientLevelPrediction")

  function(trainData, featureEngineeringSettings, done = FALSE) {
    baseFun(
      trainData = safeCopyPlpData(trainData),
      featureEngineeringSettings = featureEngineeringSettings,
      done = done
    )
  }
}

getPLPImputationMetaKey <- function(funName) {
  switch(
    funName,
    simpleImpute = "simpleImputer",
    iterativeImpute = "iterativeImputer",
    sklearnIterativeImpute = "sklearnIterativeImputer",
    stop("Unknown PLP imputation function: ", funName)
  )
}

getFittedPLPImputerSettings <- function(imputedTrainData, imputerSettings) {
  funName <- attr(imputerSettings, "fun")
  metaKey <- getPLPImputationMetaKey(funName)

  attr(imputedTrainData$covariateData, "metaData")$
    featureEngineering[[metaKey]]$settings$featureEngineeringSettings
}

runPlpImputers <- function(trainData,
                           testData = NULL,
                           imputationMethods,
                           imputationOverview = NULL) {
  if (missing(imputationMethods) || length(imputationMethods) == 0) {
    stop("Provide at least one imputation name.")
  }

  if (is.null(imputationOverview)) {
    imputationOverview <- createPLPImputationMethods()
  }

  unknownMethods <- setdiff(imputationMethods, names(imputationOverview))
  if (length(unknownMethods) > 0) {
    stop("Unknown imputation method(s): ", paste(unknownMethods, collapse = ", "))
  }

  results <- lapply(imputationMethods, function(methodName) {
    imputerSettings <- imputationOverview[[methodName]]
    funName <- attr(imputerSettings, "fun")
    imputerFun <- getPLPImputer(funName)

    trainImputed <- imputerFun(
      trainData = trainData,
      featureEngineeringSettings = imputerSettings,
      done = FALSE
    )

    fittedSettings <- getFittedPLPImputerSettings(
      imputedTrainData = trainImputed,
      imputerSettings = imputerSettings
    )

    testImputed <- NULL
    if (!is.null(testData)) {
      testImputed <- imputerFun(
        trainData = testData,
        featureEngineeringSettings = fittedSettings,
        done = TRUE
      )
    }

    list(
      methodName = methodName,
      functionName = funName,
      originalSettings = imputerSettings,
      fittedSettings = fittedSettings,
      trainImputed = trainImputed,
      testImputed = testImputed
    )
  })

  names(results) <- imputationMethods
  results
}

subsetPlpDataRows <- function(plpData, rowIds) {
  rowIds <- sort(unique(rowIds))

  labels <- as.data.frame(plpData$labels) %>%
    dplyr::filter(.data$rowId %in% rowIds)

  covariates <- plpData$covariateData$covariates %>%
    dplyr::filter(.data$rowId %in% rowIds) %>%
    collectIfNeeded()

  covariateRef <- plpData$covariateData$covariateRef %>% collectIfNeeded()
  analysisRef <- plpData$covariateData$analysisRef %>% collectIfNeeded()

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

completeCasePlp <- function(trainData, testData, targetCovariateIds) {
  targetCovariateIds <- sort(unique(targetCovariateIds))

  getObservedRows <- function(plpData, covariateIds) {
    plpData$covariateData$covariates %>%
      dplyr::filter(.data$covariateId %in% covariateIds) %>%
      dplyr::distinct(.data$rowId, .data$covariateId) %>%
      collectIfNeeded() %>%
      dplyr::count(.data$rowId, name = "nObservedTargets") %>%
      dplyr::filter(.data$nObservedTargets == length(covariateIds)) %>%
      dplyr::pull(.data$rowId)
  }

  trainRows <- getObservedRows(trainData, targetCovariateIds)
  testRows <- getObservedRows(testData, targetCovariateIds)

  list(
    methodName = "completeCase",
    trainImputed = subsetPlpDataRows(trainData, trainRows),
    testImputed = subsetPlpDataRows(testData, testRows)
  )
}
