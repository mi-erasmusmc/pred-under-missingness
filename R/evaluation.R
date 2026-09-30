#'
#' @title Extract evaluation tables from one scenario run
#' @description
#' This helper function extracts the standard PLP evaluation tables from a single evaluation object and appends simulation run, scenario,
#' imputation method and prediction model information so the tables can be combined across runs.
#' @param evaluation The evaluation object returned for one model run.
#' @param simulation The simulation run
#' @param mechanism The missingenss mechanism
#' @param ratio The missingness ratio
#' @param type The missingness score type
#' @param targetVariables The covariateIds of the target variables made missing
#' @param causeVariables The covariateIds of the variables driving missingness
#' @param imputation The imputation method name
#' @param predictionModel The prediction model name
#' @return A list containing the extracted evaluateion tables: `evaluationStatistics`, `calibrationSummary`, `thresholdSummary`,
#' `demographicSummary`, and `predictionDistribution`. Each table is returned as a dataframe with the identifier columns appended.


extractEvaluationObjects <- function(evaluation,
                                     simulation,
                                     mechanism,
                                     ratio,
                                     type,
                                     targetVariables,
                                     causeVariables,
                                     imputation,
                                     predictionModel) {
  # Return an empty set of tables when no evaluation output is available
  if (is.null(evaluation)) {
    return(list(
      evaluationStatistics = NULL,
      calibrationSummary = NULL,
      thresholdSummary = NULL,
      demographicSummary = NULL,
      predictionDistribution = NULL
    ))
  }

  # Attach scenario identifiers to a single evaluation  table
  addTable <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }

    as.data.frame(x) %>%
      dplyr::mutate(
        simulation = simulation,
        mechanism = mechanism,
        ratio = ratio,
        type = type,
        targetVariables = targetVariables,
        causeVariables = causeVariables,
        imputation = imputation,
        predictionModel = predictionModel
      )
  }

  # Return the set of evaluation tables
  list(
    evaluationStatistics = addTable(evaluation$evaluationStatistics),
    calibrationSummary = addTable(evaluation$calibrationSummary),
    thresholdSummary = addTable(evaluation$thresholdSummary),
    demographicSummary = addTable(evaluation$demographicSummary),
    predictionDistribution = addTable(evaluation$predictionDistribution)
  )
}

#' @title Collect evaluation tables across simulation runs
#' @description
#' This function iterates over the results of one simulation scenario,
#' extracts the standard PLP evaluation tables for each imputation-prediction model combination,
#' and combines the into one table per evaluation output type
#' @param modelResults The model results object containing one or more evaluated simulation runs
#' @param simulation The simulation run
#' @param mechanism The missingness mechanism
#' @param ratio The missingness ratio
#' @param type The missingness score type
#' @param targetVariables The covariateIds of the target variables made missing
#' @param causeVariables The covariateIds of the variables driving missingness
#' @return A list of data frames: `evaluationStatistics`, `calibrationSummary`,
#' `thresholdSummary`, `demographicsSummary`, `predictionDistribution`. Each table
#' combines the corresponding evaluation output across all model runs for the scenario
#'
#' @export
collectEvaluationTables <- function(modelResults,
                                    simulation,
                                    mechanism,
                                    ratio,
                                    type,
                                    targetVariables,
                                    causeVariables) {
  # Initialize tables
  allTables <- list(
    evaluationStatistics = list(),
    calibrationSummary = list(),
    thresholdSummary = list(),
    demographicSummary = list(),
    predictionDistribution = list()
  )

  # Extract the labeled tables for each scenario run
  for (i in seq_along(modelResults$modelResults)) {
    result <- modelResults$modelResults[[i]]
    tables <- extractEvaluationObjects(
      evaluation = result$evaluation,
      simulation = simulation,
      mechanism = mechanism,
      ratio = ratio,
      type = type,
      targetVariables = targetVariables,
      causeVariables = causeVariables,
      imputation = result$imputation,
      predictionModel = result$predictionModel
    )

    # Store each extracted table under its corresponding output type
    for (tableName in names(tables)) {
      allTables[[tableName]][[i]] <- tables[[tableName]]
    }
  }

  # Bind the collected tables into one dataframe per evaluation output
  lapply(allTables, dplyr::bind_rows)
}

#' @title Append a table to a csv file
#' @description This function appends a data frame to an existing csv file or creates a new
#' file if it does not yet exist.
#' @param table A dataframe to write to the csv
#' @param filePath The path to the csv file
#' @return invisible file path
appendCsvTable <- function(table, filePath) {
  # Skip writing when there is no new table or rows to append
  if (is.null(table) || nrow(table) == 0) {
    return(invisible(NULL))
  }

  # Create a new file if it does not yet exist
  dir.create(dirname(filePath), recursive = TRUE, showWarnings = FALSE)

  # Append rows to the existing files
  # Create a file with headers if the file does not yet exist
  readr::write_csv(
    table,
    file = filePath,
    append = file.exists(filePath),
    col_names = !file.exists(filePath)
  )

  invisible(filePath)
}

#' @title Parse evaluation table values
#' @description
#' This function parses evaluation metrics into a consistent format such that they can be summarised
#' across simulation runs in a later stage.
#' @param table A data frame containing the evaluation results
#' @return A data frame with parsed evaluation values for aggregation
#'
evaluationStatsHelper <- function(table) {
  # Return an empty table when no numeric metrics are available
  if (is.null(table) || nrow(table) == 0) {
    return(table)
  }

  # Parse values in the evaluation output to the correct format
  table %>%
    dplyr::mutate(
      evaluation = as.character(.data$evaluation),
      metric = as.character(.data$metric),
      value = readr::parse_double(as.character(.data$value))
    )
}

#' @title Append scenario outputs to raw result files
#' @description
#' This function appends the evaluation tables, model metrics, and scenario progress
#' information for one scenario to the raw csv files stored in the output folder.
#' Evaluation statistics are parsed into a consistent format before they are written.
#' @param evaluationTables A list of evaluation tables to append
#' @param metrics An optional metrics table to append
#' @param progress An optional scneario progress table to append
#' @param folder The output folder where the raw results are stored
#' @param ... Additional unused arguments passed through when calling the function
#' @return Invisibly return `NULL`.
#'
#' @export
appendScenarioResults <- function(evaluationTables,
                                  metrics = NULL,
                                  progress = NULL,
                                  folder,
                                  ...) {
  # Store the per-scenario retults under the raw results subfolder
  rawFolder <- file.path(folder, "raw")

  # Append each evaluation table to its matching raw csv file
  for (tableName in names(evaluationTables)) {
    table <- evaluationTables[[tableName]]

    # Parse evaluation statistics before writing to csv
    if (tableName == "evaluationStatistics" && !is.null(table) && nrow(table) > 0) {
      table <- evaluationStatsHelper(table)
    }

    appendCsvTable(table, file.path(rawFolder, paste0(tableName, ".csv")))
  }

  # Append per scenario model metrics and progress tables
  appendCsvTable(metrics, file.path(rawFolder, "metrics.csv"))
  appendCsvTable(progress, file.path(rawFolder, "scenarioProgress.csv"))
}

#' @title Load raw evaluation tables
#' @description
#' This function reads the raw evaluation tables. Missing files are returned as `NULL`
#' @param folder The simulation output folder containing the `raw` subfolder
#' @return A list containing the raw evaluatin tables:
#' `evaluationTables`, `calibrationSummary`, `thresholdSummary`,`demographicsSummary`,
#' and `predictionDistribution`
#' @export
loadRawEvaluationTables <- function(folder) {
  # Read all raw evaluation tables
  rawFolder <- file.path(folder, "raw")
  tableNames <- c(
    "evaluationStatistics",
    "calibrationSummary",
    "thresholdSummary",
    "demographicSummary",
    "predictionDistribution"
  )

  # Load each table if it exists, otherwise keep it as `NULL`
  tables <- lapply(tableNames, function(tableName) {
    filePath <- file.path(rawFolder, paste0(tableName, ".csv"))

    if (!file.exists(filePath)) {
      return(NULL)
    }

    readr::read_csv(filePath, show_col_types = FALSE)
  })

  # Return loaded tables as a named list
  names(tables) <- tableNames
  tables
}

#' @title Summarise evaluation results across simulation runs
#' @description
#' This function aggregates an evaluation table across simulation runs by
#' grouping on the scenario and result identifiers and computing mean and standard deviation for the numeric results columns
#' @param table A raw evaluation table containing evaluation statistics of multiple simulation runs
#' @return A data frame summarising numeric results across simulations
#' @export
summariseOverSims <- function(table) {
  # Return an emtpy table if there is nothing to aggregate
  if (is.null(table) || nrow(table) == 0) {
    return(table)
  }

  # Parse evaluation statistics before aggregation
  if (all(c("evaluation", "metric", "value") %in% names(table))) {
    table <- evaluationStatsHelper(table)
  }

  # Keep the identifier columns tthat define one scenario result
  groupVars <- intersect(
    c(
      "simulation",
      "mechanism",
      "ratio",
      "type",
      "imputation",
      "targetVariables",
      "causeVariables",
      "predictionModel",
      "evaluation",
      "evaluationType",
      "metric",
      "category",
      "threshold"
    ),
    names(table)
  )

  # Summarise only the numeric results columns
  numericV <- table %>%
    dplyr::select(dplyr::where(is.numeric)) %>%
    names()

  numericV <- setdiff(numericV, c("simulation", groupVars))

  # Aggregate the numeric results across simulation runs
  table %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(setdiff(groupVars, "simulation")))) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(numericV),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ stats::sd(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      nSimulation = dplyr::n_distinct(.data$simulation),
      .groups = "drop"
    )
}


#' @title Save summarised evaluation results
#' @description This function writes the aggregated evaluation tables to disk. It stores
#' one master csv per table type and also writes scenario specific subsets
#' @param aggregatedTables A list of aggregated evaluation tables
#' @param folder The output folder where the summarised results should be saved
#' @return Invisibly returns `NULL`
#' @export
saveSummarisedResults <- function(aggregatedTables, folder) {
  # Create the output folders for master tables and scenario specific tables
  mainFolder <- file.path(folder, "master")
  specificFolder <- file.path(folder, "specifications")

  dir.create(mainFolder, recursive = TRUE, showWarnings = FALSE)
  dir.create(specificFolder, recursive = TRUE, showWarnings = FALSE)

  # Write each aggregated table to the master output folder
  for (tableName in names(aggregatedTables)) {
    table <- aggregatedTables[[tableName]]

    if (is.null(table) || nrow(table) == 0) {
      next
    }

    readr::write_csv(table, file.path(mainFolder, paste0(tableName, ".csv")))

    # Skip scenario specific outputs when the identifying columns are absent
    if (!all(c("mechanism", "ratio", "type", "imputation", "predictionModel") %in% names(table))) {
      next
    }

    # Identify the unique scenario/model configurations present in the table
    configurations <- table %>%
      dplyr::distinct(
        .data$mechanism,
        .data$ratio,
        .data$type,
        .data$imputation,
        .data$predictionModel
      )

    # Write one table per configuration to a nested subfolder
    for (i in seq_len(nrow(configurations))) {
      config <- configurations[i, ]

      configTab <- table %>%
        dplyr::filter(
          .data$mechanism == config$mechanism,
          .data$ratio == config$ratio,
          .data$type == config$type,
          .data$imputation == config$imputation,
          .data$predictionModel == config$predictionModel
        )

      outputDirectory <- file.path(
        specificFolder,
        config$mechanism,
        paste0("ratio_", config$ratio),
        config$imputation,
        config$predictionModel
      )

      dir.create(outputDirectory, recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(configTab, file.path(outputDirectory, paste0(tableName, ".csv")))
    }
  }
}

#' @title Collapse identifiers into one label
#' @description
#' This function converts a vector of identifiers into a single label separated by underscores
#' @param x A vector of identifiers
#' @return A single character string containing the sorted unique identifiers, or NA_character_ when `x` is `NULL`
#'
#' @export
collapseIds <- function(x) {
  # Return a missing label when no identifiers are supplied
  if (is.null(x)) {
    return(NA_character_)
  }

  # Build one label from the sorted unique identifiers
  paste(sort(unique(x)), collapse = "_")
}

#' @title Extract identifiers from missingness patterns
#' @description
#'  This function extracts an identifier field from a list of pattern specifications
#'  and combines the values across all patterns
#'  @param patterns A list of pattern specifications
#'  @param field The field name to extract from each pattern
#'  @return A vector of extracted identifiers
#'
#' @export
extractPatternIds <- function(patterns, field) {
  # Collect identifiers across all pattern definitions
  ids <- unlist(lapply(patterns, function(pattern) pattern[[field]]), use.names = FALSE)

  # Return NULL when the requested fields are absent
  if (length(ids) == 0) {
    return(NULL)
  }

  ids
}

#' @title Resolve a mechanism specific value
#' @description This function returns the value defined for one missingness mechanism from a named list. When the
#' mechanism is absent a default value is returned instead
#' @param valuesPerMech A list of values keyed by mechanism
#' @param mech The mechanism name to resolve
#' @param defaultValue A default value to return to when `mech` is not present
#' @return The mechanism specific value
resolvePerMechanismValue <- function(valuesPerMech, mech, defaultValue = NULL) {
  # Fall back to the default value when no per-mechanism values are supplied
  if(is.null(valuesPerMech)) {
    return(defaultValue)
  }

  # Look up the value for the requested mechanism
  value <- valuesPerMech[[mech]]

  # Fall back to the default value when the mechanism has no explicit value
  if (is.null(value)) {
    return(defaultValue)
  }

  value
}

#' @title Determine the default simulation covariates
#' @description
#'  This function collects the target and cause covariate identifiers that must be
#'  present in the complete input data. It supports mechanism-specific inputs and
#'  multivariate pattern specifications
#' @param targetCovariateId The default target covariate identifier (or a vecotr of identifiers)
#' @param causeCovariateIds Optional default cause covariate identifiers
#' @param targetCovariateIdPerMech Optional mechanism-specific target covariate identifiers
#' @param causeCovariateIdsPerMech Optional mechanism-specific cause covariate identifiers
#' @param patterns Optional list of explicit pattern specificiations
#' @return A vector of required covariate identifiers
#'
defaultSimulationInputs <- function(targetCovariateId, causeCovariateIds = NULL, targetCovariateIdPerMech = NULL, causeCovariateIdsPerMech = NULL, patterns = NULL) {
  # When patterns are provided, extract the required inputs from them
  if (!is.null(patterns)) {
    return(defaultRequiredCovariates(
      targetCovariateId = extractPatternIds(patterns, "targetCovariateIds"),
      causeCovariateIds = extractPatternIds(patterns, "causeCovariateIds")
    ))
  }

  # Otherwise combine the default and mechanism-specific target/cause inputs
  perMechanismTargets <- unlist(targetCovariateIdPerMech, use.names = FALSE)
  perMechanismCauses <- unlist(causeCovariateIds, use.names = FALSE)

  defaultRequiredCovariates(
    targetCovariateId = c(targetCovariateId, perMechanismTargets),
    causeCovariateIds = c(causeCovariateIds, perMechanismCauses)
  )
}

#' @title Build a default missingness pattern
#' @description
#' This function constructs the standard pattern definition for one missingness
#' mechanism using the supplied target and cause covariates
#' @param mech The missingness mechanism: `MCAR`, `MAR`, `MNAR`
#' @param targetCovariateId The target covariate identifier(s)
#' @param causeCovariateIds Optional cause covariate identifiers
#' @param type The score transformation type used for `MAR` and `MNAR`
#' @return A pattern specification list for the requested mechanism
buildMechanismPattern <- function(mech, targetCovariateId, causeCovariateIds = NULL, type = "RIGHT") {
  # Construct the default pattern structure for the requested mechanism
  switch(
    mech,
    MCAR = list(targetCovariateIds = targetCovariateId, mechanism = "MCAR"),
    MAR = list(
      targetCovariateIds = targetCovariateId,
      causeCovariateIds = causeCovariateIds,
      mechanism = "MAR",
      type = type
    ),
    MNAR = list(
      targetCovariateIds = targetCovariateId,
      causeCovariateIds = targetCovariateId,
      mechanism = "MNAR",
      type = type
    ),
    stop("Unsupported mechanism: ", mech, call.=FALSE)
  )
}
