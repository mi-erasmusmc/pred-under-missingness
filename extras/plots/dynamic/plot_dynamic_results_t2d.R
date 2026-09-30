library(readr)
library(dplyr)
library(paletteer)
library(ggplot2)

masterSources <- list(
  list(path = "path/to/t2d_near_complete_master", scenario = "near complete"),
  list(path = "path/to/t2d_minimal_master", scenario = "minimal"),
  list(path = "path/to/t2d_very_mild_master", scenario = "very mild"),
  list(path = "path/to/t2d_mild_master", scenario = "mild"),
  list(path = "path/to/t2d_moderate_master", scenario = "moderate"),
  list(path = "path/to/t2d_realistic_master", scenario = "realistic"),
  list(path = "path/to/t2d_extreme_master", scenario = "extreme")
)

combinedOutputCsv <- NULL

resolveEvaluationFile <- function(path, summaryLevel = "master") {
  candidates <- unique(c(
    path,
    file.path(path, "evaluationStatistics.csv"),
    file.path(path, summaryLevel, "evaluationStatistics.csv"),
    file.path(path, "Missingness Simulation", summaryLevel, "evaluationStatistics.csv")
  ))

  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop("Could not find evaluationStatistics.csv for path: ", path, call. = FALSE)
  }

  existing[[1]]
}

readScenarioResults <- function(sources, summaryLevel = "master") {
  bind_rows(lapply(
    sources,
    function(source) {
      read_csv(
        resolveEvaluationFile(source$path, summaryLevel = summaryLevel),
        show_col_types = FALSE
      ) %>%
        mutate(scenario = source$scenario)
    }
  ))
}

recodeImputation <- function(x) {
  recode(
    x,
    "simpleMean_noIndicator" = "Mean",
    "simpleMean_withIndicator" = "Mean + indicator",
    "simpleMedian_noIndicator" = "Median",
    "simpleMedian_withIndicator" = "Median + indicator",
    "iterativePMM_noIndicator" = "Iterative PMM",
    "iterativePMM_withIndicator" = "Iterative PMM + indicator",
    "xgboost_noIndicator" = "Iterative XGBoost",
    "xgboost_withIndicator" = "Iterative XGBoost + indicator",
    "completeCase" = "Listwise Deletion",
    "noImputation" = "No imputation",
    .default = x
  )
}

recodePredictionModel <- function(x) {
  recode(
    x,
    "lasso" = "Lasso",
    "xgboost" = "XGBoost",
    "transformer" = "FT-Transformer",
    .default = x
  )
}

evaluationStats <- readScenarioResults(masterSources, summaryLevel = "master")

if (!is.null(combinedOutputCsv)) {
  write_csv(evaluationStats, combinedOutputCsv)
}

plotData <- evaluationStats %>%
  filter(
    evaluation == "Test",
    metric == "Eavg",
    predictionModel %in% c("lasso", "xgboost", "transformer")
  ) %>%
  mutate(
    predictionModel = recodePredictionModel(predictionModel),
    scenario = factor(
      scenario,
      levels = c("near complete", "minimal", "very mild", "mild", "moderate", "realistic", "extreme")
    ),
    imputation = recodeImputation(imputation),
    performance = if ("value_mean" %in% names(.)) as.numeric(value_mean) else as.numeric(value)
  )

plotData$predictionModel <- factor(
  plotData$predictionModel,
  levels = c("Lasso", "XGBoost", "FT-Transformer")
)
plotData$imputation <- factor(
  plotData$imputation,
  levels = c(
    "No imputation",
    "Listwise Deletion",
    "Mean",
    "Mean + indicator",
    "Median",
    "Median + indicator",
    "Iterative PMM",
    "Iterative PMM + indicator",
    "Iterative XGBoost",
    "Iterative XGBoost + indicator"
  )
)

ggplot(plotData, aes(x = scenario, y = imputation, fill = performance)) +
  geom_tile() +
  scale_fill_paletteer_c("ggthemes::Red-Green Diverging", direction = -1) +
  facet_wrap(~predictionModel) +
  labs(
    x = "Missingness scenario",
    y = "Imputation method",
    fill = "Eavg"
  ) +
  theme_minimal()
