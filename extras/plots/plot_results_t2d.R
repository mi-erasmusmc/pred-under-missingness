library(readr)
library(dplyr)
library(ggplot2)
library(paletteer)

# Configure the study result folders before running this script.
# Each entry may be either:
# - the study root directory
# - the "Missingness Simulation" directory
# - the evaluationStatistics.csv file itself
resultSources <- c(
  "path/to/t2d_results",
  "path/to/t2d_results_mnar",
  "path/to/t2d_results_mnar_xgb",
  "path/to/t2d_results_trans",
  "path/to/t2d_results_trans_mnar",
  "path/to/t2d_results_xgb_imp"
)

# Optional: set to a CSV path to save the combined results.
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

readBoundEvaluationStats <- function(paths, summaryLevel = "master") {
  if (length(paths) == 0) {
    stop("Provide at least one results path in `resultSources`.", call. = FALSE)
  }

  bind_rows(lapply(
    paths,
    function(path) {
      read_csv(resolveEvaluationFile(path, summaryLevel = summaryLevel), show_col_types = FALSE)
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

evalStats <- readBoundEvaluationStats(resultSources, summaryLevel = "master")

if (!is.null(combinedOutputCsv)) {
  write_csv(evalStats, combinedOutputCsv)
}

plotData <- evalStats %>%
  filter(
    evaluation == "Test",
    metric == "Eavg",
    predictionModel %in% c("lasso", "xgboost", "transformer")
  ) %>%
  mutate(
    mechanism = factor(mechanism, levels = c("MCAR", "MAR", "MNAR")),
    predictionModel = recodePredictionModel(predictionModel),
    imputation = recodeImputation(imputation),
    ratio = as.numeric(ratio),
    performance = if ("value_mean" %in% names(.)) as.numeric(value_mean) else as.numeric(value),
    performance_sd = if ("value_sd" %in% names(.)) as.numeric(value_sd) else NA_real_
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

createPlot <- function(modelName) {
  dataFrame <- plotData %>% filter(predictionModel == modelName)

  p <- ggplot(dataFrame, aes(x = ratio, y = performance, color = imputation, group = imputation)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    labs(
      x = "Missingness ratio",
      y = "Eavg",
      color = "Imputation method"
    ) +
    scale_x_continuous(breaks = sort(unique(dataFrame$ratio))) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  if ("mechanism" %in% names(dataFrame) && dplyr::n_distinct(dataFrame$mechanism) > 1) {
    p <- p + facet_wrap(~mechanism)
  }

  p
}

lassoPlot <- createPlot("Lasso")
xgbPlot <- createPlot("XGBoost")
transformerPlot <- createPlot("FT-Transformer")

lassoPlot
xgbPlot
transformerPlot

mnarPlot <- plotData %>% filter(mechanism == "MNAR")
ggplot(mnarPlot, aes(x = factor(ratio), y = imputation, fill = performance)) +
  geom_tile() +
  scale_fill_paletteer_c(
    "ggthemes::Red-Green Diverging",
    direction = 1,
    limits = c(0.5, 1)
  ) +
  facet_wrap(~predictionModel) +
  labs(
    x = "Missingness ratio",
    y = "Imputation method",
    fill = "Eavg"
  ) +
  theme_minimal()

baseline <- "No imputation"
scenarioCols <- intersect(
  c(
    "predictionModel",
    "mechanism",
    "type",
    "ratio",
    "targetVariables",
    "causeVariables",
    "evaluation",
    "metric"
  ),
  names(plotData)
)

diffData <- plotData %>%
  group_by(across(all_of(scenarioCols))) %>%
  mutate(
    baselinePerf = performance[match(baseline, imputation)],
    diffFromBaseline = performance - baselinePerf
  ) %>%
  ungroup() %>%
  filter(!is.na(baselinePerf))

ggplot(
  diffData %>% filter(mechanism == "MCAR"),
  aes(x = ratio, y = diffFromBaseline, color = imputation, group = imputation)
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  geom_point() +
  facet_wrap(~predictionModel) +
  labs(
    x = "Missingness ratio",
    y = paste("Eavg difference vs", baseline),
    color = "Imputation method"
  ) +
  theme_minimal(base_size = 12)
