library(scmamp)
library(stats)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

rawResultSources <- c(
  "path/to/raw_results_1",
  "path/to/raw_results_2"
)

cdPlotOutputDir <- NULL

resolveEvaluationFile <- function(path, summaryLevel = "raw") {
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

readBoundEvaluationStats <- function(paths, summaryLevel = "raw") {
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
    "simpleMean_withIndicator" = "Mean + I",
    "simpleMedian_noIndicator" = "Median",
    "simpleMedian_withIndicator" = "Median + I",
    "iterativePMM_noIndicator" = "Iterative PMM",
    "iterativePMM_withIndicator" = "Iterative PMM + I",
    "xgboost_noIndicator" = "Iterative XGBoost",
    "xgboost_withIndicator" = "Iterative XGBoost + I",
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

rawStats <- readBoundEvaluationStats(rawResultSources, summaryLevel = "raw")

aurocData <- rawStats %>%
  filter(
    evaluation == "Test",
    metric == "AUROC",
    ratio > 0
  ) %>%
  transmute(
    simulation,
    mechanism,
    ratio,
    predictionModel = recodePredictionModel(predictionModel),
    imputation = recodeImputation(imputation),
    AUROC = as.numeric(value)
  ) %>%
  distinct()

makeBlockMatrix <- function(df) {
  wide <- df %>%
    group_by(simulation, imputation) %>%
    summarise(AUROC = mean(AUROC, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = imputation, values_from = AUROC) %>%
    arrange(simulation)

  wide <- wide[stats::complete.cases(wide), , drop = FALSE]

  if (nrow(wide) < 2 || ncol(wide) < 3) {
    return(NULL)
  }

  as.matrix(wide[, -1, drop = FALSE])
}

runFriedmanOne <- function(df) {
  mat <- makeBlockMatrix(df)

  if (is.null(mat)) {
    return(tibble(
      n_blocks = NA_integer_,
      n_imputers = NA_integer_,
      statistic = NA_real_,
      df = NA_real_,
      p_value = NA_real_
    ))
  }

  ft <- stats::friedman.test(mat)

  tibble(
    n_blocks = nrow(mat),
    n_imputers = ncol(mat),
    statistic = unname(ft$statistic),
    df = unname(ft$parameter),
    p_value = ft$p.value
  )
}

runPosthocOne <- function(df) {
  mat <- makeBlockMatrix(df)

  if (is.null(mat)) {
    return(tibble())
  }

  pairs <- combn(colnames(mat), 2, simplify = FALSE)

  map_dfr(pairs, function(p) {
    wt <- wilcox.test(mat[, p[1]], mat[, p[2]], paired = TRUE, exact = FALSE)

    tibble(
      imputation_1 = p[1],
      imputation_2 = p[2],
      p_value = wt$p.value
    )
  }) %>%
    mutate(p_adj_holm = p.adjust(p_value, method = "holm"))
}

friedmanResults <- aurocData %>%
  group_by(mechanism, ratio, predictionModel) %>%
  group_modify(~ runFriedmanOne(.x)) %>%
  ungroup() %>%
  arrange(mechanism, ratio, predictionModel)

posthocResults <- aurocData %>%
  group_by(mechanism, ratio, predictionModel) %>%
  group_modify(~ {
    overall <- runFriedmanOne(.x)
    if (is.na(overall$p_value) || overall$p_value >= 0.05) {
      return(tibble())
    }
    runPosthocOne(.x)
  }) %>%
  ungroup()

plotCDAllGroups <- function(data, outDir, alpha = 0.05, cex = 1) {
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

  groups <- data %>%
    distinct(mechanism, ratio, predictionModel) %>%
    arrange(mechanism, ratio, predictionModel)

  summaryOut <- vector("list", nrow(groups))

  for (i in seq_len(nrow(groups))) {
    groupRow <- groups[i, ]
    subDf <- data %>%
      filter(
        mechanism == groupRow$mechanism,
        ratio == groupRow$ratio,
        predictionModel == groupRow$predictionModel
      )

    mat <- makeBlockMatrix(subDf)
    if (is.null(mat)) {
      summaryOut[[i]] <- tibble(
        mechanism = groupRow$mechanism,
        ratio = groupRow$ratio,
        predictionModel = groupRow$predictionModel,
        nBlocks = NA_integer_,
        nMethods = NA_integer_,
        friedmanP = NA_real_,
        plotted = FALSE
      )
      next
    }

    ft <- scmamp::friedmanTest(mat)
    fileName <- paste0("cd_", groupRow$mechanism, "_ratio_", groupRow$ratio, "_", groupRow$predictionModel, ".png")
    filePath <- file.path(outDir, fileName)

    png(filePath, width = 6300, height = 1100, res = 300)
    par(cex = 1.4)
    scmamp::plotCD(mat, alpha = alpha, cex = cex)
    dev.off()

    summaryOut[[i]] <- tibble(
      mechanism = groupRow$mechanism,
      ratio = groupRow$ratio,
      predictionModel = groupRow$predictionModel,
      nBlocks = nrow(mat),
      nMethods = ncol(mat),
      friedmanP = ft$p.value,
      plotted = TRUE
    )
  }

  bind_rows(summaryOut)
}

cdSummary <- NULL
if (!is.null(cdPlotOutputDir)) {
  cdSummary <- plotCDAllGroups(
    data = aurocData,
    outDir = cdPlotOutputDir,
    alpha = 0.05,
    cex = 1
  )
}

friedmanResults
posthocResults
cdSummary
