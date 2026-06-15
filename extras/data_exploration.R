################################################################################
# Libraries
################################################################################
library(DatabaseConnector)
library(CirceR)
library(readr)
library(SqlRender)
library(PatientLevelPrediction)
library(FeatureExtraction)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(VIM)
library(corrplot)
library(stringr)
library(paletteer)

################################################################################
# Specify paths
################################################################################
serverPath <- "C:\\Users\\maud5\\Documents\\database-1M_filtered.duckdb"
codePath <- "C:\\Users\\maud5\\OneDrive\\MSc 2025-2026\\Thesis applications\\Code\\"
figuresPath <- "C:\\Users\\maud5\\OneDrive\\MSc 2025-2026\\Thesis applications\\Plots\\"

dir.create(
  figuresPath,
  recursive = TRUE,
  showWarnings = FALSE
)

targetJson <- file.path(
  codePath,
  "target_cohort_mace_age40_79_strict.json"
)

outcomeJson <- file.path(
  codePath,
  "outcome_cohort_mace 1.json"
)

################################################################################
# Connect to database
################################################################################
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "duckdb",
  server = serverPath
)

connection <- DatabaseConnector::connect(connectionDetails)

################################################################################
# Generate target and outcome cohorts
################################################################################
buildAndExecuteCohort <- function(connection,
                                  jsonPath,
                                  cohortId,
                                  cdmSchema = "main",
                                  vocabularySchema = "main",
                                  cohortSchema = "main",
                                  cohortTable = "cohort") {
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
    targetDialect = "duckdb"
  )
  
  DatabaseConnector::executeSql(connection, translatedSQL)
}

buildAndExecuteCohort(
  connection = connection,
  jsonPath = targetJson,
  cohortId = 1
)

buildAndExecuteCohort(
  connection = connection,
  jsonPath = outcomeJson,
  cohortId = 2
)

################################################################################
# Extract PLP data
################################################################################
cdmDbSchema <- "main"
cohortsDbSchema <- "main"
cohortsDbTable <- "cohort"

covariateSettings <- FeatureExtraction::createCovariateSettings(
  useDemographicsGender = TRUE,
  useDemographicsAge = TRUE,
  useDemographicsAgeGroup = TRUE,
  useConditionOccurrenceLongTerm = TRUE,
  useDrugExposureLongTerm = TRUE,
  useMeasurementValueLongTerm = TRUE,
  longTermStartDays = -365,
  endDays = 0
)

databaseDetails <- PatientLevelPrediction::createDatabaseDetails(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDbSchema,
  cdmDatabaseName = "",
  cohortDatabaseSchema = cohortsDbSchema,
  cohortTable = cohortsDbTable,
  targetId = 1,
  outcomeDatabaseSchema = cohortsDbSchema,
  outcomeTable = cohortsDbTable,
  outcomeIds = 2,
  cdmVersion = 5
)

restrictPlpDataSettings <- PatientLevelPrediction::createRestrictPlpDataSettings()

plpData <- PatientLevelPrediction::getPlpData(
  databaseDetails = databaseDetails,
  covariateSettings = covariateSettings,
  restrictPlpDataSettings = restrictPlpDataSettings
)

################################################################################
# Create study population
################################################################################
populationSettings <- PatientLevelPrediction::createStudyPopulationSettings(
  binary = TRUE,
  includeAllOutcomes = TRUE,
  firstExposureOnly = TRUE,
  washoutPeriod = 365,
  removeSubjectsWithPriorOutcome = TRUE,
  priorOutcomeLookback = 99999,
  requireTimeAtRisk = FALSE,
  minTimeAtRisk = 0,
  riskWindowStart = 1,
  startAnchor = "cohort start",
  riskWindowEnd = 1095,
  endAnchor = "cohort start"
)

population <- PatientLevelPrediction::createStudyPopulation(
  plpData = plpData,
  outcomeId = 2,
  populationSettings = populationSettings
)

populationRowIds <- population$rowId

covariateAnalysis <- plpData$covariateData$analysisRef %>% collect()
covariateRef <- plpData$covariateData$covariateRef %>% collect()

covariatePopulation <- plpData$covariateData$covariates %>%
  filter(.data$rowId %in% populationRowIds) %>%
  collect()

################################################################################
# Keep only age, age group, gender, and the five measurements
################################################################################
measurementConceptIds <- c(
  3004249, # Systolic Blood Pressure
  3038553, # BMI
  3027114, # Total cholesterol
  3007070, # HDL cholesterol
  3009966  # LDL cholesterol
)

measurementConceptLabels <- c(
  "3004249" = "Systolic BP",
  "3038553" = "BMI",
  "3027114" = "Total Cholesterol",
  "3007070" = "HDL",
  "3009966" = "LDL"
)

bpConceptId <- 3004249

demographicCovariateIds <- c(
  8003, 9003, 10003, 11003, 12003, 13003, 14003, 15003, # Age group
  1002, # Age
  8507001, 8532001 # Gender
)

demographicCategories <- covariatePopulation %>%
  filter(covariateId %in% demographicCovariateIds) %>%
  left_join(
    covariateRef %>% select(covariateId, covariateName),
    by = "covariateId"
  ) %>%
  mutate(
    genderCat = case_when(
      covariateId == 8507001 ~ "Male",
      covariateId == 8532001 ~ "Female",
      TRUE ~ NA_character_
    ),
    ageGroupCat = case_when(
      covariateId == 8003 ~ "40-44",
      covariateId == 9003 ~ "45-49",
      covariateId == 10003 ~ "50-54",
      covariateId == 11003 ~ "55-59",
      covariateId == 12003 ~ "60-64",
      covariateId == 13003 ~ "65-69",
      covariateId == 14003 ~ "70-74",
      covariateId == 15003 ~ "75-79",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(rowId) %>%
  summarise(
    genderCat = first(na.omit(genderCat)),
    ageGroupCat = first(na.omit(ageGroupCat)),
    .groups = "drop"
  )

demographicsRef <- covariateRef %>%
  filter(.data$covariateId %in% demographicCovariateIds)

demographicsAnalysis <- covariateAnalysis %>%
  semi_join(demographicsRef, by = "analysisId")

demographicsCovariateValues <- covariatePopulation %>%
  filter(.data$covariateId %in% demographicCovariateIds) %>%
  select(rowId, covariateId, covariateValue)

measurementRef <- covariateRef %>%
  filter(
    !is.na(.data$conceptId),
    .data$conceptId %in% measurementConceptIds
  ) %>%
  distinct(covariateId, .keep_all = TRUE)

measurementAnalysis <- covariateAnalysis %>%
  semi_join(measurementRef, by = "analysisId")

measurementCovariateValues <- covariatePopulation %>%
  filter(.data$covariateId %in% measurementRef$covariateId) %>%
  select(rowId, covariateId, covariateValue)

analysisCovariateRef <- bind_rows(demographicsRef, measurementRef)
analysisCovariateAnalysis <- bind_rows(demographicsAnalysis, measurementAnalysis)

################################################################################
# Helper functions
################################################################################
plpDataHelper <- function(labels,
                          folds = NULL,
                          covariates,
                          covariateRef,
                          analysisRef,
                          templatePLPData = plpData) {
  covariateData <- Andromeda::andromeda(
    covariates = covariates,
    covariateRef = covariateRef,
    analysisRef = analysisRef
  )
  
  class(covariateData) <- "CovariateData"
  attr(class(covariateData), "package") <- "FeatureExtraction"
  attr(covariateData, "metaData") <- attr(templatePLPData$covariateData, "metaData")
  
  plpDataOut <- list(
    labels = labels,
    covariateData = covariateData,
    metaData = templatePLPData$metaData
  )
  
  if (!is.null(folds)) {
    plpDataOut$folds <- folds
  }
  
  class(plpDataOut) <- "plpData"
  plpDataOut
}

populationSubset <- function(population, selectedRowIds) {
  filteredPopulation <- population %>%
    dplyr::filter(.data$rowId %in% selectedRowIds)
  
  attr(filteredPopulation, "metaData") <- attr(population, "metaData")
  filteredPopulation
}

buildPopulationPLPData <- function(selectedRowIds, includeFolds = FALSE) {
  selectedRowIds <- sort(unique(selectedRowIds))
  
  finalPopulation <- populationSubset(population, selectedRowIds)
  
  analysisCovariates <- bind_rows(
    demographicsCovariateValues %>% filter(.data$rowId %in% selectedRowIds), 
    measurementCovariateValues %>% filter(.data$rowId %in% selectedRowIds)
  )
  
  folds <- NULL
  if (includeFolds) {
    folds <- data.frame(
      rowId = selectedRowIds,
      index=1L
    )
  }
  
  list(
    plpData = plpDataHelper(
      labels = finalPopulation, 
      folds = folds, 
      covariates = analysisCovariates, 
      covariateRef = analysisCovariateRef, 
      analysisRef = analysisCovariateAnalysis
    ), population = finalPopulation
  )
}

################################################################################
# Build full and complete-case datasets
################################################################################
buildPopulation <- buildPopulationPLPData(populationRowIds)
fullPopulationPLPData <- buildPopulation$plpData
fullPopulation <- buildPopulation$population

# Extract bp covariate id
bpCovariateId <- measurementRef %>%
  filter(.data$conceptId == bpConceptId) %>%
  pull(.data$covariateId) %>%
  unique()

bpCompleteCaseRowIds <- measurementCovariateValues %>%
  filter(.data$covariateId == bpCovariateId) %>%
  distinct(rowId) %>%
  pull(rowId) %>%
  sort()

buildPopulationCompleteBp <- buildPopulationPLPData(bpCompleteCaseRowIds)
completePopulationBpPLPData <- buildPopulationCompleteBp$plpData
completePopulationBp <- buildPopulationCompleteBp$population

# allCompleteCaseRowIds <- measurementCovariateValues %>%
#   filter(.data$covariateId %in% measurementRef$covariateId) %>%
#   distinct(rowId, covariateId) %>%
#   count(rowId, name = "nrMeasurementsObserved") %>%
#   filter(.data$nrMeasurementsObserved == length(unique(measurementRef$covariateId))) %>%
#   pull(rowId) %>%
#   sort()
# 
# buildPopulationCompleteAll <- buildPopulationPLPData(allCompleteCaseRowIds)
# completePopulationAllPLPData <- buildPopulationCompleteAll$plpData
# completePopulationAll <- buildPopulationCompleteAll$population

populationCovData <- fullPopulationPLPData$covariateData
completePopulationBpCovData <- completePopulationBpPLPData$covariateData
#completePopulationAllCovData <- completePopulationAllPLPData$covariateData

################################################################################
# Table 1
################################################################################
aggregatedTabOneCovData <- FeatureExtraction::aggregateCovariates(populationCovData)
tableOne <- FeatureExtraction::createTable1(
  aggregatedTabOneCovData,
  output = "one column",
  showCounts = TRUE,
  showPercent = TRUE
)

aggregatedTabOneCovDataCompleteBp <- FeatureExtraction::aggregateCovariates(completePopulationBpCovData)
tableOneCompleteBp <- FeatureExtraction::createTable1(
  aggregatedTabOneCovDataCompleteBp,
  output = "one column",
  showCounts = TRUE,
  showPercent = TRUE
)

# aggregatedTabOneCovDataCompleteAll <- FeatureExtraction::aggregateCovariates(completePopulationAllCovData)
# tableOneCompleteAll <- FeatureExtraction::createTable1(
#   aggregatedTabOneCovDataCompleteAll,
#   output = "one column",
#   showCounts = TRUE,
#   showPercent = TRUE
# )

################################################################################
# Measurement summary table from covariatesContinuous
################################################################################
measurementSummaryFull <- aggregatedTabOneCovData$covariatesContinuous %>%
  collect() %>%
  inner_join(
    measurementRef %>%
      select(covariateId, conceptId, covariateName),
    by = "covariateId"
  ) %>%
  mutate(
    conceptLabel = measurementConceptLabels[as.character(conceptId)]
  ) %>%
  select(
    conceptId,
    conceptLabel,
    covariateId,
    nObserved = countValue,
    mean = averageValue,
    sd = standardDeviation,
    median = medianValue,
    p25 = p25Value,
    p75 = p75Value,
    min = minValue,
    max = maxValue
  ) %>%
  arrange(match(conceptId, measurementConceptIds))

measurementSummaryCompleteBp <- aggregatedTabOneCovDataCompleteBp$covariatesContinuous %>%
  collect() %>%
  inner_join(
    measurementRef %>%
      select(covariateId, conceptId, covariateName),
    by = "covariateId"
  ) %>%
  mutate(
    conceptLabel = measurementConceptLabels[as.character(conceptId)]
  ) %>%
  select(
    conceptId,
    conceptLabel,
    covariateId,
    nObserved = countValue,
    mean = averageValue,
    sd = standardDeviation,
    median = medianValue,
    p25 = p25Value,
    p75 = p75Value,
    min = minValue,
    max = maxValue
  ) %>%
  arrange(match(conceptId, measurementConceptIds))

measurementSummaryFull
measurementSummaryCompleteBp

################################################################################
# Plotting datasets
################################################################################
makePlotData <- function(data, finalPopulation) {
  covariateLookup <- data$covariateData$covariateRef %>%
    collect() %>%
    transmute(
      covariateId,
      variableName = dplyr::coalesce(
        unname(measurementConceptLabels[as.character(conceptId)]),
        covariateName
      )
    )
  
  data$covariateData$covariates %>%
    collect() %>%
    inner_join(covariateLookup, by = "covariateId") %>%
    group_by(rowId, variableName) %>%
    summarise(
      covariateValue = mean(covariateValue, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = variableName,
      values_from = covariateValue
    ) %>%
    left_join(
      finalPopulation %>% select(rowId, subjectId),
      by = "rowId"
    ) %>%
    left_join(demographicCategories, by = "rowId")
}

plotDataFull <- makePlotData(fullPopulationPLPData, fullPopulation)
plotDataCompleteBp <- makePlotData(completePopulationBpPLPData, completePopulationBp)

################################################################################
# Plotting data
################################################################################
plotVariables <- c(
  "rowId",
  "subjectId",
  "ageGroupCat",
  "genderCat",
  "age in years",
  unname(measurementConceptLabels[as.character(measurementConceptIds)])
)

makeAnalysisPlotData <- function(data) {
  data %>%
    rename(
      `Age group` = ageGroupCat,
      Gender = genderCat,
      `Age in years` = `age in years`
    ) %>%
    select(
      rowId,
      subjectId,
      `Age group`,
      Gender,
      `Age in years`,
      all_of(unname(measurementConceptLabels[as.character(measurementConceptIds)]))
    )
}

plotDataFullClean <- makeAnalysisPlotData(plotDataFull)
plotDataCompleteBpClean <- makeAnalysisPlotData(plotDataCompleteBp)

################################################################################
# Missingness overview
################################################################################
png(
  filename = file.path(figuresPath, "aggregation_plot.png"),
  width = 8,
  height = 6,
  units = "in",
  res = 300
)

par(mar = c(12, 8, 2, 2))

VIM::aggr(
  plotDataFullClean,
  numbers = FALSE,
  prop = TRUE,
  sortVars = TRUE,
  cex.axis = 0.55,
  las = 2,
  gap = 3
)

dev.off()

measurementLabels <- c(
  "Systolic BP",
  "LDL",
  "HDL",
  "Total Cholesterol",
  "BMI"
)

measurementObserved <- rowSums(!is.na(plotDataFullClean[measurementLabels]))

allFiveMiss <- mean(measurementObserved == 0) * 100
atLeastOneObs <- mean(measurementObserved >= 1) * 100
allFiveObs <- mean(measurementObserved == 5) * 100

subsetData <- plotDataFullClean[
  measurementObserved > 0,
]


png(
  filename = file.path(figuresPath, "matrix_plot.png"),
  width = 8,
  height = 6,
  units = "in",
  res = 300
)

VIM::matrixplot(
  subsetData[measurementLabels],
  sortby = "Systolic BP",
  interactive = FALSE
)

dev.off()

################################################################################
# Demographic distributions in full population
################################################################################
ageDensity <- ggplot(plotDataFullClean, aes(x = .data[["Age in years"]])) +
  geom_density() +
  labs(
    title = "Age density",
    x = "Age",
    y = "Density"
  ) + theme_minimal()

ggsave(
  filename = file.path(figuresPath, "age_density.png"),
  plot = ageDensity, 
  width = 8,
  height = 6, 
  dpi = 300)

plotPercentBar <- function(data, varName, title) {
  ggplot(data, aes(x = .data[[varName]])) +
    geom_bar(aes(y = after_stat(count / sum(count)))) +
    scale_y_continuous(labels = percent_format()) +
    labs(
      x = varName,
      y = "Percentage",
      title = title
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_minimal()
}

ageDistr <- plotPercentBar(plotDataFullClean, "Age group", "Age group distribution")
ggsave(
  filename = file.path(figuresPath, "age_group_distribution.png"),
  plot = ageDistr, 
  width = 8,
  height = 6, 
  dpi = 300)

genderDistr <- plotPercentBar(plotDataFullClean, "Gender", "Gender distribution")
ggsave(
  filename = file.path(figuresPath, "gender_distribution.png"),
  plot =genderDistr, 
  width = 8,
  height = 6, 
  dpi = 300)
################################################################################
# Measurement distributions in full population
################################################################################
measurementPlotSpecifications <- list(
  list(var = "Systolic BP", label = "Systolic BP", x_label = "Systolic Blood Pressure (mmHg)"),
  list(var = "Total Cholesterol", label = "Total Cholesterol", x_label = "Total Cholesterol (mg/dL)"),
  list(var = "HDL", label = "HDL", x_label = "HDL (mg/dL)"),
  list(var = "LDL", label = "LDL", x_label = "LDL (mg/dL)"),
  list(var = "BMI", label = "BMI", x_label = expression(BMI~(kg/m^2)))
)

plotMeasurementDensity <- function(data, var, label, x_label) {
  ggplot(data, aes(x = .data[[var]])) +
    geom_density(na.rm = TRUE) +
    labs(
      title = paste(label, "distribution"),
      x = x_label,
      y = "Density"
    ) +
    theme_minimal()
}

plotMeasurementHistogram <- function(data, var, label, x_label, bins = 30) {
  ggplot(data, aes(x = .data[[var]])) +
    geom_histogram(bins = bins, na.rm = TRUE) +
    labs(
      title = paste(label, "distribution"),
      x = x_label,
      y = "Count"
    ) +
    theme_minimal()
}

densityPlots <- purrr::map(
  measurementPlotSpecifications,
  ~ plotMeasurementDensity(
    data = plotDataFullClean,
    var = .x$var,
    label = .x$label,
    x_label = .x$x_label
  )
)

histogramPlots <- purrr::map(
  measurementPlotSpecifications,
  ~ plotMeasurementHistogram(
    data = plotDataFullClean,
    var = .x$var,
    label = .x$label,
    x_label = .x$x_label
  )
)

purrr::walk2(
  densityPlots,
  measurementPlotSpecifications,
  ~ ggsave(
    filename = file.path(
      figuresPath,
      paste0(
        "density_",
        str_replace_all(tolower(.y$var), " ", "_"),
        ".png"
      )
    ),
    plot = .x,
    width = 7,
    height = 5,
    dpi = 300
  )
)

purrr::walk2(
  histogramPlots,
  measurementPlotSpecifications,
  ~ ggsave(
    filename = file.path(
      figuresPath,
      paste0(
        "histogram_",
        str_replace_all(tolower(.y$var), " ", "_"),
        ".png"
      )
    ),
    plot = .x,
    width = 7,
    height = 5,
    dpi = 300
  )
)

################################################################################
# Correlation plot
################################################################################
corrData <- plotDataFullClean %>%
  select(`Age in years`, all_of(measurementLabels))

corrMatrix <- cor(
  corrData,
  use = "pairwise.complete.obs"
)

png(
  filename = file.path(figuresPath, "correlation_plot.png"),
  width = 8,
  height = 6,
  units = "in",
  res = 300
)

corrplot(
  corrMatrix,
  method = "color",
  type = "upper",
  addCoef.col = "black",
  tl.cex = 0.8,
  number.cex = 0.7
)

dev.off()

################################################################################
# Margin plots
################################################################################
marginplotOutputPath <- file.path(
  figuresPath,
  "Marginplots"
)

marginplotSpecs <- list(
  list(x = "Systolic BP", y = "Age in years", file_name = "01_bp_vs_age.png"),
  list(x = "Total Cholesterol", y = "Systolic BP", file_name = "02_total_chol_vs_bp.png"),
  list(x = "BMI", y = "Systolic BP", file_name = "03_bmi_vs_bp.png"),
  list(x = "BMI", y = "Age in years", file_name = "04_bmi_vs_age.png"),
  list(x = "BMI", y = "Total Cholesterol", file_name = "05_bmi_vs_total_chol.png"),
  list(x = "Total Cholesterol", y = "Age in years", file_name = "06_total_chol_vs_age.png"),
  list(x = "Total Cholesterol", y = "LDL", file_name = "07_total_chol_vs_ldl.png"),
  list(x = "Total Cholesterol", y = "HDL", file_name = "08_total_chol_vs_hdl.png"),
  list(x = "Systolic BP", y = "HDL", file_name = "09_bp_vs_hdl.png"),
  list(x = "Systolic BP", y = "LDL", file_name = "10_bp_vs_ldl.png"),
  list(x = "BMI", y = "HDL", file_name = "11_bmi_vs_hdl.png"),
  list(x = "BMI", y = "LDL", file_name = "12_bmi_vs_ldl.png"),
  list(x = "HDL", y = "Age in years", file_name = "13_hdl_vs_age.png"),
  list(x = "LDL", y = "Age in years", file_name = "14_ldl_vs_age.png")
)

if (!is.null(marginplotOutputPath)) {
  dir.create(marginplotOutputPath, recursive = TRUE, showWarnings = FALSE)
}

runMarginplot <- function(spec,
                          data,
                          outputPath = NULL,
                          width = 1600,
                          height = 1600,
                          res = 150) {
  marginplotData <- data %>%
    select(all_of(c(spec$x, spec$y)))
  
  if (!is.null(outputPath)) {
    grDevices::png(
      filename = file.path(outputPath, spec$file_name),
      width = width,
      height = height,
      res = res
    )
    on.exit(grDevices::dev.off(), add = TRUE)
    VIM::marginplot(marginplotData)
  } else {
    VIM::marginplot(marginplotData)
  }
}

start <- Sys.time()

invisible(lapply(
  marginplotSpecs,
  runMarginplot,
  data = plotDataFullClean,
  outputPath = marginplotOutputPath
))

end <- Sys.time()
print(end - start)

################################################################################
# Comparison plots: full data vs BP complete case
################################################################################
plotDataFullClean <- plotDataFullClean %>%
  mutate(bpMissing = is.na(.data[["Systolic BP"]]))

bpAgeGroupComparison <- ggplot(
  plotDataFullClean,
  aes(x = .data[["Age group"]], fill = bpMissing)
) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(
    values = c(
      "TRUE" = "#F89C74FF",
      "FALSE" = "#9EB9F3FF"
    ),
    labels = c(
      "FALSE" = "Observed",
      "TRUE" = "Missing"
    )
  ) +
  labs(
    x = "Age group",
    y = "Percentage",
    fill = "Blood Pressure",
    title = "Blood Pressure Missingness by Age Group"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_minimal()

ggsave(
  filename = file.path(figuresPath, "bp_missingness_per_age_group.png"),
  plot =bpAgeGroupComparison, 
  width = 8,
  height = 6, 
  dpi = 300)

missingnessSummary <- plotDataFullClean %>%
  group_by(.data[["Age group"]]) %>%
  summarise(
    pctMissingBP = mean(is.na(.data[["Systolic BP"]])) * 100,
    n = n(),
    .groups = "drop"
  )

dataCombined <- bind_rows(
  plotDataFullClean %>% mutate(sample = "Full cohort"),
  plotDataCompleteBpClean %>% mutate(sample = "Complete-case cohort")
)

ageDensityPlotComp <- ggplot(
  dataCombined,
  aes(x = .data[["Age in years"]], fill = sample)
) +
  geom_density(alpha = 0.5) +
  paletteer::scale_fill_paletteer_d(
    "rcartocolor::Pastel",
    labels = c(
      "Full cohort" = "Full cohort",
      "Complete-case cohort" = "Complete-case cohort (on blood pressure)"
    )
  ) +
  labs(
    title = "Age density: full cohort vs complete-case cohort on blood pressure",
    x = "Age",
    y = "Density",
    fill = "Dataset"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(figuresPath, "age_density_comparison_missing.png"),
  plot =ageDensityPlotComp, 
  width = 8,
  height = 6, 
  dpi = 300)

ageGroupComparison <- ggplot(dataCombined, aes(x = sample, fill = .data[["Age group"]])) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_paletteer_d("rcartocolor::Pastel") +
  labs(
    x = NULL,
    y = "Percentage within sample",
    fill = "Age group",
    title = "Age group distribution: full cohort vs complete-case cohort on blood pressure"
  ) + theme_minimal()

ggsave(
  filename = file.path(figuresPath, "age_group_comparison_missing.png"),
  plot =ageGroupComparison, 
  width = 8,
  height = 6, 
  dpi = 300)

genderComparison <- ggplot(dataCombined, aes(x = sample, fill = .data[["Gender"]])) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_paletteer_d("rcartocolor::Pastel") +
  labs(
    x = NULL,
    y = "Percentage within sample",
    fill = "Gender",
    title = "Gender distribution: full cohort vs complete-case cohort on blood pressure"
  ) + theme_minimal()

ggsave(
  filename = file.path(figuresPath, "gender_comparison.png"),
  plot =genderComparison, 
  width = 8,
  height = 6, 
  dpi = 300)

