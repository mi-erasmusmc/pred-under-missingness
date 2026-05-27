# Libraries
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
################################################################################
# Connect to database
################################################################################
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "duckdb",
  server = "C:\\Users\\maud5\\Documents\\database-1M_filtered.duckdb"
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
  jsonPath = "C:\\Users\\maud5\\OneDrive\\MSc 2025-2026\\Thesis applications\\Code\\target_cohort_mace_age40_79_strict.json",
  cohortId = 1
)

buildAndExecuteCohort(
  connection = connection,
  jsonPath = "C:\\Users\\maud5\\OneDrive\\MSc 2025-2026\\Thesis applications\\Code\\outcome_cohort_mace 1.json",
  cohortId = 2
)

################################################################################
# Extract PLP data
################################################################################
cdmDbSchema <- "main"
cohortsDbSchema <- "main"
cohortsDbTable <- "cohort"

measurementConceptIds <- c(
  3004249, # Systolic Blood Pressure
  3038553, # BMI
  3027114, # Total cholesterol
  3007070, # HDL cholesterol
  3009966  # LDL cholesterol
)

baseCovariateSettings <- FeatureExtraction::createCovariateSettings(
  useDemographicsGender = TRUE,
  useDemographicsAge = TRUE,
  useDemographicsAgeGroup = TRUE,
  useDemographicsEthnicity = TRUE,
  useDemographicsRace = TRUE,
  useConditionOccurrenceAnyTimePrior = TRUE,
  useDrugExposureAnyTimePrior = TRUE,
  useDeviceExposureAnyTimePrior = TRUE,
  useProcedureOccurrenceAnyTimePrior = TRUE,
  useMeasurementAnyTimePrior = FALSE,
  useMeasurementValueAnyTimePrior = FALSE,
  useObservationAnyTimePrior = TRUE
)

measurementValueSettings <- FeatureExtraction::createCovariateSettings(
  useMeasurementValueAnyTimePrior = TRUE,
  includedCovariateConceptIds = measurementConceptIds,
  addDescendantsToInclude = FALSE
)


covariateSettings <- list(
  baseCovariateSettings,
  measurementValueSettings
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

# Can use this to restrict the sample size
restrictPlpDataSettings <- PatientLevelPrediction::createRestrictPlpDataSettings()

plpData <- PatientLevelPrediction::getPlpData(
  databaseDetails = databaseDetails,
  covariateSettings = covariateSettings,
  restrictPlpDataSettings = restrictPlpDataSettings
)

################################################################################
# Create population
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
  filter(.data$rowId %in% !!populationRowIds) %>%
  collect()

# PLP demographic covariates
demographicCovariateIds <- c(
  8003, 9003, 10003, 11003, 12003, 13003, 14003, 15003, # Age group
  1002, # Age
  8507001, 8532001, # Gender
  38003564005, 38003563005, # Ethnicity
  8516004, 8515004, 8527004 # Race
)

# Helper functions for PLP long-format data
plpDataHelper <- function(labels,
                            folds,
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
  attr(covariateData, "metaData") <- attr(templatePLPData$covariateData, "metaData")
  
  plpData <- list(
    labels = labels,
    folds = folds,
    covariateData = covariateData
  )
  
  class(plpData) <- "plpData"
  attr(plpData, "metaData") <- attr(templatePLPData, "metaData")
  
  plpData
}

# To get population subset from PLPData 
populationSubset <- function(population, selectedRowIds) {
  filteredPopulation <- population %>%
    filter(.data$rowId %in% selectedRowIds)
  
  attr(filteredPopulation, "metaData") <- attr(population, "metaData")
  filteredPopulation
}

################################################################################
# Build clean analysis PLP data
################################################################################
demographicsCovariateValues <- covariatePopulation %>%
  filter(.data$covariateId %in% demographicCovariateIds) %>%
  select("rowId", "covariateId", "covariateValue")

demographicsRef <- covariateRef %>%
  filter(.data$covariateId %in% demographicCovariateIds)

demographicsAnalysis <- covariateAnalysis %>%
  semi_join(demographicsRef, by = "analysisId")

measurementRef <- covariateRef %>%
  filter(
    .data$analysisId == 705,
    .data$conceptId %in% measurementConceptIds
  )

measurementAnalysis <- covariateAnalysis %>%
  semi_join(measurementRef, by = "analysisId")

measurementCovariateValues <- covariatePopulation %>%
  filter(.data$covariateId %in% measurementRef$covariateId) %>%
  select("rowId", "covariateId", "covariateValue")

populationRef <- bind_rows(demographicsRef, measurementRef)
populationAnalysis <- bind_rows(demographicsAnalysis, measurementAnalysis)

buildPopulationPLPData <- function(selectedRowIds) {
  selectedRowIds <- sort(unique(selectedRowIds))
  
  studyLabels <- data.frame(rowId = selectedRowIds)
  studyFolds <- data.frame(rowId = selectedRowIds, index = 1L)
  
  populationCovariates <- bind_rows(
    demographicsCovariateValues %>%
      filter(.data$rowId %in% selectedRowIds),
    measurementCovariateValues %>%
      filter(.data$rowId %in% selectedRowIds)
  )
  
  finalPopulation <- populationSubset(population, selectedRowIds)
  
  list(
    plpData = plpDataHelper(
      labels = studyLabels,
      folds = studyFolds,
      covariates = populationCovariates,
      covariateRef = populationRef,
      analysisRef = populationAnalysis
    ),
    population = finalPopulation
  )
}

# Full analysis data: all subjects in the PLP population.
buildPopulation <- buildPopulationPLPData(populationRowIds)
fullPopulationPLPData <- buildPopulation$plpData
fullPopulation <- buildPopulation$population

# Identify the systolic blood pressure covariate extracted by FeatureExtraction.
bpCovariateIds <- measurementRef %>%
  filter(.data$conceptId == 3004249) %>%
  pull(.data$covariateId) %>%
  unique()

bpCovariateId <- bpCovariateIds[[1]]

# BP complete-case analysis data: subjects with observed blood pressure at baseline.
bpCompleteCaseRowIds <- measurementCovariateValues %>%
  filter(.data$covariateId == bpCovariateId) %>%
  distinct(.data$rowId) %>%
  pull(.data$rowId) %>%
  sort()

buildPopulationCompleteBp <- buildPopulationPLPData(bpCompleteCaseRowIds)
completePopulationBpPLPData <- buildPopulationCompleteBp$plpData
completePopulationBp <- buildPopulationCompleteBp$population

# Complete case on all five measurements
allCompleteCaseRowIds <- measurementCovariateValues %>% 
  filter(.data$covariateId %in% measurementRef$covariateId) %>%
  distinct(rowId, covariateId) %>% 
  count(rowId, name="nrMeasurementsObserved") %>% 
  filter(nrMeasurementsObserved == length(unique(measurementRef$covariateId))) %>%
  pull(rowId) %>%
  sort()

buildPopulationCompleteAll <- buildPopulationPLPData(allCompleteCaseRowIds)
completePopulationAllPLPData <- buildPopulationCompleteAll$plpData
completePopulationAll <- buildPopulationCompleteAll$population

# Covariate data
populationCovData <- fullPopulationPLPData$covariateData
completePopulationBpCovData <- completePopulationBpPLPData$covariateData
completePopulationAllCovData <- completePopulationAllPLPData$covariateData

################################################################################
# Try table one from Feature Extraction package

################################################################################
# Full data
aggregatedTabOneCovData <- FeatureExtraction::aggregateCovariates(populationCovData)
tableOne <- FeatureExtraction::createTable1(
  aggregatedTabOneCovData,
  output = "one column",
  showCounts = TRUE,
  showPercent = TRUE
)

# Complete on Blood Pressure
aggregatedTabOneCovDataCompleteBp <- FeatureExtraction::aggregateCovariates(completePopulationBpCovData)
tableOneCompleteBp <- FeatureExtraction::createTable1(
  aggregatedTabOneCovDataCompleteBp,
  output="one column",
  showCounts = TRUE,
  showPercent = TRUE
)

# Complete on all five measurements
aggregatedTabOneCovDataCompleteAll <- FeatureExtraction::aggregateCovariates(completePopulationAllCovData)
tableOneCompleteAll <- FeatureExtraction::createTable1(
  aggregatedTabOneCovDataCompleteAll,
  output="one column",
  showCounts = TRUE,
  showPercent = TRUE
)

################################################################################
# Plotting datasets: full population and BP complete-case population
################################################################################

# Make wide data format to use for plotting
makePlotData <- function(studyPLPData, finalPopulation) {
  studyPLPData$covariateData$covariates %>%
    collect() %>%
    select(rowId, covariateId, covariateValue) %>%
    pivot_wider(
      names_from = covariateId,
      values_from = covariateValue
    ) %>%
    left_join(
      finalPopulation %>% select(rowId, subjectId),
      by = "rowId"
    )
}

plotDataFull <- makePlotData(
  studyPLPData = fullPopulationPLPData,
  finalPopulation = fullPopulation
)

plotDataCompleteBp <- makePlotData(
  studyPLPData = completePopulationBpPLPData,
  finalPopulation = completePopulationBp
)
################################################################################
# Create simplified demographic variables
################################################################################

addDemographicCategories <- function(data) {
  data <- data %>% mutate(
    genderCat = case_when(
      !is.na(`8507001`) ~ "Male",
      !is.na(`8532001`) ~ "Female",
      TRUE ~ NA_character_
    ),
    
    raceCat = case_when(
      !is.na(`8516004`) ~ "Black/African American",
      !is.na(`8515004`) ~ "Asian",
      !is.na(`8527004`) ~ "White",
      TRUE ~ NA_character_
    ),
    
    ethnicityCat = case_when(
      !is.na(`38003564005`) ~ "Non-Hispanic",
      !is.na(`38003563005`) ~ "Hispanic",
      TRUE ~ NA_character_
    ),
    ageGroupCat = case_when(
      !is.na(`8003`)  ~ "40-44",
      !is.na(`9003`)  ~ "45-49",
      !is.na(`10003`) ~ "50-54",
      !is.na(`11003`) ~ "55-59",
      !is.na(`12003`) ~ "60-64",
      !is.na(`13003`) ~ "65-69",
      !is.na(`14003`) ~ "70-74",
      !is.na(`15003`) ~ "75-79",
      TRUE ~ NA_character_
    )
  )
  data
}

plotDataFull <- addDemographicCategories(plotDataFull)
plotDataCompleteBp <- addDemographicCategories(plotDataCompleteBp)

################################################################################
# VIM plotting datasets
################################################################################

ageCovariateId <- "1002"
measurementIds <- as.character(measurementRef$covariateId)
covariateLabels <- c(
  "3027114840705" = "Total Cholesterol",
  "3007070840705" = "HDL", 
  "3009966840705" = "LDL",
  "3038553531705" = "BMI",
  "3004249876705" = "Systolic BP",
  "age.in.years" = "Age in years",
  "raceCat" = "Race",
  "ageGroupCat" = "Age group",
  "genderCat" = "Gender",
  "ethnicityCat" = "Ethnicity"
)

makeVIMPlotData <- function(data) {
  data %>%
    select(
      rowId,
      subjectId,
      ageGroupCat,
      raceCat,
      genderCat,
      ethnicityCat,
      any_of(measurementIds),
      any_of(ageCovariateId)
    ) %>%
    rename(age.in.years = all_of(ageCovariateId))
}

covariateRef %>%
  filter(covariateId == 1002) %>%
  select(covariateId, covariateName, analysisId)

plotDataVimFull <- makeVIMPlotData(plotDataFull)
plotDataVimCompleteBp <- makeVIMPlotData(plotDataCompleteBp)

plotDataVimFull <- plotDataVimFull %>%
  rename_with(
    ~ ifelse(.x %in% names(covariateLabels),
             covariateLabels[.x],
             .x)
  )

plotDataVimCompleteBp <- plotDataVimCompleteBp %>%
  rename_with(
    ~ ifelse(.x %in% names(covariateLabels),
             covariateLabels[.x],
             .x)
  )

# Aggregation plot to view general missingness proportions per variable + missingness occurence
par(mar = c(12, 8, 2, 2))
VIM::aggr(
  plotDataVimFull,
  numbers = FALSE,
  prop = TRUE,
  sortVars = TRUE,
  cex.axis = 0.55,
  las = 2,
  gap = 3
)

measurementLabels <- c(
  "Systolic BP",
  "LDL",
  "HDL",
  "Total Cholesterol",
  "BMI"
)
measurementObserved <- rowSums(
  !is.na(plotDataVimFull[measurementLabels])
)

# Number of observations with all five measures missing (90%)
allFiveMiss <- mean(measurementObserved == 0) * 100
# Number of observations with at least one measure (9.7%)
atLeastOneObs <- mean(measurementObserved >= 1) * 100
# Number of observations with all five measures observed (0.003%)
allFiveObs <- mean(measurementObserved == 5) * 100


# Matrix plot
subsetData <- plotDataVimFull[
  measurementObserved > 0,
]

VIM::matrixplot(
  subsetData[measurementLabels],
  sortby = "Systolic BP",
  interactive = FALSE
)

################################################################################
# Demographic distributions in the full population, without filtering

################################################################################
# Age density
ggplot(plotDataVimFull, aes(x = .data[["Age in years"]])) +
  geom_density() +
  labs(
    title = "Age density",
    x = "Age"
  )

plotPercentBar <- function(data, varName, title) {
  ggplot(data, aes(x = .data[[varName]])) +
    geom_bar(aes(y = after_stat(count / sum(count)))) +
    scale_y_continuous(labels = percent_format()) +
    labs(
      x = varName,
      y = "Percentage",
      title = title
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Age group distribution
plotPercentBar(plotDataVimFull,"Age group","Age group distribution")

# Gender distribution
plotPercentBar(plotDataVimFull,"Gender","Gender distribution")

# Race distribution
plotPercentBar(plotDataVimFull,"Race","Race distribution")

# Ethnicity distribution
plotPercentBar(plotDataVimFull,"Ethnicity","Ethnicity distribution")

################################################################################
# Measurements distributions in the full population, without filtering

################################################################################
measurementPlotSpecifications <- list(
  list(var = "Systolic BP",label = "Systolic BP",x_label = "Systolic Blood Pressure (mmHg)"),
  list(var = "Total Cholesterol", label = "Total Cholesterol", x_label = "Total Cholesterol (mg/dL)"),
  list(var = "HDL",label = "HDL", x_label = "HDL (mg/dL)"),
  list(var = "LDL",label = "LDL",x_label = "LDL (mg/dL)"),
  list(var = "BMI",label = "BMI",x_label = expression(BMI~(kg/m^2)))
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
    data = plotDataVimFull,
    var = .x$var,
    label = .x$label,
    x_label = .x$x_label
  )
)

histogramPlots <- purrr::map(
  measurementPlotSpecifications,
  ~ plotMeasurementHistogram(
    data = plotDataVimFull,
    var = .x$var,
    label = .x$label,
    x_label = .x$x_label
  )
)

densityPlots
histogramPlots


################################################################################
# Correlation plot
# Note that in the synthetic data, this probably leads to unstable results, 
# because each variable has 90% missingness and the correlation is based on
# pairwise complete observations
################################################################################

# Correlation plot
corrData <- plotDataVimFull %>%
  select("Age in years", all_of(measurementLabels))

corrMatrix <- cor(
  corrData,
  use = "pairwise.complete.obs"
)

corrplot(
  corrMatrix,
  method = "color",
  type = "upper",
  addCoef.col = "black",
  tl.cex = 0.8,
  number.cex = 0.7
)

################################################################################
# Margin plots for each variable combination
################################################################################
# Set path
marginplotOutputPath <- "C:\\Users\\maud5\\OneDrive\\MSc 2025-2026\\Thesis applications\\Data Figures\\1M data\\test plots"

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

start <- Sys.time()
runMarginplot <- function(spec,
                          data,
                          outputPath = NULL,
                          width = 1600,
                          height = 1600,
                          res = 150) {
  
  plotData <- data %>%
    select(all_of(c(spec$x, spec$y)))
  
  VIM::marginplot(plotData)
  
  if (!is.null(outputPath)) {
    grDevices::png(
      filename = file.path(outputPath, spec$file_name),
      width = width,
      height = height,
      res = res
    )
    on.exit(grDevices::dev.off(), add = TRUE)
    
    VIM::marginplot(plotData)
  }
}

invisible(lapply(
  marginplotSpecs,
  runMarginplot,
  data = plotDataVimFull,
  outputPath = marginplotOutputPath
))

end <- Sys.time()
print(end-start) # approx 11 minutes

################################################################################
# Comparison plots of complete blood pressure data and the full dataset

################################################################################
# Missing BP column
plotDataVimFull <- plotDataVimFull %>%
  mutate(bpMissing = is.na(.data[["Systolic BP"]]))

# Age group distribution comparison plot
ggplot(
  plotDataVimFull,
  aes(x = .data[["Age group"]], fill = bpMissing)
) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(
    values = c(
      "TRUE"= "#F5836E",
      "FALSE" = "#B0F0F7"
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Summarize pct of missingness per age group + total nr of observations per age group
plotDataVimFull %>%
  group_by(.data[["Age group"]]) %>%
  summarise(
    pctMissingBP = mean(is.na(.data[["Systolic BP"]])) * 100,
    n = n(),
    .groups = "drop"
  )

# BP: Full data vs complete case comparison
dataCombined <- bind_rows(
  plotDataVimFull %>%
    mutate(sample = "Full data"),
  plotDataVimCompleteBp %>%
    mutate(sample = "Complete case")
)

# Age density
ggplot(dataCombined,
       aes(x = .data[["Age in years"]],
           #color = sample,
           fill = sample)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(
    values = c(
      "Full data" = "#D699E8",
      "Complete case" = "#82EDC2"
    ),
    labels = c(
      "Full data" = "Full cohort",
      "Complete case" = "Complete case (on blood pressure)"
    )
  ) +
  labs(
    title = "Age density: full data vs complete case",
    x = "Age",
    y = "Density",
    fill = "Dataset"
  ) + guides(color="none")

# Age group density
ggplot(dataCombined, aes(x = sample, fill = .data[["Age group"]])) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(
    values = c(
      "40-44" = "#ED82AD",
      "45-49" = "#F8B591",
      "50-54" = "#E2ED82",
      "55-59" = "#82EDC2",
      "60-64" = "#82E2ED",
      "65-69" = "#8D82ED",
      "70-74" = "#EA9EE3",
      "75-79" = "#E3C5FC"
    )
  ) +
  labs(
    x = NULL,
    y = "Percentage within sample",
    fill = "Age group",
    title = "Age group distribution: full data vs complete case"
  )

# Gender
ggplot(dataCombined, aes(x = sample, fill = .data[["Gender"]])) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(
    values = c(
      "Female" = "#D699E8",
      "Male" = "#82EDC2"
    ),
    na.value = "grey70"
  ) +
  labs(
    x = NULL,
    y = "Percentage within sample",
    fill = "Gender",
    title = "Gender distribution: full data vs complete case"
  )

# Race
ggplot(dataCombined, aes(x = sample, fill = .data[["Race"]])) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(
    values = c(
      "White" = "#82E2ED",
      "Black/African American" = "#F8B591",
      "Asian" = "#E3C5FC",
      "NA" = "grey70"
    ),
    na.value = "grey70"
  ) +
  labs(
    x = NULL,
    y = "Percentage within sample",
    fill = "Race",
    title = "Race distribution: full data vs complete case"
  )

# Ethnicity
ggplot(dataCombined, aes(x = sample, fill = .data[["Ethnicity"]])) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(
    values = c(
      "Hispanic" = "#82E2ED",
      "Non-Hispanic" = "#F8B591"
    ),
    na.value = "grey70"
  ) +
  labs(
    x = NULL,
    y = "Percentage within sample",
    fill = "Ethnicity",
    title = "Ethnicity distribution: full data vs complete case"
  )

# Check complete case on BP
# plotDataVimFull %>%
#   group_by(.data[["Age group"]]) %>%
#   summarise(
#     pctObservedBP = mean(!(is.na(.data[["Systolic BP"]]))) * 100,
#     totalObservedBP = sum(!(is.na(.data[["Systolic BP"]]))),
#     n = n(),
#     .groups = "drop"
#   )
# 
# # Check complete case on all five measures
# plotDataVimFull %>%
#   group_by(.data[["Age group"]]) %>%
#   summarise(
#     pctCompleteAllFive = mean(
#       rowSums(!is.na(across(all_of(measurementLabels)))) == 5
#     ) * 100,
#     n = n(),
#     .groups = "drop"
#   )

