# =============================================================================
# META-ANALYSIS: Bird Responses to Noise Pollution
# Effect Size Calculation and Validation Script
# =============================================================================

# Load required libraries
library(metafor)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(car)
library(robustbase)

# =============================================================================
# 1. DATA LOADING AND INITIAL EXPLORATION
# =============================================================================

# Load the data
data_raw <- read.csv("DataExtraction.csv", stringsAsFactors = FALSE)

# Summary
cat("Dataset dimensions:", dim(data_raw), "\n")
cat("Number of studies:", length(unique(data_raw$Reference)), "\n")
cat("Number of effect sizes:", nrow(data_raw), "\n")

# =============================================================================
# 2. EFFECT SIZE RECALCULATION AND VALIDATION
# =============================================================================

# Using Cohen's d and sample sizes already calculated
recalculate_cohens_d <- function(row_data) {
  d <- as.numeric(row_data$`Cohens.d.value`)
  n <- as.numeric(row_data$Sample_size_n)
  variance <- NA
  method <- "student"
  if (isTRUE(!is.na(d) && !is.na(n) && n > 0)) {
    variance <- (4/n) + (d^2)/(2*n)
  }
  return(list(
    cohens_d = d,
    variance = variance,
    method = method,
    original_effect = d,
    equation_used = row_data$`Cohens.d.value_equations`
  ))
}

# Apply recalculation to all rows
cat("\nRecalculating effect sizes...\n")
recalculated_effects <- vector("list", nrow(data_raw))
suppressWarnings({
  for (i in seq_len(nrow(data_raw))) {
    recalculated_effects[[i]] <- recalculate_cohens_d(data_raw[i, ])
  }
})

# Summary of conversion results
cat("\nRecalculation Summary:\n")
cat("Total rows processed:", length(recalculated_effects), "\n")

# Extract results
data_raw$cohens_d_recalculated <- vapply(
  recalculated_effects,
  function(x) if (is.null(x$cohens_d) || length(x$cohens_d) == 0) NA_real_ else as.numeric(x$cohens_d),
  numeric(1)
)
data_raw$variance_recalculated <- vapply(
  recalculated_effects,
  function(x) if (is.null(x$variance) || length(x$variance) == 0) NA_real_ else as.numeric(x$variance),
  numeric(1)
)
data_raw$calculation_method <- vapply(
  recalculated_effects,
  function(x) if (is.null(x$method) || length(x$method) == 0) NA_character_ else as.character(x$method),
  character(1)
)
data_raw$original_effect <- vapply(
  recalculated_effects,
  function(x) if (is.null(x$original_effect) || length(x$original_effect) == 0) NA_real_ else as.numeric(x$original_effect),
  numeric(1)
)
data_raw$equation_used <- vapply(
  recalculated_effects,
  function(x) if (is.null(x$equation_used) || length(x$equation_used) == 0) NA_character_ else as.character(x$equation_used),
  character(1)
)

# Summary of extracted results
cat("\nExtraction Summary:\n")
cat("Non-NA Cohens d values:", sum(!is.na(data_raw$cohens_d_recalculated)), "\n")
cat("Non-NA variance values:", sum(!is.na(data_raw$variance_recalculated)), "\n")
cat("Non-NA original effects:", sum(!is.na(data_raw$original_effect)), "\n")

# =============================================================================
# 3 EFFECT SIZE VALIDATION AND COMPARISON
# =============================================================================

# Compare original vs recalculated effect sizes
comparison_data <- data_raw[!is.na(data_raw$cohens_d_recalculated) & 
                           !is.na(data_raw$original_effect), ]

comparison_data$cohens_d_recalculated <- as.numeric(comparison_data$cohens_d_recalculated)
comparison_data$original_effect <- as.numeric(comparison_data$original_effect)

cat("\nEffect size validation summary:\n")
cat("Number of comparable effect sizes:", nrow(comparison_data), "\n")

if (nrow(comparison_data) > 0) {
  comparison_data$difference <- comparison_data$cohens_d_recalculated - comparison_data$original_effect
  comparison_data$abs_difference <- abs(comparison_data$difference)
  
  cat("Mean absolute difference:", mean(comparison_data$abs_difference, na.rm = TRUE), "\n")
  cat("Maximum absolute difference:", max(comparison_data$abs_difference, na.rm = TRUE), "\n")
  
  # Check if differences are zero (as expected)
  if (max(comparison_data$abs_difference, na.rm = TRUE) < 0.001) {
    cat("✓ Validation passed: All differences are essentially zero\n")
  } else {
    cat("⚠ Warning: Some differences detected - investigate if needed\n")
  }
}

# =============================================================================
# 4. DATA CLEANING AND QUALITY CONTROL
# =============================================================================

# Remove rows with missing critical information
data_clean <- data_raw[!is.na(data_raw$Sample_size_n) & 
                      !is.na(data_raw$cohens_d_recalculated), ]

# Remove extreme outliers (effect sizes > 10 or < -10)
extreme_outliers <- data_clean[abs(data_clean$cohens_d_recalculated) > 10, ]
cat("\nExtreme outliers removed:", nrow(extreme_outliers), "\n")

data_clean <- data_clean[abs(data_clean$cohens_d_recalculated) <= 10, ]

# Remove studies with very small sample sizes (< 5)
small_samples <- data_clean[data_clean$Sample_size_n < 5, ]
cat("Studies with very small samples removed:", nrow(small_samples), "\n")

data_clean <- data_clean[data_clean$Sample_size_n >= 5, ]

# Final dataset summary
cat("\nFinal cleaned dataset:\n")
cat("Number of studies:", length(unique(data_clean$Reference)), "\n")
cat("Number of effect sizes:", nrow(data_clean), "\n")
cat("Effect size range:", range(data_clean$cohens_d_recalculated), "\n")

# =============================================================================
# 5. MULTI-LEVEL META-ANALYSIS
# =============================================================================

# Check for multiple effect sizes per study
effect_counts <- table(data_clean$Reference)
studies_with_multiple <- effect_counts[effect_counts > 1]

if (length(studies_with_multiple) > 0) {
  cat("\nStudies with multiple effect sizes:", length(studies_with_multiple), "\n")
  
  # Create study ID for multi-level analysis
  data_clean$study_id <- as.numeric(factor(data_clean$Reference))
  data_clean$effect_id <- 1:nrow(data_clean)
  
  # Primary multi-level meta-analysis
  primary_ma <- rma.mv(yi = cohens_d_recalculated, V = variance_recalculated,
                     random = ~ 1 | study_id / effect_id, data = data_clean)

  cat("\nPrimary Multi-Level Meta-Analysis Results:\n")
  print(primary_ma)

  # Forest plot (single-level for visualization)
  forest_ma <- rma(yi = cohens_d_recalculated, vi = variance_recalculated, 
                 data = data_clean, method = "REML")
  forest(forest_ma, main = "Forest Plot: Bird Responses to Noise Pollution")
}


# =============================================================================
#6. PUBLICATION BIAS ASSESSMENT
# =============================================================================

# Funnel plot
funnel(forest_ma, main = "Funnel Plot for Publication Bias Assessment")

# Egger's test
egger_test <- regtest(forest_ma)
cat("\nEgger's Test for Publication Bias:\n")
print(egger_test)

# =============================================================================
# 7. SENSITIVITY ANALYSES
# =============================================================================

# Single-level meta-analysis (sensitivity check)
single_level_ma <- rma(yi = cohens_d_recalculated, vi = variance_recalculated, 
                       data = data_clean, method = "REML")
cat("\nSensitivity Analysis - Single-Level Meta-Analysis:\n")
print(single_level_ma)

# Influence diagnostics
influence_results <- influence(single_level_ma) # this takes a few minutes to run

# Identify influential studies
influential_studies <- data_clean[influence_results$is.infl, ]
cat("\nInfluential studies identified:", nrow(influential_studies), "\n")

# Print influential studies details
print(influential_studies[, c("Reference", "cohens_d_recalculated", "Sample_size_n", "Habitat_Type", "Species2")])

# Sensitivity analysis without influential studies
data_sensitive <- data_clean[!influence_results$is.infl, ]
sensitive_ma <- rma(yi = cohens_d_recalculated, vi = variance_recalculated, 
                    data = data_sensitive, method = "REML")
cat("\nSensitivity Analysis - Excluding Influential Studies:\n")
print(sensitive_ma)

# =============================================================================
#8. MODERATOR ANALYSES
# =============================================================================

# Habitat type moderator
if (length(unique(data_clean$Habitat_Type)) > 1) {
  habitat_mod <- rma(yi = cohens_d_recalculated, vi = variance_recalculated,
                     mods = ~ Habitat_Type - 1, data = data_clean)
  
  cat("\nHabitat Type Moderator Analysis:\n")
  print(habitat_mod)
}

# Response type moderator
if (length(unique(data_clean$`Response_Type_Communication.Behavior.Reproduction.Physiology`)) > 1) {
  response_mod <- rma(yi = cohens_d_recalculated, vi = variance_recalculated,
                      mods = ~ `Response_Type_Communication.Behavior.Reproduction.Physiology` - 1, 
                      data = data_clean)
  
  cat("\nResponse Type Moderator Analysis:\n")
  print(response_mod)
}

# Species type moderator (using Species2)
if (length(unique(data_clean$Species2)) > 1) {
  species_mod <- rma(yi = cohens_d_recalculated, vi = variance_recalculated,
                     mods = ~ Species2 - 1, data = data_clean)
  
  cat("\nSpecies Type Moderator Analysis:\n")
  print(species_mod)
}

# Independent variable moderator
if (length(unique(data_clean$Independent_Variable_Standardized)) > 1) {
  iv_mod <- rma(yi = cohens_d_recalculated, vi = variance_recalculated,
                mods = ~ Independent_Variable_Standardized - 1, data = data_clean)
  
  cat("\nIndependent Variable Moderator Analysis:\n")
  print(iv_mod)
}

# INTERACTIVE MODERATORs

# 1. Habitat × Response Type Interaction
# Do urban birds show different communication responses than rural birds?

if (length(unique(data_clean$Habitat_Type)) > 1 && length(unique(data_clean$`Response_Type_Communication.Behavior.Reproduction.Physiology`)) > 1) {
  data_clean$habitat_response_interaction <- interaction(
    data_clean$Habitat_Type, 
    data_clean$`Response_Type_Communication.Behavior.Reproduction.Physiology`
  )
  
  habitat_response_mod <- rma(yi = cohens_d_recalculated, vi = variance_recalculated,
                             mods = ~ habitat_response_interaction - 1, data = data_clean)
  
  cat("Habitat × Response Type Interaction Analysis:\n")
  print(habitat_response_mod)
  
  # Show sample sizes for each combination
  cat("\nSample sizes for Habitat × Response Type combinations:\n")
  print(table(data_clean$Habitat_Type, data_clean$`Response_Type_Communication.Behavior.Reproduction.Physiology`))
}

# 2. Species × Habitat Interaction  
# Do songbirds respond differently in different habitats?
if (length(unique(data_clean$Species2)) > 1 && length(unique(data_clean$Habitat_Type)) > 1) {
  data_clean$species_habitat_interaction <- interaction(
    data_clean$Species2, 
    data_clean$Habitat_Type
  )
  
  species_habitat_mod <- rma(yi = cohens_d_recalculated, vi = variance_recalculated,
                             mods = ~ species_habitat_interaction - 1, data = data_clean)
  
  cat("Species × Habitat Interaction Analysis:\n")
  print(species_habitat_mod)
  
  # Show sample sizes for each combination
  cat("\nSample sizes for Species × Habitat combinations:\n")
  print(table(data_clean$Species2, data_clean$Habitat_Type))
}

# 3. Noise Type × Response Type Interaction
# Are certain noise types more disruptive to certain behaviors?
if (length(unique(data_clean$Independent_Variable_Standardized)) > 1 && length(unique(data_clean$`Response_Type_Communication.Behavior.Reproduction.Physiology`)) > 1) {
  data_clean$noise_response_interaction <- interaction(
    data_clean$Independent_Variable_Standardized, 
    data_clean$`Response_Type_Communication.Behavior.Reproduction.Physiology`
  )
  
  noise_response_mod <- rma(yi = cohens_d_recalculated, vi = variance_recalculated,
                           mods = ~ noise_response_interaction - 1, data = data_clean)
  
  cat("Noise Type × Response Type Interaction Analysis:\n")
  print(noise_response_mod)
  
  # Show sample sizes for each combination
  cat("\nSample sizes for Noise Type × Response Type combinations:\n")
  print(table(data_clean$Independent_Variable_Standardized, 
              data_clean$`Response_Type_Communication.Behavior.Reproduction.Physiology`))
}


# =============================================================================
# 9. REPORTING AND EXPORT
# =============================================================================

# Create summary report
report_summary <- list(
  original_effect_sizes = nrow(data_raw),
  recalculated_effect_sizes = sum(!is.na(data_clean$cohens_d_recalculated)),
  studies_included = length(unique(data_clean$Reference)),
  effects_per_study = length(studies_with_multiple),
  primary_effect = primary_ma$b,
  primary_se = primary_ma$se,
  primary_p = primary_ma$pval,
  single_level_effect = single_level_ma$b,
  single_level_se = single_level_ma$se,
  heterogeneity_i2 = single_level_ma$I2,
  heterogeneity_q = single_level_ma$QE,
  publication_bias_p = egger_test$pval
)

cat("\n=== META-ANALYSIS SUMMARY ===\n")
cat("Original effect sizes:", report_summary$original_effect_sizes, "\n")
cat("Recalculated effect sizes:", report_summary$recalculated_effect_sizes, "\n")
cat("Studies included:", report_summary$studies_included, "\n")
cat("Studies with multiple effects:", report_summary$effects_per_study, "\n")
cat("Primary effect size (multi-level):", round(report_summary$primary_effect, 3), "\n")
cat("Primary standard error:", round(report_summary$primary_se, 3), "\n")
cat("Primary p-value:", round(report_summary$primary_p, 10), "\n")
cat("Single-level effect size:", round(report_summary$single_level_effect, 10), "\n")
cat("Single-level standard error:", round(report_summary$single_level_se, 3), "\n")
cat("I² heterogeneity:", round(report_summary$heterogeneity_i2, 1), "%\n")
cat("Q statistic:", round(report_summary$heterogeneity_q, 2), "\n")
cat("Publication bias p-value:", round(report_summary$publication_bias_p, 4), "\n")

# Export cleaned data
write.csv(data_clean, "cleaned_meta_analysis_data.csv", row.names = FALSE)