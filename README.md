[1752715960425_README_meta_analysis (8).md](https://github.com/user-attachments/files/22355498/1752715960425_README_meta_analysis.8.md)
# Meta-Analysis: Bird Responses to Noise Pollution

## Overview
This repository contains R scripts for a meta-analysis of bird responses to noise pollution, including effect size validation, multi-level modeling, moderator analysis, and publication bias assessment.

## Files
- `meta_analysis_bird_noise.R`: Main analysis script
- `DataExtraction.csv`: Original dataset
- `cleaned_meta_analysis_data.csv`: Cleaned dataset (output)

## How to Run
1. Install all required R packages:
   ```r
   install.packages(c("metafor", "dplyr", "ggplot2", "tidyr", "stringr", "car", "robustbase"))
   ```
2. Run the analysis:
   ```r
   source("meta_analysis_bird_noise.R")
   ```

## Features
- Loads and explores the dataset
- Validates and recalculates Cohen's d effect sizes
- Removes outliers and small-sample studies
- Multi-level and single-level meta-analysis
- Moderator and interaction analyses (habitat, response type, species, noise type)
- Publication bias assessment (funnel plot, Egger's test)
- Influence diagnostics and sensitivity analysis
- Outputs cleaned data and summary statistics

## Interpreting Results
- **Effect size (Cohen's d):**
  - Small: |d| < 0.2
  - Medium: 0.2 ≤ |d| < 0.5
  - Large: |d| ≥ 0.5
- **Heterogeneity (I²):**
  - 0-25%: Low
  - 25-50%: Moderate
  - 50-75%: High
  - 75-100%: Very high
- **Publication bias:**
  - Funnel plot asymmetry and Egger's test (p < 0.05 suggests bias)

## Customization
- Change outlier threshold: `abs(data_clean$cohens_d_recalculated) > 10`
- Change minimum sample size: `data_clean$Sample_size_n < 5`
- Add moderators: see moderator analysis section in the script

## References
- Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2009). Introduction to meta-analysis. John Wiley & Sons.
- Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical Software, 36(3), 1-48.
- Konstantopoulos, S. (2011). Fixed effects and variance components estimation in three-level meta-analysis. Research Synthesis Methods, 2(1), 61-76. 
