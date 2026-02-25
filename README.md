# SEER Breast Cancer Survival Analysis

Time-to-event survival analysis of SEER breast cancer data using Kaplan–Meier estimation, log-rank tests, and Cox proportional hazards regression in SAS.

## Overview

This project evaluates the association between hormone receptor status (Estrogen and Progesterone receptors) and overall survival among female breast cancer patients diagnosed between 2006–2010.

Survival analysis techniques were applied to compare mortality risk across four hormone receptor groups:

- ER+/PR+ (Both Positive)
- ER+/PR−
- ER−/PR+
- ER−/PR− (Both Negative)

## Methods

- Kaplan–Meier survival curves
- Log-rank tests with Bonferroni adjustment
- Multivariable Cox proportional hazards models
- Proportional hazards (PH) assumption testing using Schoenfeld residuals and log-log plots

## Key Findings

- ER+/PR+ patients had the lowest hazard of death.
- ER−/PR− patients had the highest mortality risk.
- Stage, tumor size, lymph node involvement, and age were significant predictors.

## Data

The dataset is not included in this repository.

The SEER breast cancer dataset is publicly available from:
https://ieee-dataport.org/open-access/seer-breast-cancer-data

To reproduce the analysis:
1. Download the dataset.
2. Save it as `data_raw/SEER.csv`.
3. Update the file path in `survival_analysis.sas`.
