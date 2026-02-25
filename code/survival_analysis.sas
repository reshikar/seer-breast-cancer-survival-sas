# SEER Breast Cancer Survival Analysis

**Author:** Reshika Rimal  
**Methods:** Kaplan–Meier, Log-rank, Cox PHREG, Proportional Hazards diagnostics  
**Software:** SAS  
## Overview

This project evaluates the association between hormone receptor status (ER/PR) and overall survival among breast cancer patients using SEER data (2006–2010).

Kaplan–Meier survival curves and Cox proportional hazards models were used to assess mortality risk across four receptor groups.

## Statistical Methods

- Kaplan–Meier survival estimation  
- Log-rank tests  
- Cox proportional hazards regression  
- PH assumption testing using Schoenfeld residuals  

## Key Findings

- ER+/PR+ patients had the lowest hazard of death.
- ER–/PR– patients had the highest mortality risk.
- Stage, tumor size, lymph node involvement, and age were significant predictors

## Repository Structure

- `survival_analysis.sas` — main SAS script
- `docs/` — project presentation slides


## Data

The SEER dataset is not included.  
To reproduce the analysis, obtain the dataset and update the file path in the SAS script.
