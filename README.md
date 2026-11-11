# New-Jersey-Drinking-Water-PFAS-and-Prostate-Cancer-Study
Scripts for XGBoost model creation, training, and tuning to predict PFAS levels in community water systems in NJ and subsequent analysis of PFAS and mortality among prostate cancer patients in NJ
# Directory

| Directory                         | Description                                  |                         
| --------------------------------- |--------------------------------------------- |
|  `README.md` (this file)          |  Directory                                   | 
| `01 scripts model predictors`     | Predictors for XGBoost models: air pollution, tokronov, SES, industries, MIRTA, FUDS, airports, inverse distance        | 
| `02 scripts PFAS prep for XGBoost`| Prepping PFAS data for XGBoost                                         |
| `03 scripts tuning and final models` | Tuning and final XGBoost models and sensitivity analysis of testing train splits                       | 
| `04 scripts maps `                | Create map of predicted PFAS from XGBoost models               | 
| `05 scripts cleaning and prep cohort`| Add time varying variables, apply exclusions, extrapolate air pollution and NDVI, and format final analytic cohort for main analysis and lagged analysis                 | 
| `06 analysis`                     | Table 1, Kaplan Meier curves, multiple imputation, cox models, lagged analysis, splines, inverse probability censoring weights, and effect measure modication                              | 
| `07 qgcomp`                     | Quantile g-computation                            |  


