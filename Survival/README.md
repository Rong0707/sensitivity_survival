# Sensitivity Analysis of Treatment Effect to Unmeasured Confounding in Observational Studies with Survival Outcomes

## Rscripts Overview

- SimulateU_surv.R
  - Simulate the unmeasured confounder by stochastics EM method. Output is a list of simulated confounders and the corresponding posterior probabilities.
- surv_sens_var.R and surv_sens_ipw.R
  - Both are used to estimate the treatment effect assuming the existence of simulated confounders.
- emU_surv.R
  - Estimate the treatment effect by EM method. Output is the estimated regression coefficients.
- surv_final.R
  - This function is used to estimate the standard errors of estimates by EM method.

## Example
