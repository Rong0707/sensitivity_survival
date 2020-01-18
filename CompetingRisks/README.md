# Sensitivity Analysis of Treatment Effect to Unmeasured Confounding in Observational Studies with Competing Risks Outcomes

## Rscripts Overview

- SimulateU.R
  - Simulate the unmeasured confounder by stochastics EM method. Output is a list of simulated confounders and the corresponding posterior probabilities.
- competing_sens_var.R and competing_sens_ipw.R
  - Both are used to estimate the treatment effect assuming the existence of simulated confounders via regresssion adjustment and inverse probability weighting, respectively.
- emU.R
  - Estimate the treatment effect by EM method. Output is the estimated regression coefficients.
- competing_final.R
  - This function is used to estimate the standard errors of estimates by EM method.


## Example

- IBD_CD.R
  - This is how we performed sensitivity analysis on the IBD CD data.
