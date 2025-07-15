MAIVE in R: Instructions to the User
================
July 2025

This readme provides instructions on the implementation of the MAIVE
meta-analysis estimator in R: The `MAIVE` package with the
`maive.r` function.

Spurious Precision in Meta-Analysis of Observational Research by Zuzana
Irsova, Pedro R. D. Bom, Tomas Havranek, and Heiko Rachinger  
<https://meta-analysis.cz/maive>

------------------------------------------------------------------------

## Installation

Once you have installed R and RStudio Desktop, you can install MAIVE
using **devtools**:

``` r
install.packages("devtools")
devtools::install_github("meta-analysis-es/maive")
library(MAIVE)
```

For help on the `maive.r` function:

``` r
help(maive)
```

## Data Structure

The data should have the following structure:

| Column | Label | Description |
|----|----|----|
| 1 | bs | Primary estimates |
| 2 | sebs | Standard errors |
| 3 | Ns | Sample sizes |
| 4 | study id | Study identification number (needed for clustering/study fixed effects) |

## Options

The user needs to specify the following options (with defaults in
parentheses):

| Option | Description | Values |
|----|----|----|
| method | Meta-analysis method | PET=1, PEESE=2, PET-PEESE=3 (default), EK=4 |
| weighting | Weighting scheme | No weights=0 (default), Weights=1, Adjusted weights=2 |
| instrumenting | Instrument standard errors | No=0, Yes=1 (default) |
| studylevel | Study-level correlation | None=0, Fixed effects=1, Cluster=2 (default), Fixed effects and Cluster=3|
| SE | Standard errors | CR0 (Huber–White)=0, CR1 (Standard empirical correction)=1, CR2 (Bias-reduced estimator)=2, wild bootstrap=3 (default)|
| AR | Anderson-Rubin confidence interval (only for unweighted) | No=0 , Yes=1 (default)|

## Description of Default Settings

The default MAIVE meta-estimator is MAIVE-PET-PEESE with instrumented
standard errors, no weights, cluster and wild bootstrap. However, the user can adjust:

- The meta-analysis method (PET, PEESE, PET-PEESE, EK).
- The weighting (no weights, inverse-variance weights, or MAIVE-adjusted
  weights).
- Instrumentation of standard errors (yes or no).
- Accounting for study-level correlation (none, fixed effects, 
  cluster-robust methods, or fixed effects and cluster-robust methods).
- Use CR0 (Huber–White), CR1 (Standard empirical correction)=1, or CR2 (Bias-reduced estimator) for the estimation of the standard errors

## Output

The function returns:

- A MAIVE point estimate and standard error.
- A point estimate and standard error from the method chosen.
- A Hausman-type test statistic and a 5% critical value.
- If instrumenting SEs, a heteroskedasticity-robust F test of the
  first-stage regression.
- If AR option is chosen, an Anderson-Rubin confidence interval for weak
  instruments.
- Instrumented standard errors saved as `MAIVE$SE_instrumented`.
- p-value of test for publication bias / p-hacking based on instrumented FAT.

## Technical Comments

- The first stage regresses variances on a constant and inverse sample
  sizes.
- The Hausman-type test is weighted by the variance of the MAIVE
  estimator (conservative test).
- Study fixed effects are demeaned so intercept measures a grand mean.
- If no study-id column is provided, the program assumes no study-level
  correlation.
- The Anderson-Rubin confidence interval follows Keane and Neal (2023).

## Required Packages

The following R packages might be needed:

- clubSandwich
- metafor
- pracma
- sandwich
- varhandle

You can install them with:

``` r
install.packages(c("clubSandwich","metafor","pracma", "sandwich","varhandle"))

```

## References

Keane, Michael and Neal, Timothy. 2023. “Instrument strength in IV
estimation and inference: A guide to theory and practice”, Journal of
Econometrics, 235(2), 1625-1653.
<https://doi.org/10.1016/j.jeconom.2022.12.009>

Keane, Michael and Neal, Timothy. 2023. “Instrument strength in IV estimation and inference: A guide to theory and practice”, Journal of Econometrics, 235(2), 1625-1653. https://doi.org/10.1016/j.jeconom.2022.12.009

Tipton, E. (2015). Small sample adjustments for robust variance estimation with cluster-correlated data. Psychological Methods, 20(3), 375–389. https://doi.org/10.1037/met0000019
