MAIVE in R: Instructions to the User

This version: July 2024

This readme provides some instructions on the implementation of the MAIVE metaanalysis estimator in R: The meta_maive package with the maive.R function.

Spurious Precision in Meta-Analysis of Observational Research by Zuzana Irsova, Pedro R. D. Bom, Tomas Havranek, and Heiko Rachinger (https://meta-analysis.cz/maive)

0. Installation

Once you have installed R and RStudio Desktop, you can install MAIVE using devtools

  install.packages("devtools")
  
  devtools::install_github("meta-analysis-es/meta_maive")
  
  library(meta_maive)

For the help on the maive.R function 

  help(maive)
  
1. Data

The data should have the following structure:

- Column 1, to be labelled ‘bs’, should contain the primary estimates.
- Column 2, to be labelled ‘sebs’, should contain the standard errors.
- Column 3, to be labelled ‘Ns’, should contain the sample sizes.
- Column 4, to be labelled ‘study id’, should contain a study identification number.
(This column is only necessary if the user chooses clustering or study-level fixed effects. See below.)

2. Options
   
The user needs to specify the options:

- method: PET:1, PEESE:2, PET-PEESE:3, EK:4 (default 3)
  method <- 3
- weighting: default no weight: 0 ; weights: 1, adjusted weights: 2 (default 0)
  weight <- 0
- instrumenting (default 1)
  instrument <- 1 
- correlation at study level: none: 0 (default), fixed effects: 1, cluster: 2
  studylevel <-2
- Anderson-Rubin confidence interval for weak instruments (only for unweighted MAIVE -- PET, PEESE or PET-PEESE): 0 no, 1 yes
  AR <-1

default options are method=3; weight=0; instrument=1; studylevel=0; AR=0 

The default MAIVE meta-estimator is MAIVE-PET-PEESE with instrumented standard errors and no weights. However, the user can choose other options.

The code allows for the user to change the meta-analysis method (PET, PEESE, PET-PEESE, or EK), the weighting (no weights, standard inverse-variance weights, or MAIVE-
adjusted weights), the instrumentation of standard errors (yes or no), and the accounting for study-level correlation (none, study fixed effects, or cluster-robust methods).
- Method: PET=1, PEESE=2, PET-PEESE=3, EK=4. The default method option is PET-PEESE=3. 
- Weighting: no weights=0, inverse-variance weights=1, adjusted weights=2. The default weighting option is no weights=0. 
- Instrumenting the SEs: no=0, yes=1. The default instrumenting option is yes=1.
- Study-level correlation: none=0, study fixed effects=1, cluster-robust standard errors=2. 
The default study-level correlation option is none=0. 
- Anderson-Rubin confidence interval (for weak instruments) for the meta-estimate: AR=0, no (default); AR=1, yes. 
This option is available for the unweighted MAIVE – PET, PEESE and PET-PEESE versions. It is not available for study fixed effects.
Its use is recommendable if the first stage F-test statistic lies between 10 and 100 (see Keane and Neal, 2023, for further details). 
Note that obtaining the AR confidence interval takes some time.

3. Output
   
The code returns:
- A MAIVE point estimate and a standard error.
- A point estimate and standard error from the method chosen in option
- A Hausman-type test statistic and a 5% critical value. The test statistic consists of a weighted squared difference between the MAIVE point estimate and a standard
point estimate. The latter is obtained from the same method as chosen for MAIVE but with inverse-variance weights and without instrumenting the standard errors.
- When instrumenting the standard errors, the code also returns a heteroskedasticity-robust F test of the first-stage regression.
- When choosing the AR option, the code also returns a Anderson-Rubin confidence interval for weak instruments.
- The instrumented standard errors are saved as SE instrumented and can be obtained as “MAIVE$SE_instrumented”.

4. Technical comments
- In the first-stage, we regress the variances on a constant and the inverse sample sizes.
- The Hausman-type test is weighted by the variance of the MAIVE estimator and is, consequently, a conservative test. It only compares the intercepts.
- If study fixed effects are included, the dummy variables are demeaned so that the intercept measures a grand mean.
- If the user does not provide a ‘study-id’ column, the program assumes that the study-level correlation option is none.
- The Anderson-Rubin confidence interval is obtained in the standard way (see Keane and Neal, 2023, for further details).

References
Keane, Michael and Neal, Timothy. 2023. “Instrument strength in IV estimation and inference: A guide to theory and practice”, Journal of Econometrics, 235, 2, 1625-1653.
https://doi.org/10.1016/j.jeconom.2022.12.009.
