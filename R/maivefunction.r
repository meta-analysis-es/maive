#' R code for MAIVE
#'
#' R package for MAIVE: "Spurious Precision in Meta-Analysis of Observational Research" by
#' Zuzana Irsova, Pedro Bom, Tomas Havranek, and Heiko Rachinger.
#'
#' @param dat Data frame with columns bs, sebs, Ns, study_id (optional).
#' @param method 1 FAT-PET, 2 PEESE, 3 PET-PEESE, 4 EK.
#' @param weight 0 no weights, 1 standard weights, 2 adjusted weights.
#' @param instrument 1 yes, 0 no.
#' @param studylevel Correlation at study level: 0 none, 1 fixed effects, 2 cluster.
#' @param SE SE estimator: 0 CR0 (Huber–White), 1 CR1 (Standard empirical correction),
#' 2 CR2 (Bias-reduced estimator), 3 wild bootstrap.
#' @param AR Anderson Rubin corrected CI for weak instruments (only for unweighted MAIVE versions
#' of PET, PEESE, PET-PEESE, not available for fixed effects): 0 no, 1 yes.
#'
#' @details Data \code{dat} can be imported from an Excel file via:
#' \code{dat <- read_excel("inputdata.xlsx")} or from a csv file via: \code{dat <- read.csv("inputdata.csv")}
#' It should contain:
#' \itemize{
#'   \item Estimates: bs
#'   \item Standard errors: sebs
#'   \item Number of observations: Ns
#'   \item Optional: study_id
#' }
#' Default option for MAIVE: MAIVE-PET-PEESE, unweighted, instrumented, cluster SE, wild bootstrap, AR.
#'
#' @return \itemize{
#'   \item MAIVE meta-estimate and standard error
#'   \item Hausman type test: comparison between MAIVE and standard version
#'   \item When instrumenting: heteroskedastic robust F-test of the first step instrumented SEs
#'   \item p-value of test for publication bias / p-hacking based on instrumented FAT
#' }
#'
#' @export
maive <- function(dat, method, weight, instrument, studylevel, SE, AR) {

  # Manual wild cluster bootstrap function:
  # wild bootstrap
  # clustered (weights drawn per cluster, not per observation)
  # with Rademacher distribution (i.e., ±1 with 0.5 probability each)

  manual_wild_cluster_boot_se <- function(model, data, cluster_var, B = 500, seed = 123) {
    set.seed(seed)

    # Extract residuals and fitted values
    resids <- residuals(model)
    fitted_vals <- fitted(model)

    # Get cluster IDs
    clusters <- unique(as.character(data[[cluster_var]]))
    G <- length(clusters)

    # Coefficient names
    coef_names <- names(coef(model))
    k <- length(coef_names)

    # Matrix to store bootstrap coefficients
    boot_coefs <- matrix(NA, nrow = B, ncol = k)
    colnames(boot_coefs) <- coef_names

    # Loop over bootstrap replications
    for (b in 1:B) {
      # Draw Rademacher multipliers per cluster
      u_g <- sample(c(-1, 1), size = G, replace = TRUE)
      names(u_g) <- as.character(clusters)

      # Create bootstrap outcome
      data$y_boot <- fitted_vals + resids * u_g[as.character(data[[cluster_var]])]

      # Refit the same model formula on bootstrap sample
      form_boot <- update(formula(model), y_boot ~ .)
      fit_boot <- lm(form_boot, data = data)

      # Store bootstrap coefficients
      boot_coefs[b, ] <- coef(fit_boot)
    }

    # Compute bootstrap SEs
    boot_se <- apply(boot_coefs, 2, sd)

    # Compute bootstrap percentile CI
    alpha <- 0.05
    boot_ci <- t(apply(boot_coefs, 2, function(x) quantile(x, probs = c(alpha/2, 1 - alpha/2))))
    colnames(boot_ci) <- c("lower", "upper")

    # Return list
    return(list(
      boot_se = boot_se,
      boot_ci = boot_ci,
      boot_coefs = boot_coefs
    ))
  }
  # end of manual_wild_cluster_boot_se function

  methods <- c("PET","PEESE","PET-PEESE","EK")
  instrumented <- c("not instrumented","instrumented")
  weighted <- c("no weights","standardly weighted", "adjusted weights")
  studylevelcorrelation <- c("none","study level dummies", "cluster")

  if (studylevel==0){
    cluster<-0
    dummy<-0
  } else if (studylevel==1){
    cluster<-0
    dummy<-1
  } else if (studylevel==2) {
    cluster<-1
    dummy<-0
  } else if (studylevel==3) {
    cluster<-1
    dummy<-1
  }
  type_map <- c("CR0", "CR1", "CR2")
  if (SE<3){
    type_choice <- type_map[SE + 1]
  } else if (SE==3){               # not needed if bootstrap
    type_choice <- type_map[0 + 1]
  }

  # AR not available for EK, for fixed effects, and for weighted
  if (method==4|weight==1|weight==2|instrument==0) {     AR<-0
  }

  # extracting data from excel
  dat = as.data.frame(dat)
  bs<-dat[,1]
  M<-length(bs)
  sebs<-dat[,2]
  Ns<-dat[,3]

  if (dim(dat)[2]==4){
    studyid<-dat[,4]
  } else {
    studyid<-(1:M)
    dummy<-0
    cluster<-0
  }

  alpha_s <- 0.05

  # create Dummies from studyid
  df<-data.frame(studyid)
  D <- to.dummy(df,"studyid")
  D <-D-matrix(colMeans(D),nrow=M,ncol=size(D)[2], byrow = TRUE)
  D <- D[,1:(dim(D)[2]-1)]

  # g=studyid if clustered and g=(1:M)' if not clustered (gives heteroskedastic robust SE)
  if (cluster==0){
    g <- (1:M)
  } else if (cluster==1){
    g <- studyid
  }

  dat$g <- g

  # (1) Instrumenting the variances with sample size allowing for a constant and including dummies
  invNs<- 1/Ns
  sebs2<- sebs^2
  Xiv<-matrix(c(ones(M,1)[,1],invNs),nrow=M)
  varreg1 <- lm(sebs2~ 0+Xiv)
  dimiv<-2
  if (varreg1$coefficients[1]<0){
    Xiv<-invNs
    varreg1 <- lm(sebs2~ 0+Xiv)
    dimiv<-1
  }

  sebs2fit1 <- varreg1$fitted.values

  # F-statistic of first step. heteroskedasticity and autocorrelation robust variance HAC
  F_hac <- (varreg1$coefficients[dimiv]^2 /vcovCR(varreg1, cluster = g, type = type_choice)[dimiv,dimiv])

  #weight
  if (weight==0){
    w <- ones(M,1)[,1]
  } else if (weight==1){
    w <- sebs
  } else if (weight==2){
    w <- sebs2fit1^(1/2)
  }

  #instrument
  if (instrument==0){
    x <- sebs
    x2 <- sebs^2
    F_hac <-"NA"
  } else if (instrument==1){
    x <- sebs2fit1^(1/2)
    x2 <- sebs2fit1
    F_hac<-round(F_hac,3)
  }

  #choose dependent variable and regressor
  y <- bs/w
  x <- x
  x2 <- x2
  X <- matrix(c(ones(M,1)[,1], x)/w,nrow=M)
  X_d <- matrix(c(X, D/w), nrow=M)
  X2 <- matrix(c(ones(M,1)[,1], x2)/w,nrow=M)
  X2_d <- matrix(c(X2, D/w), nrow=M)

  # baseline, i.e. chosen method, with chosen options of study-level correlation
  #  but with inverse-variance weighting and without instrumenting
  y0 <- bs/sebs
  x0 <- sebs
  x20 <- sebs ^2
  X0 <- matrix(c(ones(M,1)[,1], x0)/sebs, nrow=M)
  X0_d <- matrix(c(X0, D/sebs), nrow=M)
  X20 <- matrix(c(ones(M,1)[,1], x20)/sebs, nrow=M)
  X20_d <- matrix(c(X20, D/sebs),  nrow=M)

  if (dummy==0){
    X <- X
    X0 <- X0
    X2 <- X2
    X20 <- X20
    cD <- ones(M,1)[,1]
  } else if (dummy==1){
    X <- X_d
    X0 <- X0_d
    X2 <- X2_d
    X20 <- X20_d
    cD <- matrix(c(ones(M,1)[,1],D),nrow=M)  # for EK stack constant and dummies
  }

  cD<- cD/w
  cD0<- cD/sebs

  ones_w<-ones(M,1)[,1]/w
  ones_w0<-ones(M,1)[,1]/sebs

  # Fixed effects (FE)
  wis0 <- 1/(w^2)
  fe <- sum(bs*wis0)/sum(wis0)
  varfe <- 1/sum(wis0)
  # baseline
  wis00 <- 1/(sebs^2)
  fe0 <- sum(bs*wis00)/sum(wis00)
  varfe0 <- 1/sum(wis00)

  # WLS
  wlsreg <-lm(y~ 0+ cD )
  wls <- wlsreg$coefficients[1]
  wlsse <- sqrt(vcovCR(wlsreg,cluster=g, type = type_choice)[1,1])
  # baseline
  wlsreg0 <-lm(y0~ 0+ cD0 )
  wls0 <- wlsreg0$coefficients[1]
  wlsse0 <- sqrt(vcovCR(wlsreg0, cluster=g, type = type_choice)[1,1])

  # FAT-PET - MAIVE
  fatpet <- lm(y~ 0+X)
  # FAT-PET - baseline case
  fatpet0 <- lm(y0~ 0+X0)

  # PEESE - MAIVE
  peese <- lm(y~0+X2)
  # PEESE - baseline case
  peese0 <- lm(y0~0+X20)

  # PET-PEESE - MAIVE
  if (abs(fatpet$coefficients[1]/sqrt(vcovCR(fatpet,cluster=g, type = type_choice)[1,1]))>qt(1-alpha_s/2,M-dim(X)[2]-1)){
    petpeese <- peese
  } else {
    petpeese <- fatpet
  }
  # PET-PEESE - baseline case
  if (abs(fatpet0$coefficients[1]/sqrt(vcovCR(fatpet0,cluster=g, type = type_choice)[1,1]))>qt(1-alpha_s/2,M-dim(X0)[2]-1)){
    petpeese0 <- peese0
  } else {
    petpeese0 <- fatpet0
  }

  # True effect variance - MAIVE
  Qfe0 <- sum(wlsreg$residuals*wlsreg$residuals)
  sigh2hat0 <- max(0,M*((Qfe0/(M-dim(wlsreg$model)[2]-1))-1)/sum(wis0))
  sighhat0 <- sqrt(sigh2hat0)
  #True effect variance - baseline
  Qfe00 <- sum(wlsreg0$residuals*wlsreg0$residuals)
  sigh2hat00 <- max(0,M*((Qfe00/(M-dim(wlsreg0$model)[2]-1))-1)/sum(wis00))
  sighhat00 <- sqrt(sigh2hat00)

  # Endogenous Kink (EK) Threshold- MAIVE
  if (petpeese$coefficients[1] > 1.96*sighhat0){
    a0 <- (petpeese$coefficients[1]-1.96*sighhat0)*(petpeese$coefficients[1]+1.96*sighhat0)/(2*1.96*petpeese$coefficients[1])
  } else {
    a0 <- 0
  }
  # Endogenous Kink (EK) Threshold - baseline
  if (petpeese0$coefficients[1] > 1.96*sighhat00){
    a00 <- (petpeese0$coefficients[1]-1.96*sighhat00)*(petpeese0$coefficients[1]+1.96*sighhat00)/(2*1.96*petpeese0$coefficients[1])
  } else {
    a00 <- 0
  }

  # EK - MAIVE
  if (a0>min(x)  && a0<max(x)){
    xx_w=(x-a0)*(x>a0)/w
    ekreg <- lm(y~ 0+cD+xx_w)
  } else if (a0<min(x)){
    x_w=x/w
    ekreg <- lm(y~ 0+cD+x_w)
  } else if (a0>max(x)){
    ekreg <- lm(y~ 0+cD)
  }
  ek <- ekreg$coefficients[1]

  # EK - baseline
  if (a00>min(x0)  && a00<max(x0)){
    xx0_w=(x0-a00)*(x0>a00)/sebs
    ekreg0 <- lm(y0~ 0+cD0+xx0_w)
  } else if (a00<min(x0)){
    x0_w=x0/sebs
    ekreg0 <- lm(y0~ 0+cD0+x0_w )
  } else if (a00>max(x0)){
    ekreg0 <- lm(y0~ 0+cD0 )
  }
  ek0 <- ekreg0$coefficients[1]


  # Anderson and Rubin Confidence intervals
  if (AR == 1) {

    compute_AR_CI_fast <- function(model, adjust_fun, bs, sebs, invNs, g, type_choice) {

      # Beta estimates and robust SEs
      beta0 <- model$coefficients[1]
      beta0se <- sqrt(vcovCR(model, cluster = g, type = type_choice)[1, 1])
      l0 <- beta0 - 5 * beta0se
      u0 <- beta0 + 5 * beta0se
      pr0 <- max(100, round(100 / (u0 - l0)))

      beta1 <- model$coefficients[2]
      beta1se <- sqrt(vcovCR(model, cluster = g, type = type_choice)[2, 2])
      l1 <- beta1 - 5 * beta1se
      u1 <- beta1 + 5 * beta1se
      pr1 <- max(100, round(100 / (u1 - l1)))

      b0_grid <- seq(l0, u0, by = 1/pr0)
      b1_grid <- seq(l1, u1, by = 1/pr1)

      M <- length(bs)
      ones_vec <- rep(1, M)
      Z <- cbind(ones_vec, invNs)
      PZ <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
      MZ <- diag(M) - PZ

      # Make full grid of b0 and b1
      grid <- expand.grid(b0 = b0_grid, b1 = b1_grid)

      # For each row in grid, compute AR stat
      compute_stat <- function(row) {
        b0 <- row["b0"]
        b1 <- row["b1"]

        bs_star <- adjust_fun(bs, b0, b1, sebs)

        num <- as.numeric(crossprod(bs_star, PZ %*% bs_star))
        denom <- as.numeric(crossprod(bs_star, MZ %*% bs_star))

        stat <- (M - 2) * num / denom
        stat
      }

      AR_stats <- apply(grid, 1, compute_stat)

      # Check which are accepted
      AR_accept <- matrix(AR_stats < 5.99, nrow = length(b0_grid), ncol = length(b1_grid), byrow = FALSE)

      b0_CI_all <- b0_grid[rowSums(AR_accept) > 0]
      b1_CI_all <- b1_grid[colSums(AR_accept) > 0]

      b0_CI <- c(b0_CI_all[1], b0_CI_all[length(b0_CI_all)])
      b1_CI <- c(b1_CI_all[1], b1_CI_all[length(b1_CI_all)])

      list(
        b0_CI = round(b0_CI, 3),
        b1_CI = round(b1_CI, 3)
      )
    }

    # PET adjustment
    PET_adjust <- function(bs, b0, b1, sebs) bs - b0 - b1 * sebs

    # PEESE adjustment
    PEESE_adjust <- function(bs, b0, b1, sebs) bs - b0 - b1 * sebs^2

    # Compute AR CIs for PET
    pet_res <- compute_AR_CI_fast(
      model = fatpet, adjust_fun = PET_adjust,
      bs = bs, sebs = sebs, invNs = invNs,
      g = g, type_choice = type_choice
    )
    b0_CI_AR_PET <- pet_res$b0_CI
    b1_CI_AR_PET <- pet_res$b1_CI

    # Compute AR CIs for PEESE
    peese_res <- compute_AR_CI_fast(
      model = peese, adjust_fun = PEESE_adjust,
      bs = bs, sebs = sebs, invNs = invNs,
      g = g, type_choice = type_choice
    )
    b0_CI_AR_PEESE <- peese_res$b0_CI
    b1_CI_AR_PEESE <- peese_res$b1_CI


    # PET-PEESE without weights and with AR CI
    # if PET does not contain 0, do PEESE.
    if (b0_CI_AR_PET[1]>0) {
      b0_CI_AR_PP<-b0_CI_AR_PEESE
    } else if (b0_CI_AR_PET[1]<=0) {
      b0_CI_AR_PP<-b0_CI_AR_PET
    }
  } else if (AR==0) {
    b0_CI_AR_PET= "NA"
    b0_CI_AR_PEESE= "NA"
    b0_CI_AR_PP= "NA"
  }

  "RESULTS"
  if (method==1){
    "MAIVE-FAT-PET"
    beta=fatpet$coefficients[1]
    betase=sqrt(vcovCR(fatpet,cluster=g, type = type_choice)[1,1])
    if (SE==3){
      boot_result <- manual_wild_cluster_boot_se(model = fatpet,data = dat,cluster_var = "g",B = 500)
      betase <- boot_result$boot_se[1]
    }
    "Standard FAT-PET"
    beta0=fatpet0$coefficients[1]
    beta0se=sqrt(vcovCR(fatpet0,cluster=g, type = type_choice)[1,1])
    if (SE==3){
      boot_result <- manual_wild_cluster_boot_se(model = fatpet0,data = dat,cluster_var = "g",B = 500)
      beta0se <- boot_result$boot_se[1]
    }
    "Hausman-type test"
    Hausman=(fatpet$coefficients[1]-fatpet0$coefficients[1])^2/(vcovCR(fatpet,cluster=g, type = type_choice)[1,1])
    Chi2=qchisq(p=0.05, df=1, lower.tail=FALSE)
    "AR-CI"
    b0_CI_AR=b0_CI_AR_PET
  } else if (method==2){
    "MAIVE-PEESE"
    beta=peese$coefficients[1]
    betase=sqrt(vcovCR(peese,cluster=g, type = type_choice)[1,1])
    if (SE==3){
      boot_result <- manual_wild_cluster_boot_se(model = peese,data = dat,cluster_var = "g",B = 500)
      betase <- boot_result$boot_se[1]
    }
    "Standard PEESE"
    beta0=peese0$coefficients[1]
    beta0se=sqrt(vcovCR(peese0,cluster=g, type = type_choice)[1,1])
    if (SE==3){
      boot_result <- manual_wild_cluster_boot_se(model = peese0,data = dat,cluster_var = "g",B = 500)
      beta0se <- boot_result$boot_se[1]
    }
    "Hausman-type test"
    Hausman=(peese$coefficients[1]-peese0$coefficients[1])^2/(vcovCR(peese,cluster=g, type = type_choice)[1,1])
    Chi2=qchisq(p=0.05, df=1, lower.tail=FALSE)
    "AR-CI"
    b0_CI_AR=b0_CI_AR_PEESE
  } else if (method==3){
    "MAIVE-PET-PEESE"
    beta=petpeese$coefficients[1]
    betase=sqrt(vcovCR(petpeese,cluster=g, type = type_choice)[1,1])
    if (SE==3){
      boot_result <- manual_wild_cluster_boot_se(model = petpeese,data = dat,cluster_var = "g",B = 500)
      betase <- boot_result$boot_se[1]
    }

    "Standard PET-PEESE"
    beta0=petpeese0$coefficients[1]
    beta0se=sqrt(vcovCR(petpeese0,cluster=g, type = type_choice)[1,1])
    if (SE==3){
      boot_result <- manual_wild_cluster_boot_se(model = petpeese0,data = dat,cluster_var = "g",B = 500)
      beta0se <- boot_result$boot_se[1]
    }
    "Hausman-type test"
    Hausman=(petpeese$coefficients[1]-petpeese0$coefficients[1])^2/(vcovCR(petpeese,cluster=g, type = type_choice)[1,1])
    Chi2=qchisq(p=0.05, df=1, lower.tail=FALSE)
    "AR-CI"
    b0_CI_AR=b0_CI_AR_PP
  } else if (method==4){
    "MAIVE-EK"
    beta=ekreg$coefficients[1]
    betase=sqrt(vcovCR(ekreg,cluster=g, type = type_choice)[1,1])
    if (SE==3){
      boot_result <- manual_wild_cluster_boot_se(model = ekreg,data = dat,cluster_var = "g",B = 500)
      betase <- boot_result$boot_se[1]
    }
    "Standard EK"
    beta0=ekreg0$coefficients[1]
    beta0se=sqrt(vcovCR(ekreg0,cluster=g, type = type_choice)[1,1])
    if (SE==3){
      boot_result <- manual_wild_cluster_boot_se(model = ekreg0,data = dat,cluster_var = "g",B = 500)
      beta0se <- boot_result$boot_se[1]
    }
    "Hausman-type test" # with variance of MAIVE in denominator (instead of the difference) hence is conservative
    Hausman=(ekreg$coefficients[1]-ekreg0$coefficients[1])^2/(vcovCR(ekreg,cluster=g, type = type_choice)[1,1])
    Chi2=qchisq(p=0.05, df=1, lower.tail=FALSE)
    "AR-CI"
    b0_CI_AR="NA"
  }

  if (weight==1||weight==2){
    "AR-CI"
    b0_CI_AR="NA"
  }

  "p-value of test for publication bias / p-hacking based on instrumented FAT"
  pb_p=summary(fatpet)$coefficients[2, 4]

  my_list <- list("beta"=round(beta,3), "SE"=round(betase,3),"F-test"=F_hac,"beta_standard"=round(beta0,3),"SE_standard"=round(beta0se,3),"Hausman"=round(Hausman,3), "Chi2"=round(Chi2,3), "SE_instrumented"=sebs2fit1^(1/2), "AR_CI"=b0_CI_AR, "pub bias p-value"=round(pb_p,3))
  return(my_list)
}
