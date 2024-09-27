#-------------------------------
# 0-functions.R
# shared functions
#-------------------------------

#-------------------------------
# Trachoma EU categories based on pre-selected criteria
#
# @df : a data.frame with the "study_eu" variable in it
#       study_eu is a concatenation of study_id and "EU-" and eu variables in the clean data
#
# returns the data.frame with a new variable named trachoma_cat
#-------------------------------

unpublished_EUs <- c("Niger-Ilela-2022", "Niger-Malbaza-2022", "Niger-Bagaroua-2022", 
                     "DRC-Manono-2018", "DRC-Nyemba-2018")

# 2 main categories - programmatic action:
eus_wo_PCR <- c("Sudan-El Seraif-2019", "Sudan-Saraf Omrah-2019", "Sudan-Kotom-2019",
                  "Malaysia-Sabah-2015",  "Peru-Amazonia-2020", "Papua New Guinea-Mendi-2015", 
                "Papua New Guinea-Daru-2015", "Papua New Guinea-West New Britain-2015")

# note - `eus_wo_PCR` are placed in the "Unclassified" category 

make_trachoma_cat <- function(df) {
  df <- df %>% 
    mutate(trachoma_cat = case_when(startsWith(as.character(eu_name), "Ghana") ~ "Action not needed",
                                    eu_name %in% c("Ethiopia-Alefa-2017","Niger-MORDOR/Dosso-2015","Ethiopia-Woreta Town-2017",
                                              "Togo-Anie-2017","Togo-Keran-2017","Morocco-Agdaz-2019","Morocco-Boumalne Dades-2019",
                                              "Gambia-River Regions-2014","Ethiopia-Metema-2021","Ethiopia-Woreta Town-2021",
                                              "Malawi-DHO Nkwazi-2014","Malawi-Kasisi/DHO-2014","Malawi-Luzi Kochilira-2014",
                                              "Malawi-Chapananga-2014") ~ "Action not needed",
                                    eu_name %in% c("Ethiopia-Andabet-2017","Ethiopia-Goncha-2019","Niger-PRET/Matameye-2013",
                                              "Ethiopia-Ebinat-2019","Ethiopia-TAITU/Wag Hemra-2018","Ethiopia-WUHA/Wag Hemra-2016",
                                              "Kiribati-Tarawa-2016","Kiribati-Kiritimati-2016","Tanzania-Kongwa-2018",
                                              "Tanzania-Kongwa-2013","Solomon Islands-Temotu/Rennel/Bellona-2015") ~ "Action needed",
                                    eu_name %in% eus_wo_PCR ~ "Unclassified",
                                    TRUE ~ as.character("Unclassified")),
           trachoma_cat = factor(trachoma_cat, levels = c("Action needed","Action not needed","Unclassified"))
    )
  
  return(df)
}


# 3 main categories - elimination of transmission:
eus_wo_PCR2 <- c("Sudan-El Seraif-2019","Sudan-Saraf Omrah-2019","Sudan-Kotom-2019",
                "Malaysia-Sabah-2015","Peru-Amazonia-2020","Papua New Guinea-Mendi-2015","Papua New Guinea-Daru-2015","Papua New Guinea-West New Britain-2015",
                "Togo-Anie-2017","Togo-Keran-2017","Morocco-Agdaz-2019","Morocco-Boumalne Dades-2019","Gambia-River Regions-2014")


make_trachoma_cat2 <- function(df) {
  df <- df %>% mutate(trachoma_cat = case_when(startsWith(as.character(eu_name), "Ghana") ~ "Elimination",
                                               eu_name %in% c("Ethiopia-Alefa-2017","Niger-MORDOR/Dosso-2015","Ethiopia-Woreta Town-2017","Ethiopia-Metema-2021","Ethiopia-Woreta Town-2021") ~ "Elimination",
                                               eu_name %in% c("Ethiopia-Andabet-2017","Niger-PRET/Matameye-2013","Ethiopia-TAITU/Wag Hemra-2018","Ethiopia-WUHA/Wag Hemra-2016","Kiribati-Tarawa-2016","Kiribati-Kiritimati-2016",
                                                              "Ethiopia-Ebinat-2019") ~ "Endemic",
                                               eu_name %in% eus_wo_PCR2 ~ "Unclassified",
                                  TRUE ~ as.character("Near elimination")),
                      trachoma_cat = factor(trachoma_cat, levels = c("Endemic", "Near elimination", 
                                                                     "Elimination", "Unclassified"))
  )
  
  return(df)
}

#-------------------------------
# loglik_sis()
# log-likelihood function for reversible catalytic model with fixed seroreversion rate
# akin to a susceptible-infected-susceptible (SIS) model
# 
# this LL function must then be optimized using R's optim() routine
# 
# @h    : seroconverion rate (lambda), to be estimated
# @r    : seroreversion rate (rho), fixed — should be specified
# @data : data.frame with 1 row per age group, with columns "age" "n" "k"
#         age: is age in years, 
#         n  : is the number of children measured
#         k  : is the number of children seropositive
#-------------------------------
loglik_sis <- function(h,r,data) {
  # h    : \lambda seroconversion rate (to be estimated)
  # r    : \rho seroreversion rate (fixed)
  # data  : data frame with 1 row per age group. cols = age / n / k
  h <- rep(h,nrow(data))
  r <- rep(r,nrow(data))
  t <- data[,1]
  n <- data[,2]
  k <- data[,3]
  p <- h/(r+h)*(1-exp(-(r+h)*t))
  # negative log likelihood function (to minimize with optim)
  sum( - (k)*log(p) - (n-k)*log(1-p) )
}

#-------------------------------
# estimate_scr_sis()
#
# fit a reversible catalytic model with a single seroconversion rate
# assuming a fixed seroreversion rate
#
# vectors serop, agey, and id (if specified) must be the same length!
#
# @parameters
# @serop : seropositivity indicator (1=seropos+, 0=seroneg-)
# @agey  : age in years
# @id    : ID variable for independent units (individuals or clusters)
#          if unspecified, assumes 1:length(serop) (individuals)
# @srr   : assumed sero-reversion rate (SRR), fixed
# @variance : Variance estimator: "ML", "bootstrap"
#             ML: maximum likelihood (the default), assumes i.i.d obs
#             bootstrap: non-parametric bootstrap, resampling ids w/ replacement
# @nboots    : if variance="bootstrap", number of bootstrap replicates, 
#              default is 100 (use 1000+ for publication-quality CIs)
# @seed      : set a seed for a reproducible bootstrap
#
# returns a list with the following objects
# scr       : ML estimate of the seroconversion rate (SCR)
# scr_se    : ML estimate of the SCR standard error (if using Fisher Info, NA otherwise)
# scr_min95        : lower 95% confidence interval for the SCR
# scr_max95        : upper 95% confidence interval for the SCR
# scr_ci_method    : method used to estimate the 95% CI (character) 
# minll_value      : negative log-likelihood value minimized
# sis_err_warn_msg : warning message (in the case of a failed fit) (character)
#-------------------------------
estimate_scr_sis <- function(serop,agey,id=1:length(serop),
                    srr,variance="ML",nboots=100,seed=NULL) {
  
  # set a seed (if specified)
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # confirm that serop, agey, and id are all the same length
  if((length(serop) != length(agey)) | (length(serop) != length(id))) {
    e <- simpleError(paste0("\n\nThe length of serop (",length(serop),"), agey (",length(agey), "), and id (",length(id),") vectors are not the same.\nCheck to ensure they are the same length."))
    stop(e)
  }
  
  # create a data frame and restrict
  # to non-missing observations
  df  <- data.frame(id,serop,agey) %>%
    filter(!is.na(serop) & !is.na(agey))
  
  # identify the number of unique IDs
  ids <- unique(df$id)
  
  # if bootstrap is not specified but there are fewer unique IDs
  # than observations, notify with a warning
  if(variance == "ML" & (length(ids)< length(df$serop))) {
    cat(paste0("\nAttention! The number of unique IDs (",length(ids),")\nis fewer than the number of observations in the dataset (",length(serop),"),\n but bootstrap=FALSE. \n The SEs could be wrong if the data are clustered at the level of the ID variable. \n Did you intend to specify variance='ML'?"))
  }
  
  # collapse the data to 1 obs per unique age in the dataset
  # this is just b/c the loglik_sis function takes data in that form
  df2 <- df %>%
    group_by(agey) %>%
    summarise(seropos = sum(serop), nobs = n(),
              .groups = "drop") %>% 
    select(agey, nobs, seropos)
  # Fit SIS model
  # maximum likelihood solution
  # with fixed sero-reversion rate
  minll <- tryCatch(optim(par = c(0.1),
                          fn = loglik_sis,
                          r = srr,
                          data = df2,
                          hessian = TRUE,
                          method = "Brent",
                          lower = 0,
                          upper = 1),
                    error = function(cond){return(cond$message)},
                    warning = function(cond){return(cond$message)})
  
  if(class(minll)[1] == "character") {
    minll_value <- NA
    scr_hat <- NA
    scr_se  <- NA
    scr_min95 <- NA
    scr_max95 <- NA
    scr_ci_method <- NA
    sis_err_warn_msg <- minll
  } else {
    minll_value <- minll$value
    scr_hat <- minll$par
    # estimate SE and CI from the inverse information matrix
    I <- solve(minll$hessian)
    scr_se <- sqrt(diag(I))
    scr_min95 <- scr_hat-1.96*scr_se
    scr_max95 <- scr_hat+1.96*scr_se
    scr_ci_method <- "the Maximum Likelihood SE"
    sis_err_warn_msg <- NA
  }
  
  # if the bootstrap is specified (bootstrap==TRUE)
  # fit the SIS model without the Hessian to increase speed
  # since we do not need the information matrix to estimate SEs
  if(variance=="bootstrap") {
    
    cat(paste0("\n\nFitting reversible catalytic model (SIS)\n with inference using a nonparametric bootstrap (",nboots," iterations)\n"))
    
    bsamp <- matrix(sample(ids,size=length(ids)*nboots,replace=TRUE),
                    nrow=length(ids),ncol=nboots)
    
    bootests <- foreach(brep=1:nboots,.combine=rbind) %do% {
      if(brep %% 50 == 0) {
        cat(".",brep,"\n")
      }else{
        cat(".")
        }
      
        di <- df %>%
          left_join(data.frame(id=bsamp[,brep]),di,by=c("id")) %>%
          group_by(agey) %>%
          summarise(seropos = sum(serop), nobs = n(),
                    .groups = "drop") %>% 
          select(agey, nobs, seropos)
        
        # maximum likelihood solution
        # with fixed sero-reversion rate
        minll <- tryCatch(optim(par = c(0.1),
                                fn = loglik_sis,
                                r = srr,
                                data = di,
                                hessian = FALSE,
                                method = "Brent",
                                lower = 0,
                                upper = 1),
                          error = function(cond){return(cond$message)},
                          warning = function(cond){return(cond$message)})
        
        if(class(minll)[1] == "character") {
          scr_hati <- NA
        } else {
          scr_hati <- minll$par
        }
        res <- data.frame(brep,scr_hati)
        return(res)
    }
    # from bootstrap, estimate the 95% CI
    scr_se <- NA
    scr_min95 <- quantile(bootests$scr_hati,prob=0.025)
    scr_max95 <- quantile(bootests$scr_hati,prob=0.975)
    scr_ci_method <- paste0("Bootstrap (",nboots," reps)")
  }
  
  cat("\n\n-----------------------------------------------\n",
      "SIS model fit of the seroconversion rate (SCR)\n",
      "SCR units are seroconversions per child-year \n",
      "----------------------------------------------\n",
      "Non-missing observations:   ", nrow(df),"\n",
      "Number of independent units:",length(ids),"\n",
      "95% CI estimated using", scr_ci_method,"\n",
      "Assumed (fixed) seroreversion rate:",srr,"\n",
      "Estimates:\n",
      "SCR estimate:", sprintf("%1.4f",scr_hat),"\n",
      "SCR min 95  :", sprintf("%1.4f",scr_min95),"\n",
      "SCR max 95  :", sprintf("%1.4f",scr_max95),"\n",
      "----------------------------------------------\n")
  
  return(list(scr=scr_hat,scr_se=scr_se,scr_min95=scr_min95,scr_max95=scr_max95,scr_ci_method=scr_ci_method,srr=srr,minll_value=minll_value,sis_err_warn_msg=sis_err_warn_msg))
  
}

#-------------------------------
# estimate_scr_glm()
#
# fit a catalytic model with a single seroconversion rate using a Generalized Linear Model (GLM)
# this is equivalent to an SIR model
#
# vectors serop, agey, and id (if specified) must be the same length!
#
# @serop : seropositivity indicator (1=seropos+, 0=seroneg-)
# @agey  : age in years
# @id    : ID variable for independent units (individuals or clusters)
#          if unspecified, assumes 1:length(serop) (individuals)
# @variance : Variance estimator: "ML", "robust", "bootstrap"
#             ML: maximum likelihood (the default), assumes i.i.d obs
#             robust: Huber-White robust SEs, clustered at the level of id
#             bootstrap: non-parametric bootstrap, resampling ids w/ replacement
# @nboots    : number of bootstrap replicates, 
#              default is 100 (use 1000+ for publication-quality CIs)
# 
# returns a list with the following objects
# scr       : ML estimate of the seroconversion rate (SCR)
# scr_se    : ML estimate of the SCR standard error (NA if using a bootstrap)
# scr_min95        : lower 95% confidence interval for the SCR
# scr_max95        : upper 95% confidence interval for the SCR
# scr_ci_method    : method used to estimate the 95% CI (character) 
# glm_err_warn_msg : warning message (in the case of a failed fit) (character)
#-------------------------------
estimate_scr_glm <- function(serop,agey,id=1:length(serop),
                    variance="ML",nboots=100,seed=NULL) {
  
  # set a seed (if specified)
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # confirm that serop, agey, and id are all the same length
  if((length(serop) != length(agey)) | (length(serop) != length(id))) {
    e <- simpleError(paste0("\n\nThe length of serop (",length(serop),"), agey (",length(agey), "), and id (",length(id),") vectors are not the same.\nCheck to ensure they are the same length."))
    stop(e)
  }
  
  # create a data frame and restrict
  # to non-missing observations
  df  <- data.frame(id,serop,agey) %>%
    filter(!is.na(serop) & !is.na(agey))
  
  # identify the number of unique IDs
  ids <- unique(df$id)
  
  # fit the GLM model
  glm_fit <- tryCatch(glm(serop~1, offset = log(agey),
                          data = df,
                          family = binomial(link = "cloglog")),
                      error = function(cond){return(cond$message)},
                      warning = function(cond){return(cond$message)})
  
  if(class(glm_fit)[1] == "character") {
    scr_hat <- NA
    logscr_se  <- NA
    scr_min95 <- NA
    scr_max95 <- NA
    scr_ci_method <- NA
    glm_err_warn_msg <- glm_fit
  } else {
    # retrieve the SCR and the SE of the log SCR
    scr_hat <- exp(glm_fit$coefficients)
    logscr_se <- sqrt(summary(glm_fit)$cov.unscaled)
    scr_ci_method = "Maximum likelihood SE"
    glm_err_warn_msg <- NA
    cat("\n   -------------------------------",
        "\n   GLM model fit",
        "\n   -------------------------------")
    print(summary(glm_fit))
    
    # if variance=="robust"
    # adjust the SEs using A clustered, Huber-White robust variance.
    if(variance == "robust") {
      glm_fit_rb <- coeftest(glm_fit,vcov. = vcovCL(glm_fit,cluster=df$id))
      logscr_se <- glm_fit_rb[1,2]
      scr_ci_method = "Robust, Huber-White SE"
      cat("\n   -------------------------------",
          "\n   GLM model fit, with Robust, SEs",
          "\n   -------------------------------")
      print(glm_fit_rb)
    } 
  }
  
  # estimate the 95% CIs
  scr_min95 = exp(log(scr_hat) - 1.96*logscr_se)
  scr_max95 = exp(log(scr_hat) + 1.96*logscr_se)
  
  # if variance = "bootstrap"
  # resample the ID variable with replacement and re-estimate the SCR
  if(variance=="bootstrap") {
    
    cat(paste0("\n\nFitting a catalytic model (SIR) using GLM\n with inference using a nonparametric bootstrap (",nboots," iterations)\n"))
    
    bsamp <- matrix(sample(ids,size=length(ids)*nboots,replace=TRUE),
                    nrow=length(ids),ncol=nboots)
    
    bootests <- foreach(brep=1:nboots,.combine=rbind) %do% {
      if(brep %% 50 == 0) {
        cat(".",brep,"\n")
      }else{
        cat(".")
      }
    
      di <- df %>%
        left_join(data.frame(id=bsamp[,brep]),di,by=c("id")) 
      
      # fit the GLM model
      glm_fiti <- tryCatch(glm(serop~1, offset = log(agey),
                              data = di,
                              family = binomial(link = "cloglog")),
                          error = function(cond){return(cond$message)},
                          warning = function(cond){return(cond$message)})
      
      if(class(glm_fiti)[1] == "character") {
        scr_hati <- NA
      } else {
        scr_hati <- exp(glm_fiti$coefficients)
      }
      res <- data.frame(brep,scr_hati)
      return(res)
    }
    # from bootstrap, estimate the 95% CI
    logscr_se <- NA
    scr_min95 <- quantile(bootests$scr_hati,prob=0.025)
    scr_max95 <- quantile(bootests$scr_hati,prob=0.975)
    scr_ci_method <- paste0("Bootstrap (",nboots," reps)")
    
  }
  
  cat("\n\n-----------------------------------------------\n",
      "GLM model fit of the seroconversion rate (SCR)\n",
      "SCR units are seroconversions per child-year \n",
      "----------------------------------------------\n",
      "Non-missing observations:   ", nrow(df),"\n",
      "Number of independent units:",length(ids),"\n",
      "95% CI estimated using", scr_ci_method,"\n",
      "Estimates:\n",
      "SCR estimate:", sprintf("%1.4f",scr_hat),"\n",
      "SCR min 95  :", sprintf("%1.4f",scr_min95),"\n",
      "SCR max 95  :", sprintf("%1.4f",scr_max95),"\n",
      "----------------------------------------------\n")
  
  return(list(scr=scr_hat,logscr_se=logscr_se,scr_min95=scr_min95,scr_max95=scr_max95,scr_ci_method=scr_ci_method,glm_err_warn_msg=glm_err_warn_msg))
  
}
  
#-------------------------------
# estimate_seroprev_glm()
#
# estimate population seroprevalence using a Generalized Linear Model (GLM)
# this allows for correction for clustering in the survey design
#
# vectors serop, agey, and id (if specified) must be the same length!
#
# @serop : seropositivity indicator (1=seropos+, 0=seroneg-)
# @agey  : age in years
# @id    : ID variable for independent units (individuals or clusters)
#          if unspecified, assumes 1:length(serop) (individuals)
# @variance : Variance estimator: "ML", "robust", "bootstrap"
#             ML: maximum likelihood (the default), assumes i.i.d obs
#             robust: Huber-White robust SEs, clustered at the level of id
#             bootstrap: non-parametric bootstrap, resampling ids w/ replacement
# @nboots    : number of bootstrap replicates, 
#              default is 100 (use 1000+ for publication-quality CIs)
# 
# returns a list with the following objects
# logodds               : ML estimate of the log odds of seroprevalence
# logodds_se            : ML estimate of the seroprevalence standard error on the log odds scale (NA if using a bootstrap)
# seroprev              : ML estimate of the seroprevalence
# seroprev_min95        : lower 95% confidence interval for seroprevalence
# seroprev_max95        : upper 95% confidence interval for seroprevalence
# seroprev_ci_method    : method used to estimate the 95% CI (character) 
# glm_err_warn_msg : warning message (in the case of a failed fit) (character)
#-------------------------------
estimate_seroprev_glm <- function(serop,agey,id=1:length(serop),
                             variance="ML",nboots=100,seed=NULL) {
  
  # set a seed (if specified)
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # confirm that serop, agey, and id are all the same length
  if((length(serop) != length(agey)) | (length(serop) != length(id))) {
    e <- simpleError(paste0("\n\nThe length of serop (",length(serop),"), agey (",length(agey), "), and id (",length(id),") vectors are not the same.\nCheck to ensure they are the same length."))
    stop(e)
  }
  
  # create a data frame and restrict
  # to non-missing observations
  df  <- data.frame(id,serop,agey) %>%
    filter(!is.na(serop) & !is.na(agey))
  
  # identify the number of unique IDs
  ids <- unique(df$id)
  
  # fit the GLM model
  glm_fit <- tryCatch(glm(serop~1,data=df,family=binomial(link="logit")),
                      error = function(cond){return(cond$message)},
                      warning = function(cond){return(cond$message)})
  
  if(class(glm_fit)[1] == "character") {
    logodds <- NA
    logodds_se  <- NA
    seroprev_hat <- NA
    seroprev_min95 <- NA
    seroprev_max95 <- NA
    seroprev_ci_method <- NA
    glm_err_warn_msg <- glm_fit
  } else {
    # retrieve the seroprevalence and the SE of the seroprevalence
    logodds     <- glm_fit$coefficients
    logodds_se  <- sqrt(summary(glm_fit)$cov.unscaled)
    seroprev_hat <- plogis(glm_fit$coefficients)
    seroprev_ci_method <- "Maximum likelihood SE"
    glm_err_warn_msg <- NA
    cat("\n   -------------------------------",
        "\n   GLM model fit",
        "\n   -------------------------------")
    print(summary(glm_fit))
    
    # if variance=="robust"
    # adjust the SEs using A clustered, Huber-White robust variance.
    if(variance == "robust") {
      glm_fit_rb <- coeftest(glm_fit,vcov. = vcovCL(glm_fit,cluster=df$id))
      logodds_se <- glm_fit_rb[1,2]
      seroprev_ci_method <- "Robust, Huber-White SE"
      cat("\n   -------------------------------",
          "\n   GLM model fit, with Robust, SEs",
          "\n   -------------------------------")
      print(glm_fit_rb)
    } 
  }
  
  # estimate the 95% CIs
  seroprev_min95 = plogis(logodds - 1.96*logodds_se)
  seroprev_max95 = plogis(logodds + 1.96*logodds_se)
  
  # if variance = "bootstrap"
  # resample the ID variable with replacement and re-estimate the SCR
  if(variance=="bootstrap") {
    
    cat(paste0("\n\nEstimating prevalence with a logistic model using GLM\n with inference using a nonparametric bootstrap (",nboots," iterations)\n"))
    
    bsamp <- matrix(sample(ids,size=length(ids)*nboots,replace=TRUE),
                    nrow=length(ids),ncol=nboots)
    
    bootests <- foreach(brep=1:nboots,.combine=rbind) %do% {
      if(brep %% 50 == 0) {
        cat(".",brep,"\n")
      }else{
        cat(".")
      }
      
      di <- df %>%
        left_join(data.frame(id=bsamp[,brep]),di,by=c("id")) 
      
      # fit the GLM model
      glm_fiti <- tryCatch(glm(serop~1,data=df,family=binomial(link="logit")),
                           error = function(cond){return(cond$message)},
                           warning = function(cond){return(cond$message)})
      
      if(class(glm_fiti)[1] == "character") {
        seroprev_hati <- NA
      } else {
        logodds <- glm_fiti$coefficients
        seroprev_hati <- plogis(glm_fiti$coefficients)
      }
      res <- data.frame(brep,seroprev_hati)
      return(res)
    }
    # from bootstrap, estimate the 95% CI
    logodds_se <- NA
    seroprev_min95 <- quantile(bootests$seroprev_hati,prob=0.025)
    seroprev_max95 <- quantile(bootests$seroprev_hati,prob=0.975)
    seroprev_ci_method <- paste0("Bootstrap (",nboots," reps)")
  }
 
  cat("\n\n-----------------------------------------------\n",
      "GLM estimate of seroprevalence \n",
      "Seroprevalence is in percentage points (%) \n",
      "----------------------------------------------\n",
      "Non-missing observations:   ", nrow(df),"\n",
      "Number of independent units:",length(ids),"\n",
      "95% CI estimated using", seroprev_ci_method,"\n",
      "Estimates:\n",
      "Seroprevalence estimate:", sprintf("%1.1f",seroprev_hat*100),"%\n",
      "Seroprevalence min 95  :", sprintf("%1.1f",seroprev_min95*100),"%\n",
      "Seroprevalence max 95  :", sprintf("%1.1f",seroprev_max95*100),"%\n",
      "----------------------------------------------\n")
  
  return(list(logodds = logodds,logodds_se=logodds_se,seroprev=seroprev_hat,seroprev_min95=seroprev_min95,seroprev_max95=seroprev_max95, seroprev_ci_method=seroprev_ci_method,glm_err_warn_msg=glm_err_warn_msg))
  
}
