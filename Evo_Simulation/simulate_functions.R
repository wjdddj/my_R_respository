library(survival)

################################################################################################################
## Assign responder or non-responder based on input cutoff survival date. 
## For short censored cases, RANDOMLY split them into either groups proportionally.
## R/NR are for both arms and with the same cutoff
## input:
##    df_surv. a data.frame containing 
##      1. column named 'time', numeric;
##      2. column named 'event', binary with 0 for censored cases; 1 for event
## output:
##    a vector of true labels at the same length and order of nrow of the input df_surv
################################################################################################################
# assignTrueLabel <- function(df_surv, surv_cutoff = NULL){
#   
#   df_surv$trueLabel = ifelse(df_surv$time > surv_cutoff, 'R', ifelse(df_surv$event == 0, 'unknown', 'NR'))
#   
#   n_R = sum(df_surv$trueLabel == 'R')
#   n_NR = sum(df_surv$trueLabel == 'NR')
#   n_unk = sum(df_surv$trueLabel == 'unknown')
#   r = n_R/(n_R + n_NR)
#   idx_unk = which(df_surv$trueLabel == 'unknown')
#   idx_unk_R = sample(idx_unk, round(r*n_unk, 0), replace = F)
#   idx_unk_NR = setdiff(idx_unk, idx_unk_R)
#   df_surv$trueLabel[idx_unk_R] <- 'R'
#   df_surv$trueLabel[idx_unk_NR] <- 'NR'
#   
#   trueLabel <- df_surv$trueLabel
#   trueLabel <- factor(trueLabel, levels = c('R', 'NR'))
#   z <- list(trueLabel = trueLabel)
# }

assignTrueLabel <- function(
  df_surv, evo_cutoff = NULL, plb_cutoff = NULL
){
  # df_surv <- df_OS
  df_surv$trueLabel = NA
  df_surv$trueLabel[df_surv$Group == 'TH-302'] = ifelse(
    df_surv$time[df_surv$Group == 'TH-302'] > evo_cutoff, 
    'R', 
    ifelse(
      df_surv$event[df_surv$Group == 'TH-302'] == 0, 'unknown', 'NR'
    )
  )
  df_surv$trueLabel[df_surv$Group == 'PLACEBO'] = ifelse(
    df_surv$time[df_surv$Group == 'PLACEBO'] > plb_cutoff, 
    'R', 
    ifelse(
      df_surv$event[df_surv$Group == 'PLACEBO'] == 0, 'unknown', 'NR'
    )
  )
  
  n_R_Exp = sum(df_surv$Group == 'TH-302' & df_surv$trueLabel == 'R')
  n_NR_Exp = sum(df_surv$Group == 'TH-302' & df_surv$trueLabel == 'NR')
  n_unk_Exp = sum(df_surv$Group == 'TH-302' & df_surv$trueLabel == 'unknown')
  n_R_PLB = sum(df_surv$Group == 'PLACEBO' & df_surv$trueLabel == 'R')
  n_NR_PLB = sum(df_surv$Group == 'PLACEBO' & df_surv$trueLabel == 'NR')
  n_unk_PLB = sum(df_surv$Group == 'PLACEBO' & df_surv$trueLabel == 'unknown')
  r_Exp = n_R_Exp/(n_R_Exp + n_NR_Exp)
  r_PLB = n_R_PLB/(n_R_PLB + n_NR_PLB)
  
  idx_unk_Exp = which(df_surv$Group == 'TH-302' & df_surv$trueLabel == 'unknown')
  idx_unk_PLB = which(df_surv$Group == 'PLACEBO' & df_surv$trueLabel == 'unknown')
  idx_unk <- which(df_surv$trueLabel == 'unknown')
  
  idx_unk_R = as.numeric(c(
    sample(as.character(idx_unk_Exp), round(r_Exp*n_unk_Exp, 0), replace = F),
    sample(as.character(idx_unk_PLB), round(r_PLB*n_unk_PLB, 0), replace = F)
  ))
  idx_unk_NR = setdiff(idx_unk, idx_unk_R)
  df_surv$trueLabel[idx_unk_R] <- 'R'
  df_surv$trueLabel[idx_unk_NR] <- 'NR'
  
  trueLabel <- df_surv$trueLabel
  trueLabel <- factor(trueLabel, levels = c('R', 'NR'))
  z <- list(trueLabel = trueLabel)
}

################################################################################################################
## Assign test label based on input sensitivity and specificity
## two modes: 
##    arms == 'single', sensitivity/specificity are based on experimental arm; placebo arm is split randomly;
##    arms == 'combine', sensitivity/specificity are based on both arms with the same cutoff;
##    arms == 'double', sensitivity/specificity are different between arms.
## input:
##    df_surv. a data.frame containing 
##      1. column named 'trueLabel' with factor level c('R', 'NR');
##      2. column named 'Group' with factor level c('TH-302', 'PLACEBO')
## output:
##    a vector of test labels at the same length and order of nrow of the input df_surv
################################################################################################################
assignTestLabel <- function(
  df_surv, arms = 'single', 
  se, sp, # when arms == 'single', se/sp for experimental arm only.
  sp_plb # for placebo arm, only input when arms == 'double'
){
  
  if(arms == 'combine'){
    n_R = sum(df_surv$trueLabel == 'R')
    n_NR = sum(df_surv$trueLabel == 'NR')
    pos_rate = (n_R*se + n_NR*(1-sp))/(n_NR + n_R)
    idx_test_R <- c(
      sample(which(df_surv$trueLabel == 'R'), round(n_R*se, 0), replace = F),
      sample(which(df_surv$trueLabel == 'NR'), round(n_NR*(1-sp), 0), replace = F)
    )
    df_surv$testLabel <- 'NR'
    df_surv$testLabel[idx_test_R] <- 'R'
    se_plb = sp_plb = NA
    
  }else if(arms == 'single'){
    n_R_Exp = sum(df_surv$Group == 'TH-302' & df_surv$trueLabel == 'R')
    n_NR_Exp = sum(df_surv$Group == 'TH-302' & df_surv$trueLabel == 'NR')
    n_PLB = sum(df_surv$Group == 'PLACEBO')
    n_Exp = sum(df_surv$Group == 'TH-302')
    pos_rate = (n_R_Exp*se + n_NR_Exp*(1-sp))/n_Exp # two arms are randomized, must share the same positive rate
    idx_test_R <- c(
      sample(which(df_surv$Group == 'TH-302' & df_surv$trueLabel == 'R'), round(n_R_Exp*se, 0), replace = F),
      sample(which(df_surv$Group == 'TH-302' & df_surv$trueLabel == 'NR'), round(n_NR_Exp*(1-sp), 0), replace = F),
      sample(which(df_surv$Group == 'PLACEBO'), round(n_PLB*pos_rate, 0), replace = F)
    )
    df_surv$testLabel <- 'NR'
    df_surv$testLabel[idx_test_R] <- 'R'
    se_plb = sp_plb = NA
    
  }else if(arms == 'double'){
    n_R_Exp = sum(df_surv$Group == 'TH-302' & df_surv$trueLabel == 'R')
    n_NR_Exp = sum(df_surv$Group == 'TH-302' & df_surv$trueLabel == 'NR')
    n_R_PLB = sum(df_surv$Group == 'PLACEBO' & df_surv$trueLabel == 'R')
    n_NR_PLB = sum(df_surv$Group == 'PLACEBO' & df_surv$trueLabel == 'NR')
    n_Exp = sum(df_surv$Group == 'TH-302')
    n_PLB = sum(df_surv$Group == 'PLACEBO')
    
    # two arms are randomized, must share the same positive rate
    pos_rate = (n_R_Exp*se + n_NR_Exp*(1-sp))/n_Exp 
    se_plb = (n_PLB*pos_rate - (1-sp_plb)*n_NR_PLB)/n_R_PLB
    se_plb = max(min(se_plb, 1), 0)
    if(se_plb == 1){
      sp_plb = 1 - ((n_PLB*pos_rate - se_plb*n_R_PLB)/n_NR_PLB)
    }
    if(se_plb == 0){
      sp_plb = n_PLB*pos_rate/n_NR_PLB
    }
    
    idx_test_R <- c(
      sample(which(df_surv$Group == 'TH-302' & df_surv$trueLabel == 'R'), round(n_R_Exp*se, 0), replace = F),
      sample(which(df_surv$Group == 'TH-302' & df_surv$trueLabel == 'NR'), round(n_NR_Exp*(1-sp), 0), replace = F),
      sample(which(df_surv$Group == 'PLACEBO' & df_surv$trueLabel == 'R'), round(n_R_PLB*se_plb, 0), replace = F),
      sample(which(df_surv$Group == 'PLACEBO' & df_surv$trueLabel == 'NR'), round(n_NR_PLB*(1-sp_plb), 0), replace = F)
    )
    df_surv$testLabel <- 'NR'
    df_surv$testLabel[idx_test_R] <- 'R'
  }
  
  testLabel <- df_surv$testLabel
  testLabel <- factor(testLabel, levels = c('R', 'NR'))
  z <- list(
    testLabel = testLabel,
    sesp_plb = list(
      pos_rate = pos_rate, 
      se_plb_adjusted = se_plb, # adjusted based on positive rate
      sp_plb_adjusted = sp_plb # adjusted based on positive rate
    )
  )
}


################################################################################################################
################################################################################################################
getSurvMetricsGivenTest <- function(
  df_surv
){
  # df_surv <- df_OS
  df_surv <- df_surv[, c('time', 'event', 'Group', 'testLabel')]
  df_surv1 <- df_surv[df_surv$testLabel == 'R', ]
  df_surv2 <- df_surv[df_surv$testLabel == 'NR', ]
  
  # get log-rank test p value
  fitSdiff1 <- survdiff(Surv(time, event) ~ Group, data = df_surv1, rho = 0)
  logrankR_p <- 1 - pchisq(fitSdiff1$chisq, length(fitSdiff1$n) - 1)
  fitSdiff2 <- survdiff(Surv(time, event) ~ Group, data = df_surv2, rho = 0)
  logrankNR_p <- 1 - pchisq(fitSdiff2$chisq, length(fitSdiff2$n) - 1)
  
  # fit coxph model
  fitCox1 <- coxph(Surv(time, event) ~ Group, data = df_surv1)
  coxSum1 <- summary(fitCox1)
  fitCox2 <- coxph(Surv(time, event) ~ Group, data = df_surv2)
  coxSum2 <- summary(fitCox2)
  
  # fit KM curve
  fitSurv1 <- survfit(Surv(time, event) ~ Group, data = df_surv1)
  survSum1 <- summary(fitSurv1)
  fitSurv2 <- survfit(Surv(time, event) ~ Group, data = df_surv2)
  survSum2 <- summary(fitSurv2)
  fitSurv <- survfit(Surv(time, event) ~ testLabel + Group, data = df_surv)
  survSum <- summary(fitSurv)

  z <- list(
    survR = survSum1$table,
    survNR = survSum2$table,
    survBoth = survSum$table,
    logrankR_p = logrankR_p,
    logrankNR_p = logrankNR_p,
    coxR = coxSum1$coefficients,
    coxR_scTestP = coxSum1$sctest['pvalue'], 
    coxNR = coxSum2$coefficients,
    coxNR_scTestP = coxSum2$sctest['pvalue'], 
    fitSurvR = fitSurv1,
    fitSurvNR = fitSurv2,
    fitSurvBoth = fitSurv
  )
}

################################################################################################################
################################################################################################################
generateParameters <- function(
  evo_cutoff = c(150),
  plb_cutoff = c(50), 
  # surv_cutoff = c(100, 150, 200), # cutoffs for determining true label
  arms = c('single', 'double', 'combine'), # arms = c('single', 'double', 'combine')
  se = seq(0.1, 1, 0.1),
  sp = seq(0.1, 1, 0.1), 
  # sesp = list(c(0.9, 0.5), c(0.9, 0.9), c(0.5, 0.9)), # sensitivity/specificity for experimental arm, regardless of arms parameter.
  sp_plb = seq(0.2, 0.9, 0.1) # specificity for placebo arm (rate between placebo and experimental arm), when arms == 'double'.
){
  
  df_param <- expand.grid(
    arms = arms, 
    evo_cutoff = evo_cutoff,
    plb_cutoff = plb_cutoff,
    se = se, 
    sp = sp, 
    sp_plb = sp_plb 
  )
  
  for(i in 1:nrow(df_param)){
    if(df_param[i, 'arms'] != 'double'){
      df_param[i, 'sp_plb'] = NA
    }
    if(df_param[i, 'se'] + df_param[i, 'sp'] <= 1){
      df_param[i, ] <- NA
    }
  }
  df_param <- df_param[!is.na(df_param$se), ]
  df_param <- unique(df_param)
  
  paramSets <- apply(df_param, 1, function(x){
    list(
      arms = as.character(x['arms']), 
      evo_cutoff = as.numeric(x['evo_cutoff']), 
      plb_cutoff = as.numeric(x['plb_cutoff']),
      se = as.numeric(x['se']), 
      sp = as.numeric(x['sp']),
      sp_plb = as.numeric(x['sp_plb'])
    )
  })
}

################################################################################################################
## input:
##    df_surv. a data.frame containing 
##      1. column named 'time', numeric;
##      2. column named 'event', binary with 0 for censored cases; 1 for event
##      3. column named 'Group' with factor level c('TH-302', 'PLACEBO')
##    params. a list containing:
##      list(surv_cutoff, arms, se, sp, se_plb, sp_plb, pos_rate)
################################################################################################################
coreSimulation <- function(
  df_surv, params
){
  # df_surv <- df_OS
  assignTrueLabel.obj <- assignTrueLabel(df_surv, evo_cutoff = params$evo_cutoff, plb_cutoff = params$plb_cutoff)
  df_surv$trueLabel <- assignTrueLabel.obj$trueLabel
  assignTestLabel.obj <- assignTestLabel(
    df_surv, arms = params$arms, 
    se = params$se, sp = params$sp, 
    sp_plb = params$sp_plb
  )
  df_surv$testLabel <- assignTestLabel.obj$testLabel
  
  surv_fit <- getSurvMetricsGivenTest(df_surv)
  surv_fit <- surv_fit[1:9] # only collecting numbers, ignore survfit objects. can be created later
  z <- list(
    params = c(params, assignTestLabel.obj$sesp_plb),
    Labels = df_surv[,c('trueLabel', 'testLabel')],
    surv_fit = surv_fit
  )
}

################################################################################################################
################################################################################################################
repeatCoreSimulation <- function(
  df_surv, params, r = 100
){
  result <- list()
  for(i in 1:r){
    sim <- coreSimulation(df_surv, params)
    cat(unlist(sim$params), '\n')
    result[[i]] <- list(
      surv_fit = sim$surv_fit,
      params = sim$params,
      Labels = sim$Labels
    )
  }
  
  params <- c(
    params, 
    list(
      pos_rate = round(mean(sapply(result, function(x) x$params$pos_rate)), 3),
      se_plb_adjusted = round(mean(sapply(result, function(x) x$params$se_plb_adjusted)), 3),
      sp_plb_adjusted = round(mean(sapply(result, function(x) x$params$sp_plb_adjusted)), 3)
    )
  )
  
  z <- list(
    params = params,
    result = result
  )
}

################################################################################################################
################################################################################################################
performSimulation <- function(
  df_surv,
  paramSets,
  r = 100,
  n_thread = 8
){
  library(doParallel)
  library(foreach)
  library(plyr)
  library(cvTools)
  
  n_para <- length(paramSets)
  
  set.seed(417)
  map_para <- cvFolds(n_para, n_thread)
  map_para <- dlply(data.frame(set = map_para$subsets, thread = map_para$which), .variables = 'thread')
  
  cl <- makeCluster(n_thread)
  registerDoParallel(cl)
  result <- foreach(
    i = 1:n_thread,
    .combine = 'c',
    .packages = c('survival'),
    .export = c(
      'assignTrueLabel',
      'assignTestLabel',
      'repeatCoreSimulation', 
      'coreSimulation',
      'getSurvMetricsGivenTest'
    )
  ) %dopar% {
    simByParamSet <- list()
    for(j in 1:nrow(map_para[[i]])){
      params = paramSets[[map_para[[i]]$set[j]]]
      set.seed(417)
      simByParamSet[[j]] <- repeatCoreSimulation(df_surv, params, r = r)
    }
    simByParamSet
  }
  stopCluster(cl)
  
  ## rearrange by parameters (arms, surv_cutoff, se)
  df_param <- data.frame(t(sapply(result, function(x) x$params)), stringsAsFactors = F)
  df_param <- sapply(df_param, unlist, recursive = F)
  df_param <- data.frame(df_param, stringsAsFactors = F)
  result <- result[with(df_param, order(arms, evo_cutoff, plb_cutoff, se, sp))]
  result
}







