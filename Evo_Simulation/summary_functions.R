################################################################################################################
################################################################################################################
summarizeSimulation <- function(
  performSimulation.obj
){
  # performSimulation.obj <- test
  n <- length(performSimulation.obj)
  simSum <- lapply(1:n, function(i){
    
    param <- performSimulation.obj[[i]]$params
    mean_HR_GivenR <- mean(
      sapply(performSimulation.obj[[i]]$result, function(x) x$surv_fit$coxR[1, 2])
    )
    mean_HR_GivenR <- round(mean_HR_GivenR, 3)
    CI95_HR_GivenR <- quantile(
      sapply(performSimulation.obj[[i]]$result, function(x) x$surv_fit$coxR[1, 2]), 
      c(0.025, 0.975)
    )
    CI95_HR_GivenR <- round(CI95_HR_GivenR, 3)
    CI95_HR_GivenR <- paste(CI95_HR_GivenR, collapse = ',')
    CI95_HR_GivenR <- paste0('(', CI95_HR_GivenR, ')')
    
    mean_sctestP_GivenR <- mean(
      sapply(performSimulation.obj[[i]]$result, function(x) x$surv_fit$coxR_scTestP)
    )
    mean_sctestP_GivenR <- format(signif(mean_sctestP_GivenR, 3), scientific = T)
    CI95_sctestP_GivenR <- quantile(
      sapply(performSimulation.obj[[i]]$result, function(x) x$surv_fit$coxR_scTestP), 
      c(0.025, 0.975)
    )
    CI95_sctestP_GivenR <- format(signif(CI95_sctestP_GivenR, 3), scientific = T)
    CI95_sctestP_GivenR <- paste(CI95_sctestP_GivenR, collapse = ',')
    CI95_sctestP_GivenR <- paste0('(', CI95_sctestP_GivenR, ')')
    
    mean_logrankP_GivenR <- mean(
      sapply(performSimulation.obj[[i]]$result, function(x) x$surv_fit$logrankR_p)
    )
    mean_logrankP_GivenR <- format(signif(mean_logrankP_GivenR, 3), scientific = T)
    CI95_logrankP_GivenR <- quantile(
      sapply(performSimulation.obj[[i]]$result, function(x) x$surv_fit$logrankR_p), 
      c(0.025, 0.975)
    )
    CI95_logrankP_GivenR <- format(signif(CI95_logrankP_GivenR, 3), scientific = T)
    CI95_logrankP_GivenR <- paste(CI95_logrankP_GivenR, collapse = ',')
    CI95_logrankP_GivenR <- paste0('(', CI95_logrankP_GivenR, ')')
    
    medianSurv_GivenR <- sapply(performSimulation.obj[[i]]$result, function(x) x$surv_fit$survR[, 'median'])
    mean_medianSurv_GivenR <- apply(medianSurv_GivenR, 1, mean, na.rm = T) 
    mean_medianSurv_GivenR <- round(mean_medianSurv_GivenR, 0)
    
    out <- data.frame(
      t(unlist(param)),
      mean_HR_GivenR = mean_HR_GivenR,
      CI95_HR_GivenR = CI95_HR_GivenR,
      mean_sctestP_GivenR = mean_sctestP_GivenR,
      CI95_sctestP_GivenR = CI95_sctestP_GivenR,
      Group = names(mean_medianSurv_GivenR), 
      mean_medianSurv_GivenR = mean_medianSurv_GivenR,
      ratio_medianSurv_GivenR = mean_medianSurv_GivenR[2]/mean_medianSurv_GivenR[1],
      mean_logrankP_GivenR = mean_logrankP_GivenR,
      CI95_logrankP_GivenR = CI95_logrankP_GivenR, 
      stringsAsFactors = F
    )
    rownames(out) <- NULL
    out
  })
  
  simSum <- do.call(rbind, simSum)
  rownames(simSum) <- NULL
  colnames(simSum) <- c(
    'Analysis Scheme',
    'Cutoff Time for R/NR Evo Arm',
    'Cutoff Time for R/NR Placebo Arm',
    'Sensitivity for Evo Arm',
    'Specificity for Evo Arm',
    'input Specificity for Placebo Arm',
    'Positive Rate of the Test',
    'adjusted Sensitivity for Placebo Arm',
    'adjusted Specificity for Placebo Arm',
    'mean(HR)|Responder',
    '95%CI HR|Responder',
    'mean(score test pvalue)|Responder',
    '95%CI(score test pvalue)|Responder',
    'Treatment',
    'mean(median survival)|Responder',
    'ratio(median survival)|Responder',
    'mean(log-rank pvalue)|Responder',
    '95%CI(log-rank pvalue)|Responder'
  )
  simSum
}


################################################################################################################
################################################################################################################
plotSurvResponder <- function(
  getSurvMetricsGivenTest.obj, main
){
  plot(
    getSurvMetricsGivenTest.obj$fitSurvR,
    mark.time = T, 
    col = c('red', 'blue'), lwd = 2,
    main = main,
    xlab = 'survival time (days)', ylab = 'probability of survival', 
    cex.main = 1.5, cex = 1.5, cex.lab = 1.5
  )
  
  legend(
    'topright', 
    legend = c(
      paste0(
        rownames(getSurvMetricsGivenTest.obj$survR), 
        ', n=', getSurvMetricsGivenTest.obj$survR[, 'records'], 
        ', event=', getSurvMetricsGivenTest.obj$survR[, 'events']
      ),
      paste0(
        '(TH-302 | Responder): HR = ', round(getSurvMetricsGivenTest.obj$coxR[1, 2], 3), 
        '; p = ', round(getSurvMetricsGivenTest.obj$logrankR_p, 3)
      )
    ),
    col = c('red', 'blue', 'white'),
    horiz = FALSE,
    cex = 1, lty = 1, lwd = 2, 
    bty = 'n'
  )
}


################################################################################################################
################################################################################################################
plotSurvNonResponder <- function(
  getSurvMetricsGivenTest.obj, main
){
  plot(
    getSurvMetricsGivenTest.obj$fitSurvNR,
    mark.time = T, 
    col = c('red', 'blue'), lwd = 2,
    main = main,
    xlab = 'survival time (days)', ylab = 'probability of survival', 
    cex.main = 1.5, cex = 1.5, cex.lab = 1.5
  )
  
  legend(
    'topright', 
    legend = c(
      paste0(
        rownames(getSurvMetricsGivenTest.obj$survNR), 
        ', n=', getSurvMetricsGivenTest.obj$survNR[, 'records'], 
        ', event=', getSurvMetricsGivenTest.obj$survNR[, 'events']
      ),
      paste0(
        '(TH-302 | Non-Responder): HR = ', round(getSurvMetricsGivenTest.obj$coxNR[1, 2], 3), 
        '; p = ', round(getSurvMetricsGivenTest.obj$logrankNR_p, 3)
      )
    ),
    col = c('red', 'blue', 'white'),
    horiz = FALSE,
    cex = 1, lty = 1, lwd = 2, 
    bty = 'n'
  )
}

################################################################################################################
################################################################################################################
plotSurvBoth <- function(
  getSurvMetricsGivenTest.obj, main
){
  # getSurvMetricsGivenTest <- surv_fit
  plot(
    getSurvMetricsGivenTest.obj$fitSurvBoth,
    mark.time = T, 
    col = c('red', 'blue', 'black', 'green'), lwd = 2,
    main = main,
    xlab = 'survival time (days)', ylab = 'probability of survival', 
    cex.main = 1.5, cex = 1.5, cex.lab = 1.5
  )
  legend(
    'topright', 
    legend = c(
      paste0(
        rownames(getSurvMetricsGivenTest.obj$survBoth), 
        ', n=', getSurvMetricsGivenTest.obj$survBoth[, 'records'], 
        ', event=', getSurvMetricsGivenTest.obj$survBoth[, 'events']
      ),
      paste0(
        '(TH-302 | Responder): HR = ', round(getSurvMetricsGivenTest.obj$coxR[1, 2], 3), 
        '; p = ', round(getSurvMetricsGivenTest.obj$logrankR_p, 3)
      ),
      paste0(
        '(TH-302 | Non-Responder): HR = ', round(getSurvMetricsGivenTest.obj$coxNR[1, 2], 3), 
        '; p = ', round(getSurvMetricsGivenTest.obj$logrankNR_p, 3)
      )
    ), 
    col = c('red', 'blue', 'black', 'green', 'white', 'white'),
    horiz = FALSE,
    cex = 1, lty = 1, lwd = 2, 
    bty = 'n'
  )
}


################################################################################################################
################################################################################################################
plotSurv <- function(
  getSurvMetricsGivenTest.obj, main, 
  plot_type = 'both'
){
  if(plot_type == 'R'){
    plotSurvResponder(getSurvMetricsGivenTest.obj, main)
  }else if(plot_type == 'NR'){
    plotSurvNonResponder(getSurvMetricsGivenTest.obj, main)
  }else if(plot_type == 'both'){
    plotSurvBoth(getSurvMetricsGivenTest.obj, main)
  }else{
    stop("wrong parameter for plot. Please choose from ['R'|'NR'|'both']!\n")
  }
}


################################################################################################################
################################################################################################################
plotSampleSurv <- function(
  df_surv, performSimulation.obj, plot_type = 'R'
){
  
  today <- format(Sys.Date(), '%Y%m%d')
  
  n <- length(performSimulation.obj)
  
  lapply(1:n, function(i){
    param <- performSimulation.obj[[i]]$params
    header = paste(c(plot_type, param), collapse = '_')
    df_surv <- data.frame(
      df_surv,
      performSimulation.obj[[i]]$result[[1]]$Labels
    )
    pdf(paste0(today, '_(', header, ').pdf'), 10, 10)
    plotSurv(
      getSurvMetricsGivenTest(df_surv),
      plot_type = plot_type,
      main = header
    )
    dev.off()
  })
  
}




################################################################################################################
################################################################################################################
plotHREvo <- function(
  dfHM, evo_cutoff = 100, plb_cutoff = 100, scheme = 'double', addition = NULL
){
  ggplot(
    aes(x = spec_evo, y = sens_evo), 
    data = dfHM[dfHM$scheme == scheme&dfHM$evo_cutoff == evo_cutoff&dfHM$plb_cutoff == plb_cutoff,]
  ) +
    geom_tile(aes(fill = hr), colour = 'white') +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0.6, limits = c(0,1.5)) + 
    scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
    scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
    labs(
      x = 'Specificity for Evo Arm', 
      y = 'Sensitivity for Evo Arm', 
      title = paste0(
        'HR given performance in EVO Arm', 
        '\n(evo cutoff = ', evo_cutoff, '; plb cutoff = ', plb_cutoff, '; scheme = "', scheme, '"', addition, ')',
        '\nwhite = 0.6'
      )
    ) + 
    theme(
      title = element_text(size = 15), 
      axis.title = element_text(size = 15), 
      axis.text = element_text(size = 15)
    )
}

################################################################################################################
################################################################################################################
plotHRPlacebo <- function(
  dfHM, evo_cutoff = 100, plb_cutoff = 100, scheme = 'double', addition = NULL
){
  ggplot(
    aes(x = spec_plb, y = sens_plb), 
    data = dfHM[dfHM$scheme == scheme&dfHM$evo_cutoff == evo_cutoff&dfHM$plb_cutoff == plb_cutoff,]
  ) +
    geom_tile(aes(fill = hr), colour = 'white') +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0.6, limits = c(0,1.5)) + 
    scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
    scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
    labs(
      x = 'Specificity for Placebo Arm', 
      y = 'Sensitivity for Placebo Arm', 
      title = paste0(
        'HR given performance in PLACEBO Arm', 
        '\n(evo cutoff = ', evo_cutoff, '; plb cutoff = ', plb_cutoff, '; scheme = "', scheme, '"', addition, ')',
        '\nwhite = 0.6'
      )
    ) + 
    theme(
      title = element_text(size = 15), 
      axis.title = element_text(size = 15), 
      axis.text = element_text(size = 15)
    )
}

################################################################################################################
################################################################################################################
plotRatioEvo <- function(
  dfHM, evo_cutoff = 100, plb_cutoff = 100, scheme = 'double', addition = NULL
){
  ggplot(
    aes(x = spec_evo, y = sens_evo), 
    data = dfHM[dfHM$scheme == scheme&dfHM$evo_cutoff == evo_cutoff&dfHM$plb_cutoff == plb_cutoff,]
  ) +
    geom_tile(aes(fill = surRatio), colour = 'white') +
    scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue', midpoint = 1.3, limits = c(0,5)) + 
    scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
    scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
    labs(
      x = 'Specificity for Evo Arm', 
      y = 'Sensitivity for Evo Arm', 
      title = paste0(
        'median survival ratio given performance in EVO Arm',
        '\n(evo cutoff = ', evo_cutoff, '; plb cutoff = ', plb_cutoff, '; scheme = "', scheme, '"', addition, ')',
        '\nwhite = 1.3'
      )
    ) + 
    theme(
      title = element_text(size = 15), 
      axis.title = element_text(size = 15), 
      axis.text = element_text(size = 15)
    )
}

################################################################################################################
################################################################################################################
plotRatioPlacebo <- function(
  dfHM, evo_cutoff = 100, plb_cutoff = 100, scheme = 'double', addition = NULL
){
  ggplot(
    aes(x = spec_plb, y = sens_plb), 
    data = dfHM[dfHM$scheme == scheme&dfHM$evo_cutoff == evo_cutoff&dfHM$plb_cutoff == plb_cutoff,]
  ) +
    geom_tile(aes(fill = surRatio), colour = 'white') +
    scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue', midpoint = 1.3, limits = c(0,5)) + 
    scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
    scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
    labs(
      x = 'Specificity for Placebo Arm', 
      y = 'Sensitivity for Placebo Arm', 
      title = paste0(
        'median survival ratio given performance in PLACEBO Arm',
        '\n(evo cutoff = ', evo_cutoff, '; plb cutoff = ', plb_cutoff, '; scheme = "', scheme, '"', addition, ')',
        '\nwhite = 1.3'
      )
    ) + 
    theme(
      title = element_text(size = 15), 
      axis.title = element_text(size = 15), 
      axis.text = element_text(size = 15)
    )
}

################################################################################################################
################################################################################################################
plotLogRankPEvo <- function(
  dfHM, evo_cutoff = 100, plb_cutoff = 100, scheme = 'double', addition = NULL
){
  ggplot(
    aes(x = spec_evo, y = sens_evo), 
    data = dfHM[dfHM$scheme == scheme&dfHM$evo_cutoff == evo_cutoff&dfHM$plb_cutoff == plb_cutoff,]
  ) +
    geom_tile(aes(fill = lrP), colour = 'white') +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0.05, limits = c(0, 0.5)) + 
    scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
    scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
    labs(
      x = 'Specificity for Evo Arm', 
      y = 'Sensitivity for Evo Arm', 
      title = paste0(
        'log-rank p given performance in EVO Arm',
        '\n(evo cutoff = ', evo_cutoff, '; plb cutoff = ', plb_cutoff, '; scheme = "', scheme, '"', addition, ')',
        '\nwhite = 0.05'
      )
    ) + 
    theme(
      title = element_text(size = 15), 
      axis.title = element_text(size = 15), 
      axis.text = element_text(size = 15)
    )
}

################################################################################################################
################################################################################################################
plotLogRankPPlacebo <- function(
  dfHM, evo_cutoff = 100, plb_cutoff = 100, scheme = 'double', addition = NULL
){
  ggplot(
    aes(x = spec_plb, y = sens_plb), 
    data = dfHM[dfHM$scheme == scheme&dfHM$evo_cutoff == evo_cutoff&dfHM$plb_cutoff == plb_cutoff,]
  ) +
    geom_tile(aes(fill = lrP), colour = 'white') +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0.05, limits = c(0, 0.5)) + 
    scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
    scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
    labs(
      x = 'Specificity for Placebo Arm', 
      y = 'Sensitivity for Placebo Arm', 
      title = paste0(
        'log-rank p given performance in PLACEBO Arm',
        '\n(evo cutoff = ', evo_cutoff, '; plb cutoff = ', plb_cutoff, '; scheme = "', scheme, '"', addition, ')',
        '\nwhite = 0.05'
      )
    ) + 
    theme(
      title = element_text(size = 15), 
      axis.title = element_text(size = 15), 
      axis.text = element_text(size = 15)
    )
}

################################################################################################################
################################################################################################################
plotPosRate <- function(
  dfHM, evo_cutoff = 100, plb_cutoff = 100, scheme = 'double', addition = NULL
){
  ggplot(
    aes(x = sp, y = se), 
    data = dfHM[dfHM$scheme == scheme&dfHM$evo_cutoff == evo_cutoff&dfHM$plb_cutoff == plb_cutoff,]
  ) +
    geom_tile(aes(fill = pos_rate), colour = 'white') +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0.5) + 
    scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
    scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
    labs(
      x = 'Specificity for Placebo Arm', 
      y = 'Sensitivity for Placebo Arm', 
      title = paste0(
        'Test Positive Rate',
        '\n(evo cutoff = ', evo_cutoff, '; plb cutoff = ', plb_cutoff, '; scheme = "', scheme, '"', addition,')',
        '\nwhite = 0.5'
      )
    ) + 
    theme(
      title = element_text(size = 15), 
      axis.title = element_text(size = 15), 
      axis.text = element_text(size = 15)
    )
}