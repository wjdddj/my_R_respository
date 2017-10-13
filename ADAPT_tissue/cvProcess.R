generateParameters <- function(
  Methods = c('cv', 'LOOCV'),
  Ks = c(5, 10, 20)
){
  df_params <- expand.grid(method = Methods, K = Ks)
  df_params$K[df_params$method == 'LOOCV'] <- NA
  df_params <- unique(df_params)
  df_params <- df_params[order(df_params$method), ]
  rownames(df_params) <- NULL
  paramSets <- apply(df_params, 1, function(x){
    list(
      method = as.character(x[1]), 
      K = as.numeric(x[2])
    )
  })
}

cvProcess <- function(
  df_train,
  paramSets,
  R = 100,
  doCIAUC = T, 
  Permute = F,
  nPerm = 1000,
  nThread = 8
){
  
  results <- list()
  for(i in 1:length(paramSets)){
    method = paramSets[[i]]$method
    K = paramSets[[i]]$K
    cat('method =', method, '; K =', K, '\n')
    
    ## perform cross validation
    trueAUCs <- c()
    cvfits <- runCV(
      df_train, 
      method = method, 
      K = K, 
      R = R, 
      Parallel = T, 
      nThread = nThread
    )
    trueAUCs <- sapply(cvfits, function(x){
      auc(x$pred$obs, x$pred$R, direction = '<')
    })
    trueAUCs <- round(trueAUCs, 3)
    
    ## save represented ROC curve
    cvfit <- cvfits[[which.min(abs(trueAUCs - mean(trueAUCs, na.rm = T)))]]
    
    ## average cross validation score
    preds <- lapply(cvfits, function(x) x$pred[,c('R', 'rowIndex')])
    pred <- Reduce(
      function(x, y){
        merge(
          x, y, 
          by = 'rowIndex'
        )
      },
      preds
    )
    mean_pred <- data.frame(
      rowIndex = pred$rowIndex,
      R = apply(pred[,2:ncol(pred)], 1, mean, na.rm = T)
    )
    mean_pred <- merge(
      mean_pred,
      cvfit$pred[,c('rowIndex', 'obs')],
      by = 'rowIndex'
    )
    
    ci_AUC <- round(ci.auc(
      auc(cvfit$finalModel$data$.outcome, cvfit$finalModel$fitted.values, direction = '<')
    ), 3)
    ## calculate CI of cvAUC
    ci_cvAUC <- NULL
    if(doCIAUC){
      AUC <- auc(mean_pred$obs, mean_pred$R, direction = '<')
      ci_cvAUC <- round(ci.auc(AUC), 3)
    }
    
    ## perform permutation on cv
    PermCV = permAUCs <- NULL
    if(Permute){
      PermCV <- permuteCV(
        df_train, method = method,
        K = K, nPerm = nPerm, nThread = nThread
      )
      permAUCs <- sapply(PermCV, function(x){
        AUC <- auc(x$permLabel, x$R, direction = '<')
      })
    }

    result <- list(
      method = as.character(method),
      K = K,
      pred = mean_pred,
      trueAUCs = trueAUCs,
      PermCV = PermCV,
      permAUCs = permAUCs,
      cvfits = cvfits,
      cvfit = cvfit,
      ci_AUC = ci_AUC,
      ci_cvAUC = ci_cvAUC
    )
    results[[i]] <- result
  }
  results
}


plotCVROC <- function(cvfit, main){
  #cvfit <- results[[1]]$cvfit
  cvfit$pred$obs <- factor(cvfit$pred$obs, levels = c('NR', 'R'))
  auc <- auc(cvfit$pred$obs, cvfit$pred$R, direction = '<')
  auc <- round(auc, 3)
  roc(
    cvfit$pred$obs, cvfit$pred$R, direction = '<',
    plot = T,
    legacy.axes = T, cex = 2, 
    main = paste0(main, '\nAUC = ', auc)
  )
}

plotROC <- function(finalModel, main){
  finalModel$data$.outcome <- factor(finalModel$data$.outcome, levels = c('NR', 'R'))
  auc <- auc(finalModel$data$.outcome, finalModel$fitted.values, direction = '<')
  auc <- round(auc, 3)
  roc(
    finalModel$data$.outcome, finalModel$fitted.values, direction = '<',
    plot = T,
    legacy.axes = T, cex = 2, 
    main = paste0(main, '\nAUC = ', auc)
  )
}

plotPermAUC <- function(result, main){
  meanTrueAUC <- mean(result$trueAUCs, na.rm = T)
  pvalue <- sum(result$permAUCs > meanTrueAUC)/length(result$permAUCs)
  hist(
    result$permAUCs, breaks = seq(0, 1, 0.05),
    xlim = c(0, 1),
    xlab = 'AUCs with permuted labels',
    ylab = 'Frequency',
    cex = 2, 
    main = paste0(main, '\nP = ', pvalue)
  )
  abline(v = meanTrueAUC, col = 'red', lwd = 2)
}




