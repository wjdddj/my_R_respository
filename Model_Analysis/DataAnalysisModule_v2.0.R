##################################################################################################################################
#### MODULE 1: Feature filtering
##################################################################################################################################

##################################################################################################################################
# base functions for univariate analysis
# Calculate AUC
pauc = function(data,min.spec=0){
  intensity = as.numeric(data[,2])
  cc = data[,1]
  w.cases = which(cc==1&intensity!="NA")
  w.controls = which(cc==0&intensity!="NA")
  x = intensity[w.cases]
  y = intensity[w.controls]
  nx = length(x)
  ny = length(y)  
  u = 0
  for(i in 1:nx){
    area = (length(which(y<x[i]))+length(which(y==x[i]))/2)/ny
    area = area-min.spec
    if(area<0){area=0} 
    u = area*ny + u 
  }
  u = u/(nx*ny)
  z <- list(pAUC=u, diag.area = 0.5*(1-min.spec)^2, specificity=min.spec)
}

# Calculate sensitivity
sens <- function(data, spec=0.95){
  intensity = as.numeric(data[,2])
  cc = data[,1]
  w.cases = which(cc==1)
  w.controls = which(cc==0)
  x = intensity[w.cases]
  y = intensity[w.controls]
  nx = length(x)
  ny = length(y)
  q = quantile(y,spec)
  #q = max(y)
  sensitivity = length(intensity[which(cc==1&intensity>=q)])/nx
  specificity = length(intensity[which(cc==0&intensity<q)])/ny
  z <- list(sensitivity=sensitivity,specificity=specificity,cutoff=q)
}

# Calculate sensitivity
sens_low <- function(data, spec=0.95){
  intensity = as.numeric(data[,2])
  cc = data[,1]
  w.cases = which(cc==1)
  w.controls = which(cc==0)
  x = intensity[w.cases]
  y = intensity[w.controls]
  nx = length(x)
  ny = length(y)
  q = quantile(y,1-spec)
  #q = max(y)
  sensitivity = length(intensity[which(cc==1&intensity<=q)])/nx
  specificity = length(intensity[which(cc==0&intensity>q)])/ny
  z <- list(sensitivity=sensitivity,specificity=specificity,cutoff=q)
}


# Calculate specificity
spec <- function(data, sens=0.95){
  intensity = as.numeric(data[,2])
  cc = data[,1]
  w.cases = which(cc==1&intensity!="NA")
  w.controls = which(cc==0&intensity!="NA")
  x = intensity[w.cases]
  y = intensity[w.controls]
  nx = length(x)
  ny = length(y)
  q = quantile(x,sens)
  sensitivity = length(intensity[which(cc==1&intensity<=q)])/nx
  specificity = length(intensity[which(cc==0&intensity>q)])/ny
  z <- list(sensitivity=sensitivity,specificity=specificity,cutoff=q)
}

# Calculate specificity
spec_low <- function(data, sens=0.95){
  intensity = as.numeric(data[,2])
  cc = data[,1]
  w.cases = which(cc==1&intensity!="NA")
  w.controls = which(cc==0&intensity!="NA")
  x = intensity[w.cases]
  y = intensity[w.controls]
  nx = length(x)
  ny = length(y)
  q = quantile(x,1-sens)
  sensitivity = length(intensity[which(cc==1&intensity>=q)])/nx
  specificity = length(intensity[which(cc==0&intensity<q)])/ny
  z <- list(sensitivity=sensitivity,specificity=specificity,cutoff=q)
}

# Calculate sensitivity using preset cutoff
ss_pre <- function(data, cutoff=0.95){
  intensity = as.numeric(data[,2])
  cc = data[,1]
  w.cases = which(cc==1&intensity!="NA")
  w.controls = which(cc==0&intensity!="NA")
  x = intensity[w.cases]
  y = intensity[w.controls]
  nx = length(x)
  ny = length(y)
  sensitivity = length(intensity[which(cc==1&intensity>=cutoff)])/nx
  specificity = length(intensity[which(cc==0&intensity<cutoff)])/ny
  accuracy = (sensitivity * nx + specificity * ny)/(nx + ny)
  z <- list(sensitivity=sensitivity,specificity=specificity, accuracy = accuracy)
}

# K statistics
k.stat <- function(data,x1=0.975,x2=0.8,y1=0.975,y2=0.025){
  intensity = as.numeric(data[,2])
  cc = data[,1]
  w.cases = which(cc==1&intensity!="NA")
  w.controls = which(cc==0&intensity!="NA")
  x = intensity[w.cases]
  y = intensity[w.controls]
  nx = length(x)
  ny = length(y)
  k = (quantile(x,c(x1,x2))[1]-quantile(x,c(x1,x2))[2])/(quantile(y,c(y1,y2))[1]-quantile(y,c(y1,y2))[2])
  z <- list(k=k, case=c(x1,x2), control=c(y1,y2))
}

##################################################################################################################################
# main function for feature filtering
# Initial filter features using simple univariate stats
# input dataFrame: last column is the response
##################################################################################################################################
feature_filtering <- function(
  dataFrame, na.rm = T, 
  ks_test = T, t_test = T, wilcox_test = T, cutoff_pValue = 0.05, 
  FC = T, cutoff_FC = 1.2, 
  spec = F, sens = F, spec_low = F, sens_low = F, cutoffCalculation = 0.95, cutoffSelection = 0.1,
  pAUC_spec = F, pAUC_sens = F, cutoffPAUC = 0.95, selCoef = 2,
  k_stat = F
){
  
  library(pROC)
  n <- ncol(dataFrame)
  response <- dataFrame[,n]
  metrics <- data.frame(featureID = colnames(dataFrame)[-n])
  idxList <- list()
  
  ##################################################################################################################################
  ## group stats
  if(ks_test|t_test|wilcox_test){
    # cutoff_pValue = 0.05
    cat(paste('p value cutoff is ', cutoff_pValue, '\n', sep = ''))
  }
  
  if(ks_test){
    pVecKStest <- apply(dataFrame[,-n], 2, function(x){
      response1 <- response
      if(na.rm){
        response1 = response1[!is.na(x)]
        x = x[!is.na(x)]
      }
      test <- try(ks.test(x[response1==0|response1=='control'], x[response1==1|response1=='case']), silent = T)
      result <- ifelse(is(test, 'try-error'), NA, test$p.value)
    })
    idx_ks <- which(pVecKStest <= cutoff_pValue)
    idxList[[length(idxList)+1]] <- idx_ks
    metrics <- data.frame(metrics, pKStest = pVecKStest)
  }
  
  if(wilcox_test){
    pVecWiltest <- apply(dataFrame[,-n], 2, function(x){
      response1 <- response
      if(na.rm){
        response1 = response1[!is.na(x)]
        x = x[!is.na(x)]
      }
      test <- try(wilcox.test(x[response1==0|response1=='control'], x[response1==1|response1=='case']), silent = T)
      result <- ifelse(is(test, 'try-error'), NA, test$p.value)
    })
    idx_wtest <- which(pVecWiltest <= cutoff_pValue)
    idxList[[length(idxList)+1]] <- idx_wtest
    metrics <- data.frame(metrics, pWtest = pVecWiltest)
  }
  
  if(t_test){
    pVecTtest <- apply(dataFrame[,-n], 2, function(x){
      response1 <- response
      if(na.rm){
        response1 = response1[!is.na(x)]
        x = x[!is.na(x)]
      }
      test <- try(t.test(x[response1==0|response1=='control'], x[response1==1|response1=='case'], var.equal = FALSE), silent = T)
      result <- ifelse(is(test, 'try-error'), NA, test$p.value)
    })
    idx_ttest <- which(pVecTtest <= cutoff_pValue)
    idxList[[length(idxList)+1]] <- idx_ttest
    metrics <- data.frame(metrics, pTtest = pVecTtest)
  }
  
  ##################################################################################################################################
  ## fold change
  if(FC){
    # cutoff_FC = 1.2
    cat(paste('fold change cutoff is ', cutoff_FC, '\n', sep =''))
    
    foldChange_log <- apply(dataFrame[,-n], 2, function(x){
      response1 <- response
      if(na.rm){
        response1 = response1[!is.na(x)]
        x = x[!is.na(x)]
      }
      log2(mean(x[response1==1|response1=='case'])/mean(x[response1==0|response1=='control']))
    })
    idx_fc <- which(abs(foldChange_log)>log2(cutoff_FC))
    idxList[[length(idxList)+1]] <- idx_fc
    metrics <- data.frame(metrics, FClog = foldChange_log)
  }
  
  ##################################################################################################################################
  ## sensitivity and specificity
  if(spec|sens|spec_low|sens_low){
    # cutoffCalculation = 0.95
    # cutoffSelection = 0.09
    cat(paste('sensitivity and specificity cutoff is ', cutoffCalculation, '\n', sep = ''))
    cat(paste('selection cutoff is ', cutoffSelection, '\n', sep = ''))
  }
  
  if(spec){
    spec_100sens <- apply(dataFrame[,-n], 2, function(x){
      spec(data.frame(response, x), sens = cutoffCalculation)$specificity
    })
    idx_spec <- which(spec_100sens > cutoffSelection)
    if(length(idx_spec)==0){
      cat('not enough feature\n')
    }
    idxList[[length(idxList)+1]] <- idx_spec
    metrics <- data.frame(metrics, specAT100sens = spec_100sens)
  }
  
  if(spec_low){
    specLow_100sens <- apply(dataFrame[,-n], 2, function(x){
      spec_low(data.frame(response, x), sens = cutoffCalculation)$specificity
    })
    idx_specLow <- which(specLow_100sens > cutoffSelection)
    if(length(idx_specLow)==0){
      cat('not enough feature\n')
    }
    idxList[[length(idxList)+1]] <- idx_specLow
    metrics <- data.frame(metrics, specLowAT100sens = specLow_100sens)
  }
  
  if(sens){
    sens_100spec <- apply(dataFrame[,-n], 2, function(x){
      sens(data.frame(response, x), spec = cutoffCalculation)$sensitivity
    })
    idx_sens <- which(sens_100spec > cutoffSelection)
    if(length(idx_sens)==0){
      cat('not enough feature\n')
    }
    idxList[[length(idxList)+1]] <- idx_sens
    metrics <- data.frame(metrics, sensAT100spec = sens_100spec)
  }
  
  if(sens_low){
    sensLow_100spec <- apply(dataFrame[,-n], 2, function(x){
      sens_low(data.frame(response, x), spec = cutoffCalculation)$sensitivity
    })
    idx_sensLow <- which(sensLow_100spec > cutoffSelection)
    if(length(idx_sensLow)==0){
      cat('not enough feature\n')
    }
    idxList[[length(idxList)+1]] <- idx_sensLow
    metrics <- data.frame(metrics, sensLowAT100spec = sensLow_100spec)
  }
  
  ##################################################################################################################################
  ## pAUC
  if(pAUC_spec|pAUC_sens){
    # cutoffPAUC = 0.95
    areaRandom = 0.5*(1-cutoffPAUC)^2
    # selCoef = 2
    areaDiagnos = selCoef * areaRandom
    cat(paste('pAUC cutoff is ', cutoffPAUC, '\n', sep = ''))
    cat(paste('random diagnostic AUC is ', areaRandom, '\n', sep = ''))
    cat(paste('diagnostic AUC cutoff is ', areaDiagnos, '\n', sep = ''))
  }
  
  if(pAUC_spec){
    pAUC_highspec <- apply(dataFrame[,-n], 2, function(x){
      auc(response, x, partial.auc = c(1, cutoffPAUC), partial.auc.focus = 'specificity')
    })
    idx_pAUC_highspec <- which(pAUC_highspec > areaDiagnos)
    # idx_pAUC_highspec <- order(pAUC_highspec, decreasing = T)[1:100]
    idxList[[length(idxList)+1]] <- idx_pAUC_highspec
    metrics <- data.frame(metrics, pAUC_highspec = pAUC_highspec)
  }
  
  if(pAUC_sens){
    pAUC_highsens <- apply(dataFrame[,-n], 2, function(x){
      auc(response, x, partial.auc = c(1, cutoffPAUC), partial.auc.focus = 'sensitivity')
    })
    idx_pAUC_highsens <- which(pAUC_highsens > areaDiagnos)
    # idx_pAUC_highsens <- order(pAUC_highsens, decreasing = T)[1:100]
    idxList[[length(idxList)+1]] <- idx_pAUC_highsens
    metrics <- data.frame(metrics, pAUC_highsens = pAUC_highsens)
  }
  
  ##################################################################################################################################
  ## k stat
  if(k_stat){
    caseTop = 0.975
    caseBot = 0.8
    controlTop = 0.975
    controlBot = 0.025
    kStatCutoff = 1.0
    cat(paste('case span is calculated between ', caseBot, ' quantile and ', caseTop, ' quantile\n', sep = ''))
    cat(paste('control span is calculated between ', controlBot, ' quantile and ', controlTop, ' quantile\n', sep = ''))
    cat(paste('selection cutoff is ', kStatCutoff, '\n', sep = ''))
    
    kStat <- apply(dataFrame[,-n], 2, function(x){
      k.stat(data.frame(dataFrame[,n], x), x1 = caseTop, x2 = caseBot, y1 = controlTop, y2 = controlBot)$k
    })
    idx_kStat <- which(kStat > kStatCutoff)
    idxList[[length(idxList)+1]] <- idx_kStat
    metrics <- data.frame(metrics, kStat = kStat)
  }
  
  names(idxList) <- colnames(metrics)[-1]
  rownames(metrics) <- NULL
  idx_sig <- unique(unlist(idxList))
  cat(paste('total ', length(idx_sig), ' features were selected\n', sep = ''))
  
  z <- list(idxList = idxList,
            idx_sig = idx_sig,
            metrics = metrics)
}

##################################################################################################################################
#### MODULE 2: Modeling by Random Forest
#### feature selection is implemented by either statistic filter or variable importance
#### functions with cross-validation, done in parallel
##################################################################################################################################

##################################################################################################################################
# random forest model building using feature_filtering()
##################################################################################################################################
randomF_featureFilter <- function(dataFrame, ntree = 1000, mtry = NULL, replace = F, featureFilter = T){
  
  library(randomForest)
  n <- ncol(dataFrame)
  response <- as.character(dataFrame[,n])
  response <- as.factor(response)
  dataFrame <- data.frame(dataFrame[,-n], response=response)  
  n.sample <- nrow(dataFrame)
  
  if(featureFilter){
    cat('filtering features\n')
    selFeature <- feature_filtering(dataFrame, ks_test = T, t_test = T, wilcox_test = T, cutoff_pValue = 0.05, 
                                    FC = T, cutoff_FC = 1.2, 
                                    spec = F, sens = F, spec_low = F, sens_low = F, cutoffCalculation = 0.95, cutoffSelection = 0.1,
                                    pAUC_spec = F, pAUC_sens = F, cutoffPAUC = 0.95, selCoef = 2,
                                    k_stat = F)
    idx_sig <- selFeature$idx_sig
    idx_newDF <- c(idx_sig, n)
    dataFrame_fs <- dataFrame[,idx_newDF]
  }
  else{
    dataFrame_fs <- dataFrame
  }
  
  cat('fit random forest model\n')
  
  if(is.null(mtry)){
    fitRF <- randomForest(response~., data = dataFrame_fs, ntree = ntree, replace = replace)
  }else{
    fitRF <- randomForest(response~., data = dataFrame_fs, ntree = ntree, mtry = mtry, replace = replace)
  }
  
  predOut <- list(cbind(c(1:n.sample),fitRF$votes))
  
  if(featureFilter){
    z <- list(predOut = predOut, 
              sigFeatureIndex = idx_sig, 
              sigFeature = colnames(dataFrame)[idx_sig], 
              selFeature = selFeature, 
              RFobj = fitRF)  
  }else{
    z <- list(predOut = predOut, RFobj = fitRF)
  }
  z
}


##################################################################################################################################
# original cross validation with random forest model
# feature selection: Statistic filtering
# with repeat cross validation
##################################################################################################################################
randomF_cv1 <- function(wholeDataframe, featureFilter = T, K=10, R=1, ntree=1000, mtry = NULL, replace = F){
  
  library(cvTools)
  
  n <- ncol(wholeDataframe)
  response <- wholeDataframe[,n]
  response <- as.factor(response)
  wholeDataframe <- data.frame(wholeDataframe[,-n], response=response)
  
  #prepare folds for cross-validation
  n.sample <- nrow(wholeDataframe)
  predOut <- list()
  sigFeatureByRepeat <- list()
  
  for(j in 1:R){
    
    cat(paste('\nround', j, 'of', R, '\n'))
    subset <- cvFolds(n.sample, K)
    predAll <- c()
    sigFeatureByFold <- list()
    
    for(i in 1:K){
      cat(paste('fold', i, 'of', K, '\n'))
      #select features by filtering the training set
      idx_test <- subset$subsets[subset$which==i]
      train_df <- wholeDataframe[-idx_test,]
      fit <- randomF_featureFilter(train_df, ntree = 1000, mtry = NULL, replace = F, featureFilter = featureFilter)
      
      test_df <- wholeDataframe[idx_test,c(fit$sigFeatureIndex, n)]
      pred <- predict(fit$RFobj, newdata=test_df, type='prob')
      
      #record predicted probabilities and sample index
      outTab <- cbind(idx_test, pred)
      predAll <- rbind(predAll, outTab)
      sigFeatureByFold[[length(sigFeatureByFold)+1]] <- fit$sigFeatureIndex
    }
    
    predOut[[length(predOut)+1]] <- predAll
    sigFeatureByRepeat[[length(sigFeatureByRepeat)+1]] <- sigFeatureByFold
  }
  
  z <- list(predOut = predOut, 
            sigFeatureByRepeat = sigFeatureByRepeat)
}


##################################################################################################################################
# cross validation using random forest model, done in parallel
# number of threads is determined by number of folds
# feature selection: Statistic filtering, related parameters are set in randomF_featureFilter()
# with repeat cross validation
##################################################################################################################################
randomF_cv1_para <- function(wholeDataframe, featureFilter = T, K=10, R=1, ntree=1000, mtry = NULL, replace = F){
  
  library(foreach)
  library(doParallel)
  library(cvTools)
  cl <- makeCluster(K)
  registerDoParallel(cl)
  
  n <- ncol(wholeDataframe)
  response <- wholeDataframe[,n]
  response <- as.factor(response)
  wholeDataframe <- data.frame(wholeDataframe[,-n], response=response)
  
  #prepare folds for cross-validation
  n.sample <- nrow(wholeDataframe)
  predOut <- list()
  sigFeatureByRepeat <- list()
  
  cvCombine1 <- function(x, y){
    Out1 <- rbind(x[[1]], y[[1]])     # x and y must be list
    Out2 <- c(x[[2]], y[[2]])     # x and y must be list
    z <- list(predAll = Out1, sigFeatureByFold = Out2)
  }
  
  for(j in 1:R){
    
    cat(paste('\nround', j, 'of', R, '\n'))
    subset <- cvFolds(n.sample, K)
    predAll <- c()
    sigFeatureByFold <- list()
    
    crossVal <- foreach(i = 1:K, 
                        .combine = cvCombine1, 
                        .export = c('feature_filtering', 'randomF_featureFilter')) %dopar% {
                          
                          idx_test <- subset$subsets[subset$which==i]
                          train_df <- wholeDataframe[-idx_test,]
                          #fit model using selected features 
                          
                          fit <- randomF_featureFilter(train_df, ntree=ntree, mtry=mtry, replace = replace, featureFilter = featureFilter)
                          
                          # test on holdout set
                          test_df <- wholeDataframe[idx_test,c(fit$sigFeatureIndex, n)]
                          pred <- predict(fit$RFobj, newdata=test_df, type='prob')
                          
                          # record predicted probabilities and sample index
                          outTab <- cbind(idx_test, pred)
                          sigFeatureByFold[[length(sigFeatureByFold)+1]] <- fit$sigFeatureIndex
                          
                          paraOut <- list(predAll = outTab, sigFeatureByFold = sigFeatureByFold)
                        }
    
    predOut[[length(predOut)+1]] <- crossVal$predAll
    sigFeatureByRepeat[[length(sigFeatureByRepeat)+1]] <- crossVal$sigFeatureByFold
  }
  stopCluster(cl)
  z <- list(predOut = predOut, 
            sigFeatureByRepeat = sigFeatureByRepeat)
}


##################################################################################################################################
# random forest model using varSelRF() - variable importance for feature selection
##################################################################################################################################
randomRF_varSelRF <- function(runDataFrame, nInitialTree = 5000, ntreeIterat = 1000, vars.drop.frac = 0.2, recompute.var.imp = F){
  
  library(varSelRF)
  
  n <- ncol(runDataFrame)
  response <- as.factor(runDataFrame[,n])
  cat('performing variable selection...\n')
  varSelRF1 <- varSelRF(runDataFrame[,-n], response, 
                        c.sd=1, mtryFactor = 2, 
                        ntree = nInitialTree, ntreeIterat = ntreeIterat, 
                        vars.drop.num = NULL, vars.drop.frac = vars.drop.frac, whole.range = TRUE, 
                        recompute.var.imp = recompute.var.imp, verbose = TRUE,
                        returnFirstForest = TRUE, fitted.rf = NULL, 
                        keep.forest = FALSE)
  selVar <- varSelRF1$selected.var
  sfIdx <- which(colnames(runDataFrame)%in%selVar)
  runDFSF <- data.frame(runDataFrame[,sfIdx], response = response)
  
  cat('building final model...\n')
  fitRF <- randomForest(response~., data = runDFSF, ntree = ntreeIterat, mtry = ceiling(sqrt(length(sfIdx))), replace = F)
  z <- list(variableSel = varSelRF1, 
            fittedModel = fitRF, 
            predOut = fitRF$votes[,2], 
            sfIdx = sfIdx, 
            selectedVar = selVar)
}

##################################################################################################################################
# cross validation with random forest model
# feature selection: varSelRF-importance
# with repeat cross validation
##################################################################################################################################
randomF_cv2 <- function(wholeDataframe, 
                        K = 10, R = 1, 
                        ntree = 1000, nInitialTree = 5000, ntreeIterat = 2000, 
                        vars.drop.frac = 0.3, recompute.var.imp = T, replace = F){
  
  library(cvTools)
  library(varSelRF)
  
  n <- ncol(wholeDataframe)
  response <- as.character(wholeDataframe[,n])
  response[response==1] <- 'case'
  response[response==0] <- 'control'
  response <- as.factor(response)
  wholeDataframe <- cbind(wholeDataframe[,-n],response=response)  
  n.sample <- nrow(wholeDataframe)
  predOut <- list()
  selectedVar <- list()
  for(j in 1:R){
    cat(paste('round 1 of', R, '\n'))
    #prepare folds for cross-validation
    subset <- cvFolds(n.sample, K)
    predAll <- c()
    selectedVarIter <- list()
    for(i in 1:K){
      #select features by filtering the training set
      idx_test <- subset$subsets[subset$which==i]
      train_df <- wholeDataframe[-idx_test,]
      #idx_selectedFeature <- c(1,feature_filtering(train_df)$idx_sig)
      #train_rf_df <- train_df[,idx_selectedFeature]
      #test_rf_df <- wholeDataframe[idx_test,idx_selectedFeature]
      
      selVar <- varSelRF(train_df[,-n], as.factor(train_df[,n]), 
                         c.sd=1, mtryFactor = 2, 
                         ntree = nInitialTree, ntreeIterat = ntreeIterat, 
                         vars.drop.num = NULL, vars.drop.frac = vars.drop.frac, whole.range = TRUE, 
                         recompute.var.imp = recompute.var.imp, verbose = TRUE,
                         returnFirstForest = TRUE, fitted.rf = NULL, 
                         keep.forest = FALSE)$selected.var
      
      idx_selectedFeature <- c(which(colnames(train_df)%in%selVar), n)
      train_rf_df <- train_df[,idx_selectedFeature]
      test_rf_df <- wholeDataframe[idx_test,idx_selectedFeature]
      
      #fit model using selected features, and test on holdout set 
      
      fit <- randomForest(response~., data=train_rf_df, ntree=ntree, mtry=ceiling(sqrt(length(idx_selectedFeature)-1)), replace = replace)
      pred <- predict(fit, newdata=test_rf_df, type='prob')
      
      #record predicted probabilities and sample index
      outTab <- cbind(idx_test, pred)
      predAll <- rbind(predAll, outTab)
      selectedVarIter[[length(selectedVarIter)+1]] <- selVar
      cat(paste(i, ',', sep=''))
    }
    selectedVar[[length(selectedVar)+1]] <- selectedVarIter
    predOut[[length(predOut)+1]] <- predAll
    
  }
  
  z <- list(predOut = predOut, 
            selectedVar = selectedVar)
}


##################################################################################################################################
# cross validation with random forest model
# feature selection: varSelRF-importance
# with repeat cross validation
# with parallel computing
# number of threads is determined by number of folds
##################################################################################################################################
randomF_cv2_para <- function(wholeDataframe, 
                             K = 10, R = 1, 
                             ntree = 1000, nInitialTree = 5000, ntreeIterat = 1000, 
                             vars.drop.frac = 0.2, recompute.var.imp = F, replace = F){
  library(randomForest)
  library(foreach)
  library(doParallel)
  library(cvTools)
  
  n <- ncol(wholeDataframe)
  response <- as.factor(wholeDataframe[,n])
  wholeDataframe <- data.frame(wholeDataframe[,-n], response=response)  
  n.sample <- nrow(wholeDataframe)
  predOut <- list()
  
  # A combine function to combine the result from parallel computing
  cvCombine <- function(x, y){
    Out1 <- rbind(x[[1]], y[[1]])     # x and y must be list
    Out2 <- c(x[[2]], y[[2]])     # x and y must be list
    Out3 <- c(x[[3]], y[[3]])     # x and y must be list
    z <- list(predOut = Out1, selectedVar = Out2, sfIdx = Out3)
  }
  
  cat('start cross-validation...\n\n')
  for(j in 1:R){
    
    cat(paste('Repeat 1 of', R, '\n', sep = ' '))
    cat(paste(Sys.time(), '\n', sep=''))
    # prepare folds for cross-validation
    subset <- cvFolds(n.sample, K)
    predAll <- c()
    selectedVarIter <- list()
    sfIdxIter <- list()
    cl <- makeCluster(K)
    registerDoParallel(cl)
    crossVal <- foreach(i = 1:K, 
                        .combine = cvCombine,
                        .export = c('randomRF_varSelRF')) %dopar% {
                          
                          library(varSelRF)
                          
                          idx_test <- subset$subsets[subset$which==i]
                          train_df <- wholeDataframe[-idx_test,] 
                          trainResponse <- as.factor(train_df[,n])
                          
                          # select feature and fit model using selected features  
                          fitVarSelRF <- randomRF_varSelRF(train_df, nInitialTree = nInitialTree, ntreeIterat = ntreeIterat, vars.drop.frac = 0.2, recompute.var.imp = F)
                          
                          sfIdx <- fitVarSelRF$selectedVar
                          train_rf_df <- data.frame(train_df[,sfIdx], trainResponse)
                          test_rf_df <- wholeDataframe[idx_test,sfIdx]
                          
                          # test on holdout set
                          pred <- predict(fitVarSelRF$fittedModel, newdata=test_rf_df, type='prob')
                          
                          #record predicted probabilities and sample index
                          outTab <- cbind(idx_test, pred)
                          selectedVarIter[[length(selectedVarIter)+1]] <- fitVarSelRF$selectedVar
                          sfIdxIter[[length(sfIdxIter)+1]] <- fitVarSelRF$sfIdx
                          z1 <- list(outTab = outTab, selectedVarIter = selectedVarIter, sfIdxIter = sfIdxIter)
                        }
    
    cat(paste(Sys.time(), '\n', sep=' '))
    predOut[[length(predOut)+1]] <- crossVal
    
  }
  stopCluster(cl)
  predOut
}


##################################################################################################################################
# combine data from each cross validation during parallel computing using randomF_cv2_para()
##################################################################################################################################
fittedRF_cv <- function(fitRF_cv){
  R <- length(fitRF_cv)
  K <- length(fitRF_cv[[1]]$predOut)
  fittedAllRepeats <- list()
  for(i in 1:R){
    fittedFold <- c()
    for(j in 1:K){
      fittedFold <- rbind(fittedFold, fitRF_cv[[i]]$predOut[[j]])
    }
    fittedAllRepeats[[length(fittedAllRepeats)+1]] <- fittedFold
    names(fittedAllRepeats)[i] <- paste('Repeat', i)
  }
  fittedAllRepeats
}


##################################################################################################################################
# plot ROC for each cross-validation
# fittedRF is an object out put from 
#   randomF_featureFilter()$predOut, randomF_cv1()$predOut, randomF_cv2()$predOut, randomF_cv1_para()$predOut
##################################################################################################################################
fitRFplot <- function(fittedRF, runDataFrame, R=1, mfrow = c(1,1), cex.main = 3, cex.lab = 2, main = NULL, col = 'black'){
  n <- ncol(runDataFrame)
  library(pROC)
  response <- runDataFrame[,n]
  if(!is.null(mfrow))par(mfrow = mfrow)
  for(i in 1:R){
    if(is.null(main)){
      roc(response[fittedRF[[i]][,1]], fittedRF[[i]][,2], plot = T,
          main = round(roc(response[fittedRF[[i]][,1]], fittedRF[[i]][,2], direction = '<')$auc, 3),
          cex.main = cex.main, cex.lab = cex.lab, col = col, direction = '<')
    }
    else{
      roc(response[fittedRF[[i]][,1]], fittedRF[[i]][,2], plot = T,
          main = main,
          cex.main = cex.main, cex.lab = cex.lab, direction = '<')
    }
  }
}


##################################################################################################################################
#### MODULE 3: Permutation Analysis of randomForest AUCs
## @@@@@@ PROBLEM @@@@@@@
## @@@@@@ PROBLEM @@@@@@@ Package pROC automatically adjust the direction when calculating auc. Need to change to another package
## @@@@@@ PROBLEM @@@@@@@
###################################################################################################################################


##################################################################################################################################
# permutation test on OOB prediction generated AUC
# random forest, with feature selection by statistical filtering, without tuning
##################################################################################################################################
rf_permute <- function(wholeDataframe, featureFilter = T, r=100){
  n.sample <- nrow(wholeDataframe)
  y <- as.factor(wholeDataframe[,ncol(wholeDataframe)])
  fit_orig <- randomF_featureFilter(wholeDataframe, featureFilter = featureFilter, ntree = 1000, mtry = NULL, replace = F)
  auc_orig <- roc(y, fit_orig$predOut[[1]][,2], direction = '<')$auc
  aucAll <- auc_orig
  
  for(i in 1:r){
    y_perm <- sample(y)
    data_perm <- data.frame(wholeDataframe[,-ncol(wholeDataframe)], y_perm)
    fit_perm <- randomF_featureFilter(data_perm, featureFilter = featureFilter, ntree = 1000, mtry = NULL, replace = F)
    auc_perm <- roc(y_perm, fit_perm$predOut[[1]][,2], direction = '<')$auc
    aucAll <- c(aucAll, auc_perm)
    cat(i)
    cat(' ')
  }
  aucAll
}

##################################################################################################################################
# permutation test on OOB prediction generated AUC in Parallel using doParallel & foreach
# random forest, with feature selection by statistical filtering, without tuning
##################################################################################################################################
rf_permute_para <- function(wholeDataframe, featureFilter = T, r=100, thread = 10){
  
  library(doParallel)
  library(foreach)
  
  cl <- makeCluster(thread)
  registerDoParallel(cl)
  aucTab <- foreach(j = 1:thread, 
                    .combine = rbind,
                    .export = c('rf_permute', 'randomF_featureFilter', 'feature_filtering')) %dopar% {
                      library(randomForest)
                      library(pROC)                  
                      permAUC <- rf_permute(wholeDataframe, featureFilter = featureFilter, r = ceiling(r/thread))
                      
                    }
  stopCluster(cl)
  colnames(aucTab) <- c('original', c(1:ceiling(r/thread)))
  aucTab
}

##################################################################################################################################
# permutation test on 10-fold cross-validation generated AUC
# random forest, with feature selection by statistical filtering, without tuning
##################################################################################################################################
rf_permute_cv <- function(wholeDataframe, featureFilter = T, r=100){
  
  n.sample <- nrow(wholeDataframe)
  y <- as.factor(wholeDataframe[,ncol(wholeDataframe)])
  fitRFCV_ori <- randomF_cv1(wholeDataframe, featureFilter = featureFilter, K=10, R=1, ntree=1000, mtry = NULL, replace = F)
  auc_ori <- roc(y[fitRFCV_ori$predOut[[1]][,1]], fitRFCV_ori$predOut[[1]][,2], plot = F, direction = '<')$auc
  
  cat('cross-validation using original dataset: done!')
  cat('\n')
  
  aucAll <- auc_ori
  for(i in 1:r){
    y_perm <- sample(y)
    data_perm <- data.frame(wholeDataframe[,-ncol(wholeDataframe)], y_perm)
    fit_perm <- randomF_cv1(data_perm, featureFilter = featureFilter, K=10, R=1, ntree=1000, mtry = NULL, replace = F)
    auc_perm <- roc(y_perm[fit_perm$predOut[[1]][,1]], fit_perm$predOut[[1]][,2], plot = F, direction = '<')$auc
    aucAll <- c(aucAll, auc_perm)
    cat(paste(i, sep = ','))
    cat('\n')
  }
  aucAll
}


##################################################################################################################################
# permutation test on 10-fold cross-validation generated AUC in Parallel using doParallel & foreach
# random forest, with feature selection by statistical filtering, without tuning
##################################################################################################################################
# permutation test with cross validation in parallel
rf_permute_cv_para <- function(wholeDataframe, featureFilter = T, r=1000, thread = 20){
  
  library(doParallel)
  library(foreach)
  
  cl <- makeCluster(thread)
  registerDoParallel(cl)
  aucTab <- foreach(j = 1:thread, 
                    .combine = rbind, 
                    .export = c('rf_permute_cv', 'randomF_featureFilter', 'feature_filtering', 'randomF_cv1'))  %dopar% {
                      
                      library(randomForest)
                      library(pROC) 
                      permClust1_cv1 <- rf_permute_cv(wholeDataframe, featureFilter = featureFilter, r = ceiling(r/thread))
                    }
  stopCluster(cl)
  colnames(aucTab) <- c('original', c(1:ceiling(r/thread)))
  aucTab
}


##################################################################################################################################
#### MODULE 4: Permutation Analysis of Univariate test 
##################################################################################################################################

##################################################################################################################################
# Perform permutation test on number of significant aptamers after each univariate test
# method could be either 'ttest' or 'ANOVA'
##################################################################################################################################
UnivariatePermTest <- function(runDataFrame, method = 'ttest', filename, r = 1000, thread = 10){
  library(doParallel)
  library(foreach)
  cl <- makeCluster(thread)
  registerDoParallel(cl)
  
  n <- ncol(runDataFrame)
  response = runDataFrame[, n]
  
  if(method == 'ttest'){
    cat('performing permutation analysis of T test selected significant features...\n')
    pVal1 <- apply(runDataFrame[,-n], 2, function(x){
      vec1 <- x[response=='Cancer'|response==1]
      vec2 <- x[!(response=='Cancer'|response==1)]
      t.test(vec1, vec2, alternative = 'two.sided')$p.value
    })
    sigSumOri1 <- sum(pVal1<0.05)
    pperm1 <- foreach(i = 1:thread, .combine = cbind) %dopar% {
      pperm1_1 <- c()
      for(i in 1:ceiling(r/thread)){
        yperm1 <- sample(response)
        pVal1_tmp <- apply(runDataFrame[,-n], 2, function(x){
          vec1 <- x[yperm1=='Cancer'|yperm1==1]
          vec2 <- x[!(yperm1=='Cancer'|yperm1==1)]
          t.test(vec1, vec2, alternative = 'two.sided', var.equal = TRUE)$p.value
        })
        pperm1_1 <- cbind(pperm1_1, pVal1_tmp)
      }
      pperm1_1
    }
    sigSum1 <- apply(pperm1, 2, function(x){
      sum(x<0.05)
    })
    pvalue1 <- sum(sigSum1>sigSumOri1)/length(sigSum1)
    cat(paste('p value = ', pvalue1, '\n'))
    
    pdf(paste(filename, '.pdf', sep = ''), 10, 6)
    hist(sigSum1, nclass = 100, 
         main = 'Distribution of the number of significant aptamers\nafter 1000 permutation (T test)',
         xlab = 'number of significant aptamers (p<0.05)',
         cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
    abline(v = sigSumOri1, col = 'red')
    dev.off()
    z <- c(sigSumOri1, sigSum1)
    write.csv(z, file = paste(filename, '.csv', sep = ''))
    write.csv(pvalue1, file = paste(filename, '_pvalue.csv', sep = ''))
  }
  
  if(method == 'ANOVA'){
    cat('performing permutation analysis of ANOVA F test selected significant features...\n')
    pVal3 <- apply(runDataFrame[,-n], 2, function(x){
      summary(aov(x~response))[[1]][1,5]
    })
    sigSumOri3 <- sum(pVal3<0.05)
    pperm3 <- foreach(i = 1:thread, .combine = cbind) %dopar% {
      pperm3_1 <- c()
      for(i in 1:ceiling(r/thread)){
        yperm3 <- sample(response)
        pVal3_tmp <- apply(runDataFrame[,-n], 2, function(x){
          summary(aov(x~yperm3))[[1]][1,5]
        })
        pperm3_1 <- cbind(pperm3_1, pVal3_tmp)
      }
      pperm3_1
    }
    sigSum3 <- apply(pperm3, 2, function(x){
      sum(x<0.05)
    })
    pvalue2 = sum(sigSum3>sigSumOri3)/length(sigSum3)
    cat(paste('p value = ', pvalue2, '\n'))
    
    pdf(paste(filename, '.pdf', sep = ''), 10, 6)
    hist(sigSum3, nclass = 100, 
         main = 'Distribution of the number of significant aptamers\nafter 1000 permutation (ANOVA F test)',
         xlab = 'number of significant aptamers (p<0.05)',
         cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
    abline(v = sigSumOri3, col = 'red')
    dev.off()
    z <- c(sigSumOri3, sigSum3)
    write.csv(z, file = paste(filename, '.csv', sep = ''))
    write.csv(pvalue2, file = paste(filename, '_pvalue.csv', sep = ''))
  }
  stopCluster(cl)
  z
}






















