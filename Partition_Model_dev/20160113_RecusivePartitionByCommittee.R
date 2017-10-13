options(stringsAsFactors = F)
library(pROC)

#Calculate sensitivity when cases above cutoff are positive
sens <- function(data, spec = 0.95){
  intensity = as.numeric(data[,2])
  cc = data[,1]
  w.cases <- which(cc == 1 & intensity != "NA")
  w.controls <- which(cc == 0 & intensity != "NA")
  x = intensity[w.cases]
  y = intensity[w.controls]
  nx = length(x)
  ny = length(y)
  q = quantile(y, spec, na.rm = T)
  #q = max(y)
  sensitivity = length(intensity[which(cc == 1&intensity > q)])/nx
  specificity = length(intensity[which(cc == 0&intensity < q)])/ny
  z <- list(sensitivity = sensitivity, specificity = specificity, cutoff = q)
}

#Calculate sensitivity when cases below cutoff are positive
sens_low <- function(data, spec = 0.95){
  intensity = as.numeric(data[,2])
  cc = data[,1]
  w.cases <- which(cc == 1 & intensity != "NA")
  w.controls <- which(cc == 0 & intensity != "NA")
  x = intensity[w.cases]
  y = intensity[w.controls]
  nx = length(x)
  ny = length(y)
  q = quantile(y, 1-spec, na.rm = T)
  #q = max(y)
  sensitivity = length(intensity[which(cc == 1&intensity < q)])/nx
  specificity = length(intensity[which(cc == 0&intensity > q)])/ny
  z <- list(sensitivity = sensitivity, specificity = specificity, cutoff = q)
}


#Calculate specificity when controls above cutoff are positive
spec <- function(data, sens = 0.95){
  intensity = as.numeric(data[,2])
  cc = data[,1]
  w.cases <- which(cc == 1 & intensity != "NA")
  w.controls <- which(cc == 0 & intensity != "NA")
  x = intensity[w.cases]
  y = intensity[w.controls]
  nx = length(x)
  ny = length(y)
  q = quantile(x, sens, na.rm = T)
  sensitivity = length(intensity[which(cc == 1&intensity < q)])/nx
  specificity = length(intensity[which(cc == 0&intensity > q)])/ny
  z <- list(sensitivity = sensitivity, specificity = specificity, cutoff = q)
}

#Calculate specificity when controls below cutoff are positive
spec_low <- function(data, sens = 0.95){
  intensity = as.numeric(data[,2])
  cc = data[,1]
  w.cases <- which(cc == 1 & intensity != "NA")
  w.controls <- which(cc == 0 & intensity != "NA")
  x = intensity[w.cases]
  y = intensity[w.controls]
  nx = length(x)
  ny = length(y)
  q = quantile(x, 1-sens, na.rm = T)
  sensitivity = length(intensity[which(cc == 1&intensity > q)])/nx
  specificity = length(intensity[which(cc == 0&intensity < q)])/ny
  z <- list(sensitivity = sensitivity, specificity = specificity, cutoff = q)
}

# Calculate sensitivity and specificity using preset cutoff
ss_pre <- function(data, cutoff = 0.95, testCasePositive = 'aboveCutoff', testLabel = F){
  intensity = as.numeric(data[,2])
  cc = data[,1]
  w.cases <- which(cc == 1 & intensity != "NA")
  w.controls <- which(cc == 0 & intensity != "NA")
  x = intensity[w.cases]
  y = intensity[w.controls]
  nx = length(x)
  ny = length(y)
  if(testCasePositive == 'aboveCutoff'){
    sensitivity = length(intensity[which(cc == 1&intensity > cutoff)])/nx
    specificity = length(intensity[which(cc == 0&intensity < cutoff)])/ny
    z <- list(sensitivity = sensitivity, specificity = specificity)
    if(testLabel){
      testLabel <- intensity
      testLabel[intensity > cutoff] <- 1
      testLabel[intensity <= cutoff] <- 0
      z <- list(sensitivity = sensitivity, specificity = specificity, testLabel = testLabel)
    }
  }
  else if(testCasePositive == 'belowCutoff'){
    sensitivity = length(intensity[which(cc == 1&intensity < cutoff)])/nx
    specificity = length(intensity[which(cc == 0&intensity > cutoff)])/ny
    z <- list(sensitivity = sensitivity, specificity = specificity)
    if(testLabel){
      testLabel <- intensity
      testLabel[intensity < cutoff] <- 1
      testLabel[intensity >= cutoff] <- 0
      z <- list(sensitivity = sensitivity, specificity = specificity, testLabel = testLabel)
    }
  }
  else{
    cat('Please choose \'aboveCutoff\' or \'belowCutoff\' for testPositive\n')
  }
  z
}


####################################################################################################################################
## 
surveyCandidates <- function(dataFrame, CaseAboveCutoff = 1, CaseBelowCutoff = 1, ControlAboveCutoff = 1, ControlBelowCutoff = 1){
  nColumn <- ncol(dataFrame)
  response <- dataFrame[, nColumn]
  feature <- colnames(dataFrame[, -nColumn])
  # three columns: feature, directionType, performance
  performanceDF <- c()
  if(!is.null(CaseAboveCutoff)){
    performance1 <- apply(dataFrame[, -nColumn], 2, function(x){
      sens1 <- sens(data.frame(response, x), spec = CaseAboveCutoff)
      out <- c(sens1$sensitivity, sens1$cutoff)
    })
    performance1 <- data.frame(t(performance1))
    colnames(performance1) <- c('performance', 'cutoff')
    performanceMatrix1 <- data.frame(feature = feature, directionType = 'CaseAboveCutoff', performance = performance1$performance, cutoff = performance1$cutoff, cutoffPercent = CaseAboveCutoff)
    performanceDF <- rbind(performanceDF, performanceMatrix1)
  }
  
  if(!is.null(CaseBelowCutoff)){
    performance2 <- apply(dataFrame[, -nColumn], 2, function(x){
      sens_low1 <- sens_low(data.frame(response, x), spec = CaseBelowCutoff)
      out <- c(sens_low1$sensitivity, sens_low1$cutoff)
    })
    performance2 <- data.frame(t(performance2))
    colnames(performance2) <- c('performance', 'cutoff')
    performanceMatrix2 <- data.frame(feature = feature, directionType = 'CaseBelowCutoff', performance = performance2$performance, cutoff = performance2$cutoff, cutoffPercent = CaseBelowCutoff)
    performanceDF <- rbind(performanceDF, performanceMatrix2)
  }
  
  if(!is.null(ControlAboveCutoff)){
    performance3 <- apply(dataFrame[, -nColumn], 2, function(x){
      spec1 <- spec(data.frame(response, x), sens = ControlAboveCutoff)
      out <- c(spec1$specificity, spec1$cutoff)
    })
    performance3 <- data.frame(t(performance3))
    colnames(performance3) <- c('performance', 'cutoff')
    performanceMatrix3 <- data.frame(feature = feature, directionType = 'ControlAboveCutoff', performance = performance3$performance, cutoff = performance3$cutoff, cutoffPercent = ControlAboveCutoff)
    performanceDF <- rbind(performanceDF, performanceMatrix3)
  }
  
  if(!is.null(ControlBelowCutoff)){
    performance4 <- apply(dataFrame[, -nColumn], 2, function(x){
      spec_low1 <- spec_low(data.frame(response, x), sens = ControlBelowCutoff)
      out <- c(spec_low1$specificity, spec_low1$cutoff)
    })
    performance4 <- data.frame(t(performance4))
    colnames(performance4) <- c('performance', 'cutoff')
    performanceMatrix4 <- data.frame(feature = feature, directionType = 'ControlBelowCutoff', performance = performance4$performance, cutoff = performance4$cutoff, cutoffPercent = ControlBelowCutoff)
    performanceDF <- rbind(performanceDF, performanceMatrix4)
  }
  
  rownames(performanceDF) <- NULL
  z <- list(performanceDF = performanceDF)
  
}

## rule based model
combineRule1 <- function(dataFrame, method = 'case-specific', topCaseFeatures, atLeastCase = 2, topControlFeatures, atLeastControl = 2){

  if(method == 'case-specific' | method == 'both'){
    trueCaseTable <- sweep(dataFrame, 2, topCaseFeatures$cutoff, '>')
    caseCall <- data.frame(sampleID = rownames(dataFrame), 
                           predict = as.numeric(apply(trueCaseTable, 1, sum) >= atLeastCase))
    z <- list(caseCall = caseCall)
  }
  else if(method == 'control-specific' | method == 'both'){
    trueControlTable <- sweep(dataFrame, 2, topControlFeatures$cutoff, '>')
    controlCall <- data.frame(sampleID = rownames(dataFrame), 
                              predict = as.numeric(apply(trueControlTable, 1, sum) < atLeastControl))
    z <- list(controlCall = controlCall)
  }
  
  if(method == 'both'){
    
    z <- list(caseCall = caseCall, controlCall = controlCall)
  }
  z
}


## recursive rules
recursiveRules_v1 <- function(dataFrame, method = 'case-specific', 
                              CaseAboveCutoff = 0.99, atLeastCase = 2, 
                              ControlAboveCutoff = 0.99, atLeastControl = 2,
                              numberOfTopFeatures = 50, maxRound = 100, minSampleRemain = 10){
  # dataFrame = runDF
  nColumn <- ncol(dataFrame)
  response <- dataFrame[,nColumn]
  
  # initialize loop
  dataFrameInLoop <- dataFrame
  responseInLoop <- response
  stopCriteria = F
  nRound <- 1 
  
  # output parameters
  predictTable <- c()
  parameters <- c()
  roundSummary <- c()
  
  while(stopCriteria == F){
    cat(paste('node ', nRound, ', ', sep = ''))
    if(method == 'case-specific'){
      grandSurvey1 <- surveyCandidates(dataFrameInLoop, CaseAboveCutoff = CaseAboveCutoff, CaseBelowCutoff = NULL, ControlAboveCutoff = NULL, ControlBelowCutoff = NULL)
      topCaseFeaturesIDX <- order(grandSurvey1$performanceDF$performance, decreasing = T)[1:numberOfTopFeatures]
      topCaseFeatures <- grandSurvey1$performanceDF[topCaseFeaturesIDX, ]
      topCaseFeatures <- topCaseFeatures[(topCaseFeatures$performance + CaseAboveCutoff) > 1, ]
      numberOfFeaturesPerRound <- nrow(topCaseFeatures)
      cat(paste('number of features: ', numberOfFeaturesPerRound, sep = ''))
      # print(topCaseFeatures)
      
      # stopCriteria checkpoint 1
      if(numberOfFeaturesPerRound < atLeastCase | numberOfFeaturesPerRound == 1){
        cat('not enough features! stop now!\n')
        stopCriteria = T
        lastPredict <- data.frame(sampleID = rownames(dataFrameInLoop), predict = 0, nRound = nRound)
        predictTable <- rbind(predictTable, lastPredict)
      }
      else{
        selectedFeaturesPerRoundIDX <- match(topCaseFeatures$feature, colnames(dataFrameInLoop))
        diagCaseMatrix <- dataFrameInLoop[, selectedFeaturesPerRoundIDX]
        # cat(paste('; ', nrow(diagCaseMatrix), ', ', ncol(diagCaseMatrix), sep = ''))
        predictCase <- combineRule1(diagCaseMatrix, method = 'case-specific', topCaseFeatures = topCaseFeatures, atLeastCase = atLeastCase)
        predictCaseIDX <- which(predictCase$caseCall$predict == 1)

        # stopCriteria checkpoint 1
        nSamplePartitioned <- length(predictCaseIDX)
        if(nSamplePartitioned == 0){
          cat('no sample is above cutoff! stop now!\n')
          stopCriteria = T
          lastPredict <- data.frame(sampleID = rownames(dataFrameInLoop), predict = 0, nRound = nRound)
          predictTable <- rbind(predictTable, lastPredict)
        }
        else{
          # update dataFrameInLoop, responseInLoop
          dataFrameInLoop <- dataFrameInLoop[-predictCaseIDX, -selectedFeaturesPerRoundIDX]
          responseInLoop <- responseInLoop[-predictCaseIDX]
          nSampleRemain = nrow(dataFrameInLoop)
          nFeatureRemain = ncol(dataFrameInLoop) - 1
          nClass = length(unique(responseInLoop))
          
          predictTable <- rbind(predictTable, data.frame(predictCase$caseCall[predictCaseIDX, ], nRound = nRound))
          parameters <- rbind(parameters, data.frame(nRound, topCaseFeatures, atLeastCase))
          roundSummary <- rbind(roundSummary, data.frame(nRound, numberOfFeaturesPerRound, nFeatureRemain, nSamplePartitioned, nSampleRemain))
          
          cat(paste('; remain feature number: ', nFeatureRemain,
                    '; nSamplePartitioned: ', nSamplePartitioned, 
                    '; remain sample size: ', nSampleRemain, '\n', sep = ''))
          
          # stopCriteria checkpoint 2
          if(nRound == maxRound | nSampleRemain < minSampleRemain | nClass == 1){
            cat('finished loop, stop now! \n')
            stopCriteria = T
            if(nSampleRemain != 0){
              lastPredict <- data.frame(sampleID = rownames(dataFrameInLoop), predict = 0, nRound = nRound + 1)
              predictTable <- rbind(predictTable, lastPredict)
            }
          }
          else{
            nRound = nRound + 1
          }
        }
      }
    }
    
    else if(method == 'both'){
      grandSurvey1 <- surveyCandidates(dataFrameInLoop, CaseAboveCutoff = CaseAboveCutoff, CaseBelowCutoff = NULL, ControlAboveCutoff = NULL, ControlBelowCutoff = NULL)
      topCaseFeaturesIDX <- which(grandSurvey1$performanceDF$performance + CaseAboveCutoff > 1.02)
      topCaseFeatures <- grandSurvey1$performanceDF[topCaseFeaturesIDX, ]
      numberOfTopCaseFeatures <- nrow(topCaseFeatures)
      
      grandSurvey2 <- surveyCandidates(dataFrameInLoop, CaseAboveCutoff = NULL, CaseBelowCutoff = NULL, ControlAboveCutoff = ControlAboveCutoff, ControlBelowCutoff = NULL)
      topControlFeaturesIDX <- which(grandSurvey2$performanceDF$performance + ControlAboveCutoff > 1.02)
      topControlFeatures <- grandSurvey2$performanceDF[topControlFeaturesIDX, ]
      numberOfTopControlFeatures <- nrow(topControlFeatures)
      
      if(numberOfTopControlFeatures > numberOfTopCaseFeatures){
        cat('ControlAboveCutoff is selected, ')
      }
      else{
        cat('CaseAboveCutoff is selected, ')
        
      }
      
      
      numberOfFeaturesPerRound <- nrow(topCaseFeatures)
      cat(paste('number of features: ', numberOfFeaturesPerRound, sep = ''))
      
    }
    
  }
  
  z <- list(predictTable = predictTable, parameters = data.frame(parameters), roundSummary = roundSummary)
}

## predict recursive rules
predictRecurRules <- function(recursiveRules, method = 'case-specific', newData){
  
  # initiate loop to predict newData
  nRounds <- length(recursiveRules$roundSummary$nRound)
  newDataInLoop <- newData
  nRound = 1
  stopCriteria = F
  nSampleRemain = nrow(newData)
  
  predictTable <- c()
  if(method == 'case-specific'){
    while(stopCriteria == F){
      cat(paste('node ', nRound, ', ', sep = ''))
      topCaseFeatures <- recursiveRules$parameters[recursiveRules$parameters$nRound == nRound, ]
      numberOfFeaturesPerRound <- nrow(topCaseFeatures)
      selectedFeaturesPerRoundIDX <- match(topCaseFeatures$feature, colnames(newDataInLoop))
      diagCaseMatrix <- newDataInLoop[, selectedFeaturesPerRoundIDX]
      predictCase <- combineRule1(diagCaseMatrix, method = 'case-specific', topCaseFeatures = topCaseFeatures, atLeastCase = unique(topCaseFeatures$atLeastCase))
      predictCaseIDX <- which(predictCase$caseCall$predict == 1)
      nSamplePartitioned <- length(predictCaseIDX)
      if(nSamplePartitioned == 0){
        cat('no sample above cutoff, skip this round! \n')
        nRound = nRound + 1
        
        # stopCriteria Checkpoint
        if(nRound >= nRounds){
          cat('stop now! \n')
          stopCriteria = T
          if(nSampleRemain != 0){
            lastPredict <- data.frame(sampleID = rownames(newDataInLoop), predict = 0)
            predictTable <- rbind(predictTable, lastPredict)
          }
        }
      }
      else{
        predictTable <- rbind(predictTable, predictCase$caseCall[predictCaseIDX, ])
        newDataInLoop <- newDataInLoop[-predictCaseIDX, -selectedFeaturesPerRoundIDX]
        # responseInLoop <- responseInLoop[-predictCaseIDX]
        nSampleRemain = nrow(newDataInLoop)
        
        cat(paste('number of features: ', numberOfFeaturesPerRound, 
                  '; nSamplePartitioned: ', nSamplePartitioned, 
                  '; remain sample size: ', nSampleRemain, '\n', sep = ''))
        
        # stopCriteria Checkpoint
        if(nRound >= nRounds | nSampleRemain == 0){
          cat('stop now! \n')
          stopCriteria = T
          if(nSampleRemain != 0){
            lastPredict <- data.frame(sampleID = rownames(newDataInLoop), predict = 0)
            predictTable <- rbind(predictTable, lastPredict)
          }
        }
        else{
          nRound = nRound + 1
        }
      }
    }
  }
  else if(method == 'control-specific'){
    
  }
  else if(method == 'both'){
    
  }  
  z <- list(predictTable = predictTable)
}


## leave one out cross-validation
recursiveRules_cv <- function(dataFrame, method = 'case-specific',
                              K = 10, R = 1,
                              CaseAboveCutoff = 0.99, atLeastCase = 2, 
                              # ControlAboveCutoff = 0.99, atLeastControl = 2,
                              numberOfTopFeatures = 50, 
                              maxRound = 10, minSampleRemain = 20){
  library(cvTools)
  # library(foreach)
  # library(doParallel)
  # numberOfThreads <- min(K, 40)
  # cl <- makeCluster(numberOfThreads)
  # registerDoParallel(cl)
  
  nSample <- nrow(dataFrame)
  n <- ncol(dataFrame)
  response <- dataFrame[,n]
  response <- as.factor(response)
  dataFrame <- data.frame(dataFrame[,-n], response=response)
  predOut <- list()
  for(j in 1:R){
    cat(paste('\nround', j, 'of', R, '\n'))
    subset <- cvFolds(nSample, K)
    predAll <- c()
    for(i in 1:K){
      cat(paste('\nfold', i, 'of', K, '\n'))
      
      #select features by filtering the training set
      idx_test <- subset$subsets[subset$which==i]
      trainDF <- dataFrame[-idx_test, ]
      testDF <- dataFrame[idx_test, -n]
      
      fitRecRules <- recursiveRules_v1(trainDF, method = method, 
                                       CaseAboveCutoff = CaseAboveCutoff, atLeastCase = atLeastCase, 
                                       # ControlAboveCutoff = ControlAboveCutoff, atLeastControl = atLeastControl,
                                       numberOfTopFeatures = numberOfTopFeatures, 
                                       maxRound = maxRound, minSampleRemain = minSampleRemain)
      cat('predicting...\n')
      pred <- predictRecurRules(fitRecRules, newData = testDF)
      cat(paste('test sample size: ', length(idx_test), '\n', sep =''))
      outTab <- cbind(idx_test, pred$predictTable)
      predAll <- rbind(predAll, outTab)
      
    }
    
    colnames(predAll) <- c('testIDX', 'SampleID', 'ProbCase')
    predAll <- data.frame(predAll)
    predOut[[length(predOut)+1]] <- predAll
  }
  z <- list(predOut = predOut)
}


##
rocRecurRules <- function(dataFrame, method = 'case-specific', performanceType = 'training',
                          K = 10, R = 1,
                          atLeastCase = 2, 
                          # ControlAboveCutoff = ControlAboveCutoff, atLeastControl = atLeastControl,
                          numberOfTopFeatures = 50, 
                          maxRound = 10, minSampleRemain = 20,
                          seqCutoff = seq(0.95, 1, 0.01), plot = T){

  library(pracma)
  n <- ncol(dataFrame)
  response <- dataFrame[,n]
  sensTable <- data.frame(sensitivity = 1, specificity = 0)
  nRound = 1
  
  for(eachCutoff in seqCutoff){
    cat(paste('\nroc round ', nRound, ', cutoff is ', eachCutoff, '\n', sep = ''))
    
    if(performanceType == 'training'){
      predictInLoop <- recursiveRules_v1(dataFrame, method = 'case-specific', 
                                         CaseAboveCutoff = eachCutoff, atLeastCase = atLeastCase, 
                                         # ControlAboveCutoff = 0.99, atLeastControl = 2,
                                         numberOfTopFeatures = numberOfTopFeatures, 
                                         maxRound = maxRound, minSampleRemain = minSampleRemain)
      
      predictDF <- data.frame(response = response[match(predictInLoop$predictTable$sampleID, rownames(dataFrame))],
                              prediction = predictInLoop$predictTable$predict)
    }
    else if(performanceType == 'cross-validation'){
      predictInLoop <- recursiveRules_cv(dataFrame, method = method,
                                         K = K, R = R,
                                         CaseAboveCutoff = eachCutoff, atLeastCase = atLeastCase, 
                                         # ControlAboveCutoff = 0.99, atLeastControl = 2,
                                         numberOfTopFeatures = numberOfTopFeatures, 
                                         maxRound = maxRound, minSampleRemain = minSampleRemain)
      
      predictDF <- data.frame(response = response[match(predictInLoop$predOut[[1]]$SampleID, rownames(dataFrame))],
                              prediction = predictInLoop$predOut[[1]]$ProbCase)
    }
    
    sensInLoop <- ss_pre(predictDF, cutoff = 0.5)
    sensTable <- rbind(sensTable, data.frame(sensitivity = sensInLoop$sensitivity, specificity = sensInLoop$specificity))
    nRound = nRound + 1
  }
  
  plotTable <- data.frame(TPR = sensTable$sensitivity,
                          FPR = (1 - sensTable$specificity))
  plotTable <- plotTable[order(plotTable$FPR), ]
  AUC <- round(abs(trapz(plotTable$FPR, plotTable$TPR)), 3)
  
  if(plot){
    par(mar = c(5,5,5,1))
    plot(plotTable$TPR~plotTable$FPR, type = 'l', xlim = c(0,1), ylim = c(0,1), lwd = 1.5, 
         xlab = 'FPR', ylab = 'TPR', main = AUC, cex.main = 2, cex.lab = 1.5)
    abline(a = 0, b = 1, lwd = 1.5)
  }
  
  z <- list(sensTable = sensTable, plotTable = plotTable, AUC = AUC)
}


# try different parameter combination
# CaseAboveCutoff: seq(0.91, 1, 0.01)
# numberOfTopFeatures: seq(10, 100, 10)
# atLeastCase: c(1, 2, 3, 4)
## tuning parameters
tune_recurRules <- function(dataFrame, method = 'case-specific',
                            performanceType = 'training',
                            K = 10, R = 1,
                            numberOfTopFeatures = c(),
                            CaseAboveCutoff = c(),
                            atLeastCase = c(),
                            maxRound = 1, 
                            minSampleRemain = 20){
  
  library(foreach)
  library(doParallel)
  
  n <- ncol(dataFrame)
  response <- dataFrame[,n]
  
  nTrials_numberOfTopFeatures <- length(numberOfTopFeatures)
  nTrials_CaseAboveCutoff <- length(CaseAboveCutoff)
  nTrials_atLeastCase <- length(atLeastCase)
  nTrials_maxRound <- length(maxRound)

  # numberOfThreads <- max(nTrials_numberOfTopFeatures, nTrials_CaseAboveCutoff, nTrials_atLeastCase, nTrials_maxRound, nTrials_minSampleRemain)
  numberOfThreads <- nTrials_numberOfTopFeatures
  cl <- makeCluster(numberOfThreads)
  registerDoParallel(cl)
  
  performOut <- foreach(i = 1:numberOfThreads, 
                        .export = c('sens', 'sens_low', 'spec', 'spec_low', 'ss_pre', 'combineRule1', 'surveyCandidates', 'recursiveRules_v1', 'predictRecurRules', 'recursiveRules_cv'), 
                        .combine = 'rbind') %dopar% {
    
    each_numberOfTopFeatures <- numberOfTopFeatures[i]
    performanceTable <- data.frame()
    
    for(each_atLeastCase in atLeastCase){
      for(each_CaseAboveCutoff in CaseAboveCutoff){
        for(each_maxRound in maxRound){
          if(performanceType == 'training'){
            predictInLoop <- recursiveRules_v1(dataFrame, method = 'case-specific', 
                                               CaseAboveCutoff = each_CaseAboveCutoff, atLeastCase = each_atLeastCase, 
                                               # ControlAboveCutoff = 0.99, atLeastControl = 2,
                                               numberOfTopFeatures = each_numberOfTopFeatures, maxRound = each_maxRound, minSampleRemain = 10)
            predictDF <- data.frame(response = response[match(predictInLoop$predictTable$sampleID, rownames(dataFrame))],
                                    prediction = predictInLoop$predictTable$predict)
            
          }
          else if(performanceType == 'cross-validation'){
            predictInLoop <- recursiveRules_cv(dataFrame, method = 'case-specific',
                                            K = K, R = R,
                                            CaseAboveCutoff = each_CaseAboveCutoff, atLeastCase = each_atLeastCase, 
                                            # ControlAboveCutoff = 0.99, atLeastControl = 2,
                                            numberOfTopFeatures = each_numberOfTopFeatures, 
                                            maxRound = each_maxRound, minSampleRemain = 20)
            predictDF <- data.frame(response = response[match(predictInLoop$predOut[[1]]$SampleID, rownames(dataFrame))],
                                    prediction = predictInLoop$predOut[[1]]$ProbCase)
          }
          
          sensInLoop <- ss_pre(predictDF, cutoff = 0.5)
          performanceInLoop <- data.frame(numberOfTopFeatures = each_numberOfTopFeatures,
                                          atLeastCase = each_atLeastCase,
                                          CaseAboveCutoff = each_CaseAboveCutoff,
                                          maxRound = each_maxRound, 
                                          sensitivity = sensInLoop$sensitivity, 
                                          specificity = sensInLoop$specificity)
          performanceTable <- rbind(performanceTable, performanceInLoop)
        }
      }
    }
    performanceTable
  }
  
  performOut
  
}



# try different parameter combination
# CaseAboveCutoff: seq(0.91, 1, 0.01)
# numberOfTopFeatures: seq(10, 100, 10)
# atLeastCase: c(1, 2, 3, 4)
## tuning parameters
tune_recurRules_simple <- function(dataFrame, method = 'case-specific',
                                   performanceType = 'training',
                                   K = 10, R = 1,
                                   numberOfTopFeatures = c(),
                                   CaseAboveCutoff = c(),
                                   atLeastCase = c(),
                                   maxRound = 1, 
                                   minSampleRemain = 20){
  
  n <- ncol(dataFrame)
  response <- dataFrame[,n]
  
  nRound = 1
  performanceTable <- data.frame()
  for(each_numberOfTopFeatures in numberOfTopFeatures){
    for(each_atLeastCase in atLeastCase){
      for(each_CaseAboveCutoff in CaseAboveCutoff){
        for(each_maxRound in maxRound){
          cat(paste('\ntune round ', nRound, 
                    ', numberOfTopFeatures: ', each_numberOfTopFeatures, 
                    ', atLeastCase: ', each_atLeastCase,
                    ', CaseAboveCutoff: ', each_CaseAboveCutoff,
                    ', maxRound: ', each_maxRound, '\n', sep = ''))
          
          if(performanceType == 'training'){
            predictInLoop <- recursiveRules_v1(dataFrame, method = 'case-specific', 
                                               CaseAboveCutoff = each_CaseAboveCutoff, atLeastCase = each_atLeastCase, 
                                               # ControlAboveCutoff = 0.99, atLeastControl = 2,
                                               numberOfTopFeatures = each_numberOfTopFeatures, maxRound = each_maxRound, minSampleRemain = 10)
            predictDF <- data.frame(response = response[match(predictInLoop$predictTable$sampleID, rownames(dataFrame))],
                                    prediction = predictInLoop$predictTable$predict)
            
          }
          else if(performanceType == 'cross-validation'){
            predictInLoop <- recursiveRules_cv(dataFrame, method = 'case-specific',
                                               K = K, R = R,
                                               CaseAboveCutoff = each_CaseAboveCutoff, atLeastCase = each_atLeastCase, 
                                               # ControlAboveCutoff = 0.99, atLeastControl = 2,
                                               numberOfTopFeatures = each_numberOfTopFeatures, 
                                               maxRound = each_maxRound, minSampleRemain = 20)
            predictDF <- data.frame(response = response[match(predictInLoop$predOut[[1]]$SampleID, rownames(dataFrame))],
                                    prediction = predictInLoop$predOut[[1]]$ProbCase)
          }
          
          sensInLoop <- ss_pre(predictDF, cutoff = 0.5)
          performanceInLoop <- data.frame(numberOfTopFeatures = each_numberOfTopFeatures,
                                          atLeastCase = each_atLeastCase,
                                          CaseAboveCutoff = each_CaseAboveCutoff,
                                          maxRound = each_maxRound, 
                                          sensitivity = sensInLoop$sensitivity, 
                                          specificity = sensInLoop$specificity)
          performanceTable <- rbind(performanceTable, performanceInLoop)
          nRound = nRound + 1
        }
      }
    }
  }
  
  performOut <- performanceTable
  
}





















