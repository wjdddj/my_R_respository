options(stringsAsFactors = F)
library(pROC)

############################################################################################################################################################################
# simple non-parametric statistics calculation
# data: first column is binary response, 
############################################################################################################################################################################
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
    sensitivity = length(intensity[which(cc == 1 & intensity > cutoff)])/nx
    specificity = length(intensity[which(cc == 0 & intensity < cutoff)])/ny
    z <- list(sensitivity = sensitivity, specificity = specificity)
    if(testLabel){
      testLabel <- intensity
      testLabel[intensity > cutoff] <- 1
      testLabel[intensity <= cutoff] <- 0
      z <- list(sensitivity = sensitivity, specificity = specificity, testLabel = testLabel)
    }
  }
  else if(testCasePositive == 'belowCutoff'){
    sensitivity = length(intensity[which(cc == 1 & intensity < cutoff)])/nx
    specificity = length(intensity[which(cc == 0 & intensity > cutoff)])/ny
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


############################################################################################################################################################################
## compute sensitivity and/or specficity at different directions for all candidate features 
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

############################################################################################################################################################################
# rule based model
# topCaseFeatures: selected features to build the rule-based model, it shares the same structure as the surveyCandidates()$performanceDF
############################################################################################################################################################################
combineRule1 <- function(dataFrame, method = 'case-specific', topCaseFeatures, atLeastCase = 2, topControlFeatures, atLeastControl = 2){
  
  if(method == 'case-specific'){
    trueCaseTable <- sweep(dataFrame, 2, topCaseFeatures$cutoff, '>')
    caseCall <- data.frame(sampleID = rownames(dataFrame), 
                           predict = as.numeric(apply(trueCaseTable, 1, sum) >= atLeastCase))
    z <- caseCall
  }
  else if(method == 'control-specific'){
    trueControlTable <- sweep(dataFrame, 2, topControlFeatures$cutoff, '>')
    controlCall <- data.frame(sampleID = rownames(dataFrame), 
                              predict = as.numeric(apply(trueControlTable, 1, sum) < atLeastControl))
    z <- controlCall
  }
  z
}


############################################################################################################################################################################
# recursive rules model
# Parameters:
#    method = c('both', 'case-specific', 'control-specific'). Type of analysis, the direction of partition.
#    numberOfTopFeatures: positive integer or 'max'. 
#                         positive integer: The number of top ranked features to feed into the rule-based model.
#                         ‘max’: no limits, all features that match the sensitivity + specificity criteria will be fed into the rule-based model.
#    specPLUSsensCutoff: numeric. For each feature, sensitivity + specificity needs to meet the cutoff to be included in the rule-based model.
#    CaseAboveCutoff / ControlAboveCutoff: numeric (0,1]. The cutoff of specificity or sensitivity.
#    atLeastCase / atLeaseControl: positive integer. The number of features a sample need to pass in order to be considered positive in the rule-based model.
#    maxRound: positive integer. The maxmum number of rounds of partition allowed.
#    minSampleRemain: integer [0, infinite). The minimum number of remaining samples allowed after recursive partitions.
############################################################################################################################################################################
recursiveRules_v1 <- function(dataFrame, 
                              method = 'case-specific', 
                              numberOfTopFeatures = 50, specPLUSsensCutoff = 1.02, 
                              CaseAboveCutoff = 0.99, atLeastCase = 2, 
                              ControlAboveCutoff = 0.99, atLeastControl = 2,
                              maxRound = 100, 
                              minSampleRemain = 10){
  # dataFrame = runDF
  nColumn <- ncol(dataFrame)
  response <- dataFrame[,nColumn]
  casePrevalence <- sum(response==1)/length(response)
  
  # initialize loop
  dataFrameInLoop <- dataFrame
  responseInLoop <- response
  stopCriteria = F
  nRound <- 1 
  nSampleRemain <- nrow(dataFrameInLoop)
  
  # output parameters
  predictTable <- c()
  parameters <- c()
  roundSummary <- c()
  
  while(stopCriteria == F){
    cat(paste('node ', nRound, ', ', sep = ''))
    # compute sensitivity or specificity according to method input argument, decide modeOfAnlaysis at current round
    if(method == 'case-specific' | method == 'both'){
      grandSurvey1 <- surveyCandidates(dataFrameInLoop, CaseAboveCutoff = CaseAboveCutoff, CaseBelowCutoff = NULL, ControlAboveCutoff = NULL, ControlBelowCutoff = NULL)
      if(numberOfTopFeatures == 'max'){
        topCaseFeaturesIDX <- which(grandSurvey1$performanceDF$performance + CaseAboveCutoff > specPLUSsensCutoff)
        topCaseFeatures <- grandSurvey1$performanceDF[topCaseFeaturesIDX, ]
      }
      else{
        numberOfTopFeatures <- as.numeric(numberOfTopFeatures)
        topCaseFeaturesIDX <- order(grandSurvey1$performanceDF$performance, decreasing = T)[1:numberOfTopFeatures]
        topCaseFeatures <- grandSurvey1$performanceDF[topCaseFeaturesIDX, ]
        topCaseFeatures <- topCaseFeatures[(topCaseFeatures$performance + CaseAboveCutoff) > specPLUSsensCutoff, ]
      }
      numberOfTopCaseFeatures <- nrow(topCaseFeatures)
      modeOfAnalysis <- 'CaseAboveCutoff'
    }
    
    if(method == 'control-specific' | method == 'both'){
      grandSurvey2 <- surveyCandidates(dataFrameInLoop, CaseAboveCutoff = NULL, CaseBelowCutoff = NULL, ControlAboveCutoff = ControlAboveCutoff, ControlBelowCutoff = NULL)
      if(numberOfTopFeatures == 'max'){
        topControlFeaturesIDX <- which(grandSurvey2$performanceDF$performance + ControlAboveCutoff > specPLUSsensCutoff)
        topControlFeatures <- grandSurvey2$performanceDF[topControlFeaturesIDX, ]
      }
      else{
        topControlFeaturesIDX <- order(grandSurvey2$performanceDF$performance, decreasing = T)[1:numberOfTopFeatures]
        topControlFeatures <- grandSurvey2$performanceDF[topControlFeaturesIDX, ]
        topControlFeatures <- topControlFeatures[(topControlFeatures$performance + CaseAboveCutoff) > specPLUSsensCutoff, ]
      }
      numberOfTopControlFeatures <- nrow(topControlFeatures)
      modeOfAnalysis <- 'ControlAboveCutoff'
    }
    
    if(method == 'both'){
      if(numberOfTopCaseFeatures > numberOfTopControlFeatures | method == 'case-specific'){
        cat('CaseAboveCutoff is selected; ')
        modeOfAnalysis <- 'CaseAboveCutoff'
      }
      else{
        cat('ControlAboveCutoff is selected; ')
        modeOfAnalysis <- 'ControlAboveCutoff'
      }
    }
    
    # start prediction mode
    if(modeOfAnalysis == 'CaseAboveCutoff'){
      topFeatures <- topCaseFeatures
      numberOfFeaturesPerRound <- numberOfTopCaseFeatures
      atLeastNumber <- atLeastCase
    }
    else if(modeOfAnalysis == 'ControlAboveCutoff'){
      topFeatures <- topControlFeatures
      numberOfFeaturesPerRound <- numberOfTopControlFeatures
      atLeastNumber <- atLeastControl
    }
    cat(paste('number of features: ', numberOfFeaturesPerRound, 
              '; sensitivity + specificity > ', specPLUSsensCutoff, sep = ''))
    
    # stopCriteria checkpoint 1
    if(numberOfFeaturesPerRound < atLeastNumber | numberOfFeaturesPerRound == 1){
      cat(', not enough features! stop now!\n')
      stopCriteria = T
      lastCaseProb <- sum(responseInLoop == 1)/nSampleRemain
      if(lastCaseProb > casePrevalence){
        lastProb <- 1
      }
      else{
        lastProb <- 0
      }
      lastPredict <- data.frame(sampleID = rownames(dataFrameInLoop), predict = lastProb, nRound = nRound, modeOfAnalysis = 'LastRound')
      predictTable <- rbind(predictTable, lastPredict)
    }
    else{
      if(modeOfAnalysis == 'CaseAboveCutoff'){
        selectedFeaturesPerRoundIDX <- match(topCaseFeatures$feature, colnames(dataFrameInLoop))
        diagMatrix <- dataFrameInLoop[, selectedFeaturesPerRoundIDX]
        # cat(paste('; ', nrow(diagMatrix), ', ', ncol(diagMatrix), sep = ''))
        predictDF <- combineRule1(diagMatrix, method = 'case-specific', topCaseFeatures = topCaseFeatures, atLeastCase = atLeastCase)
        predictIDX <- which(predictDF$predict == 1)
      }
      else if(modeOfAnalysis == 'ControlAboveCutoff'){
        selectedFeaturesPerRoundIDX <- match(topControlFeatures$feature, colnames(dataFrameInLoop))
        diagMatrix <- dataFrameInLoop[, selectedFeaturesPerRoundIDX]
        # cat(paste('; ', nrow(diagMatrix), ', ', ncol(diagMatrix), sep = ''))
        predictDF <- combineRule1(diagMatrix, method = 'control-specific', topControlFeatures = topControlFeatures, atLeastControl = atLeastControl)
        predictIDX <- which(predictDF$predict == 0)
      }
      
      # stopCriteria checkpoint 2
      nSamplePartitioned <- length(predictIDX)
      if(nSamplePartitioned == 0){
        cat(', no sample is above cutoff! stop now!\n')
        stopCriteria = T
        lastCaseProb <- sum(responseInLoop == 1)/nSampleRemain
        if(lastCaseProb > casePrevalence){
          lastProb <- 1
        }
        else{
          lastProb <- 0
        }
        lastPredict <- data.frame(sampleID = rownames(dataFrameInLoop), predict = lastProb, nRound = nRound, modeOfAnalysis = 'LastRound')
        predictTable <- rbind(predictTable, lastPredict)
      }
      else{
        # update dataFrameInLoop, responseInLoop
        # dataFrameInLoop <- dataFrameInLoop[-predictCaseIDX, -selectedFeaturesPerRoundIDX]
        dataFrameInLoop <- dataFrameInLoop[-predictIDX, ]
        responseInLoop <- responseInLoop[-predictIDX]
        nSampleRemain = nrow(dataFrameInLoop)
        # nFeatureRemain = ncol(dataFrameInLoop) - 1
        nClass = length(unique(responseInLoop))
        
        if(modeOfAnalysis == 'CaseAboveCutoff'){
          predictTable <- rbind(predictTable, data.frame(predictDF[predictIDX, ], nRound = nRound, modeOfAnalysis = modeOfAnalysis))
        }
        else if(modeOfAnalysis == 'ControlAboveCutoff'){
          predictTable <- rbind(predictTable, data.frame(predictDF[predictIDX, ], nRound = nRound, modeOfAnalysis = modeOfAnalysis))
        }
        parameters <- rbind(parameters, data.frame(nRound, modeOfAnalysis, topFeatures, atLeastNumber))
        roundSummary <- rbind(roundSummary, data.frame(nRound, modeOfAnalysis, numberOfFeaturesPerRound, specPLUSsensCutoff, nSamplePartitioned, nSampleRemain))
        
        cat(paste('; nSamplePartitioned: ', nSamplePartitioned, 
                  '; remain sample size: ', nSampleRemain, '\n', sep = ''))
        
        # stopCriteria checkpoint 3
        if(nRound == maxRound | nSampleRemain <= minSampleRemain | nClass == 1){
          if(nRound == maxRound){
            cat('maxRound reached, stop now!\n')
          }
          if(nSampleRemain < minSampleRemain){
            cat('no more samples remain, stop now!\n')
          }
          if(nClass == 1){
            cat('all remaining samples are pure, stop now!\n')
          }
          cat('loop finished!\n')
          stopCriteria = T
          if(nSampleRemain != 0){
            lastCaseProb <- sum(responseInLoop == 1)/nSampleRemain
            if(lastCaseProb > casePrevalence){
              lastProb <- 1
            }
            else{
              lastProb <- 0
            }
            lastPredict <- data.frame(sampleID = rownames(dataFrameInLoop), predict = lastProb, nRound = nRound + 1, modeOfAnalysis = 'LastRound')
            predictTable <- rbind(predictTable, lastPredict)
          }
        }
        else{
          nRound = nRound + 1
        }
      }
    }
  }
  if(is.null(roundSummary)){
    modelType = 'failed'
  }
  else{
    modelType = 'successful'
  }
  
  z <- list(predictTable = predictTable, parameters = data.frame(parameters), roundSummary = roundSummary, modelType = modelType, method = method)
}

############################################################################################################################################################################
# predict recursive rules
############################################################################################################################################################################
predictRecurRules <- function(recursiveRules, newData){
  
  method = recursiveRules$method
  
  # initiate loop to predict newData
  nRounds <- length(recursiveRules$roundSummary$nRound)
  newDataInLoop <- newData
  nRound = 1
  stopCriteria = F
  nSampleRemain = nrow(newData)
  predictTable <- c()
  
  # check the training modelType
  if(recursiveRules$modelType == 'failed'){
    cat('No model detected.\nWill not perform prediction.\nAll samples are assigned as unpredictable.\n')
    if(method == 'case-specific'){
      lastPredict <- data.frame(sampleID = rownames(newDataInLoop), predict = 0)
      predictTable <- rbind(predictTable, lastPredict)
    }
    else if(method == 'control-specific'){
      lastPredict <- data.frame(sampleID = rownames(newDataInLoop), predict = 1)
      predictTable <- rbind(predictTable, lastPredict)
    }
    else if(method == 'both'){
      lastPredict <- data.frame(sampleID = rownames(newDataInLoop), predict = 'unpredictable')
      predictTable <- rbind(predictTable, lastPredict)
    }
  }
  else{
    while(stopCriteria == F){
      cat(paste('node ', nRound, ', ', sep = ''))
      if(recursiveRules$roundSummary$modeOfAnalysis[nRound] == 'CaseAboveCutoff'){
        cat('CaseAboveCutoff is selected; ')
        topCaseFeatures <- recursiveRules$parameters[recursiveRules$parameters$nRound == nRound, ]
        numberOfFeaturesPerRound <- nrow(topCaseFeatures)
        selectedFeaturesPerRoundIDX <- match(topCaseFeatures$feature, colnames(newDataInLoop))
        diagMatrix <- newDataInLoop[, selectedFeaturesPerRoundIDX]
        predictDF <- combineRule1(diagMatrix, method = 'case-specific', topCaseFeatures = topCaseFeatures, atLeastCase = unique(topCaseFeatures$atLeastNumber))
        predictIDX <- which(predictDF$predict == 1)
      }
      
      if(recursiveRules$roundSummary$modeOfAnalysis[nRound] == 'ControlAboveCutoff'){
        cat('ControlAboveCutoff is selected; ')
        topControlFeatures <- recursiveRules$parameters[recursiveRules$parameters$nRound == nRound, ]
        numberOfFeaturesPerRound <- nrow(topControlFeatures)
        selectedFeaturesPerRoundIDX <- match(topControlFeatures$feature, colnames(newDataInLoop))
        diagMatrix <- newDataInLoop[, selectedFeaturesPerRoundIDX]
        predictDF <- combineRule1(diagMatrix, method = 'control-specific', topControlFeatures = topControlFeatures, atLeastControl = unique(topControlFeatures$atLeastNumber))
        predictIDX <- which(predictDF$predict == 0)
      }
      cat(paste('number of features: ', numberOfFeaturesPerRound, sep = ''))
      
      nSamplePartitioned <- length(predictIDX)
      if(nSamplePartitioned == 0){
        cat('; no sample above cutoff, skip this round! \n')
        nRound = nRound + 1
        
        # stopCriteria Checkpoint
        if(nRound > nRounds){
          cat('\nmax round reached, stop now!\nloop finished!\n')
          stopCriteria = T
          if(nSampleRemain != 0){
            if(method == 'case-specific'){
              lastPredict <- data.frame(sampleID = rownames(newDataInLoop), predict = 0)
              predictTable <- rbind(predictTable, lastPredict)
            }
            else if(method == 'control-specific'){
              lastPredict <- data.frame(sampleID = rownames(newDataInLoop), predict = 1)
              predictTable <- rbind(predictTable, lastPredict)
            }
            else if(method == 'both'){
              lastPredict <- data.frame(sampleID = rownames(newDataInLoop), predict = 'unpredictable')
              predictTable <- rbind(predictTable, lastPredict)
            }
          }
        }
      }
      else{
        predictTable <- rbind(predictTable, predictDF[predictIDX, ])
        newDataInLoop <- newDataInLoop[-predictIDX, ]
        # responseInLoop <- responseInLoop[-predictIDX]
        nSampleRemain = nrow(newDataInLoop)
        
        cat(paste('; nSamplePartitioned: ', nSamplePartitioned, 
                  '; remain sample size: ', nSampleRemain, '\n', sep = ''))
        
        # stopCriteria Checkpoint
        if(nRound >= nRounds | nSampleRemain == 0){
          if(nRound >= nRounds){
            cat('max round reached, stop now!\n')
          }
          if(nSampleRemain == 0){
            cat('no more samples remain, stop now!\n')
          }
          cat('loop finished!\n')
          stopCriteria = T
          if(nSampleRemain != 0){
            if(method == 'case-specific'){
              lastPredict <- data.frame(sampleID = rownames(newDataInLoop), predict = 0)
              predictTable <- rbind(predictTable, lastPredict)
            }
            else if(method == 'control-specific'){
              lastPredict <- data.frame(sampleID = rownames(newDataInLoop), predict = 1)
              predictTable <- rbind(predictTable, lastPredict)
            }
            else if(method == 'both'){
              lastPredict <- data.frame(sampleID = rownames(newDataInLoop), predict = 'unpredictable')
              predictTable <- rbind(predictTable, lastPredict)
            }
          }
        }
        else{
          nRound = nRound + 1
        } 
      }
    }
  }
  z <- list(predictTable = predictTable)
}


############################################################################################################################################################################
# cross-validation
# no parallel computing
############################################################################################################################################################################
recursiveRules_cv <- function(dataFrame, 
                              K = 10, R = 1,
                              method = 'case-specific',
                              numberOfTopFeatures = 50, specPLUSsensCutoff = 1.02, 
                              CaseAboveCutoff = 0.99, atLeastCase = 2, 
                              ControlAboveCutoff = 0.99, atLeastControl = 2,
                              maxRound = 100, minSampleRemain = 10){
  library(cvTools)
  
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
                                       numberOfTopFeatures = numberOfTopFeatures, specPLUSsensCutoff = specPLUSsensCutoff, 
                                       CaseAboveCutoff = CaseAboveCutoff, atLeastCase = atLeastCase, 
                                       ControlAboveCutoff = ControlAboveCutoff, atLeastControl = atLeastControl,
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

############################################################################################################################################################################
# cross-validation in parallel
# parallel computing implemented by folds
############################################################################################################################################################################
recursiveRules_cv_para <- function(dataFrame, verbose = F, maxProcesses = 40, 
                                   K = 10, R = 1, 
                                   method = 'case-specific', 
                                   numberOfTopFeatures = 50, specPLUSsensCutoff = 1.02, 
                                   CaseAboveCutoff = 0.99, atLeastCase = 2, 
                                   ControlAboveCutoff = 0.99, atLeastControl = 2,
                                   maxRound = 100, minSampleRemain = 10){
  library(cvTools)
  library(foreach)
  library(doParallel)
  numberOfThreads <- min(K, maxProcesses)
  
  nSample <- nrow(dataFrame)
  n <- ncol(dataFrame)
  response <- dataFrame[,n]
  response <- as.factor(response)
  dataFrame <- data.frame(dataFrame[,-n], response=response)
  predOut <- list()
  
  for(j in 1:R){
    
    # whether output to console
    if(verbose){
      cl <- makeCluster(numberOfThreads, outfile = '')
    }
    else{
      cl <- makeCluster(numberOfThreads)
    }
    registerDoParallel(cl)
    
    cat(paste('\nround', j, 'of', R, '\n'))
    
    # identify which samples belong to which fold and thread
    subset <- cvFolds(nSample, K)
    if(K > maxProcesses){
      foldsPerThread <- cvFolds(K, maxProcesses)
      idxAll <- c()
      for(i in 1:numberOfThreads){
        foldInfo <- foldsPerThread$subsets[foldsPerThread$which == i]
        for(each_fold in foldInfo){
          idxByFold <- data.frame(IDX = subset$subsets[subset$which == each_fold], 
                                  fold = each_fold,
                                  thread = i)
          idxAll <- rbind(idxAll, idxByFold)
        }
      }
    }
    else{
      idxAll <- c()
      for(i in 1:numberOfThreads){
        idxByFold <- data.frame(IDX = subset$subsets[subset$which == i], 
                                fold = i,
                                thread = i)
        idxAll <- rbind(idxAll, idxByFold)
      }
    }
    

    # start parallel computing based on fold
    predAll <- foreach(i = 1:numberOfThreads, 
                       .export = c('sens', 'sens_low', 'spec', 'spec_low', 'ss_pre', 'combineRule1', 'surveyCandidates', 'recursiveRules_v1', 'predictRecurRules'), 
                       .combine = 'rbind') %dopar% {
                         
                         # cat(paste('thread: ', i, '\n', sep = ''))
                         # foldInfo <- foldsPerThread$subsets[foldsPerThread$which == i]
                         foldInfo <- idxAll[idxAll$thread == i,]
                         uniqueFolds <- unique(foldInfo$fold)
                         nFolds <- length(uniqueFolds)
                         predictByThread <- c()
                         for(each_foldIDX in 1:nFolds){
                           # idx_test <- subset$subsets[subset$which == each_fold]
                           idx_test <- foldInfo$IDX[foldInfo$fold == uniqueFolds[each_foldIDX]]
                           trainDF <- dataFrame[-idx_test, ]
                           testDF <- dataFrame[idx_test, -n]
                           
                           fitRecRules <- recursiveRules_v1(trainDF, method = method, 
                                                            CaseAboveCutoff = CaseAboveCutoff, atLeastCase = atLeastCase, 
                                                            ControlAboveCutoff = ControlAboveCutoff, atLeastControl = atLeastControl,
                                                            numberOfTopFeatures = numberOfTopFeatures, 
                                                            maxRound = maxRound, minSampleRemain = minSampleRemain)
                           cat('predicting...\n')
                           pred <- predictRecurRules(fitRecRules, newData = testDF)
                           cat(paste('test sample size: ', length(idx_test), '\n', sep =''))
                           outTab <- cbind(idx_test, pred$predictTable)
                           predictByThread <- rbind(predictByThread, outTab)
                         }
                         predictByThread
                       }
    
    colnames(predAll) <- c('testIDX', 'SampleID', 'ProbCase')
    predAll <- data.frame(predAll)
    predOut[[length(predOut)+1]] <- predAll
    stopCluster(cl)
  }
  z <- list(predOut = predOut)
}


############################################################################################################################################################################
# need to re-write
############################################################################################################################################################################
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


############################################################################################################################################################################
# tuning parameters
# parallel computed implemented by AboveCutoff
############################################################################################################################################################################
tune_recurRules <- function(dataFrame, 
                            performanceType = 'training',
                            K = 10, R = 1,
                            methods = c(),
                            numberOfTopFeatures = c(),
                            specPLUSsensCutoff = c(), 
                            AboveCutoff = c(),
                            atLeastNumber = c(),
                            maxRound = 1, 
                            minSampleRemain = 10){
  
  library(foreach)
  library(doParallel)
  
  n <- ncol(dataFrame)
  response <- dataFrame[,n]
  
  nTrials_numberOfTopFeatures <- length(numberOfTopFeatures)
  nTrials_specPLUSsensCutoff <- length(specPLUSsensCutoff)
  nTrials_AboveCutoff <- length(AboveCutoff)
  nTrials_atLeastNumber <- length(atLeastNumber)
  nTrials_maxRound <- length(maxRound)
  nTrials_methods <- length(methods)
  
  # numberOfThreads <- max(nTrials_numberOfTopFeatures, nTrials_CaseAboveCutoff, nTrials_atLeastCase, nTrials_maxRound, nTrials_minSampleRemain)
  numberOfThreads <- nTrials_AboveCutoff
  cl <- makeCluster(numberOfThreads, outfile = '')
  registerDoParallel(cl)
  
  nRound = 1
  performOut <- foreach(i = 1:numberOfThreads, 
                        .export = c('sens', 'sens_low', 'spec', 'spec_low', 'ss_pre', 'combineRule1', 'surveyCandidates', 'recursiveRules_v1', 'predictRecurRules', 'recursiveRules_cv'), 
                        .combine = 'rbind') %dopar% {
                          
    # output console to a file for each thread                    
    filename <- paste(Sys.time(), '_AboveCutoff_', AboveCutoff[i], '.txt', sep = '')
    sink(filename, append = T, type = 'output')
    
    each_AboveCutoff <- AboveCutoff[i]
    performanceTable <- data.frame()
    
    for(each_method in methods){
      for(each_numberOfTopFeatures in numberOfTopFeatures){
        for(each_specPLUSsensCutoff in specPLUSsensCutoff){
          for(each_atLeastNumber in atLeastNumber){
            for(each_maxRound in maxRound){
              
              cat(paste('\ntune round ', nRound, 
                        ', method: ', each_method,
                        ', numberOfTopFeatures: ', each_numberOfTopFeatures, 
                        ', specPLUSsensCutoff: ', each_specPLUSsensCutoff, 
                        ', atLeastNumber: ', each_atLeastNumber,
                        ', AboveCutoff: ', each_AboveCutoff,
                        ', maxRound: ', each_maxRound, '\n', sep = ''))
              
              if(performanceType == 'training'){
                predictInLoop <- recursiveRules_v1(dataFrame, method = each_method, 
                                                   CaseAboveCutoff = each_AboveCutoff, atLeastCase = each_atLeastNumber, 
                                                   ControlAboveCutoff = each_AboveCutoff, atLeastControl = each_atLeastNumber,
                                                   numberOfTopFeatures = each_numberOfTopFeatures, specPLUSsensCutoff = each_specPLUSsensCutoff, 
                                                   maxRound = each_maxRound, minSampleRemain = 10)
                predictDF <- data.frame(response = response[match(predictInLoop$predictTable$sampleID, rownames(dataFrame))],
                                        prediction = predictInLoop$predictTable$predict)
                numberOfPredicatable <- 500
              }
              else if(performanceType == 'cross-validation'){
                predictInLoop <- recursiveRules_cv(dataFrame, method = each_method,
                                                   K = K, R = R,
                                                   CaseAboveCutoff = each_AboveCutoff, atLeastCase = each_atLeastNumber, 
                                                   ControlAboveCutoff = each_AboveCutoff, atLeastControl = each_atLeastNumber,
                                                   numberOfTopFeatures = each_numberOfTopFeatures, specPLUSsensCutoff = each_specPLUSsensCutoff, 
                                                   maxRound = each_maxRound, minSampleRemain = 10)
                predictableIDX <- which(predictInLoop$predOut[[1]]$ProbCase != 'unpredictable')
                predictDF <- data.frame(response = response[predictInLoop$predOut[[1]]$testIDX[predictableIDX]],
                                        prediction = predictInLoop$predOut[[1]]$ProbCase[predictableIDX])
                numberOfPredicatable <- length(predictableIDX)
              }
              
              if(numberOfPredicatable == 0){
                sensInLoop <- list(sensitivity = NA, specificity = NA)
              }
              else{
                sensInLoop <- ss_pre(predictDF, cutoff = 0.5)
              }
              performanceInLoop <- data.frame(methods = each_method,
                                              numberOfTopFeatures = each_numberOfTopFeatures,
                                              specPLUSsensCutoff = each_specPLUSsensCutoff, 
                                              atLeastNumber = each_atLeastNumber,
                                              AboveCutoff = each_AboveCutoff,
                                              maxRound = each_maxRound, 
                                              numberOfPredicatable = numberOfPredicatable,
                                              sensitivity = sensInLoop$sensitivity, 
                                              specificity = sensInLoop$specificity)
              performanceTable <- rbind(performanceTable, performanceInLoop)
              
              nRound = nRound + 1 
            }
          }
        }
      }
    }
    performanceTable
  }
  stopCluster(cl)
  
  performOut
  
}



############################################################################################################################################################################
# tuning parameters
# simple version without parallel computing
############################################################################################################################################################################
tune_recurRules_simple <- function(dataFrame, filename,
                                   maxProcesses, 
                                   performanceType = 'training',
                                   K = 10, R = 1,
                                   methods = c(),
                                   numberOfTopFeatures = c(),
                                   specPLUSsensCutoff = c(), 
                                   AboveCutoff = c(),
                                   atLeastNumber = c(),
                                   maxRound = c(), 
                                   minSampleRemain = 10){
  
  n <- ncol(dataFrame)
  response <- dataFrame[,n]
  
  nRound = 1
  performanceTable <- data.frame()
  fileName <- paste(Sys.Date(), '_', filename, '.txt', sep = '')
  write.csv(performanceTable, file = fileName)
  for(each_method in methods){
    for(each_AboveCutoff in AboveCutoff){
      for(each_maxRound in maxRound){
        for(each_atLeastNumber in atLeastNumber){
          for(each_specPLUSsensCutoff in specPLUSsensCutoff){
            for(each_numberOfTopFeatures in numberOfTopFeatures){
              cat(paste('\ntune round ', nRound, 
                        ', method: ', each_method,
                        ', numberOfTopFeatures: ', each_numberOfTopFeatures, 
                        ', specPLUSsensCutoff: ', each_specPLUSsensCutoff, 
                        ', atLeastNumber: ', each_atLeastNumber,
                        ', AboveCutoff: ', each_AboveCutoff,
                        ', maxRound: ', each_maxRound, '\n', sep = ''))
              
              if(performanceType == 'training'){
                predictInLoop <- recursiveRules_v1(dataFrame, method = each_method, 
                                                   CaseAboveCutoff = each_AboveCutoff, atLeastCase = each_atLeastNumber, 
                                                   ControlAboveCutoff = each_AboveCutoff, atLeastControl = each_atLeastNumber,
                                                   numberOfTopFeatures = each_numberOfTopFeatures, 
                                                   maxRound = each_maxRound, minSampleRemain = minSampleRemain)
                predictDF <- data.frame(response = response[match(predictInLoop$predictTable$sampleID, rownames(dataFrame))],
                                        prediction = predictInLoop$predictTable$predict)
                numberOfPredicatable <- 0
              }
              else if(performanceType == 'cross-validation'){
                predictInLoop <- recursiveRules_cv_para(dataFrame, method = each_method, maxProcesses = maxProcesses, verbose = F,
                                                        K = K, R = R,
                                                        CaseAboveCutoff = each_AboveCutoff, atLeastCase = each_atLeastNumber, 
                                                        ControlAboveCutoff = each_AboveCutoff, atLeastControl = each_atLeastNumber,
                                                        numberOfTopFeatures = each_numberOfTopFeatures, specPLUSsensCutoff = each_specPLUSsensCutoff, 
                                                        maxRound = each_maxRound, minSampleRemain = minSampleRemain)
                predictableIDX <- which(predictInLoop$predOut[[1]]$ProbCase != 'unpredictable')
                predictDF <- data.frame(response = response[predictInLoop$predOut[[1]]$testIDX[predictableIDX]],
                                        prediction = predictInLoop$predOut[[1]]$ProbCase[predictableIDX])
                numberOfPredicatable <- length(predictableIDX)
              }
              
              if(numberOfPredicatable == 0){
                sensInLoop <- list(sensitivity = NA, specificity = NA)
              }
              else{
                sensInLoop <- ss_pre(predictDF, cutoff = 0.5)
              }
              
              performanceInLoop <- data.frame(method = each_method,
                                              numberOfTopFeatures = each_numberOfTopFeatures,
                                              specPLUSsensCutoff = each_specPLUSsensCutoff, 
                                              atLeastNumber = each_atLeastNumber,
                                              AboveCutoff = each_AboveCutoff,
                                              maxRound = each_maxRound, 
                                              numberOfPredicatable = numberOfPredicatable, 
                                              sensitivity = sensInLoop$sensitivity, 
                                              specificity = sensInLoop$specificity)
              
              
              output <- paste('\ntune round ', nRound, 
                              ', method: ', each_method,
                              ', numberOfTopFeatures: ', each_numberOfTopFeatures, 
                              ', specPLUSsensCutoff: ', each_specPLUSsensCutoff, 
                              ', atLeastNumber: ', each_atLeastNumber,
                              ', AboveCutoff: ', each_AboveCutoff,
                              ', maxRound: ', each_maxRound,
                              ', numberOfPredicatable: ', numberOfPredicatable, 
                              ', sensitivity: ', sensInLoop$sensitivity, 
                              ', specificity: ', sensInLoop$specificity,'\n', sep = '')
              write(output, file = fileName, append = TRUE)
              
              performanceTable <- rbind(performanceTable, performanceInLoop)
              nRound = nRound + 1
            }
          }
        }
      }
    }
  }
  
  performOut <- performanceTable
  
}