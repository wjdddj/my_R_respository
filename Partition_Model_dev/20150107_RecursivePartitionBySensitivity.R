options(stringsAsFactors = F)
library(pROC)

#Calculate sensitivity when cases above cutoff are positive
sens <- function(data, spec = 0.95){
  intensity = as.numeric(data[,2])
  cc = data[,1]
  w.cases = which(cc == 1)
  w.controls = which(cc == 0)
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
  w.cases = which(cc == 1)
  w.controls = which(cc == 0)
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
  w.cases = which(cc == 1&intensity != "NA")
  w.controls = which(cc == 0&intensity != "NA")
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
  w.cases = which(cc == 1&intensity != "NA")
  w.controls = which(cc == 0&intensity != "NA")
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
  w.cases = which(cc == 1&intensity != "NA")
  w.controls = which(cc == 0&intensity != "NA")
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

##
rPartition_v1 <- function(dataFrame, 
                          # minimum number of samples partitioned during each round
                          minSamplePerRoundCoef = 0.1, 
                          # maximum number of rounds
                          maxRound = 30, 
                          # minimum number of samples remained
                          minRemain = 20, 
                          # performance surveyCandidates() parameters
                          CaseAboveCutoff = 1, 
                          CaseBelowCutoff = 1, 
                          ControlAboveCutoff = 1, 
                          ControlBelowCutoff = 1){
  
  n <- ncol(dataFrame)
  response <- dataFrame[,n]
  casePrevalence <- sum(response==1)/length(response)
  
  # initiate while loop
  dataFrameInLoop <- dataFrame
  responseInLoop <- dataFrame[,n]
  nRound <- 1
  stopCriteria = F
  
  # output parameters
  predictTable <- c()
  # selectedNodes <- list()
  parameters <- c()
  
  while(stopCriteria == F){
    cat(paste('node ', nRound, ', ', sep = ''))
    grandSurvey <- surveyCandidates(dataFrameInLoop, 
                                    CaseAboveCutoff = CaseAboveCutoff, CaseBelowCutoff = CaseBelowCutoff, 
                                    ControlAboveCutoff = ControlAboveCutoff, ControlBelowCutoff = ControlBelowCutoff)
    
    # stopCriteria check point 1
    nAvailableCandidates <- sum(apply(grandSurvey$performanceDF, 1, function(x){
      as.numeric(x[3]) + as.numeric(x[5]) > 1
    }))
    if(nAvailableCandidates == 0){
      cat('no more candidates! \n')
      cat('stop now!\n')
      stopCritera = T
      lastCasePrevalence <- sum(responseInLoop == 1)/length(responseInLoop)
      if(lastCasePrevalence > casePrevalence){
        casePrediction = 1
      }
      else{casePrediction = 0}
      lastPredict <- data.frame(sampleID = rownames(dataFrameInLoop), 
                                casePrediction = casePrediction, controlPrediction = 1 - casePrediction, nRound = nRound)
      predictTable <- rbind(predictTable, lastPredict)
      colnames(predictTable) <- c('sampleID', 'casePrediction', 'controlPrediction', 'nRound')
    }
    else{
      # selected the best node
      selectedNodeIDX <- which.max(grandSurvey$performanceDF$performance)
      selectedNode <- grandSurvey$performanceDF[selectedNodeIDX,]
      selectedNodeDF <- data.frame(response = responseInLoop, dataFrameInLoop[, colnames(dataFrameInLoop) == selectedNode$feature])
      
      # collect prediction
      if(selectedNode$directionType == 'CaseAboveCutoff'){
        ss <- ss_pre(selectedNodeDF, cutoff = selectedNode$cutoff, testCasePositive = 'aboveCutoff', testLabel = T)
        predictedSampleIDX <- which(ss$testLabel == '1')
        predictedInLoop <- data.frame(sampleID = rownames(dataFrameInLoop)[predictedSampleIDX],
                                     casePrediction = 1, controlPrediction = 0, nRound = nRound)
      }
      else if(selectedNode$directionType == 'CaseBelowCutoff'){
        ss <- ss_pre(selectedNodeDF, cutoff = selectedNode$cutoff, testCasePositive = 'belowCutoff', testLabel = T)
        predictedSampleIDX <- which(ss$testLabel == '1')
        predictedInLoop <- data.frame(sampleID = rownames(dataFrameInLoop)[predictedSampleIDX],
                                      casePrediction = 1, controlPrediction = 0, nRound = nRound)
      }
      else if(selectedNode$directionType == 'ControlAboveCutoff'){
        ss <- ss_pre(selectedNodeDF, cutoff = selectedNode$cutoff, testCasePositive = 'belowCutoff', testLabel = T)
        predictedSampleIDX <- which(ss$testLabel == '0')
        predictedInLoop <- data.frame(sampleID = rownames(dataFrameInLoop)[predictedSampleIDX],
                                      casePrediction = 0, controlPrediction = 1, nRound = nRound)
      }
      else if(selectedNode$directionType == 'ControlBelowCutoff'){
        ss <- ss_pre(selectedNodeDF, cutoff = selectedNode$cutoff, testCasePositive = 'aboveCutoff', testLabel = T)
        predictedSampleIDX <- which(ss$testLabel == '0')
        predictedInLoop <- data.frame(sampleID = rownames(dataFrameInLoop)[predictedSampleIDX],
                                      casePrediction = 0, controlPrediction = 1, nRound = nRound)
      }
      
      # stopCriteria check point 2
      nSamplePartitioned <- length(predictedSampleIDX)
      nCaseInLoop <- sum(responseInLoop == 1)
      nControlInLoop <- sum(responseInLoop == 0)
      if(selectedNode$directionType == 'CaseAboveCutoff' | selectedNode$directionType == 'CaseBelowCutoff' ){
        minSamplePerRound <- max ( ceiling(minSamplePerRoundCoef * nCaseInLoop), 3 )
      }
      else{
        minSamplePerRound <- max ( ceiling(minSamplePerRoundCoef * nControlInLoop), 3 )
      }
      if(nSamplePartitioned < minSamplePerRound){
        cat('selected samples are less than expected! \n')
        cat('stop now!\n')
        stopCriteria = T
        
        lastCasePrevalence <- sum(responseInLoop == 1)/length(responseInLoop)
        if(lastCasePrevalence > casePrevalence){
          casePrediction = 1
        }
        else{casePrediction = 0}
        
        lastPredict <- data.frame(sampleID = rownames(dataFrameInLoop), 
                                  casePrediction = casePrediction, controlPrediction = 1 - casePrediction, nRound = nRound)
        predictTable <- rbind(predictTable, lastPredict)
        colnames(predictTable) <- c('sampleID', 'casePrediction', 'controlPrediction', 'nRound')
      }
      else{
        # update dataFrameInLoop, responseInLoop
        dataFrameInLoop <- dataFrameInLoop[-predictedSampleIDX, ]
        responseInLoop <- responseInLoop[-predictedSampleIDX]
        nRemain = nrow(dataFrameInLoop)
        nClass = length(unique(responseInLoop))
        
        # collect parameters after each node selection
        predictTable <- rbind(predictTable, predictedInLoop)
        # selectedNodes[[nRound]] <- selectedNode
        parameters <- rbind(parameters, c(nRound, selectedNode$feature, selectedNode$directionType, selectedNode$performance, selectedNode$cutoff, selectedNode$cutoffPercent, nSamplePartitioned, nRemain))
        colnames(parameters) <- c('nRound', 'selectedNodeName', 'directionType', 'performance', 'cutoff', 'cutoffPercent', 'nSamplePartitioned', 'nRemain')
        
        cat(paste('selected node: ', selectedNode$feature, 
                  '; directionType: ', selectedNode$directionType, 
                  '; performance: ', round(selectedNode$performance, 3), 
                  '; cutoff: ', round(selectedNode$cutoff, 3), 
                  '; cutoff Percent: ', selectedNode$cutoffPercent, 
                  '; nSamplePartitioned: ', nSamplePartitioned, 
                  '; remain sample size: ', nRemain, '\n', sep = ''))
        
        # stopCriteria check point 3
        if(nRound == maxRound | nRemain < minRemain | nClass == 1){
          cat('stop now!\n')
          stopCriteria = T
          lastCasePrevalence <- sum(responseInLoop == 1)/length(responseInLoop)
          if(lastCasePrevalence > casePrevalence){
            casePrediction = 1
          }
          else{casePrediction = 0}
          lastPredict <- data.frame(sampleID = rownames(dataFrameInLoop), 
                                    casePrediction = casePrediction, controlPrediction = 1 - casePrediction, nRound = nRound + 1)
          predictTable <- rbind(predictTable, lastPredict)
          colnames(predictTable) <- c('sampleID', 'casePrediction', 'controlPrediction', 'nRound')
        }
        else{nRound = nRound + 1}
      }

    }
    
  }
  
  z <- list(predictTable = predictTable,
            parameters = data.frame(parameters))
  
}





























