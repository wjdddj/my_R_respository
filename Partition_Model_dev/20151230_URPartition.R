options(stringsAsFactors = F)

################################################################################################################################################
# label outliers and determine cutoff point
################################################################################################################################################
DataSplit <- function(dataVec, method = 'outlier', 
                      # if method == 'outlier', coefficient is coef*mid50percentile
                      # if method == 'quantile', coefficient is quantile cutoff
                      coef = 2){
  if(method == 'outlier'){
    splitPoint = quantile(dataVec, 0.75) + coef * (quantile(dataVec, 0.75) - quantile(dataVec, 0.25))
    splitLabel <- dataVec > splitPoint
    nHighSample <- sum(splitLabel)
  }
  
  if(method == 'quantile'){
    splitPoint <- quantile(dataVec, coef)
    splitLabel <- dataVec >= splitPoint
    nHighSample <- sum(splitLabel)
  }
  
  z <- list(splitPoint = splitPoint, splitLabel = splitLabel, nHighSample = nHighSample)
}


################################################################################################################################################
# build submodels by either rpart() 
# or assign probability to 1 when already pure
# dataFrame: the outlier dataFrame
################################################################################################################################################
rPartLoop <- function(dataFrame){
  library(rpart)
  library(pROC)
  n <- ncol(dataFrame)
  response = dataFrame[,n]
  
  if(length(unique(dataFrame$response))==1){
    predLabel <- unique(dataFrame$response)
    if(predLabel == 1){
      predictLoop <- cbind(rep(0, nrow(dataFrame)), rep(1, nrow(dataFrame)))
    }else{
      predictLoop <- cbind(rep(1, nrow(dataFrame)), rep(0, nrow(dataFrame)))
    }
    rownames(predictLoop) <- rownames(dataFrame)
    colnames(predictLoop) <- c(0,1)
    AUC = 1
    z <- list(predictLoop = predictLoop, AUC = AUC)
  }
  else{
    fitRP <- rpart(response~., data = dataFrame)
    predictLoop <- predict(fitRP, newdata = dataFrame)
    AUC <- roc(response, predictLoop[,2])$auc
    z <- list(predictLoop = predictLoop, rPartObej = fitRP, AUC = AUC)
  }
  z
}

################################################################################################################################################
# scan through all possible splits in parallel
# dataFrame: the input dataFrame at the begining of each round
# DataSplits: the output by lapply()::DataSplit()
################################################################################################################################################
rPartPara <- function(dataFrame, DataSplits, 
                      # parameters to define DfExcludeThreshold, when number of outliers is less than DfExcludeThreshold, the split will be excluded from analysis
                      # DfExcludeThreshold <- max(ceiling(DfExcludeCoef * numberOfSample), minSampleSizeToSplit)
                      DfExcludeCoef = 0.1, 
                      minSampleSizeToSplit = 5, 
                      # number of threads to run in parallel
                      nThreads = 10){
  library(doParallel)
  library(foreach)
  cl <- makeCluster(nThreads)
  registerDoParallel(cl)
  
  # combine function for parallel processes
  rPartCombine1 <- function(x,y){
    rankTable <- rbind(x$rankTable, y$rankTable)
    rPartModels <- c(x$rPartModels, y$rPartModels)
    rPartDFs <- c(x$rPartDFs, y$rPartDFs)
    dataSplitsPara <- c(x$dataSplitsPara, y$dataSplitsPara)
    #splitLabelPara <- c(x$splitLabelPara, y$splitLabelPara)
    #splitPointPara <- c(x$splitPointPara, y$splitPointPara)
    z <- list(rankTable = rankTable, rPartModels = rPartModels, rPartDFs = rPartDFs, dataSplitsPara = dataSplitsPara)
  }
  
  numberOfCol <- ncol(dataFrame)
  response = dataFrame[, numberOfCol]
  numberOfFeature <- numberOfCol - 1
  numberOfSample <- nrow(dataFrame)
  
  splitLabelsInput <- lapply(DataSplits, function(x)x$splitLabel)

  # split features into batches for parallel processes
  numberOfSeqPerBatch <- ceiling(numberOfFeature / nThreads)
  batchLabel <- c()
  for(j in 1:(nThreads-1)){
    batchLabel <- c(batchLabel, rep(j, numberOfSeqPerBatch))
  }
  batchLabel <- c(batchLabel, rep(j+1, numberOfFeature - (nThreads-1) * numberOfSeqPerBatch))

  # start parallel processes
  splitRpart <- foreach (i = 1:nThreads, .combine = rPartCombine1, .export = c('rPartLoop', 'DataSplit')) %dopar% {
    library(rpart)
    library(pROC)

    dataFramePara <- dataFrame[, -numberOfCol][, batchLabel == i]
    responsePara <- response
    dataSplitsPara <- DataSplits[batchLabel == i]
    splitLabelPara <- splitLabelsInput[batchLabel == i]

    
    # remove splits when splitted sample size is less than allowed
    DfExcludeThreshold <- max(ceiling(DfExcludeCoef * numberOfSample), minSampleSizeToSplit)
    dataSplitsPara[sapply(splitLabelPara, function(x)sum(x)) < DfExcludeThreshold] <- NULL
    splitLabelPara[sapply(splitLabelPara, function(x)sum(x)) < DfExcludeThreshold] <- NULL

    # prepare dataframes for rpart() trees
    rPartDFs <- lapply(splitLabelPara, function(x){
      data.frame(dataFramePara[x,], response = responsePara[x])
    })
    
    # fit and predict training data 
    rPartModels <- lapply(rPartDFs, function(x){
      rPartLoop(x)
    })
    
    # collect AUCs and sizes
    AUCs <- unlist(sapply(rPartModels, function(x)x$AUC))
    rankTable <- data.frame(AUCs = AUCs, size = unlist(sapply(rPartDFs, nrow)))
    
    outList <- list(rankTable = rankTable, rPartModels = rPartModels, rPartDFs = rPartDFs, dataSplitsPara = dataSplitsPara)
  }
  stopCluster(cl)
  splitRpart
}


################################################################################################################################################
# unbiased recursive partition model, requires rpart package
################################################################################################################################################
urPartition_v1 <- function(dataFrame, 
                           # stop recursive partitioning when number of round exceed maxRound
                           maxRound = 50, 
                           # stop recursive partitioning when number of remaining samples below minRemain
                           minRemain = 20, 
                           # stop recursive partitioning when submodel AUC is less than minAUC
                           minAUC = 0.6){
  
  n <- ncol(dataFrame)
  response <- dataFrame[,n]
  casePrevalence <- sum(response==1)/length(response)
  
  # split samples based on each feature
  dataSplits <- lapply(dataFrame[,-n], function(x){
    split <- DataSplit(x, method = 'outlier', coef = 2)
  })
  
  # initiate loop
  dataFrameInLoop <- dataFrame
  responseInLoop <- dataFrame[,n]
  dataSplitsInLoop <- dataSplits
  nRound <- 1
  stopCritera = F
  
  # output parameters
  predictTable <- c()
  selectedNodes <- list()
  parameters <- c()
  partitionModels <- list()
  while(stopCritera == F){
    cat(paste('node ', nRound, ', ', sep = ''))
    
    # survey all possible splits in parallel, compute AUC for each possible split
    grandFitRPart <- rPartPara(dataFrameInLoop, DataSplits = dataSplitsInLoop, DfExcludeCoef = 0.05, minSampleSizeToSplit = 10, nThreads = 10)
    
    # examine number of possible splits
    # stop critera check 1, before update input: nSplitDFs
    nSplitDFs <- length(grandFitRPart$rPartDFs)
    if(nSplitDFs == 0){
      cat('no more outlier samples! \n')
      cat('stop now!\n')
      stopCritera = T
      
      lastPredict <- cbind(rep(1-casePrevalence, nRemain), rep(casePrevalence, nRemain), rep(nRound, nRemain))
      rownames(lastPredict) <- rownames(dataFrameInLoop)
      predictTable <- rbind(predictTable, lastPredict)
      colnames(predictTable) <- c('Prob0', 'Prob1', 'nRound')
    }
    else{
      # select the best node by ranking the performance after each possible split
      selectedNodeIDX <- order(grandFitRPart$rankTable$AUCs, grandFitRPart$rankTable$size, decreasing = T)[1]
      selectedNode <- grandFitRPart$dataSplitsPara[[selectedNodeIDX]]
      selectedNodeName <- names(grandFitRPart$dataSplitsPara)[selectedNodeIDX]
      splitCutoff <- selectedNode$splitPoint
      selectedNodes[[nRound]] <- selectedNode
      nSampleSelected <- selectedNode$nHighSample
      selectedAUC <- grandFitRPart$rankTable$AUCs[selectedNodeIDX]
      
      # stop critera check 2, before update input: AUC
      if(selectedAUC < minAUC){
        cat('AUC below cutoff! \n')
        cat('stop now!\n')
        stopCritera = T
        lastPredict <- cbind(rep(1-casePrevalence, nRemain), rep(casePrevalence, nRemain), rep(nRound, nRemain))
        rownames(lastPredict) <- rownames(dataFrameInLoop)
        predictTable <- rbind(predictTable, lastPredict)
        colnames(predictTable) <- c('Prob0', 'Prob1', 'nRound')
      }
      else{
        # update input data
        sampleToExclude <- which(selectedNode$splitLabel)
        dataFrameInLoop <- dataFrameInLoop[-sampleToExclude,]
        dataSplitsInLoop <- lapply(dataSplitsInLoop, function(x){
          x$splitLabel <- x$splitLabel[-sampleToExclude]
          x$nHighSample <- sum(x$splitLabel)
          x
        })
        nRemain <- nrow(dataFrameInLoop)
        nClass <- length(unique(dataFrameInLoop$response))
        
        # collect the submodels and parameters after each selected node
        partitionModels[[nRound]] <- grandFitRPart$rPartModels[[selectedNodeIDX]]
        predictTable <- rbind(predictTable, cbind(grandFitRPart$rPartModels[[selectedNodeIDX]]$predictLoop, nRound))
        parameters <- rbind(parameters, c(nRound, selectedNodeName, splitCutoff, nSampleSelected, selectedAUC, nRemain))
        colnames(parameters) <- c('nRound', 'selectedNodeName', 'splitCutoff', 'nSampleSelected', 'selectedAUC', 'nRemain')
        
        cat(paste('selected node: ', selectedNodeName, 
                  '; splitCutoff: ', splitCutoff, 
                  '; splitDF sample size: ', nSampleSelected, 
                  '; AUC: ', round(selectedAUC, 3), 
                  '; remain sample size: ', nRemain, '\n', sep = ''))
        
        # stop critera check 3, after update input: nRound, nRemain, nClass
        if(nRound > (maxRound - 1) | nRemain < minRemain | nClass == 1){
          cat('stop now!\n')
          stopCritera = T
          lastPredict <- cbind(rep(1-casePrevalence, nRemain), rep(casePrevalence, nRemain), rep(nRound, nRemain))
          rownames(lastPredict) <- rownames(dataFrameInLoop)
          predictTable <- rbind(predictTable, lastPredict)
          colnames(predictTable) <- c('Prob0', 'Prob1', 'nRound')
        }
        else{nRound = nRound + 1}
      }
    }
  }
  
  z <- list(selectedNodes = selectedNodes, 
            partitionModels = partitionModels, 
            predictTable = data.frame(predictTable), 
            parameters = data.frame(parameters), 
            lastPredict = lastPredict)
}

################################################################################################################################################
# unbiased partition model, predict function
################################################################################################################################################
predictUP <- function(urPartition, newData){
  
  # initiate loop to identify proper submodel for each sample
  nRounds <- length(urPartition$parameters$nRound)
  newDataInLoop <- newData
  nRound = 1
  stopCriteria = F
  SampleModelLookUp <- c()
  nRemain = nrow(newData)
  
  cat('assigning submodels...\n')
  while(stopCriteria == F){
    selectedNodeName <- urPartition$parameters$selectedNodeName[nRound]
    cat(paste('node ', nRound, ', nodeName: ', selectedNodeName, '; ', sep = ''))
    splitSampleIDX <- which(newDataInLoop[,colnames(newDataInLoop)==selectedNodeName] > as.numeric(urPartition$parameters$splitCutoff[nRound]))
    if(length(splitSampleIDX)==0){
      nRound = nRound + 1
      cat('no sample above cutoff, skip this node\n')
    }
    else{
      splitSample <- rownames(newDataInLoop)[splitSampleIDX]
      # cat(paste(length(splitSampleIDX), '\n'))
      SampleModelLookUp <- rbind(SampleModelLookUp, cbind(splitSample, nRound))
      newDataInLoop <- newDataInLoop[-splitSampleIDX,]
      nRemain <- nrow(newDataInLoop)
      nRound = nRound + 1
      cat(paste('split sample size: ', length(splitSampleIDX), '; remain sample size: ', nRemain, '\n', sep = ''))
    }
    
    # stop criteria check
    if(nRound > nRounds | nRemain == 0){
      stopCriteria = T
      cat('stop now!\n')
    }
  }
  SampleModelLookUp <- data.frame(SampleModelLookUp)
  SampleModelLookUp$nRound <- as.numeric(SampleModelLookUp$nRound)
  
  # predict predictable samples using their corresponding submodel
  predictableIDX <- match(SampleModelLookUp$splitSample, rownames(newData))
  newDataPredictable <- newData[predictableIDX, ]
  ModelToUse <- unique(SampleModelLookUp$nRound)
  nModelToUse <- length(ModelToUse)
  cat(paste('total number of nodes: ', nRounds, '\n', sep = ''))
  cat(paste('number of submodels to predict: ', nModelToUse, '\n', sep = ''))
  predOut <- c()
  
  cat('start prediction...\n')
  for(i in 1:nModelToUse){
    cat(paste('node ', i, ', ', sep = ''))
    ModelLoop <- urPartition$partitionModels[[ModelToUse[i]]]$rPartObej
    newDataPredictLoop <- newDataPredictable[SampleModelLookUp$nRound==ModelToUse[i],]
    
    # examine if splitSamples are already pure or requires rpart
    if(is.null(ModelLoop)){
      cat('submodel is null, assigning predict probability...\n')
      ClassType <- c(0,1)[urPartition$partitionModels[[ModelToUse[i]]]$predictLoop[1,]==1]
      if(ClassType==1){
        predictLoop <- cbind(rep(0, nrow(newDataPredictLoop)), rep(1, nrow(newDataPredictLoop)), ModelToUse[i])
        rownames(predictLoop) <- rownames(newDataPredictLoop)
        predOut <- rbind(predOut, predictLoop)
      }
      else{
        predictLoop <- cbind(rep(1, nrow(newDataPredictLoop)), rep(0, nrow(newDataPredictLoop)), ModelToUse[i])
        rownames(predictLoop) <- rownames(newDataPredictLoop)
        predOut <- rbind(predOut, predictLoop)
      }
    }
    else{
      cat('fitting rpart...\n')
      predictLoop <- cbind(predict(ModelLoop, newdata = newDataPredictLoop), ModelToUse[i])
      predOut <- rbind(predOut, predictLoop)
    }
  }
  cat(paste('predictable sample size: ', nrow(predOut),'\n', sep = ''))
  predOut <- cbind(predOut, TRUE)
  
  # assign 50% probability for remaining unpredictable samples
  if(nRemain != 0){
    unPredictable <- cbind(rep(0.5, nRemain), rep(0.5, nRemain), rep(ModelToUse[i]+1, nRemain), rep(FALSE, nRemain))
    rownames(unPredictable) <- rownames(newDataInLoop)
    cat(paste('unpredictable sample size: ', nrow(unPredictable),'\n', sep = ''))
    predOut <- rbind(predOut, unPredictable)
  }
  
  colnames(predOut) <- c('Prob0', 'Prob1', 'ModelToUse', 'Predictable')
  matchIDX <- match(rownames(newData), rownames(predOut))
  predOut <- data.frame(predOut[matchIDX, ])
  predOut
}

# roc(response[match(rownames(predOut), rownames(runDF))], predOut[,2], plot = T)


################################################################################################################################################
# unbiased partition model, cross-validation function
################################################################################################################################################
urPartition_cv <- function(dataFrame, K = 10, R = 10){
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
      
      fitUP <- urPartition_v1(trainDF, maxRound = 50, minRemain = 20, minAUC = 0.6)
      cat('predicting...\n')
      pred <- predictUP(fitUP, newData = testDF)
      cat(paste('test sample size: ', length(idx_test), '\n', sep =''))
      outTab <- cbind(idx_test, pred)
      predAll <- rbind(predAll, outTab)
      
    }
    
    colnames(predAll) <- c('testIDX', 'Prob0', 'Prob1', 'ModelToUse', 'Predictable')
    predAll <- data.frame(predAll)
    predOut[[length(predOut)+1]] <- predAll
  }
  z <- list(predOut = predOut)
}



########################################################################################################################################################################
load('~/Documents/Project1_ADAPT/894Probing/20151208_Day1ToDay4/QC/DataTables/2015-12-08-DataTable_Norm_Average.RData')
response <- AvgSampleInfo$CaseStatus
response[response!='Cancer'] <- 0
response[response=='Cancer'] <- 1
runDF <- data.frame(t(AvgNormDF$AvgDF), response = as.factor(response))

fiturPartition <- urPartition_v1(runDF, maxRound = 20, minRemain = 20, minAUC = 0.6)
fiturPartition$predictTable
pdf('20151230_TrainAUC_ExpID4574_urPartition.pdf', 8, 8)
AUC <- round(roc(response[match(rownames(fiturPartition$predictTable), rownames(runDF))], fiturPartition$predictTable[,2])$auc, 4)
roc(response[match(rownames(fiturPartition$predictTable), rownames(runDF))], fiturPartition$predictTable[,2], plot = T, main = paste('AUC = ', AUC, sep = ''), cex.main = 2)
dev.off()

predictableSamples <- fiturPartition$predictTable[fiturPartition$predictTable$nRound!=19,]
pdf('20151230_TrainAUC_ExpID4574_urPartition_onlyPredictableSamples.pdf', 8, 8)
AUC_predictable <- round(roc(response[match(rownames(predictableSamples), rownames(runDF))], predictableSamples[,2])$auc, 4)
roc(response[match(rownames(predictableSamples), rownames(runDF))], predictableSamples[,2], plot = T, main = paste('AUC = ', AUC_predictable, sep = ''), cex.main = 2)
dev.off()


fitUPCV <- urPartition_cv(runDF, K = 10, R = 1)
rownames(fitUPCV$predOut[[1]])
pdf('20151230_10FoldCVAUC_ExpID4574_urPartition.pdf', 8, 8)
AUC_cv <- round(roc(response[fitUPCV$predOut[[1]][,1]], fitUPCV$predOut[[1]][,3])$auc, 4)
roc(response[fitUPCV$predOut[[1]][,1]], fitUPCV$predOut[[1]][,3], plot = T, main = paste('AUC = ', AUC_cv, sep = ''), cex.main = 2)
dev.off()

predictableSamples_CV <- fitUPCV$predOut[[1]][fitUPCV$predOut[[1]]$Predictable==1, ]
pdf('20151230_10FoldCVAUC_ExpID4574_urPartition_onlyPredictableSamples.pdf', 8, 8)
AUC_cv_predictable <- round(roc(response[predictableSamples_CV[,1]], predictableSamples_CV[,3])$auc, 4)
roc(response[predictableSamples_CV[,1]], predictableSamples_CV[,3], plot = T, main = paste('AUC = ', AUC_cv_predictable, sep = ''), cex.main = 2)
dev.off()
























