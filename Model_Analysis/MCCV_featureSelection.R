#####################################################################################################################################################################################
## Date: 2016-04-26
## Author: Ryan Wang
##
## Description:
#####################################################################################################################################################################################

##################################################################################################################################################################
# generate a panel of plots from a list of ggplot objects
##################################################################################################################################################################
grid_plot_ggplot2 <- function(ggplotList, ncol, nrow){
  library(gridExtra)
  library(methods)
  plotList <- c(ggplotList, nrow = nrow, ncol = ncol)
  do.call(grid.arrange, plotList)
}

#####################################################################################################################################################################################
## normalization
#####################################################################################################################################################################################
normalization <- function(dataFrame, margin){
  if(margin == 'column'){
    colmeans <- colMeans(dataFrame)
    colstds <- apply(dataFrame, 2, sd)
    dataFrame <- sweep(dataFrame, 2, colmeans, '-')
    dataFrame <- sweep(dataFrame, 2, colstds, '/')
  }else if(margin == 'row'){
    rowmeans <- rowMeans(dataFrame)
    rowstds <- apply(dataFrame, 1, sd)
    dataFrame <- sweep(dataFrame, 1, rowmeans, '-')
    dataFrame <- sweep(dataFrame, 1, rowstds, '/')
  }
  dataFrame
}

#####################################################################################################################################################################################
## make jitter plot
## input: dataFrame, must contain a column called 'response' with the true labels, as well as the attribute of interest
#####################################################################################################################################################################################
jitplot <- function(dataFrame, attr_name){
  library(ggplot2)
  testDF <- data.frame(response = dataFrame$response,
                       data = dataFrame[, attr_name])
  colnames(testDF) = c('response', 'data')
  p = ggplot(testDF, aes(x = factor(response), y = data)) +
    geom_boxplot(notch = FALSE) + 
    # geom_violin() + 
    geom_dotplot(binaxis = 'y', stackdir = 'center', position = 'dodge', binwidth = max(testDF$data)/60) +
    stat_summary(fun.y = median, geom = "point", shape = 18, size = 8, color = "red") +
    labs(x = 'Case Status', title = attr_name) + 
    # scale_fill_grey() + 
    theme(title = element_text(size = 15), axis.title = element_text(size = 20), axis.text = element_text(size = 15),
          panel.grid = element_line(size = 1), legend.position = "none")
  p
}

#####################################################################################################################################################################################
## Feature selection by logistic regression coupled with random subsetting features and repeated monte carlo cross validation.
##
##    During each iteration, a certain number of features are randomly selected (numberOfFeaturePerIter), then MCCV are performed for certain 
## times (numberOfMCCV). After MCCV, mean predicted values are used to compute the performance metrics (pAUCs). This process is repeated
## iteratively for certain times (numberOfFeatureSets). Eventually, the performance metrics for each feature set during the iteration
## are collected as output. Parallel computing is implemented for the feature sets iteration.
##
##    Currently, only the following pAUCs are calculated:
##    paucInLoop_90sens, paucInLoop_95sens, paucInLoop_98sens, paucInLoop_90spec, paucInLoop_95spec, paucInLoop_98spec, aucInLoop
## 
## input: dataFrame, columns are attributes, rows are samples. Last column is the response variable.
## output: performSummary, the performance metrics table for all feature sets evaluated.
#####################################################################################################################################################################################
MCCV_evaluation <- function(
  dataFrame, 
  # number of features selected for each feature set per iteration
  numberOfFeaturePerIter = 20, 
  # number of feature sets to evaluate
  numberOfFeatureSets = 100, 
  numberOfThread = 10, 
  # number of MCCV for each feature set
  numberOfMCCV = 100,
  # fraction of cases to holdout during MCCV
  holdout = 0.2
){

  library(doParallel)
  library(foreach)
  library(cvTools)
  
  nSample <- nrow(dataFrame)
  nFeature <- ncol(dataFrame) - 1
  featureVector <- colnames(dataFrame[,1:nFeature])
  response <- dataFrame[, nFeature + 1]
  
  # log10 transform and normalize by feature (column)
  dataFrame <- data.frame(normalization(log10(dataFrame[, -(nFeature + 1)]), margin = 'column'),
                          response = response)
  
  # if(numberOfFeatureSets%%numberOfThread != 0){
  #   stop('numberOfFeatureSets cannot be exactly divided by numberOfThread.\n
  #        Please input another set of numberOfFeatureSets and numberOfThread!\n')
  # }
  foldsPerThread <- ceiling(numberOfFeatureSets/numberOfThread)
  
  cl <- makeCluster(numberOfThread)
  registerDoParallel(cl)
  system('mkdir logs')
  # subset <- cvFolds(numberOfFeatureSets, numberOfThread)
  performSummary <- foreach( i = 1:numberOfThread, .combine = rbind, .export = c()) %dopar% {
    library(pROC)
    performTab <- c()
    log_file <- paste('logs/', Sys.Date(), '_log_thread', i, '.txt', sep = '')
    sink(log_file, append = T, type = 'output')
    
    for( j in 1:foldsPerThread ){
      cat('current iteration:', j, ';', foldsPerThread - j, 'to go\n')
      # sampling features
      fs_FeatureIDX <- sample(nFeature, numberOfFeaturePerIter, replace = F)
      
      # repeat MCCV numberOfMCCV times
      predictRepeat <- data.frame(names = rownames(dataFrame), stringsAsFactors = F)
      for(fold in 1:numberOfMCCV){
        # cat(fold, ';')
        fs_trainIDX <- sample(nSample, ceiling((1 - holdout) * nSample), replace = F)
        trainDF <- dataFrame[fs_trainIDX, c(fs_FeatureIDX, (nFeature + 1))]
        testDF <- dataFrame[-fs_trainIDX, c(fs_FeatureIDX, (nFeature + 1))]
        
        fitGLMInLoop <- glm(response~., data = trainDF, family = binomial())
        predictInLoop <- predict(fitGLMInLoop, newdata = testDF, type = 'response')
        predictInLoop <- data.frame(names = names(predictInLoop), test_probability = predictInLoop)
        predictRepeat <- suppressWarnings(merge(predictRepeat, predictInLoop, by.x = 'names', by.y = 'names', all = TRUE))
      }
      # the mean predicted value during the MCCV is collected
      predictGLM <- apply(predictRepeat[,-1], 1, function(x)mean(as.numeric(x), na.rm = T))
      names(predictGLM) <- predictRepeat$names
      
      # compute pAUCs using the mean cross validation prediction
      paucInLoop_90sens <- auc(response, predictGLM, partial.auc.focus = 'sensitivity', partial.auc = c(0.9, 1))
      paucInLoop_95sens <- auc(response, predictGLM, partial.auc.focus = 'sensitivity', partial.auc = c(0.95, 1))
      paucInLoop_98sens <- auc(response, predictGLM, partial.auc.focus = 'sensitivity', partial.auc = c(0.98, 1))
      paucInLoop_90spec <- auc(response, predictGLM, partial.auc.focus = 'specificity', partial.auc = c(0.9, 1))
      paucInLoop_95spec <- auc(response, predictGLM, partial.auc.focus = 'specificity', partial.auc = c(0.95, 1))
      paucInLoop_98spec <- auc(response, predictGLM, partial.auc.focus = 'specificity', partial.auc = c(0.98, 1))
      aucInLoop <- auc(response, predictGLM)
      
      performInLoop <- data.frame(t(featureVector[fs_FeatureIDX]), 
                                  paucInLoop_90sens,
                                  paucInLoop_95sens,
                                  paucInLoop_98sens,
                                  paucInLoop_90spec,
                                  paucInLoop_95spec,
                                  paucInLoop_98spec,
                                  aucInLoop)
      performTab <- rbind(performTab, performInLoop)
    }
    performTab
  }
  stopCluster(cl)
  # collect parameters
  parameters <- list(numberOfFeaturePerIter = numberOfFeaturePerIter,
                     numberOfFeatureSets = foldsPerThread * numberOfThread,
                     numberOfMCCV = numberOfMCCV,
                     holdout = holdout)
  
  z <- list(performSummary = performSummary, parameters = parameters)
}

#####################################################################################################################################################################################
## Summarize the frequencies of each feature appeared in a certain number of top panels (numberOfTopPanelsToAnalyze). 
## input: MCCV_evaluation, output from MCCV_evaluation() function, the performance metrics and feature sets to analyze.
##        focus, the specific metric to rank
## 
#####################################################################################################################################################################################
featureScore <- function(MCCV_evaluation, dataFrame, 
                           focus = 'aucInLoop', 
                           numberOfTopPanelsToAnalyze = 100){
  
  nFeatures <- ncol(dataFrame) - 1
  numberOfFeaturePerIter <- MCCV_evaluation$parameters$numberOfFeaturePerIter
  performSummary <- MCCV_evaluation$performSummary
  
  # check if the focus is properly specified
  if(!focus%in%colnames(performSummary)){
    stop('Parameter focus is not properly specified!\nPlease only choose from the following options:\npaucInLoop_90sens, paucInLoop_95sens, paucInLoop_98sens, paucInLoop_90spec, paucInLoop_95spec, paucInLoop_98spec, aucInLoop\n')
  }
  
  # compute the number of expected frequencies after random selection
  expRandomFreq <- numberOfFeaturePerIter / nFeatures * numberOfTopPanelsToAnalyze
  # compute the feature frequencies in the top panels
  Features <- colnames(dataFrame)[1:nFeatures]
  FeaturesInTopPanels <- performSummary[order(performSummary[, focus], decreasing = T)[1:numberOfTopPanelsToAnalyze], 1:numberOfFeaturePerIter]
  FeatureVec <- as.vector(t(FeaturesInTopPanels))
  FeatureFreq <- sapply(1:nFeatures, function(x){
    sum(FeatureVec == Features[x])
  })
  FreqSummary <- data.frame(Features = Features, FeatureFreq = FeatureFreq, score = FeatureFreq/expRandomFreq)
  FreqSummary <- FreqSummary[order(FreqSummary$FeatureFreq, decreasing = T), ]
  
  z <- list(FreqSummary = FreqSummary,
            expRandomFreq = expRandomFreq)
}


#####################################################################################################################################################################################
## stepwise fitting logistic regression model removing features with unstable coefficient
#####################################################################################################################################################################################
stepFeatureAnalysis <- function(
  dataFrame, 
  # number of iterations during MCCV process
  numberOfIter = 30, 
  # number of threads for parallel processing, 
  # the actual numberOfIter is ceiling(numberOfIter/numberOfThread) * numberOfThread
  numberOfThread = 3,
  holdout = 0.2
){
  library(doParallel)
  library(foreach)
  library(pROC)
  
  nSample <- nrow(dataFrame)
  numberOfFeatures <- ncol(dataFrame) - 1
  response <- dataFrame[, (numberOfFeatures + 1)]
  
  # collect remaining feature during each step
  feature_step <- list()
  # collect the logistic regression model during each step
  fitGLM_step <- list()
  # collect the order of feature that gets removed
  FeatureRemoveOrder <- c()
  # collect the standard deviation of each features during each step
  coef_std_step <- list()
  # collect the predicted values of MCCV during each step
  predictGLM_step <- list()
  # collect the performance metrics during each step
  performTab <- c()
  
  dataFrameInLoop <- data.frame(normalization(log10(dataFrame[, -(numberOfFeatures + 1)]), margin = 'column'),
                                response = response)
  # start the step wise process. During each round, the least stable feature gets removed from the dataset (dataFrameInLoop)
  for(i in 1: (numberOfFeatures - 1)){
    cat('round', i)
    # build model and collect information using all samples and remaining features at the first step
    fitGLM_step[[i]] <- glm(response~., data = dataFrameInLoop, family = binomial())
    feature_step[[i]] <- colnames(dataFrameInLoop)[1:(ncol(dataFrameInLoop) - 1)]
    
    # assess feature coefficient stability by MCCV
    ## multithread processing
    foldsPerThread <- ceiling( numberOfIter/numberOfThread )
    cl <- makeCluster(numberOfThread)
    registerDoParallel(cl)
    
    paraCombine <- function(x,y){
      # predictRepeat is a list, canbe combined by c() during foreach()
      predictRepeatList <- c(x$predictRepeat, y$predictRepeat)
      # coefs is a matrix, canbe combined by rbind() during foreach()
      coefs <- rbind(x$coefs, y$coefs)
      z <- list(predictRepeatList = predictRepeatList, coefs = coefs)
    }
    
    ## MCCV on glm models in parallel
    MCCV_analysis <- foreach( j = 1:numberOfThread, .combine = 'paraCombine') %dopar% {
      coefGLM <- c()
      predictRepeat <- data.frame(names = rownames(dataFrame), stringsAsFactors = F)
      # within each thread, loop to perform MCCV
      for(fold in 1:foldsPerThread){
        trainIDX <- sample(nSample, ceiling((1 - holdout) * nSample), replace = F)
        trainDF <- dataFrameInLoop[trainIDX, ]
        testDF <- dataFrameInLoop[-trainIDX, ]
        fitGLMInLoop <- glm(response~., data = trainDF, family = binomial())
        coefGLM <- rbind(coefGLM, coef(fitGLMInLoop))
        
        predictInLoop <- predict(fitGLMInLoop, newdata = testDF, type = 'response')
        predictInLoop <- data.frame(names = names(predictInLoop), test_probability = predictInLoop)
        predictRepeat <- suppressWarnings(merge(predictRepeat, predictInLoop, by.x = 'names', by.y = 'names', all = TRUE))
      }
      # very important step, to make predictRepeat a list so that it can be combined by the paraCombine() function using c()
      predictRepeat = list(predictRepeat)
      # collect predicted values and coefficients 
      z <- list(predictRepeat = predictRepeat, coefs = coefGLM)
    }
    stopCluster(cl)
    
    # compute performance after a MCCV process
    predictRepeat <- MCCV_analysis$predictRepeatList[[1]]
    # merge predicted values among threads after parallel computing
    for(i_thread in 1:(numberOfThread - 1)){
      predictRepeat <- suppressWarnings(merge(predictRepeat, MCCV_analysis$predictRepeatList[[i_thread+1]], by.x = 'names', by.y = 'names', all = TRUE))
    }
    colnames(predictRepeat) <- c('names', paste('iter', c(1:numberOfIter), sep = '_'))
    # compute the average predicted value for all MCCV iterations
    predictGLM <- apply(predictRepeat[,-1], 1, function(x)mean(as.numeric(x), na.rm = T))
    names(predictGLM) <- predictRepeat$names
    predictGLM_step[[i]] <- predictGLM
    
    # compute the performance metrics
    paucInLoop_90sens <- auc(response, predictGLM, partial.auc.focus = 'sensitivity', partial.auc = c(0.9, 1))
    paucInLoop_95sens <- auc(response, predictGLM, partial.auc.focus = 'sensitivity', partial.auc = c(0.95, 1))
    paucInLoop_98sens <- auc(response, predictGLM, partial.auc.focus = 'sensitivity', partial.auc = c(0.98, 1))
    paucInLoop_90spec <- auc(response, predictGLM, partial.auc.focus = 'specificity', partial.auc = c(0.9, 1))
    paucInLoop_95spec <- auc(response, predictGLM, partial.auc.focus = 'specificity', partial.auc = c(0.95, 1))
    paucInLoop_98spec <- auc(response, predictGLM, partial.auc.focus = 'specificity', partial.auc = c(0.98, 1))
    aucInLoop <- auc(response, predictGLM)
    
    performInLoop <- data.frame(round = i, 
                                paucInLoop_90sens,
                                paucInLoop_95sens,
                                paucInLoop_98sens,
                                paucInLoop_90spec,
                                paucInLoop_95spec,
                                paucInLoop_98spec,
                                aucInLoop)
    performTab <- rbind(performTab, performInLoop)
    
    # evaluate the stability of coefficients for each feature
    coef_std <- apply(MCCV_analysis$coefs, 2, sd, na.rm = T)
    coef_std_step[[i]] <- coef_std
    # remove feature with the least stability at each step
    rm_idx <- which.max(coef_std[-1])
    # collect features by the order of removal
    FeatureRemoveOrder <- c(FeatureRemoveOrder, colnames(dataFrameInLoop)[rm_idx])
    cat(': removed feature is', colnames(dataFrameInLoop)[rm_idx], '\n')
    # update the data frame by removing the feature
    dataFrameInLoop <- dataFrameInLoop[, -rm_idx]
  }
  
  z <- list(fitGLM_step = fitGLM_step, 
            feature_step = feature_step,
            coef_std_step = coef_std_step,
            FeatureRemoveOrder = FeatureRemoveOrder,
            performTab = performTab,
            predictGLM_step = predictGLM_step,
            numberOfIter = numberOfIter)
}

#####################################################################################################################################################################################
#####################################################################################################################################################################################
getMetrics <- function(  
  dataFrame,
  prev = 0.2
){
  colnames(dataFrame) <- c('response', 'predicted_value')
  library(ROCR)
  library(data.table)
  library(ggplot2)
  ROCR_prediction <- prediction(dataFrame$predicted_value, dataFrame$response)
  sensitivity <- performance(ROCR_prediction, measure = 'sens')@y.values[[1]]
  specificity <- performance(ROCR_prediction, measure = 'spec')@y.values[[1]]
  cutoffs <- performance(ROCR_prediction, measure = 'tpr')@x.values[[1]]
  auc <- round(performance(ROCR_prediction, measure = 'auc')@y.values[[1]], 3)
  npv <- (specificity * (1 - prev)) / (specificity * (1 - prev) + (1 - sensitivity) * prev)
  ppv <- (sensitivity * prev) / (sensitivity * prev + (1 - specificity) * (1 - prev))
  
  df_metrics <- data.table(cutoffs = rep(cutoffs, 4),
                           values = c(sensitivity, specificity, ppv, npv),
                           metrics = c(rep('sensitivity', length(sensitivity)), 
                                       rep('specificity', length(specificity)), 
                                       rep('ppv', length(ppv)), 
                                       rep('npv', length(npv))))
  z <- list(df_metrics = df_metrics, prev = prev)
}

#####################################################################################################################################################################################
## plot sensitivity, specificity, ppv and/or npv against cutoffs for each model
## input: dataFrame, first column is true label, second column is predicted label.
##        prev, assumed prevalence for npv, ppv calculation
##        metric_to_plot, select the desired metric to plot
## output: a ggplot object
#####################################################################################################################################################################################
metric_plot <- function(
  df_metrics, 
  prev,
  metric_to_plot = c('sensitivity', 'specificity', 'npv', 'ppv'), 
  add_line = 100
){
  
  df_metrics_to_plot <- df_metrics[df_metrics$metrics %in% metric_to_plot, ]
  metric_plot <- ggplot(df_metrics_to_plot, aes(x = cutoffs, y = values, color = metrics)) +
    geom_line(size = 2) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    geom_vline(xintercept = add_line, color = 'black') + 
    labs(x = 'cutoffs', y = 'percentatge', title = paste('Performance by cutoffs (Prevalence=', prev, ')', sep = '')) + 
    theme(title = element_text(size = 20, face = 'bold'), axis.title = element_text(size = 20), axis.text = element_text(size = 20),
          panel.background = element_rect(fill = "white"),
          axis.line.x = element_line(colour = "black", size = 1, linetype = "solid"),
          axis.line.y = element_line(colour = "black", size = 1, linetype = "solid"),
          legend.title = element_text(size = 20), legend.text = element_text(size = 10), 
          legend.background = element_rect(fill="white"), 
          legend.key.size = unit(15, units = 'mm'))
  metric_plot
}


#####################################################################################################################################################################################
## plot roc curve
## input: dataFrame, first column is true label, second column is predicted label.
## output: a ggplot object
#####################################################################################################################################################################################
roc_plot <- function(
  dataFrame
){
  colnames(dataFrame) <- c('response', 'predicted_value')
  library(ROCR)
  library(data.table)
  library(ggplot2)
  
  ROCR_prediction <- prediction(dataFrame$predicted_value, dataFrame$response)
  sensitivity <- performance(ROCR_prediction, measure = 'sens')@y.values[[1]]
  specificity <- performance(ROCR_prediction, measure = 'spec')@y.values[[1]]
  cutoffs <- performance(ROCR_prediction, measure = 'tpr')@x.values[[1]]
  auc <- round(performance(ROCR_prediction, measure = 'auc')@y.values[[1]], 3)
  df_metrics_to_plot <- data.table(TPR = sensitivity, 
                                   FPR = 1 - specificity)
  roc_plot <- ggplot(df_metrics_to_plot, aes(x = FPR, y = TPR)) +
    geom_line(size = 2) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    geom_abline(intercept = 0, slope = 1) + 
    labs(x = 'FPR', y = 'TPR', title = paste('AUC = ', auc, sep = '')) + 
    theme(title = element_text(size = 25, face = 'bold'), axis.title = element_text(size = 20), axis.text = element_text(size = 20),
          panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA, color = 'black', size = 1.5), 
          legend.title = element_text(size = 20), legend.text = element_text(size = 10), 
          legend.background = element_rect(fill="white"), 
          legend.key.size = unit(15, units = 'mm'))
  roc_plot
}























