#########################################################################################################################################################
#### MODULE 5: Main functions for data analysis
#########################################################################################################################################################

#########################################################################################################################################################
# Perform univariate permutation analysis
#########################################################################################################################################################
Main_UnivariatePerm <- function(AvgNormDF, AvgSampleInfo, ANOVA = T, CancerVsRest = T, CancerVsNonCancer = T, CancerVsNormal = T, r = numberOfPermutation, thread = numberOfThread){
  cat('\nPerforming univariate permutation analysis...\n')
  system('mkdir UnivariateAnalysis')
  setwd('./UnivariateAnalysis')
  naInspection <- apply(AvgNormDF$AvgDF, 1, function(x)sum(is.na(x)))
  AvgDF <- AvgNormDF$AvgDF[naInspection==0,]
  
  runDF <- data.frame(t(AvgDF), response = AvgSampleInfo$CaseStatus)
  if(ANOVA) {
    ANOVAperm <- UnivariatePermTest(runDF, method = 'ANOVA', filename = paste(Sys.Date(), '_NumSignificantApt_permute_ANOVA', sep = ''), r = r, thread = thread)
  }
  
  if(CancerVsRest){
    cat('\nComparing Cancer vs Rest\n')
    response1 = AvgSampleInfo$CaseStatus
    response1[response1!='Cancer'] <- 0
    response1[response1=='Cancer'] <- 1
    runDF1 <- data.frame(t(AvgDF), response = as.factor(response1))
    ttestPerm1 <- UnivariatePermTest(runDF1, method = 'ttest', filename = paste(Sys.Date(), '_NumSignificantApt_permute_ttest_CancerVsRest', sep = ''), r = r, thread = thread)
  }
  
  if(CancerVsNonCancer){
    cat('\nComparing Cancer vs Non-Cancer\n')
    response2 = AvgSampleInfo$CaseStatus[AvgSampleInfo$CaseStatus!='Normal_Female']
    response2[response2!='Cancer'] <- 0
    response2[response2=='Cancer'] <- 1
    runDF2 <- data.frame(t(AvgDF[,AvgSampleInfo$CaseStatus!='Normal_Female']), response = as.factor(response2))
    ttestPerm2 <- UnivariatePermTest(runDF2, method = 'ttest', filename = paste(Sys.Date(), '_NumSignificantApt_permute_ttest_CancerVsNonCancer', sep = ''), r = r, thread = thread)
  }
  
  if(CancerVsNormal){
    cat('\nComparing Cancer vs Normal\n')
    response3 = AvgSampleInfo$CaseStatus[AvgSampleInfo$CaseStatus!='Non-Cancer']
    response3[response3!='Cancer'] <- 0
    response3[response3=='Cancer'] <- 1
    runDF3 <- data.frame(t(AvgDF[,AvgSampleInfo$CaseStatus!='Non-Cancer']), response = as.factor(response3))
    ttestPerm3 <- UnivariatePermTest(runDF3, method = 'ttest', filename = paste(Sys.Date(), '_NumSignificantApt_permute_ttest_CancerVsNormal', sep = ''), r = r, thread = thread)
  } 
  
  cat('\nUnivariate permutation analysis done!\n')
  #z <- list(ANOVAperm = list(ANOVAperm, contrast = 'Cancer vs Non-Cancer vs Normal'),
  #          ttestPerm1 = list(ttestPerm1, contrast = 'Cancer vs Rest'),
  #          ttestPerm2 = list(ttestPerm2, contrast = 'Cancer vs Non-Cancer'),
  #          ttestPerm3 = list(ttestPerm3, contrast = 'Cancer vs Normal'))
}


#########################################################################################################################################################
# build random forest models with feature filetering and cross-validation
#########################################################################################################################################################
Main_RFModels <- function(AvgNormDF, AvgSampleInfo, CancerVsRest = T, CancerVsNonCancer = T, CancerVsNormal = T){
  cat('\nBuilding random forest models...\n')
  system('mkdir GeneralModels')
  setwd('./GeneralModels')
  naInspection <- apply(AvgNormDF$AvgDF, 1, function(x)sum(is.na(x)))
  AvgDF <- AvgNormDF$AvgDF[naInspection==0,]
  if(CancerVsRest){
    cat('\nComparing Cancer vs Rest\n')
    response1 = AvgSampleInfo$CaseStatus
    response1[response1!='Cancer'] <- 0
    response1[response1=='Cancer'] <- 1
    runDF1 <- data.frame(t(AvgDF), response = as.factor(response1))
    
    cat('fitting random forest model...\n')
    fitRF1 <- randomF_featureFilter(runDF1, ntree = 1000, mtry = NULL, replace = F, featureFilter = T)
    cat('10-fold cross-validation...\n')
    fitCVRF1 <- randomF_cv1_para(runDF1, featureFilter = T, K=10, R=1, ntree=1000, mtry = NULL, replace = F)
    
    cat('plotting ROC curves\n')
    pdf(paste(Sys.Date(), '_CancerVsRest_OOB.pdf', sep = ''), 6, 6)
    fitRFplot(fitRF1$predOut, runDF1)
    dev.off()
    
    pdf(paste(Sys.Date(), '_CancerVsRest_10FoldCV.pdf', sep = ''), 6, 6)
    fitRFplot(fitCVRF1$predOut, runDF1, R=1, mfrow = c(1,1), cex.main = 3, cex.lab = 2)
    dev.off()
  }
  
  if(CancerVsNonCancer){
    cat('\nComparing Cancer vs Non-Cancer\n')
    response2 = AvgSampleInfo$CaseStatus[AvgSampleInfo$CaseStatus!='Normal_Female']
    response2[response2!='Cancer'] <- 0
    response2[response2=='Cancer'] <- 1
    runDF2 <- data.frame(t(AvgDF[,AvgSampleInfo$CaseStatus!='Normal_Female']), response = as.factor(response2))
    
    cat('fitting random forest model...\n')
    fitRF2 <- randomF_featureFilter(runDF2, featureFilter = T, ntree=1000, mtry = NULL, replace = F)
    cat('10-fold cross-validation...\n')
    fitCVRF2 <- randomF_cv1_para(runDF2, featureFilter = T, K=10, R=1, ntree=1000, mtry = NULL, replace = F)
    
    cat('plotting ROC curves\n')
    pdf(paste(Sys.Date(), '_CancerVsNonCancer_OOB.pdf', sep = ''), 6, 6)
    fitRFplot(fitRF2$predOut, runDF2)
    dev.off()
    
    pdf(paste(Sys.Date(), '_CancerVsNonCancer_10FoldCV.pdf', sep = ''), 6, 6)
    fitRFplot(fitCVRF2$predOut, runDF2, R=1, mfrow = c(1,1), cex.main = 3, cex.lab = 2)
    dev.off()
  }
  
  if(CancerVsNormal){
    cat('\nComparing Cancer vs Normal\n')
    response3 = AvgSampleInfo$CaseStatus[AvgSampleInfo$CaseStatus!='Non-Cancer']
    response3[response3!='Cancer'] <- 0
    response3[response3=='Cancer'] <- 1
    runDF3 <- data.frame(t(AvgDF[,AvgSampleInfo$CaseStatus!='Non-Cancer']), response = as.factor(response3))
    
    cat('fitting random forest model...\n')
    fitRF3 <- randomF_featureFilter(runDF3, ntree = 1000, mtry = NULL, replace = F, featureFilter = T)
    cat('10-fold cross-validation...\n')
    fitCVRF3 <- randomF_cv1_para(runDF3, featureFilter = T, K=10, R=1, ntree=1000, mtry = NULL, replace = F)
    
    cat('plotting ROC curves\n')
    pdf(paste(Sys.Date(), '_CancerVsNormal_OOB.pdf', sep = ''), 6, 6)
    fitRFplot(fitRF3$predOut, runDF3)
    dev.off()
    
    pdf(paste(Sys.Date(), '_CancerVsNormal_10FoldCV.pdf', sep = ''), 6, 6)
    fitRFplot(fitCVRF3$predOut, runDF3, R=1, mfrow = c(1,1), cex.main = 3, cex.lab = 2)
    dev.off()
  }
  cat('\nFinished building model!\n')
  
  # z <- list(fitRF1 = list(fitRF1, contrast = 'cancer vs rest', model.type = 'RandomForest OOB'),
  #          fitCVRF1 = list(fitCVRF1, contrast = 'cancer vs rest', model.type = 'RandomForest 10-fold cross-validation'),
  #          fitRF2 = list(fitRF2, contrast = 'cancer vs non-cancer', model.type = 'RandomForest OOB'),
  #          fitCVRF2 = list(fitCVRF2, contrast = 'cancer vs non-cancer', model.type = 'RandomForest 10-fold cross-validation'),
  #          fitRF3 = list(fitRF3, contrast = 'cancer vs normal', model.type = 'RandomForest OOB'),
  #          fitCVRF3 = list(fitCVRF3, contrast = 'cancer vs normal', model.type = 'RandomForest 10-fold cross-validation'))
}


#########################################################################################################################################################
# Perform permutation analysis on random forest models
# three options for predMethod: 'OOB', '10-fold CV' or 'both'
#########################################################################################################################################################
Main_RFPerm <- function(AvgNormDF, AvgSampleInfo, CancerVsRest = T, CancerVsNonCancer = T, CancerVsNormal = T, predMethod = 'OOB', r = 1000, thread = 20){
  cat('\nPerforming permutation analysis on RandomForest AUCs...\n')
  system('mkdir PermutationAUCs')
  setwd('./PermutationAUCs')
  naInspection <- apply(AvgNormDF$AvgDF, 1, function(x)sum(is.na(x)))
  AvgDF <- AvgNormDF$AvgDF[naInspection==0,]
  
  if(CancerVsRest){
    cat('\nComparing Cancer vs Rest\n')
    response1 = AvgSampleInfo$CaseStatus
    response1[response1!='Cancer'] <- 0
    response1[response1=='Cancer'] <- 1
    runDF1 <- data.frame(t(AvgDF), response = as.factor(response1))
    if(predMethod == 'both' | predMethod =='OOB'){
      cat('Permutation analysis of out of bag AUCs\n')
      cat(paste('start time: ', Sys.time(), '\n', sep = ''))
      perm1_OOB <- rf_permute_para(runDF1, featureFilter = T, r = r, thread = thread)
      cat(paste('end time: ', Sys.time(), '\n', sep = ''))
      pvalue1_OOB <- sum(perm1_OOB[,-1]>mean(perm1_OOB[,1]))/(thread*ceiling(r/thread))
      cat(paste('p value = ', pvalue1_OOB, '\n'))
      
      pdf(paste(Sys.Date(), '_OOB_AUC_permute_CancerVsRest.pdf', sep = ''), 10, 6)
      hist(perm1_OOB[,-1], nclass = 100,
           main = 'OOB AUC after 1000 Permutation (Cancer vs Rest)',
           xlab = 'AUC', ylab = 'Frequency',
           cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
      abline(v = mean(perm1_OOB[,1]), col = 'red')
      dev.off()
      write.csv(perm1_OOB, file = paste(Sys.Date(), '_OOB_AUC_permute_CancerVsRest.csv', sep = ''))
      write.csv(pvalue1_OOB, file = paste(Sys.Date(), '_OOB_AUC_permute_CancerVsRest_pvalue.csv', sep = ''))
    }
    if(predMethod == 'both' | predMethod =='10-fold CV'){
      cat('Permutation analysis of 10-fold CV AUCs\n')
      cat(paste('start time: ', Sys.time(), '\n', sep = ''))
      perm1_CV <- rf_permute_cv_para(runDF1, featureFilter = T, r = r, thread = thread)
      cat(paste('end time: ', Sys.time(), '\n', sep = ''))
      pvalue1_CV <- sum(perm1_CV[,-1]>mean(perm1_CV[,1]))/(thread*ceiling(r/thread))
      cat(paste('p value = ', pvalue1_CV, '\n'))
      
      pdf(paste(Sys.Date(), '_10foldCV_AUC_permute_CancerVsRest.pdf', sep = ''), 10, 6)
      hist(perm1_CV[,-1], nclass = 100,
           main = '10-fold cross-validation AUC\nafter 1000 Permutation (Cancer vs Rest)',
           xlab = 'AUC', ylab = 'Frequency',
           cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
      abline(v = mean(perm1_CV[,1]), col = 'red')
      dev.off()
      write.csv(perm1_CV, file = paste(Sys.Date(), '_10foldCV_AUC_permute_CancerVsRest.csv', sep = ''))
      write.csv(pvalue1_CV, file = paste(Sys.Date(), '_10foldCV_AUC_permute_CancerVsRest_pvalue.csv', sep = ''))
    }
  }
  
  if(CancerVsNonCancer){
    cat('\nComparing Cancer vs Non-Cancer\n')
    response2 = AvgSampleInfo$CaseStatus[AvgSampleInfo$CaseStatus!='Normal_Female']
    response2[response2!='Cancer'] <- 0
    response2[response2=='Cancer'] <- 1
    runDF2 <- data.frame(t(AvgDF[,AvgSampleInfo$CaseStatus!='Normal_Female']), response = as.factor(response2))
    if(predMethod == 'both' | predMethod =='OOB'){
      cat('Permutation analysis of out of bag AUCs\n')
      cat(paste('start time: ', Sys.time(), '\n', sep = ''))
      perm2_OOB <- rf_permute_para(runDF2, featureFilter = T, r = r, thread = thread)
      cat(paste('end time: ', Sys.time(), '\n', sep = ''))
      pvalue2_OOB <- sum(perm2_OOB[,-1]>mean(perm2_OOB[,1]))/(thread*ceiling(r/thread))
      cat(paste('p value = ', pvalue2_OOB, '\n'))
      
      pdf(paste(Sys.Date(), '_OOB_AUC_permute_CancerVsNonCancer.pdf', sep = ''), 10, 6)
      hist(perm2_OOB[,-1], nclass = 100,
           main = 'OOB AUC after 1000 Permutation (Cancer vs Non-Cancer)',
           xlab = 'AUC', ylab = 'Frequency',
           cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
      abline(v = mean(perm2_OOB[,1]), col = 'red')
      dev.off()
      write.csv(perm2_OOB, file = paste(Sys.Date(), '_OOB_AUC_permute_CancerVsNonCancer.csv', sep = ''))
      write.csv(pvalue2_OOB, file = paste(Sys.Date(), '_OOB_AUC_permute_CancerVsNonCancer_pvalue.csv', sep = ''))
    }
    if(predMethod == 'both' | predMethod =='10-fold CV'){
      cat('Permutation analysis of 10-fold CV AUCs\n')
      cat(paste('start time: ', Sys.time(), '\n', sep = ''))
      perm2_CV <- rf_permute_cv_para(runDF2, featureFilter = T, r = r, thread = thread)
      cat(paste('end time: ', Sys.time(), '\n', sep = ''))
      pvalue2_CV <- sum(perm2_CV[,-1]>mean(perm2_CV[,1]))/(thread*ceiling(r/thread))
      cat(paste('p value = ', pvalue2_CV, '\n'))
      
      pdf(paste(Sys.Date(), '_10foldCV_AUC_permute_CancerVsNonCancer.pdf', sep = ''), 10, 6)
      hist(perm2_CV[,-1], nclass = 100,
           main = '10-fold cross-validation AUC\nafter 1000 Permutation (Cancer vs Non-Cancer)',
           xlab = 'AUC', ylab = 'Frequency',
           cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
      abline(v = mean(perm2_CV[,1]), col = 'red')
      dev.off()
      write.csv(perm2_CV, file = paste(Sys.Date(), '_10foldCV_AUC_permute_CancerVsNonCancer.csv', sep = ''))
      write.csv(pvalue2_CV, file = paste(Sys.Date(), '_10foldCV_AUC_permute_CancerVsNonCancer_pvalue.csv', sep = ''))
    }
  }
  
  if(CancerVsNormal){
    cat('\nComparing Cancer vs Normal\n')
    response3 = AvgSampleInfo$CaseStatus[AvgSampleInfo$CaseStatus!='Non-Cancer']
    response3[response3!='Cancer'] <- 0
    response3[response3=='Cancer'] <- 1
    runDF3 <- data.frame(t(AvgDF[,AvgSampleInfo$CaseStatus!='Non-Cancer']), response = as.factor(response3))
    if(predMethod == 'both' | predMethod =='OOB'){
      cat('Permutation analysis of out of bag AUCs\n')
      cat(paste('start time: ', Sys.time(), '\n', sep = ''))
      perm3_OOB <- rf_permute_para(runDF3, featureFilter = T, r = r, thread = thread)
      cat(paste('end time: ', Sys.time(), '\n', sep = ''))
      pvalue3_OOB <- sum(perm3_OOB[,-1]>mean(perm3_OOB[,1]))/(thread*ceiling(r/thread))
      cat(paste('p value = ', pvalue3_OOB, '\n'))
      
      pdf(paste(Sys.Date(), '_OOB_AUC_permute_CancerVsNormal.pdf', sep = ''), 10, 6)
      hist(perm3_OOB[,-1], nclass = 100,
           main = 'OOB AUC after 1000 Permutation (Cancer vs Normal)',
           xlab = 'AUC', ylab = 'Frequency',
           cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
      abline(v = mean(perm3_OOB[,1]), col = 'red')
      dev.off()
      write.csv(perm3_OOB, file = paste(Sys.Date(), '_OOB_AUC_permute_CancerVsNormal.csv', sep = ''))
      write.csv(pvalue3_OOB, file = paste(Sys.Date(), '_OOB_AUC_permute_CancerVsNormal_pvalue.csv', sep = ''))
    }
    if(predMethod == 'both' | predMethod =='10-fold CV'){
      cat('Permutation analysis of 10-fold CV AUCs\n')
      cat(paste('start time: ', Sys.time(), '\n', sep = ''))
      perm3_CV <- rf_permute_cv_para(runDF3, featureFilter =T, r = r, thread = thread)
      cat(paste('end time: ', Sys.time(), '\n', sep = ''))
      pvalue3_CV <- sum(perm3_CV[,-1]>mean(perm3_CV[,1]))/(thread*ceiling(r/thread))
      cat(paste('p value = ', pvalue3_CV, '\n'))
      
      pdf(paste(Sys.Date(), '_10foldCV_AUC_permute_CancerVsNormal.pdf', sep = ''), 10, 6)
      hist(perm3_CV[,-1], nclass = 100,
           main = '10-fold cross-validation AUC\nafter 1000 Permutation (Cancer vs Normal)',
           xlab = 'AUC', ylab = 'Frequency',
           cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
      abline(v = mean(perm3_CV[,1]), col = 'red')
      dev.off()
      write.csv(perm3_CV, file = paste(Sys.Date(), '_10foldCV_AUC_permute_CancerVsNormal.csv', sep = ''))
      write.csv(pvalue3_CV, file = paste(Sys.Date(), '_10foldCV_AUC_permute_CancerVsNormal_pvalue.csv', sep = ''))
    }
  }
  cat('\nPermutation analysis of AUC are done!\n')
}

#########################################################################################################################################################
# use: Rscript MainProcess_DataAnalysis_v1.R XXX.RData [numberOfPermutation] [numberOfThread]
# RData file should consist of:
#   1. AvgDF, rows are features, columns are instances
#   2. AvgSampleInfo, that contains CaseStatus column
# This script depends on the following script in the same working directory
#   1. 'DataAnalysisModule_v1.0.R'
#########################################################################################################################################################

#########################################################################################################################################################
# Main Processing
#########################################################################################################################################################
#source('~/Documents/Project1_ADAPT/894Probing/DateAnalysisPipelineDev_Model/randomForest_models_v2.4.R')
#source('~/Documents/Project1_ADAPT/894Probing/DateAnalysisPipelineDev_Model/permutation_v2.1.R')
workingDIR <- getwd()
source(paste(workingDIR, '/DataAnalysisModule_v1.0.R', sep = ''))

inputArguments <- commandArgs(trailingOnly = T)
RDataFile <- inputArguments[1]
numberOfPermutation <- as.numeric(inputArguments[2])
numberOfThread <- as.numeric(inputArguments[3])
#numberOfPermutation = 100
#numberOfThread = 20
load(RDataFile)
cat('\nloading datafiles...\n')

setwd(workingDIR)
Process1 <- Main_UnivariatePerm(AvgNormDF, AvgSampleInfo, ANOVA = T, CancerVsRest = T, CancerVsNonCancer = T, CancerVsNormal = T, r = numberOfPermutation, thread = numberOfThread)
setwd(workingDIR)
Process2 <- Main_RFModels(AvgNormDF, AvgSampleInfo, CancerVsRest = T, CancerVsNonCancer = T, CancerVsNormal = T)
setwd(workingDIR)
Process3 <- Main_RFPerm(AvgNormDF, AvgSampleInfo, CancerVsRest = T, CancerVsNonCancer = T, CancerVsNormal = T, predMethod = 'both', r = numberOfPermutation, thread = numberOfThread)













































