#source('~/R_modules/clinical_data/clinical_data_pipeline.R')
source('~/R_modules/ADAPT_tissue/ADAPT_tissue.R')
library(pROC)
library(dplyr)
library(caret)
library(cvAUC)
library(doParallel)
library(foreach)

##########################################################################################
##########################################################################################
preprocess <- function(
  df_train
){
  n <- ncol(df_train)
  Group <- factor(df_train[,n], levels = c('NR', 'R'))
  df_train_lr <- df_train[, -n]
  scale_train <- Scale(log(1 + df_train_lr), scale = T, center = T)
  df_train_lr <- data.frame(
    scale_train$df_scaled,
    Group = df_train$Group
  )
  rownames(df_train_lr) <- NULL
  z <- list(
    df_train_lr = df_train_lr,
    scale_train = scale_train
  )
}

##########################################################################################
##########################################################################################
coreCV <- function(
  df_train,
  method = 'cv', K = 10
){
  n <- ncol(df_train)
  tc <- trainControl(
    method = method, number = K, repeats = 1, 
    savePredictions = T, classProbs = T, returnData = T,
    seeds = NULL
  )
  cvfit <- train(
    Group ~., 
    data = df_train,
    method = 'glm',
    family = binomial,
    trControl = tc
  )
  z <- cvfit
}

##########################################################################################
##########################################################################################
runCV <- function(
  df_train,
  method = 'cv',
  K = 10, 
  R = 100,
  Parallel = F,
  nThread = 4
){
  cvfits <- list()
  if(method == 'LOOCV'){
    cat("method is set to 'LOOCV', no repeat needed.\n")
    cvfit <- coreCV(df_train, method = method, K = K)
    cvfits[[1]] <- cvfit
  }else{
    if(Parallel){
      R <- ceiling(R/nThread)
      cl <- makeCluster(nThread)
      registerDoParallel(cl)
      cvfits <- foreach(
        i = 1:nThread,
        .combine = 'c',
        .export = c('R', 'coreCV'),
        .packages = c('caret')
      ) %dopar% {
        cvfitsInLoop <- list()
        for(i in 1:R){
          cvfit <- coreCV(df_train, method = method, K = K)
          cvfitsInLoop[[i]] <- cvfit
        }
        cvfitsInLoop
      }
      stopCluster(cl)
      registerDoSEQ()
    }else{
      for(i in 1:R){
        cvfit <- coreCV(df_train, method = method, K = K)
        cvfits[[i]] <- cvfit
      }
    }
  }
  cvfits
}

##########################################################################################
##########################################################################################
corePermuteCV <- function(
  df_train, 
  method = 'cv', K = 10,
  nPerm = 100
){
  n <- ncol(df_train)
  y <- df_train[, n]

  predList <- list()
  for(i in 1:nPerm){
    # cat(i, '\n')
    y_perm <- sample(y, replace = F)
    df_train_perm <- data.frame(
      df_train[, -n],
      Group = y_perm
    )
    cvfit <- coreCV(df_train_perm, method = method, K = K)
    # roc(cvfit$pred$obs, cvfit$pred$R, direction = '<', plot = T)
    if(method == 'LOOCV'){
      predList[[i]] <- cvfit$pred[, c(
        'rowIndex', 'obs', 'R'
      )]
    }else{
      predList[[i]] <- cvfit$pred[, c(
        'rowIndex', 'obs', 'R', 'Resample'
      )]
    }
    colnames(predList[[i]])[2] <- 'permLabel' 
  }
  predList
}

##########################################################################################
##########################################################################################
permuteCV <- function(
  df_train,
  method = 'cv', K = 10,
  nPerm = 1000,
  nThread = 4
){

  cl <- makeCluster(nThread)
  registerDoParallel(cl)
  predList <- foreach(
    j = 1:nThread, 
    .combine = 'c',
    .export = c('coreCV', 'corePermuteCV'),
    .packages = c('caret', 'pROC')
  ) %dopar% {
    predByThread <- corePermuteCV(
      df_train, 
      method = method, K = K,
      nPerm = ceiling(nPerm/nThread)
    )
  }
  stopCluster(cl)
  registerDoSEQ()
  predList
}

##########################################################################################
##########################################################################################
# coreBoots <- function(
#   df_train,
#   method = 'cv', K = 10,
#   nBoots = 100
# ){
#   n <- ncol(df_train)
#   y <- df_train[, n]
#   
#   predList <- list()
#   for(i in 1:nPerm){
#     # cat(i, '\n')
#     y_perm <- sample(y, replace = F)
#     df_train_perm <- data.frame(
#       df_train[, -n],
#       Group = y_perm
#     )
#     cvfit <- coreCV(df_train_perm, method = method, K = K)
#     # roc(cvfit$pred$obs, cvfit$pred$R, direction = '<', plot = T)
#     if(method == 'LOOCV'){
#       predList[[i]] <- cvfit$pred[, c(
#         'rowIndex', 'obs', 'R'
#       )]
#     }else{
#       predList[[i]] <- cvfit$pred[, c(
#         'rowIndex', 'obs', 'R', 'Resample'
#       )]
#     }
#     colnames(predList[[i]])[2] <- 'permLabel' 
#   }
#   predList
# }
# 
# 
# ##########################################################################################
# ##########################################################################################
# bootsCV <- function(
#   df_train,
#   method = 'cv', K = 10,
#   nBoots = 1000,
#   nThread = 4
# ){
#   
# }
# 
# 
