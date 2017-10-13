library(dplyr)

#########################################################################################################
# format score data
#########################################################################################################
formatScoreData <- function(scoreDF){
  
  #scoreDF <- dataset5DF
  scoreDF <- scoreDF %>%
    filter(
      toExclude == 'n'
    ) %>%
    select(
      AccessionNumber, Library, N.H, Sum
    )
  
  scoreDF <- merge(
    data.frame(
      AccessionNumber = scoreDF$AccessionNumber[scoreDF$Library == 'A'],
      libA_N.H = scoreDF$N.H[scoreDF$Library == 'A']
    ),
    data.frame(
      AccessionNumber = scoreDF$AccessionNumber[scoreDF$Library == 'M'],
      libM_Sum = scoreDF$Sum[scoreDF$Library == 'M']
    ),
    by = 'AccessionNumber'
  )
  
  scoreDF <- scoreDF[complete.cases(scoreDF),]
  scoreDF
}



#########################################################################################################
# center and normalize by attributes
# return a Scale object
#########################################################################################################
Scale <- function(dataFrame, center = T, scale = T){
  means <- sapply(dataFrame, mean, na.rm = T)
  stds <- sapply(dataFrame, sd, na.rm = T)
  if(center)
    df_scaled <- sweep(dataFrame, 2, means, '-')
  if(scale)
    df_scaled <- sweep(df_scaled, 2, stds, '/')
  z <- list(df_scaled = df_scaled, 
            means = means, stds = stds, 
            para = list(center = center, scale = scale))
}

#########################################################################################################
# Apply center and normalize to a new data set
# return a Scale object
#########################################################################################################
Scale_apply <- function(dataFrame, Scale){
  center <- Scale$para$center
  scale <- Scale$para$scale
  means <- Scale$means
  stds <- Scale$stds
  if(center)
    df_scaled <- sweep(dataFrame, 2, means, '-')
  if(scale)
    df_scaled <- sweep(df_scaled, 2, stds, '/')
  z <- list(df_scaled = df_scaled, 
            means = means, stds = stds, 
            para = list(center = center, scale = scale))
}

#########################################################################################################
#########################################################################################################
grid_plot_ggplot2 <- function(ggplotList, ncol, nrow){
  library(gridExtra)
  library(methods)
  plotList <- c(ggplotList, nrow = nrow, ncol = ncol)
  do.call(grid.arrange, plotList)
}

#########################################################################################################
#########################################################################################################
plot_jitter_AHC <- function(df, y_column, x_column, color, title){
  library(ggplot2)
  ggplot(df, aes(x = get(x_column), y = get(y_column), color = get(color))) +
    geom_jitter(position=position_jitter(0.2), size = 3) + 
    labs(x = x_column, y = y_column, title = title) + 
    theme(title = element_text(size = 20), 
          #legend.position="none", 
          legend.title = element_blank(),
          axis.title = element_text(size = 20), 
          axis.text = element_text(size = 20),
          line = element_line(linetype = 'solid', size = 4),
          panel.grid = element_line(size = 1),
          panel.background = element_rect(fill = "white"))
}

#########################################################################################################
#########################################################################################################
plot_scatter_AHC <- function(df, y_column, x_column, color, title){
  library(ggplot2)
  ggplot(df, aes(x = get(x_column), y = get(y_column), color = get(color))) +
    geom_point(size = 3) + 
    labs(x = x_column, y = y_column, title = title) + 
    theme(title = element_text(size = 20), 
          #legend.position="none", 
          legend.title = element_blank(),
          axis.title = element_text(size = 20), 
          axis.text = element_text(size = 20),
          line = element_line(linetype = 'solid', size = 4),
          panel.grid = element_line(size = 1),
          panel.background = element_rect(fill = "white"))
}

#########################################################################################################
# Obtain Average data
# AllDF must be a dataframe, replicateLabels must match the columns of AllDF
#########################################################################################################
DF_average <- function(
  AllDF, replicateLabels, .parallel = F,
  removeHighCV = F, CVcutoff = 0.2, 
  save.RData = F, filename
){
  #AllDF <- NormDataDF_SampleOnly$AllDF
  #replicateLabels <- NormDataDF_SampleOnly$sampleInfo$SampleID
  
  library(plyr)
  uniqID <- unique(replicateLabels)
  AllDF_split <- DF_split(AllDF, splitBy = replicateLabels)
  n_sample <- length(AllDF_split$splitDFs)
  
  if(.parallel){
    library(doParallel)
    cl <- makeCluster(8)
    registerDoParallel(cl)
    #clusterExport(cl = cl, varlist = list("AllDF_split"), envir = environment())
    avg_out <- llply(1:n_sample, .parallel = T, function(i){
      replDF <- AllDF_split$splitDFs[[i]]
      if(class(replDF) != 'data.frame'){
        cat(paste('Not enough replicate for ', uniqID[i], '. CV will not be computed.\n', sep = ''))
        avg <- replDF
        cv <- NA
      }else{
        avg <- apply(replDF, 1, mean, na.rm = T)
        std <- apply(replDF, 1, sd, na.rm = T)
        cv <- std/abs(avg)
      }
      z <- list(avg = avg, cv = cv)
    })
    stopCluster(cl)
  }
  else{
    avg_out <- llply(1:n_sample, function(i){
      replDF <- AllDF_split$splitDFs[[i]]
      if(class(replDF) != 'data.frame'){
        cat(paste('Not enough replicate for ', uniqID[i], '. CV will not be computed.\n', sep = ''))
        avg <- replDF
        cv <- NA
      }else{
        avg <- apply(replDF, 1, mean, na.rm = T)
        std <- apply(replDF, 1, sd, na.rm = T)
        cv <- std/abs(avg)
      }
      z <- list(avg = avg, cv = cv)
    })
  }
  AvgDF <- llply(avg_out, function(x)x$avg)
  AvgDF <- data.frame(do.call(cbind, AvgDF))
  CVDF <- llply(avg_out, function(x)x$cv)
  CVDF <- data.frame(do.call(cbind, CVDF))
  colnames(AvgDF) <- names(AllDF_split$splitDFs)
  colnames(CVDF) <- names(AllDF_split$splitDFs)
  
  if(removeHighCV){
    for(j in 1:nrow(AvgDF)){
      AvgDF[j, which(CVDF[j,] > CVcutoff)] <- NA
    }
  }
  datainfo_avg <- list(AvgDF=AvgDF, CVDF = CVDF, type = 'average')
  if(save.RData){
    save(datainfo_avg, file = filename)
  }
  
  z <- datainfo_avg
}

#########################################################################################################
# split data frame based on grouping factor (e.g. sampleGroup, tecanBatch, Fraction, or FlowCell etc)
# input DataFrame could be DF_prep()$AllDF, DF_combine()$AllDF, DF_average()$AvgDF
# input sample_table should be the description covariates of the input DataFrame with the same order
#########################################################################################################
DF_split <- function(DataFrame, sample_table = NULL, splitBy = NULL){
  FactorTab <- table(splitBy)
  n <- length(FactorTab)
  FactorNames <- names(FactorTab)
  splitDFs <- list()
  splitInfos <- list()
  for(i in 1:n){
    idxModule <- which(splitBy==FactorNames[i])
    splitDFs[[length(splitDFs)+1]] <- DataFrame[,idxModule]
    if(!is.null(sample_table)){
      splitInfos[[length(splitInfos)+1]] <- sample_table[idxModule,]
    }
  }
  names(splitDFs) <- FactorNames
  
  if(!is.null(sample_table)){
    names(splitInfos) <- FactorNames
    z <- list(splitDFs = splitDFs, splitInfos = splitInfos)
  }
  else{
    z <- list(splitDFs = splitDFs)
  }
  z
}

#########################################################################################################
## Add group of NR, R.
## input must have column tnt
#########################################################################################################
addResponseGroup <- function(
  dataFrame, cutoff = 180, timeColumn = 'tnt', censorColumn = 'censored'
){
  if(sum(colnames(dataFrame) == timeColumn) == 0)
    stop("no 'time' column for calculation!")
  
  Group <- ifelse(
    dataFrame[, censorColumn] == 'Y' & dataFrame[, timeColumn] <= cutoff,
    'shortCensored',
    ifelse(
      dataFrame[, timeColumn] > cutoff,
      'R', 'NR'
    )
  )
  
  #   Group <- factor(
  #     ifelse(
  #       dataFrame[, timeColumn] > cutoff,
  #       'R', 'NR'
  #     ),
  #     levels = c('NR', 'R')
  #   )
  
  if(sum(Group == 'shortCensored')>0){
    cat('following cases are short censored, please remove them from analysis.\n')
    print(dataFrame[which(Group == 'shortCensored'), ])
  }
  
  Group <- factor(Group, levels = c('NR', 'R'))
  Group
}

#########################################################################################################
## plot ROC for Tissue ADAPT
#########################################################################################################
myROC <- function(
  response, predictor, direction, main
){
  AUC <- auc(response = response, predictor = predictor, direction = direction)
  AUC <- round(AUC, 3)
  ROC <- roc(
    response = response, predictor = predictor, 
    direction = direction, plot = T, percent = T,
    legacy.axes = T, cex = 2, 
    main = paste0(main, '\nAUC = ', AUC)
  )
  z <- list(ROC = ROC, AUC = AUC)
}

#########################################################################################################
## plot library A score on multiple profile from Austria samples
#########################################################################################################
plotAScoreMultiProfile <- function(
  eachAusMultiDF
){
  eachAusMultiDF <- eachAusMultiDF[eachAusMultiDF$Library == 'A', ]
  barplot(
    eachAusMultiDF$N.H,
    names.arg = paste(
      eachAusMultiDF$SampleID,
      eachAusMultiDF$CollectionDate,
      sep = '\n'
    ),
    main = paste0(
      'Library A\n',
      'PatientID: ',
      unique(eachAusMultiDF$PatientID)
    ),
    ylab = 'libA_N.H',
    col = as.factor(eachAusMultiDF$start_firstline_after_collection)
  )
  par(new = T)
  plot(
    eachAusMultiDF$tnt, type = 'l', 
    pch = 16, xaxt = 'n', yaxt = 'n', 
    xlab = 'n', ylab = 'n',
    col = 'blue',
    cex = 2
  )
  axis(4)
  mtext('tnt', side = 4)
  abline(h = 180, col = 'purple', lty = 3)
}

#########################################################################################################
## plot library M score on multiple profile from Austria samples
#########################################################################################################
plotMScoreMultiProfile <- function(
  eachAusMultiDF
){
  eachAusMultiDF <- eachAusMultiDF[eachAusMultiDF$Library == 'M', ]
  barplot(
    eachAusMultiDF$Sum,
    names.arg = paste(
      eachAusMultiDF$SampleID,
      eachAusMultiDF$CollectionDate,
      sep = '\n'
    ),
    main = paste0(
      'Library M\n',
      'PatientID: ',
      unique(eachAusMultiDF$PatientID)
    ),
    ylab = 'libM_Sum',
    col = as.factor(eachAusMultiDF$start_firstline_after_collection)
  )
  par(new = T)
  plot(
    eachAusMultiDF$tnt, type = 'l', 
    pch = 16, xaxt = 'n', yaxt = 'n', 
    xlab = 'n', ylab = 'n',
    col = 'blue',
    cex = 2
  )
  axis(4)
  mtext('tnt', side = 4)
  abline(h = 180, col = 'purple', lty = 3)
}

#########################################################################################################
## generate steps from min to max with constant interval
#########################################################################################################
genCutoffs <- function(dataVec, steps = 100){
  span <- max(dataVec) - min(dataVec)
  span <- min(dataVec) + span/steps * (0:steps)
}

################################################################################################################
################################################################################################################
getCoxMetrics <- function(df_surv){
  library(survival)
  df_surv <- df_surv[, c('time', 'event', 'testLabel')]
  if(length(table(df_surv$testLabel)) < 2){
    z <- data.frame(
      HR = NA,
      p_value = NA
    )
  }else{
    fitCox <- coxph(Surv(time, event) ~ testLabel, data = df_surv)
    coxSum <- summary(fitCox)
    HR <- round(coxSum$coef[1,2], 3)
    p_value <- round(coxSum$coef[1,5], 3)
    z <- data.frame(
      HR = HR,
      p_value = p_value
    )
  }
  z
}

################################################################################################################
## input dataFrame must have a predictor column
################################################################################################################
scanCoxByCutoffs <- function(
  dataFrame,
  cutoffs
){
  out <- lapply(
    cutoffs,
    function(cutoff){
      df_surv <- data.frame(
        time = dataFrame$tnt,
        event = ifelse(dataFrame$censored == 'Y', 0, 1),
        testLabel = ifelse(
          dataFrame$predictor > cutoff, 
          'Test +', 'Test -'
        )
      )
      z <- getCoxMetrics(df_surv)
    }
  )
  out <- do.call(rbind, out)
  z <- data.frame(
    cutoffs = cutoffs,
    out
  )
}

################################################################################################################
## clean up unknowns
################################################################################################################
convert_unkowns <- function(x) gsub('NA|NULL|(^[Uu]nknown.*$)|(^[Nn]ot.*$)', 'Unknown', x)

################################################################################################################
## clean up stage
## consolidate I, II, and IV; only remove trailing numbers for III
################################################################################################################
convertStage <- function(x, collapseIII = T) {
  x <- gsub('^I[A-D1-9]+$', 'I', x)
  x <- gsub('^II[A-D1-9]+$', 'II', x)
  if(collapseIII){
    x <- gsub('^III[A-D1-9]+$', 'III', x)
  }else if(any(grepl('^III$', x))){
    cat('detected \'III\', set collapseIII to TRUE\n')
    x <- gsub('^III[A-D1-9]+$', 'III', x)
  }else{
    x <- ifelse(grepl('^III[A-D1-9]+$', x), gsub('[1-9]+$', '', x), x)
  }
  x <- gsub('^IV[A-D1-9]+$', 'IV', x)
  x <- convert_unkowns(x)
  x
}

################################################################################################################
## clean up grade
################################################################################################################
convertGrade <- function(x) {
  x <- gsub('\\s*/\\s*', '/', trim(x))
  x <- convert_unkowns(x)
  x
}







