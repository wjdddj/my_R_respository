#!/usr/bin/Rscript

rm(list = ls())
suppressWarnings(library(gsubfn))
suppressWarnings(library(data.table))
# library(dplyr)

####################################################################################################################################
# data prep function
# readDF is the output file from AptamerGroupAnalysis_v4 from Jason, done on hpcbioserv
# output an AptamerDF class object, which consists a list of a datamatrix and a sample information table.
####################################################################################################################################
DF_prep <- function(readDF, save.RData = F, filename){
  #readDF <- RawDF_C
  sampleInfo <- data.frame(t(readDF[1:2,2:ncol(readDF)]), stringsAsFactors = F)
  names(sampleInfo) <- c('SampleName', 'ReplicateID')
  AllDF <- readDF[3:nrow(readDF), 2:ncol(readDF)]
  AllDF <- apply(AllDF, 2, as.numeric)
  AllDF <- data.frame(AllDF)
  rownames(AllDF) <- readDF[3:nrow(readDF), 1]
  colnames(AllDF) <- sampleInfo[, 'SampleName']
  datainfo_all <- list(AllDF=AllDF, sampleInfo=sampleInfo, type = 'raw', class = 'AptamerDF')
  if(save.RData){
    save(datainfo_all, file = filename)
  }
  z <- datainfo_all
}

####################################################################################################################################
# Normalize to frequency (normalized counts)
####################################################################################################################################
DF_norm <- function(DF_prep, NormalizationFactor = NULL, save.RData = F, filename){
  colsum <- apply(DF_prep$AllDF, 2, sum, na.rm = TRUE)
  if(is.null(NormalizationFactor)){
    NormalizationFactor <- mean(colsum)
  }
  AllDF_norm <- sweep(DF_prep$AllDF, 2, colsum, '/')
  AllDF_norm <- sweep(AllDF_norm, 2, NormalizationFactor,'*')
  datainfo_norm <- list(AllDF = AllDF_norm, sampleInfo = DF_prep$sampleInfo, type = 'normalized', class = 'AptamerDF')
  if(save.RData){
    save(datainfo_norm, file = filename)
  }
  z <- datainfo_norm
}

####################################################################################################################################
####################################################################################################################################
DF_norm_simple <- function(AllDF, colsum = NULL, NormalizationFactor = NULL){
  if(is.null(colsum)){
    colsum <- apply(AllDF, 2, sum, na.rm = TRUE)
  }
  if(is.null(NormalizationFactor)){
    NormalizationFactor <- mean(colsum)
  }
  AllDF_norm <- sweep(AllDF, 2, colsum, '/')
  AllDF_norm <- sweep(AllDF_norm, 2, NormalizationFactor,'*')
  z <- AllDF_norm
}

####################################################################################################################################
# obtain total valid reads from gs data
####################################################################################################################################
DF_totalValidReads <- function(
  gs_file
){
  gs_table <- read.csv(
    gs_file,
    header = F, stringsAsFactors = F, check.names = F
  )
  gs_table <- data.frame(
    SampleName = paste(
      gs_table[3:nrow(gs_table), 1],
      gs_table[3:nrow(gs_table), 2],
      sep = '-'
    ),
    ReplicateID = gs_table[3:nrow(gs_table), 1],
    TotalValidReads = as.numeric(gs_table[3:nrow(gs_table), 3]),
    stringsAsFactors = F
  )
  rownames(gs_table) <- NULL
  gs_table
}

####################################################################################################################################
# Obtain Average data
# AllDF must be a dataframe, replicateLabels must match the columns of AllDF
####################################################################################################################################
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

####################################################################################################################################
# obtain the highest ranked top(default is 10000) features from each column, and then merge them
####################################################################################################################################
DF_top <- function(dataFrame, method = 'top', top = 10000, cutoff = 20){
  n <- ncol(dataFrame)
  if(method == 'top'){
    topIDX <- apply(dataFrame, 2, function(x){
      order(x, decreasing = T)[1:top]
    })
  }
  
  else if(method == 'cutoff'){
    topIDX <- apply(dataFrame, 2, function(x){
      which(x >= cutoff)
    })  
  }
  
  keptIDX <- unique(as.vector(unlist(topIDX)))
  newDF <- dataFrame[keptIDX, ]
  z <- list(newDF = newDF,
            keptIDX = keptIDX)
  
}


####################################################################################################################################
# split data frame based on grouping factor (e.g. sampleGroup, tecanBatch, Fraction, or FlowCell etc)
# input DataFrame could be DF_prep()$AllDF, DF_combine()$AllDF, DF_average()$AvgDF
# input sample_table should be the description covariates of the input DataFrame with the same order
####################################################################################################################################
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

####################################################################################################################################
## write process log file
####################################################################################################################################
writeProcessLog <- function(out_dir, line){
  log_file <- paste0(out_dir, '/', 'process.log')
  if(!file.exists(log_file))
    file.create(log_file)
  write(line, log_file, append = T)
}

####################################################################################################################################
## help function
####################################################################################################################################
printHelp <- function(){
  cat(paste0(
    'This script is to help process data output from the aptamer pipeline, ', 
    'after running ExtractionInformationFromNGSPrimaryAnalysis and AptamerGroupAnalysisWithSampleTypes_v4.1 ', 
    'or after running aptamer_process.sh [gs|dt] on hpcbioserv/devhpcmgmt01.\n'
  ))
  cat('USAGE:\n')
  cat('AptamerBasic.R [method] arguments\n')
  cat('AptamerBasic.R prep data_dir out_dir [normalizationFactor]. normalizationFactor is optional with default 1e7.\n')
  cat('AptamerBasic.R scatter prep_dt_file st_file out_dir out_file [label_var_name]. label_var_name is optional with default SampleID.\n')
}


####################################################################################################################################
## Prepare data (raw, normalization, average) for each flowcell
####################################################################################################################################
PrepareData <- function(
  data_dir, out_dir, 
  normalizationFactor = 1e7
){
  
  now <- format(Sys.Date(), '%Y%m%d')
  
  dir.create(out_dir, showWarnings = F)
  
  filenames <- dir(data_dir)
  dt_file <- paste0(data_dir, '/', filenames[grep('DataTable', filenames)])
  gs_file <- paste0(data_dir, '/', filenames[grep('GeneralStat', filenames)])
  st_file <- paste0(data_dir, '/', filenames[grep('sample_table|sampletable', filenames)])
  
  ## read files and process data
  RawDF <- read.csv(dt_file, header = F, stringsAsFactors = F, check.names = F)
  RawDF <- DF_prep(RawDF)
  write.csv(
    RawDF$AllDF, file = paste0(out_dir, '/', now, '_raw.csv')
  )
  
  NormDF <- DF_norm(RawDF, NormalizationFactor = normalizationFactor)
  write.csv(
    NormDF$AllDF, file = paste0(out_dir, '/', now, '_norm.csv')
  )
  
  ## read sample table and match rows to data matrix
  stDF <- read.csv(st_file, header = T, stringsAsFactors = F, check.names = F)
  if(
    sum(!RawDF$sampleInfo$ReplicateID %in% stDF$ReplicateID) > 0 | 
    nrow(RawDF$sampleInfo) != nrow(stDF)
  )
    stop('Sample table ReplicateID column must match ReplicateID in the DataTable file. Please check.')
  
  stDF <- stDF[match(RawDF$sampleInfo$ReplicateID, stDF$ReplicateID), ]
  gsDF <- DF_totalValidReads(gs_file)
  stDF <- merge(
    stDF, gsDF, by = 'SampleName'
  )
  
#   stDF$Case_Barcode <- strapply(stDF$SpecimenBarcode, '[0-9]{9}', simplify = T)
  write.csv(
    stDF, file = paste0(out_dir, '/', now, '_sample_table.csv'),
    row.names = F
  )
  
  ## prepare averaged data based on SampleID column
  avg_columns <- c('SampleID', 'CaseType', 'ER', 'PR', 'HER2', 'Library') 
  stAvgDF <- unique(stDF[, which(colnames(stDF) %in% avg_columns)])
  
  if(length(unique(stDF$SampleID)) != nrow(stAvgDF))
    stop('Columns select for averaged dataset are not unique. Please check avg_columns vector.')
  else
    write.csv(
      stAvgDF, file = paste0(out_dir, '/', now, '_sample_table_avg.csv'),
      row.names = F
    )
  
  RawAvgDF <- DF_average(RawDF$AllDF, replicateLabels = stDF$SampleID)
  write.csv(
    RawAvgDF$AvgDF, file = paste0(out_dir, '/', now, '_raw_avg.csv')
  )
  write.csv(
    RawAvgDF$CVDF, file = paste0(out_dir, '/', now, '_raw_replicateCV.csv')
  )
  NormAvgDF <- DF_average(NormDF$AllDF, replicateLabels = stDF$SampleID)
  write.csv(
    NormAvgDF$AvgDF, file = paste0(out_dir, '/', now, '_norm_avg.csv')
  )
  write.csv(
    NormAvgDF$CVDF, file = paste0(out_dir, '/', now, '_norm_replicateCV.csv')
  )
  
}


####################################################################################################################################
## plot scatterplots for each flowcell
####################################################################################################################################
scatterplot <- function(
  prep_dt_file, # output datatable file from PrepareData()
  st_file, # output sampletable file from PrepareData()
  out_dir,
  out_file, 
  label_var_name # variable name in sampletable file to label
){
  now <- format(Sys.Date(), '%Y%m%d')
  dir.create(out_dir, showWarnings = F)
  
  dtDF <- fread(prep_dt_file, header = T, stringsAsFactors = F, data.table = F)
  rownames(dtDF) <- dtDF$V1
  dtDF$V1 <- NULL
  stDF <- fread(st_file, header = T, stringsAsFactors = F, data.table = F)
  rownames(stDF) <- stDF$V1
  stDF$V1 <- NULL
  
  # check dtDF and stDF dimension match
  if(nrow(stDF) != ncol(dtDF))
    stop('Number of rows in sample table file must match number of columns in datatable file. Please check.')
  
  # check label_var_name in sample table file
  if(!label_var_name %in% colnames(stDF))
    stop('label_var_name must be one of the columns in sample table. Please check')
  
  # check label_var_name uniqueness
  if(length(unique(stDF[, label_var_name])) != nrow(stDF))
    stop('label_var_name column must contain unique values. Please check.')
  
  plotDF <- dtDF
  colnames(plotDF) <- stDF[, label_var_name]
  plotDF <- plotDF[, order(colnames(plotDF))]
  n_col <- ncol(plotDF)
  maxlim <- max(sapply(plotDF, function(x) quantile(x, 0.999, na.rm = T)), na.rm = T)
  png(paste0(out_dir, '/', now, '_', out_file, '.png'), n_col*200, n_col*200, res = 100)
  pairs(
    plotDF, pch = 16, cex = 0.5, 
    xlim = c(0, maxlim),
    ylim = c(0, maxlim)
  )
  dev.off()
}

####################################################################################################################################
## main function
####################################################################################################################################
main <- function(inputArgs){
  
  method_options <- c('prep', 'scatter')
  method = inputArgs[1]
  
  ## check method validity
  if(!method %in% method_options){
    printHelp()
    stop('method must be in ', paste(method_options, collapse = '|'), '. Please check.')
  }
  
  # prepare data
  if(method == 'prep'){
    if(is.null(inputArgs[2])){
      printHelp()
      stop('Please input argument to proceed.')
    }
    
    ## to prepare data files
    ## command line run: Rscript AptamerBasic.R prep data_dir out_dir [normalizationFactor]
    data_dir = as.character(inputArgs[2])
    out_dir = as.character(inputArgs[3])
    if(is.na(inputArgs[4]))
      normalizationFactor = 1e7
    else
      normalizationFactor = as.numeric(inputArgs[4])
    
    PrepareData(data_dir, out_dir, normalizationFactor)
    
    writeProcessLog(
      out_dir, 
      paste0(Sys.time(), '. Prepare Data. Normalization Factor: ', normalizationFactor)
    )
  }
  
  # make scatterplot
  if(method == 'scatter'){
    ## to make scatterplots
    ## command line run: Rscript AptamerBasic.R scatter prep_dt_file st_file out_dir [label_var_name]
    prep_dt_file = inputArgs[2]
    st_file = inputArgs[3]
    out_dir = inputArgs[4]
    out_file = inputArgs[5]
    if(is.na(inputArgs[6]))
      label_var_name = 'SampleID'
    else
      label_var_name = inputArgs[6]
    
    scatterplot(prep_dt_file, st_file, out_dir, out_file, label_var_name)
    
    writeProcessLog(
      out_dir, 
      paste0(
        Sys.time(), '. Make scatterplots from ', prep_dt_file, 
        ', with column ', label_var_name, ' from ', st_file
      )
    )
  }
}

#########
inputArgs = commandArgs(trailingOnly = T)
main(inputArgs)



