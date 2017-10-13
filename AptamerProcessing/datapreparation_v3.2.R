####################################################################################################################################
# Data Manipulation Module
####################################################################################################################################
# read raw data csv files
readFiles <- function(directory){
  setwd(directory)
  filenames <- dir()
  filelist <- lapply(filenames, function(x){z <- read.csv(x, header = F, stringsAsFactors = F)})
}

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
# combine datamatrix and sampleInfo from DF_prep() output
# combine raw count data from all flowcells
####################################################################################################################################
DF_combine <- function(DF_prep1, DF_prep2, sortDF = F, type = 'raw'){
  finalDF <- data.frame(DF_prep1$AllDF, DF_prep2$AllDF)
  finalSampleInfo <- rbind(DF_prep1$sampleInfo, DF_prep2$sampleInfo)
  datainfo_all <- list(AllDF = finalDF, sampleInfo = finalSampleInfo, type = type, class = 'AptamerDF')
  if(sortDF){
    sortIDX <- order(datainfo_all$sampleInfo$SampleType)
    finalDF <- datainfo_all$AllDF[,sortIDX]
    finalSampleInfo <- datainfo_all$sampleInfo[sortIDX,]
    datainfo_all <- list(AllDF = finalDF, sampleInfo = finalSampleInfo, type = type, class = 'AptamerDF')
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
# Data Quality Control Analysis Module
####################################################################################################################################
# split data frame based on grouping factor (e.g. sampleGroup, tecanBatch, Fraction, or FlowCell etc)
# input DataFrame could be DF_prep()$AllDF, DF_combine()$AllDF, DF_average()$AvgDF
# input sample_table should be the description covariates of the input DataFrame with the same order
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
# examine total count, 2000 count, ratio
####################################################################################################################################
CountInspection <- function(RawDataDF, general_stat, sample_table, cutoffCount = 500000, 
                            plot = T, plotBy,
                            writeTab = F){
  cat('\nAnalyzing total count\n')
  totalcount_allSpecies <- as.numeric(general_stat[3:nrow(general_stat),3])
  # cat(totalcount_allSpecies)
  
  names(totalcount_allSpecies) <- paste(general_stat[3:nrow(general_stat),1], general_stat[3:nrow(general_stat),2], sep = '.')
  matchIDX <- match(RawDataDF$sampleInfo$ReplicateID, general_stat[3:nrow(general_stat),1])
  # cat(matchIDX)
  totalcount_allSpecies <- totalcount_allSpecies[matchIDX]
  totalcount_2000 <- apply(RawDataDF$AllDF, 2, sum, na.rm = TRUE)
  # cat(totalcount_2000)
  ratio <- totalcount_2000/totalcount_allSpecies
  # cat(ratio)
  LessThanCutoff <- totalcount_allSpecies < cutoffCount
  countSummary <- data.frame(totalcount_allSpecies = totalcount_allSpecies, 
                             totalcount_2000 = totalcount_2000,
                             FractionOf2000 = ratio,
                             LessThanCutoff = LessThanCutoff)
  
  if(plot|writeTab){
    system('mkdir total_count_inspection_plots')
    setwd('./total_count_inspection_plots')
  }
  
  if(plot){
    pdf(paste(Sys.Date(), '_DotPlot_totalCount_2000Count.pdf', sep = ''), 10, 10)
    par(mfrow = c(3,1))
    plot(totalcount_allSpecies, pch = 16, main = 'total count for all species', col = as.factor(plotBy))
    plot(totalcount_2000, pch = 16, main = 'total count for 2000 aptamers', col = as.factor(plotBy))
    plot(ratio, pch = 16, main = 'Fraction of 2000 aptamers in total', col = as.factor(plotBy))
    dev.off()
    
    pdf(paste(Sys.Date(), '_BoxplotCovariates_totalCount_2000Count.pdf', sep = ''), 10, 35)
    par(mfrow = c(7,2))
    boxplot(log10(countSummary$totalcount_allSpecies)~sample_table$SampleFraction, main = 'totalCount vs Fraction')
    boxplot(log10(countSummary$totalcount_2000)~sample_table$SampleFraction, main = '2000Count vs Fraction')
    boxplot(log10(countSummary$totalcount_allSpecies)~sample_table$tecanBatch, main = 'totalCount vs Tecan Batch (plate)')
    boxplot(log10(countSummary$totalcount_2000)~sample_table$tecanBatch, main = '2000Count vs Tecan Batch (plate)')   
    boxplot(log10(countSummary$totalcount_allSpecies)~sample_table$FlowCell, main = 'totalCount vs Flowcell')
    boxplot(log10(countSummary$totalcount_2000)~sample_table$FlowCell, main = '2000Count vs Flowcell')
    boxplot(log10(countSummary$totalcount_allSpecies)~sample_table$ExperimentDate, main = 'totalCount vs Experiment Date')
    boxplot(log10(countSummary$totalcount_2000)~sample_table$ExperimentDate, main = '2000Count vs Experiment Date')
    boxplot(log10(countSummary$totalcount_allSpecies)~sample_table$PCRindex, main = 'totalCount vs PCRindex')
    boxplot(log10(countSummary$totalcount_2000)~sample_table$PCRindex, main = '2000Count vs PCRindex')
    plot(log10(countSummary$totalcount_allSpecies)~sample_table$DNAcon, pch = 16, main = 'totalCount vs DNA concentration')
    plot(log10(countSummary$totalcount_2000)~sample_table$DNAcon, pch = 16, main = '2000Count vs DNA concentration')
    boxplot(log10(countSummary$totalcount_allSpecies)~sample_table$replicate, main = 'totalCount vs replicate')
    boxplot(log10(countSummary$totalcount_2000)~sample_table$replicate, main = '2000Count vs replicate')
    dev.off()
    
    pdf(paste(Sys.Date(), '_DotPlot_totalCount_2000Count_logScale.pdf', sep = ''), 10, 10)
    par(mfrow = c(3,1))
    plot(log10(totalcount_allSpecies), pch = 16, main = 'total count for all species', col = as.factor(plotBy))
    plot(log10(totalcount_2000), pch = 16, main = 'total count for 2000 aptamers', col = as.factor(plotBy))
    plot(ratio, pch = 16, main = 'Fraction of 2000 aptamers in total', col = as.factor(plotBy))
    dev.off()
    
    pdf(paste(Sys.Date(), '_Histogram_totalCount_2000Count.pdf', sep = ''), 10, 10)
    par(mfrow = c(2,1))
    hist(totalcount_allSpecies, nclass = 30, main = 'total count for all species')
    hist(totalcount_2000, nclass = 30, main = 'total count for 2000 aptamers')
    dev.off()
  }
  
  if(writeTab){
    write.csv(countSummary, file = paste(Sys.Date(), '_totalCount_2000Count.csv', sep = ''))
  }
  
  z <- countSummary
}

####################################################################################################################################
# plot input library vs all samples by plate
# input should be combined dataframe of RawData or NormData by plate
####################################################################################################################################
InputVsSample <- function(SplitDF_byBatch, plot = T, plotLimit = NULL){
  cat('\nmaking scatterplot of input library vs samples\n')
  n <- length(SplitDF_byBatch$splitDFs)
  batch <- names(SplitDF_byBatch$splitDFs)
  CorInputLibA = CorInputLibB <- c()
  for(i in 1:n){
    inputLibA <- SplitDF_byBatch$splitDFs[[i]][,SplitDF_byBatch$splitInfos[[i]]$sampleGroup=='InputLibA']
    inputLibB <- SplitDF_byBatch$splitDFs[[i]][,SplitDF_byBatch$splitInfos[[i]]$sampleGroup=='InputLibB']
    SampleDF <- SplitDF_byBatch$splitDFs[[i]][,SplitDF_byBatch$splitInfos[[i]]$sampleGroup!='InputLibA'&
                                                SplitDF_byBatch$splitInfos[[i]]$sampleGroup!='InputLibB']
    m <- ncol(SampleDF)
    InputLibA.tmp <- cor(inputLibA, SampleDF, use = 'pairwise.complete.obs', method = 'spearman')
    InputLibA.tmp <- data.frame(t(InputLibA.tmp), tecanBatch = rep(i, m))
    InputLibB.tmp <- cor(inputLibB, SampleDF, use = 'pairwise.complete.obs', method = 'spearman')
    InputLibB.tmp <- data.frame(t(InputLibB.tmp), tecanBatch = rep(i, m))
    CorInputLibA <- rbind(CorInputLibA, InputLibA.tmp)
    CorInputLibB <- rbind(CorInputLibB, InputLibB.tmp)
    
    if(plot){
      cat(paste('scatterplot of plate ', i, '\n', sep = ''))
      pdf(paste(Sys.Date(),'_scatterPlot_inputA_','plate_', batch[i], '.pdf', sep = ''), 30, 30)
      par(mfrow = c(10,10))
      for(j in 1:m){
        if(is.null(plotLimit)){
          plot(inputLibA~SampleDF[,j], pch = 16, cex = 0.8, main = paste(colnames(SampleDF)[j], round(InputLibA.tmp[j,1], 3), sep = ':'), ylab = 'InputLibA (12.5pg)')
        }
        else{
          plot(inputLibA~SampleDF[,j], pch = 16, cex = 0.8, main = paste(colnames(SampleDF)[j], round(InputLibA.tmp[j,1], 3), sep = ':'), ylab = 'InputLibA (12.5pg)',
               xlim = c(0, plotLimit), ylim = c(0, plotLimit))
        }
        
      }
      dev.off()
      pdf(paste(Sys.Date(),'_scatterPlot_inputB_','plate_', batch[i], '.pdf', sep = ''), 30, 30)
      par(mfrow = c(10,10))
      for(j in 1:m){
        if(is.null(plotLimit)){
          plot(inputLibB~SampleDF[,j], pch = 16, cex = 0.8, main = paste(colnames(SampleDF)[j], round(InputLibB.tmp[j,1], 3), sep = ':'), ylab = 'InputLibB (1.25pg)')
        }
        else{
          plot(inputLibB~SampleDF[,j], pch = 16, cex = 0.8, main = paste(colnames(SampleDF)[j], round(InputLibB.tmp[j,1], 3), sep = ':'), ylab = 'InputLibB (1.25pg)', 
               xlim = c(0, plotLimit), ylim = c(0, plotLimit))
        }
      }
      dev.off()
    }
    
  }
  colnames(CorInputLibA) <- c('InputLibA','tecanBatch')
  colnames(CorInputLibB) <- c('InputLibB','tecanBatch')
  
  pdf(paste(Sys.Date(),'_CorrCoefToInput.pdf', sep = ''), 20, 8)
  plot(CorInputLibA$InputLibA, pch = 16, cex = 0.8, col = as.factor(CorInputLibA$tecanBatch))
  dev.off()
  
  write.csv(CorInputLibA, file = paste(Sys.Date(),'_scatterPlot_inputA.csv', sep = ''))
  write.csv(CorInputLibB, file = paste(Sys.Date(),'_scatterPlot_inputB.csv', sep = ''))
  z <- list(CorInputLibA = CorInputLibA, CorInputLibB = CorInputLibB)
}

# plot NoPlasma vs all samples by plate
# input should be combined dataframe of RawData or NormData by plate
NoPlamsaVsSample <- function(SplitDF_byBatch, plot = T, plotLimit = NULL){
  cat('\nmaking scatterplot of noPlasma control vs samples\n')
  n <- length(SplitDF_byBatch$splitDFs)
  batch <- names(SplitDF_byBatch$splitDFs)
  CorNoPlasma <- c()
  for(i in 1:n){
    noPlasma <- SplitDF_byBatch$splitDFs[[i]][,SplitDF_byBatch$splitInfos[[i]]$sampleGroup=='NoPlasma']
    SampleDF <- SplitDF_byBatch$splitDFs[[i]][,SplitDF_byBatch$splitInfos[[i]]$sampleGroup!='NoPlasma']
    m <- ncol(SampleDF)
    noPlasma_tmp <- cor(noPlasma, SampleDF, use = 'pairwise.complete.obs', method = 'spearman')
    noPlasma_tmp <- data.frame(t(noPlasma_tmp), tecanBatch = rep(i, m))
    CorNoPlasma <- rbind(CorNoPlasma, noPlasma_tmp)

    if(plot){
      cat(paste('scatterplot of plate ', i, '\n', sep = ''))
      pdf(paste(Sys.Date(),'_scatterPlot_noPlasma_','plate_', batch[i], '.pdf', sep = ''), 30, 30)
      par(mfrow = c(10,10))
      for(j in 1:m){
        if(is.null(plotLimit)){
          plot(noPlasma~SampleDF[,j], pch = 16, cex = 0.8, main = paste(colnames(SampleDF)[j], round(noPlasma_tmp[j,1], 3), sep = ':'), ylab = 'NoPlasma')
        }
        else{
          plot(noPlasma~SampleDF[,j], pch = 16, cex = 0.8, main = paste(colnames(SampleDF)[j], round(noPlasma_tmp[j,1], 3), sep = ':'), ylab = 'NoPlasma',
               xlim = c(0, plotLimit), ylim = c(0, plotLimit))
        }
      }
      dev.off()
    }
  }
  colnames(CorNoPlasma) <- c('NoPlasma','tecanBatch')
  write.csv(CorNoPlasma, file = paste(Sys.Date(),'_scatterPlot_noPlasma.csv', sep = ''))
  z <- list(CorNoPlasma = CorNoPlasma)
}

####################################################################################################################################
# DNA concentration vs correlation to Input
####################################################################################################################################
plotDNAvsInputCorr <- function(InputVsSample, sample_table){
  cat('plotting DNA concentration vs Correlation Coefficients to InputLibrary\n')
  DNAcon <- sample_table$DNAcon[match(rownames(InputVsSample$CorInputLibA), sample_table$ReplicateID)]
  pdf('DNAconVsCorcoeffToInput.pdf', 8, 8)
  plot(DNAcon~InputVsSample$CorInputLibA$InputLibA, pch = 16)
  dev.off()
  z <- data.frame(DNAcon, CorInputLibA = InputVsSample$CorInputLibA$InputLibA)
}


# make pairwise scatterplot for each group of controls (QC, input, NoPlasma, or Samples)
pairwise <- function(SplitDF_byGroup, sample_wise = F, multiCore = TRUE){
  cat('\nmaking pairwise scatterplot\n')
  # plot QC sample
  cat('plotting QC sample\n')
  nQCs <- ncol(SplitDF_byGroup$splitDFs$QC)
  maxLimitQC <- quantile(as.matrix(SplitDF_byGroup$splitDFs$QC), 0.995)
  png(paste(Sys.Date(),'scatterPlot_QC.png', sep = '_'), 200*nQCs, 200*nQCs, res = 200)
  pairs(SplitDF_byGroup$splitDFs$QC, pch = 16, cex = 0.4, xlim = c(0, maxLimitQC), ylim = c(0, maxLimitQC))
  dev.off()
  CorQC <- cor(SplitDF_byGroup$splitDFs$QC, use = 'pairwise.complete.obs', method = 'spearman')
  write.csv(CorQC, file = paste(Sys.Date(),'scatterPlot_QC.csv', sep = '_'))
  
  # plot input library
  cat('plotting input library\n')
  nInputs <- 2*ncol(SplitDF_byGroup$splitDFs$InputLibA)
  inputLib <- data.frame(SplitDF_byGroup$splitDFs$InputLibA, SplitDF_byGroup$splitDFs$InputLibB)
  maxLimitInput <- quantile(as.matrix(inputLib), 0.995)
  png(paste(Sys.Date(),'scatterPlot_input.png', sep = '_'), 200*nInputs, 200*nInputs, res = 200)
  pairs(inputLib, pch = 16, cex = 0.4, xlim = c(0, maxLimitInput), ylim = c(0, maxLimitInput))
  dev.off()
  CorInput <- cor(inputLib, use = 'pairwise.complete.obs', method = 'spearman')
  write.csv(CorInput, file = paste(Sys.Date(),'scatterPlot_input.csv', sep = '_'))
  
  # plot NoPlasmaControl
  cat('plotting NoPlasma Control\n')
  nNoPs <- ncol(SplitDF_byGroup$splitDFs$NoPlasma)
  maxLimitNoP <- quantile(as.matrix(SplitDF_byGroup$splitDFs$NoPlasma), 0.995)
  png(paste(Sys.Date(),'scatterPlot_NoPlasma.png', sep = '_'), 200*nNoPs, 200*nNoPs, res = 200)
  pairs(SplitDF_byGroup$splitDFs$NoPlasma, pch = 16, cex = 0.4, xlim = c(0, maxLimitNoP), ylim = c(0, maxLimitNoP))
  dev.off()
  CorNoPlasma <- cor(SplitDF_byGroup$splitDFs$NoPlasma, use = 'pairwise.complete.obs', method = 'spearman')
  write.csv(CorNoPlasma, file = paste(Sys.Date(),'scatterPlot_NoPlasma.csv', sep = '_'))
  
  z <- list(CorQC = CorQC, CorInput = CorInput, CorNoPlasma = CorNoPlasma)
  
  # plot Sample-Wise
  if(sample_wise){
    cat('plotting all sample replicates\n')
    sampleDFByBatch <- DF_split(SplitDF_byGroup$splitDFs$Sample, SplitDF_byGroup$splitInfos$Sample, splitBy = SplitDF_byGroup$splitInfos$Sample$tecanBatch)
    nBatches <- length(sampleDFByBatch$splitDFs)
    batch <- names(sampleDFByBatch$splitDFs)
    
    if(multiCore){
      library(doParallel)
      library(foreach)
      library(cvTools)
      totalcores <- detectCores()
      numberOfThreads <- min(nBatches, totalcores)
      cl <- makeCluster(numberOfThreads, outfile = '')
      registerDoParallel(cl)
      figuresPerThread <- cvFolds(nBatches, numberOfThreads)
      
      foreach(i = 1:numberOfThreads) %dopar% {
        figureInfo <- figuresPerThread$subsets[figuresPerThread$which == i]
        for(each_figure in figureInfo){
          cat(paste('samplewise scatterplot of plate ', each_figure, '\n', sep = ''))
          maxLimitSample <- quantile(as.matrix(SplitDF_byGroup$splitDFs$Sample), 0.995)
          png(paste(Sys.Date(),'scatterPlot_Sample_plate_', batch[each_figure],'.png', sep = ''), 10000, 10000, res = 200)
          pairs(sampleDFByBatch$splitDFs[[each_figure]], pch = 16, cex = 0.2, xlim = c(0, maxLimitSample), ylim = c(0, maxLimitSample))
          dev.off()
        }
      }
      stopCluster(cl)
    }

    else{
      for(i in 1:nBatches){
        cat(paste('samplewise scatterplot of plate ', batch[i], '\n', sep = ''))
        maxLimitSample <- quantile(as.matrix(SplitDF_byGroup$splitDFs$Sample), 0.995)
        png(paste(Sys.Date(),'scatterPlot_Sample_plate_', batch[i],'.png', sep = ''), 10000, 10000, res = 200)
        pairs(sampleDFByBatch$splitDFs[[i]], pch = 16, cex = 0.2, xlim = c(0, maxLimitSample), ylim = c(0, maxLimitSample))
        dev.off()
      }
    }

    CorSample <- cor(SplitDF_byGroup$splitDFs$Sample, use = 'pairwise.complete.obs', method = 'spearman')
    z <- list(CorQC = CorQC, CorInput = CorInput, CorNoPlasma = CorNoPlasma, CorSample = CorSample)
  }
  z
}

####################################################################################################################################
# intersample cv inspection
####################################################################################################################################
CV_analysis <- function(SplitDF, filename){
  n <- length(SplitDF$splitDFs)
  CV_sample <- c()
  for(i in 1:n){
    CV.tmp <- apply(SplitDF$splitDFs[[i]], 1, function(x){
      sd(x, na.rm = T)/mean(x, na.rm = T)
    })
    CV_sample <- cbind(CV_sample, CV.tmp)
  }
  pdf(paste(Sys.Date(), '_', filename, '.pdf', sep = ''), 4*ceiling(sqrt(n)), 4*ceiling(sqrt(n)))
  par(mfrow = c(ceiling(sqrt(n)), ceiling(sqrt(n))))
  for(i in 1:n){
    hist(CV_sample[,i], breaks = seq(0, 100, 0.01), xlim = c(0,1), main = names(SplitDF$splitDFs)[i])
  }
  dev.off()
}






