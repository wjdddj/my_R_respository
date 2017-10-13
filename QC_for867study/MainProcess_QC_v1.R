####################################################################################################################################
# Main Function to generate QC results
####################################################################################################################################
Main_QC <- function(filename_rawdata, filename_generalStat, filename_sampleTable, plot_data = 'both'){
  seqToExclude <- c('ATCGAGTGTCTCCCGTTTTATCTTTGATTGTTTGC', 'GAGGGCTACTGAAGTTAATGGCATTCTTCTCTATC')

  RawDataRead <- read.csv(filename_rawdata, header = F, stringsAsFactors = F)
  sample_table <- read.csv(filename_sampleTable, header = T, stringsAsFactors = T)
  general_stat <- read.csv(filename_generalStat, header = F, stringsAsFactors = F)
  
  system('mkdir QC_Analysis')
  setwd('./QC_Analysis')
  workingDIR <- getwd()
  RawDataDF <- DF_prep(RawDataRead)
  matchIDX <- match(sample_table$ReplicateID, RawDataDF$sampleInfo$ReplicateID)
  RawDataDF <- list(AllDF = RawDataDF$AllDF[,matchIDX], sampleInfo = RawDataDF$sampleInfo[matchIDX,])
  NormDataDF <- DF_norm(RawDataDF, NormalizationFactor = 1.7e6)
  
  RawDataDF$AllDF <- RawDataDF$AllDF[!rownames(RawDataDF$AllDF) %in% seqToExclude, ]
  SplitDF_byGroup_raw <- DF_split(RawDataDF$AllDF, sample_table, splitBy = sample_table$sampleGroup)
  SplitDF_byBatch_raw <- DF_split(RawDataDF$AllDF, sample_table, splitBy = sample_table$tecanBatch)
  
  # examine total count, 2000 count, ratio
  setwd(workingDIR)
  CountAnlaysis <- CountInspection(RawDataDF, general_stat, sample_table, plot = T, plotBy = sample_table$tecanBatch, cutoffCount = 500000, writeTab = T)
  
  if(plot_data=='both'|plot_data=='raw'){
    cat('\nQC Analysis using RAW data\n')
    # plot input vs samples
    setwd(workingDIR)
    system('mkdir input_scatter_Raw')
    setwd('./input_scatter_Raw')
    inputRaw <- InputVsSample(SplitDF_byBatch_raw, plot = T, plotLimit = 6000)
    plotDNAvsInputCorr(inputRaw, sample_table)
    # plot noPlasma vs samples
    setwd(workingDIR)
    system('mkdir noPlasma_scatter_Raw')
    setwd('./noPlasma_scatter_Raw')
    noPlasmaRaw <- NoPlamsaVsSample(SplitDF_byBatch_raw, plot = T, plotLimit = 6000)
    
  }
  
  NormDataDF$AllDF <- NormDataDF$AllDF[!rownames(NormDataDF$AllDF) %in% seqToExclude, ]
  SplitDF_byGroup_norm <- DF_split(NormDataDF$AllDF, sample_table, splitBy = sample_table$sampleGroup)
  SplitDF_byBatch_norm <- DF_split(NormDataDF$AllDF, sample_table, splitBy = sample_table$tecanBatch)
  if(plot_data=='both'|plot_data=='norm'){
    cat('\nQC Analysis using Normalized data\n')
    # plot input vs samples
    setwd(workingDIR)
    system('mkdir input_scatter_Norm')
    setwd('./input_scatter_Norm')
    inputNorm <- InputVsSample(SplitDF_byBatch_norm, plot = T, plotLimit = 6000)
    plotDNAvsInputCorr(inputNorm, sample_table)
    # plot noPlasma vs samples
    setwd(workingDIR)
    system('mkdir noPlasma_scatter_Norm')
    setwd('./noPlasma_scatter_Norm')
    noPlasmaNorm <- NoPlamsaVsSample(SplitDF_byBatch_norm, plot = T, plotLimit = 6000)
  }
  
  
  if(plot_data=='both'|plot_data=='raw'){
    # make pairwise scatterplots
    setwd(workingDIR)
    system('mkdir pairwise_scatterplot_Raw')
    setwd('./pairwise_scatterplot_Raw')
    pairwiseRaw <- pairwise(SplitDF_byGroup_raw, sample_wise = T, multiCore = T)
  }
  if(plot_data=='both'|plot_data=='norm'){
    # make pairwise scatterplots
    setwd(workingDIR)
    system('mkdir pairwise_scatterplot_Norm')
    setwd('./pairwise_scatterplot_Norm')
    pairwiseNorm <- pairwise(SplitDF_byGroup_norm, sample_wise = T, multiCore = T)
  }
}

####################################################################################################################################
# Main Function to generate Averaged Data, CV Data, CV filtered Data, Imputed Data, merged experimental and clinical covariates
####################################################################################################################################
Main_DataExport <- function(filename_rawdata, filename_sampleTable, filename_covariates, filename_generalStat, 
                            NormalizationFactor = 1.7e6, 
                            removeLowCount = T, cutoffCount = 500000, 
                            removeInputAffected = F, cutoffInput = 0.3){
  
  RawDataRead <- read.csv(filename_rawdata, header = F, stringsAsFactors = F)
  sample_table <- read.csv(filename_sampleTable, header = T, stringsAsFactors = T)
  ClinicalCov <- read.csv(filename_covariates, header = T, stringsAsFactors = F)
  general_stat <- read.csv(filename_generalStat, header = F, stringsAsFactors = F)
  
  system('mkdir DataTables')
  setwd('./DataTables')
  cat('\nMerging covariates tables\n')
  sampleInfo <- merge(sample_table, ClinicalCov, by.x = 'CASE_Barcode', by.y = 'CASE_Barcode')
  
  cat('\nExtracting samples from all data\n')
  RawDataDF <- DF_prep(RawDataRead)
  matchIDX <- match(sample_table$ReplicateID, RawDataDF$sampleInfo$ReplicateID)
  RawDataDF <- list(AllDF = RawDataDF$AllDF[,matchIDX], sampleInfo = RawDataDF$sampleInfo[matchIDX,])
  
  
  # Flag the low count replicates
  # and remove them from the analysis
  if(removeLowCount){
    CountAnalysis <- CountInspection(RawDataDF, general_stat, sample_table, plot = F, plotBy = sample_table$tecanBatch, cutoffCount = cutoffCount, writeTab = F)
    RawDataDF$sampleInfo <- data.frame(RawDataDF$sampleInfo, LessThan500K = CountAnalysis$LessThanCutoff)
    RawDataDF$AllDF <- RawDataDF$AllDF[, !CountAnalysis$LessThanCutoff]
    RawDataDF$sampleInfo <- data.frame(RawDataDF$sampleInfo[!CountAnalysis$LessThanCutoff, ], sample_table[!CountAnalysis$LessThanCutoff, ])
  }
  
  # split data frame by different sampleGroups: Sample, QC, InputLib, NoPlasma, NTC
  SplitDF_byGroup_raw <- DF_split(RawDataDF$AllDF, sample_table = RawDataDF$sampleInfo, splitBy = RawDataDF$sampleInfo$sampleGroup)
  NormDataDF <- DF_norm(RawDataDF, NormalizationFactor = NormalizationFactor)
  SplitDF_byGroup_norm <- DF_split(NormDataDF$AllDF, sample_table = NormDataDF$sampleInfo, splitBy = NormDataDF$sampleInfo$sampleGroup)
  
  matchIDX2 <- match(colnames(SplitDF_byGroup_raw$splitDFs$Sample), sampleInfo$ReplicateID)
  RawDataDF_SampleOnly <- list(AllDF = SplitDF_byGroup_raw$splitDFs$Sample, sampleInfo = sampleInfo[matchIDX2,])
  NormDataDF_SampleOnly <- list(AllDF = SplitDF_byGroup_norm$splitDFs$Sample, sampleInfo = sampleInfo[matchIDX2,])
  RawDataDF_QCs <- list(QC = SplitDF_byGroup_raw$splitDFs$QC, 
                        InputLibA = SplitDF_byGroup_raw$splitDFs$InputLibA, 
                        InputLibB = SplitDF_byGroup_raw$splitDFs$InputLibB,
                        NoPlasma = SplitDF_byGroup_raw$splitDFs$NoPlasma,
                        NTC = SplitDF_byGroup_raw$splitDFs$NTC)
  NormDataDF_QCs <- list(QC = SplitDF_byGroup_norm$splitDFs$QC, 
                         InputLibA = SplitDF_byGroup_norm$splitDFs$InputLibA, 
                         InputLibB = SplitDF_byGroup_norm$splitDFs$InputLibB,
                         NoPlasma = SplitDF_byGroup_norm$splitDFs$NoPlasma,
                         NTC = SplitDF_byGroup_norm$splitDFs$NTC)
  
  # compute correlation to input library only on normalized counts
  # remove sample replicates when their correlation to input is higher than cutoff
  if(removeInputAffected){
    cat('\nAnalyzing correlation to input...')
    InputLib <- apply(SplitDF_byGroup_norm$splitDFs$InputLibA, 1, mean, na.rm = T)
    SampleDF <- NormDataDF_SampleOnly$AllDF
    CorToInput <- cor(InputLib, SampleDF, use = 'pairwise.complete.obs', method = 'spearman')
    
    pdf(paste(Sys.Date(),'_CorrCoefToInput.pdf', sep = ''), 20, 8)
    plot(CorToInput[1,], pch = 16, cex = 0.8, col = as.factor(NormDataDF_SampleOnly$sampleInfo$tecanBatch))
    dev.off()
    
    cat('remove sample replicates with high correlation to input\n')
    InputAboveCutoff <- CorToInput > cutoffInput
    NormDataDF_SampleOnly$AllDF <- NormDataDF_SampleOnly$AllDF[, !InputAboveCutoff]
    NormDataDF_SampleOnly$sampleInfo <- NormDataDF_SampleOnly$sampleInfo[!InputAboveCutoff, ]
  }
  
  cat('\nSaving raw and normalized data\n')
  save(RawDataDF_SampleOnly, file = paste(Sys.Date(), '-DataTable_Raw.RData', sep =''))
  write.csv(RawDataDF_SampleOnly$AllDF, file = paste(Sys.Date(), '-DataTable_Raw.csv', sep =''))
  save(NormDataDF_SampleOnly, file = paste(Sys.Date(), '-DataTable_Norm.RData', sep =''))
  write.csv(NormDataDF_SampleOnly$AllDF, file = paste(Sys.Date(), '-DataTable_Norm.csv', sep =''))
  save(RawDataDF_QCs, file = paste(Sys.Date(), '-QCs_Raw.RData', sep =''))
  write.csv(RawDataDF_QCs, file = paste(Sys.Date(), '-QCs_Raw.csv', sep =''))
  save(NormDataDF_QCs, file = paste(Sys.Date(), '-QCs_Norm.RData', sep =''))
  write.csv(NormDataDF_QCs, file = paste(Sys.Date(), '-QCs_Norm.csv', sep =''))
  
  cat('\nComputing average and CV\n')
  #  AvgRawDF <- DF_average(RawDataDF_SampleOnly$AllDF, RawDataDF_SampleOnly$sampleInfo$SampleID)
  AvgNormDF <- DF_average(NormDataDF_SampleOnly$AllDF, NormDataDF_SampleOnly$sampleInfo$SampleID, removeHighCV = F, CVcutoff = 0.25)
  AvgNormDF_CVfilter <- DF_average(NormDataDF_SampleOnly$AllDF, NormDataDF_SampleOnly$sampleInfo$SampleID, removeHighCV = T, CVcutoff = 0.25)
  matchIDX3 <- match(colnames(AvgNormDF$AvgDF), sampleInfo$SampleID)
  AvgSampleInfo <- sampleInfo[matchIDX3,]
  rmIDX <- match(c('specimen_barcode', 'replicate', 'ReplicateID', 'DNAcon', 'PCRindex', 'plateOrder', 'X'), colnames(AvgSampleInfo))
  AvgSampleInfo <- AvgSampleInfo[,-rmIDX]
  
  cat('\nSaving averaged data and cv table\n')
  save(AvgNormDF, AvgSampleInfo, file = paste(Sys.Date(), '-DataTable_Norm_Average.RData', sep =''))
  save(AvgNormDF_CVfilter, AvgSampleInfo, file = paste(Sys.Date(), '-DataTable_Norm_Average_CVfiltered.RData', sep =''))
  write.csv(AvgNormDF$AvgDF, file = paste(Sys.Date(), '-DataTable_Norm_Average.csv', sep =''))
  write.csv(AvgNormDF_CVfilter$AvgDF, file = paste(Sys.Date(), '-DataTable_Norm_Average_CVfiltered.csv', sep =''))
  write.csv(AvgNormDF$CVDF, file = paste(Sys.Date(), '-DataTable_Norm_CV.csv', sep =''))
  
  z <- list(AvgNormDF = AvgNormDF, AvgSampleInfo = AvgSampleInfo)
}

####################################################################################################################################
# Main Function to analyzing intersample CV, replicate CV
####################################################################################################################################
Main_CoefVariance <- function(Main_DataExport){
  setwd(workingDIR)
  system('mkdir CV_analysis')
  setwd('./CV_analysis')
  cat('plotting CV distributions...\n')
  # plot CV distributions
  CV_sample <- apply(Main_DataExport$AvgNormDF$AvgDF, 1, function(x){
    sd(x, na.rm = T)/mean(x, na.rm = T)
  })
  pdf(paste(Sys.Date(), 'CV_distributions_Rep&Intersample.pdf', sep = ''), 6, 8)
  par(mfrow = c(2,1))
  CVDF <- Main_DataExport$AvgNormDF$CVDF
  CVDF <- CVDF[!is.na(CVDF)]
  hist(CVDF, breaks = seq(0, 100, 0.01), xlim = c(0,1), main = 'technical CV distribution')
  hist(CV_sample, breaks = seq(0, 100, 0.01), xlim = c(0,1), main = 'intersample CV distribution')
  dev.off()
  # plot intersample CV distribution by covariates
  SplitDF_byFC_avg <- DF_split(Main_DataExport$AvgNormDF$AvgDF, Main_DataExport$AvgSampleInfo, splitBy = Main_DataExport$AvgSampleInfo$FlowCell)
  SplitDF_byDay_avg <- DF_split(Main_DataExport$AvgNormDF$AvgDF, Main_DataExport$AvgSampleInfo, splitBy = Main_DataExport$AvgSampleInfo$ExperimentDate)
  SplitDF_byBatch_avg <- DF_split(Main_DataExport$AvgNormDF$AvgDF, Main_DataExport$AvgSampleInfo, splitBy = Main_DataExport$AvgSampleInfo$tecanBatch)
  CV_analysis(SplitDF_byFC_avg, filename = 'byFlowCell')
  CV_analysis(SplitDF_byDay_avg, filename = 'byDay')
  CV_analysis(SplitDF_byBatch_avg, filename = 'byBatch')
  cat('plotting CV distributions...done!\n')
}


####################################################################################################################################
# Main processing
# USE: Rscript QC_functions_v1.1_test.R filename_rawdata filename_generalStat filename_sampleTable filename_covariates
####################################################################################################################################
workingDIR <- getwd()
source(paste(workingDIR, '/datapreparation_v3.1.R', sep = ''))
# source(paste(workingDIR, '/QC_functions_v1.7.R', sep = ''))

inputArguments <- commandArgs(trailingOnly = T)
filename_rawdata <- inputArguments[1]
filename_generalStat <- inputArguments[2]
filename_sampleTable <- inputArguments[3]
filename_covariates <- inputArguments[4]

setwd(workingDIR)
Data <- Main_DataExport(filename_rawdata, filename_sampleTable, filename_covariates, filename_generalStat, removeLowCount = T)
setwd(workingDIR)
Main_CoefVariance(Data)
setwd(workingDIR)
Main_QC(filename_rawdata, filename_generalStat, filename_sampleTable, plot_data = 'both')



