####################################################################################################################################
# Main Function to generate QC results
####################################################################################################################################
Main_QC <- function(filename_rawdata, filename_generalStat, filename_sampleTable, plot_data = 'both'){
  
  RawDataRead <- read.csv(filename_rawdata, header = F, stringsAsFactors = F)
  sample_table <- read.csv(filename_sampleTable, header = T, stringsAsFactors = T)
  general_stat <- read.csv(filename_generalStat, header = F, stringsAsFactors = F)
  
  system('mkdir QC_Analysis')
  setwd('./QC_Analysis')
  workingDIR <- getwd()
  RawDataDF <- DF_prep(RawDataRead)
  matchIDX <- match(sample_table$ReplicateID, RawDataDF$sampleInfo$ReplicateID)
  RawDataDF <- list(AllDF = RawDataDF$AllDF[,matchIDX], sampleInfo = RawDataDF$sampleInfo[matchIDX,])
  SplitDF_byGroup_raw <- DF_split(RawDataDF$AllDF, sample_table, splitBy = sample_table$sampleGroup)
  SplitDF_byBatch_raw <- DF_split(RawDataDF$AllDF, sample_table, splitBy = sample_table$tecanBatch)
  
  # examine total count, 2000 count, ratio
  setwd(workingDIR)
  CountAnlaysis <- CountInspection(RawDataDF, general_stat, sample_table, plot = T, plotBy = sample_table$tecanBatch, cutoffCount = 500000, writeTab = T)
  
  if(plot_data=='both'|plot_data=='raw'){
    cat('\nQC Analysis using RAW data\n')
    # make pairwise scatterplots
    setwd(workingDIR)
    system('mkdir pairwise_scatterplot_Raw')
    setwd('./pairwise_scatterplot_Raw')
    pairwiseRaw <- pairwise(SplitDF_byGroup_raw, sample_wise = F)
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
  
  NormDataDF <- DF_norm(RawDataDF, NormalizationFactor = NULL)
  SplitDF_byGroup_norm <- DF_split(NormDataDF$AllDF, sample_table, splitBy = sample_table$sampleGroup)
  SplitDF_byBatch_norm <- DF_split(NormDataDF$AllDF, sample_table, splitBy = sample_table$tecanBatch)
  if(plot_data=='both'|plot_data=='norm'){
    cat('\nQC Analysis using Normalized data\n')
    # make pairwise scatterplots
    setwd(workingDIR)
    system('mkdir pairwise_scatterplot_Norm')
    setwd('./pairwise_scatterplot_Norm')
    pairwiseNorm <- pairwise(SplitDF_byGroup_norm, sample_wise = F)
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
}

####################################################################################################################################
# Main Function to generate Averaged Data, CV Data, CV filtered Data, Imputed Data, merged experimental and clinical covariates
####################################################################################################################################
Main_DataExport <- function(filename_rawdata, filename_sampleTable, filename_covariates, filename_generalStat, removeLowCount = T){
  RawDataRead <- fread(filename_rawdata, header = F, data.table = F)
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
  CountAnlaysis <- CountInspection(RawDataDF, general_stat, sample_table, plot = F, plotBy = sample_table$tecanBatch, cutoffCount = 500000, writeTab = F)
  RawDataDF$sampleInfo <- data.frame(RawDataDF$sampleInfo, LessThan500K = CountAnlaysis$LessThanCutoff)
  
  # and remove them from the analysis
  if(removeLowCount){
    RawDataDF$AllDF <- RawDataDF$AllDF[, !CountAnlaysis$LessThanCutoff]
    RawDataDF$sampleInfo <- data.frame(RawDataDF$sampleInfo[!CountAnlaysis$LessThanCutoff, ], sample_table[!CountAnlaysis$LessThanCutoff, ])
  }
  
  # split data frame by different sampleGroups: Sample, QC, InputLib, NoPlasma, NTC
  SplitDF_byGroup_raw <- DF_split(RawDataDF$AllDF, sample_table = RawDataDF$sampleInfo, splitBy = RawDataDF$sampleInfo$sampleGroup)
  NormDataDF <- DF_norm(RawDataDF, NormalizationFactor = 1.6e6)
  SplitDF_byGroup_norm <- DF_split(NormDataDF$AllDF, sample_table = NormDataDF$sampleInfo, splitBy = NormDataDF$sampleInfo$sampleGroup)
  
  matchIDX2 <- match(colnames(SplitDF_byGroup_raw$splitDFs$Sample), sampleInfo$ReplicateID)
  RawDataDF_SampleOnly <- list(AllDF = SplitDF_byGroup_raw$splitDFs$Sample, sampleInfo = sampleInfo[matchIDX2,])
  NormDataDF_SampleOnly <- list(AllDF = SplitDF_byGroup_norm$splitDFs$Sample, sampleInfo = sampleInfo[matchIDX2,])
  RawDataDF_QCs <- data.frame(QC = SplitDF_byGroup_raw$splitDFs$QC, 
                              InputLibA = SplitDF_byGroup_raw$splitDFs$InputLibA, 
                              InputLibB = SplitDF_byGroup_raw$splitDFs$InputLibB,
                              NoPlasma = SplitDF_byGroup_raw$splitDFs$NoPlasma,
                              NTC = SplitDF_byGroup_raw$splitDFs$NTC)
  NormDataDF_QCs <- data.frame(QC = SplitDF_byGroup_norm$splitDFs$QC, 
                               InputLibA = SplitDF_byGroup_norm$splitDFs$InputLibA, 
                               InputLibB = SplitDF_byGroup_norm$splitDFs$InputLibB,
                               NoPlasma = SplitDF_byGroup_norm$splitDFs$NoPlasma,
                               NTC = SplitDF_byGroup_norm$splitDFs$NTC)
  
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
  AvgNormDF <- DF_average(NormDataDF_SampleOnly$AllDF, NormDataDF_SampleOnly$sampleInfo$SampleID)
  matchIDX3 <- match(colnames(AvgNormDF$AvgDF), sampleInfo$SampleID)
  AvgSampleInfo <- sampleInfo[matchIDX3,]
  rmIDX <- match(c('specimen_barcode', 'replicate', 'ReplicateID', 'DNAcon', 'PCRindex', 'plateOrder', 'X'), colnames(AvgSampleInfo))
  AvgSampleInfo <- AvgSampleInfo[,-rmIDX]
  
  cat('\nSaving averaged data and cv table\n')
  save(AvgNormDF, AvgSampleInfo, file = paste(Sys.Date(), '-DataTable_Norm_Average.RData', sep =''))
  write.csv(AvgNormDF$AvgDF, file = paste(Sys.Date(), '-DataTable_Norm_Average.csv', sep =''))
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
}

####################################################################################################################################
# Main Data preparation for Analysis
####################################################################################################################################






####################################################################################################################################
# Main processing
# USE: Rscript QC_functions_v1.1_test.R filename_rawdata filename_generalStat filename_sampleTable filename_covariates
####################################################################################################################################
workingDIR <- getwd()
source(paste(workingDIR, '/datapreparation_v3.0.R', sep = ''))
# source(paste(workingDIR, '/QC_functions_v1.7.R', sep = ''))

inputArguments <- commandArgs(trailingOnly = T)
filename_rawdata <- inputArguments[1]
filename_generalStat <- inputArguments[2]
filename_sampleTable <- inputArguments[3]
filename_covariates <- inputArguments[4]

aptamerToExclude <- c('')

workingDIR <- getwd()
Main_QC(filename_rawdata, filename_generalStat, filename_sampleTable, plot_data = 'both')
setwd(workingDIR)
Data <- Main_DataExport(filename_rawdata, filename_sampleTable, filename_covariates, filename_generalStat, removeLowCount = T)
setwd(workingDIR)
Main_CoefVariance(Data)



