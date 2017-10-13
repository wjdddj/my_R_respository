####################################################################################################################################
# Main Function to generate Averaged Data, CV Data, CV filtered Data, Imputed Data, merged experimental and clinical covariates
####################################################################################################################################
Main_DataExport <- function(
  filename_rawdata = NULL, filename_sampleTable = NULL, filename_covariates = NULL, filename_generalStat = NULL, 
  NormalizationFactor = 1.7e6,
  removeLowCount = T, cutoffCount = 500000,
  removeInputAffected = F, cutoffInput = 0.3
){
  library(data.table)
  
  now <- format(Sys.time(), '%Y-%m-%d-%H-%M-%S')
  sample_table <- read.csv(filename_sampleTable, head = T, stringsAsFactors = F)
  general_stat <- read.csv(filename_generalStat, header = F, stringsAsFactors = F)
  if(!is.null(filename_covariates)){
    clinicalCov <- read.csv(filename_covariates, header = T, stringsAsFactors = F)
    cat('\nMerging covariates tables\n')
    sample_table <- merge(sample_table, clinicalCov, by.x = 'CASE_Barcode', by.y = 'CASE_Barcode')
  }
  RawDataDF <- fread(filename_rawdata, header = F, data.table = F)
  
  workingDIR <- getwd()
  out_folder <- file.path(workingDIR, paste0(now, '_Data'))
  dir.create(out_folder)
  setwd(out_folder)
  
  cat('\nExtracting samples from all data\n')
  RawDataDF <- DF_prep(RawDataDF)
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
  
  matchIDX2 <- match(colnames(SplitDF_byGroup_raw$splitDFs$Sample), sample_table$ReplicateID)
  RawDataDF_SampleOnly <- list(AllDF = SplitDF_byGroup_raw$splitDFs$Sample, sampleInfo = sample_table[matchIDX2,])
  NormDataDF_SampleOnly <- list(AllDF = SplitDF_byGroup_norm$splitDFs$Sample, sampleInfo = sample_table[matchIDX2,])

  
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
  save(RawDataDF_SampleOnly, file = paste(now, '-DataTable_Raw.RData', sep =''))
  write.csv(RawDataDF_SampleOnly$AllDF, file = paste(now, '-DataTable_Raw.csv', sep =''))
  save(NormDataDF_SampleOnly, file = paste(now, '-DataTable_Norm.RData', sep =''))
  write.csv(NormDataDF_SampleOnly$AllDF, file = paste(now, '-DataTable_Norm.csv', sep =''))
  save(RawDataDF_QCs, file = paste(now, '-QCs_Raw.RData', sep =''))
  write.csv(RawDataDF_QCs, file = paste(now, '-QCs_Raw.csv', sep =''))
  save(NormDataDF_QCs, file = paste(now, '-QCs_Norm.RData', sep =''))
  write.csv(NormDataDF_QCs, file = paste(now, '-QCs_Norm.csv', sep =''))
  
  cat('\nComputing average and CV\n')
  #  AvgRawDF <- DF_average(RawDataDF_SampleOnly$AllDF, RawDataDF_SampleOnly$sampleInfo$SampleID)
  AvgNormDF <- DF_average(NormDataDF_SampleOnly$AllDF, NormDataDF_SampleOnly$sampleInfo$SampleID, .parallel = T, removeHighCV = F, CVcutoff = 0.25)
  AvgNormDF_CVfilter <- DF_average(NormDataDF_SampleOnly$AllDF, NormDataDF_SampleOnly$sampleInfo$SampleID, .parallel = T, removeHighCV = T, CVcutoff = 0.25)
  matchIDX3 <- match(colnames(AvgNormDF$AvgDF), sampleInfo$SampleID)
  AvgSampleInfo <- sampleInfo[matchIDX3,]
  rmIDX <- match(c('specimen_barcode', 'replicate', 'ReplicateID', 'DNAcon', 'PCRindex', 'plateOrder', 'X'), colnames(AvgSampleInfo))
  AvgSampleInfo <- AvgSampleInfo[,-rmIDX]
  
  cat('\nSaving averaged data and cv table\n')
  save(AvgNormDF, AvgSampleInfo, file = paste(now, '-DataTable_Norm_Average.RData', sep =''))
  save(AvgNormDF_CVfilter, AvgSampleInfo, file = paste(now, '-DataTable_Norm_Average_CVfiltered.RData', sep =''))
  write.csv(AvgNormDF$AvgDF, file = paste(now, '-DataTable_Norm_Average.csv', sep =''))
  write.csv(AvgNormDF_CVfilter$AvgDF, file = paste(now, '-DataTable_Norm_Average_CVfiltered.csv', sep =''))
  write.csv(AvgNormDF$CVDF, file = paste(now, '-DataTable_Norm_CV.csv', sep =''))
  
  z <- list(AvgNormDF = AvgNormDF, AvgSampleInfo = AvgSampleInfo)
}
