##################################################################################################################################################################
## to use:
##   Rscript libAnalysis enSeq_R8Lib-HVNJHBCXX_S00_L00_R1_none.V.44046.SUM.csv R8Lib [lib|kmer] 
##################################################################################################################################################################
rm(list = ls())
#source('~/R_modules/Enrichment_Analysis/enrichment_analysis_v1.0.R')
#source('~/R_modules/QC_ADAPT_Plasma/datapreparation_v3.1.R')
source('enrichment_analysis_v1.0.R')
#source('datapreparation_v3.2.R')
source('individual_library_analysis.R')

inputArguments <- commandArgs(trailingOnly = T)
file_name <- inputArguments[1]
sample_name <- inputArguments[2]
mode <- inputArguments[3]

mainDir <- getwd()
dir.create(file.path(mainDir, sample_name), showWarnings = FALSE)
setwd(file.path(mainDir, sample_name))

if(mode == 'lib'){
  file_name <- file.path(mainDir, file_name)
  lib_analysis <- individual_library_analysis(
    file_name, sample_name,
    # export top sequences into .csv files
    ToDo_topSequence = T, 
    numberOfTopSpecies = 100000,
    # plot cumulative count vs cutoffs
    ToDo_cumCountAtCutoffs = T, ToPlot_cumCountAtCutoffs = T,
    # plot species vs count histograms, number of rows: grid_panel[1], number of columns: grid_panel[2]
    ToDo_speciesVsCount = T, ToPlot_speciesVsCount = T,
    spvscount_cutoff = 1, scales = list(c(1000, 1000), c(10000, 500), c(50000, 500), c(10000, 5000), c(50000, 5000)), 
    # plot lorenz curves for individual library
    ToDo_lorenzCurve = T, ToPlot_lorenzCurve = T, 
    cutoffs_lorenz = c(5, 10, 20, 50, 100)
  )
}

if(mode == 'kmer'){
  file_name <- file.path(mainDir, file_name)
  #cutoff_kmer <- inputArguments[4]
  cutoff_kmer <- 5
  kmer_analysis <- individual_kmer_analysis(
    file_name, sample_name,
    cutoff_kmer = cutoff_kmer, k = c(1, 2, 3)
  )
}





