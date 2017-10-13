# source('~/Documents/Project1_ADAPT/Developed_Modules/Implemented/Enrichment_modules_v1/enrichment_analysis_v1.0.R')
# source('~/Documents/Project1_ADAPT/Developed_Modules/Implemented/QC_modules_v1/datapreparation_v3.1.R')
source('enrichment_analysis_v1.0.R')
source('datapreparation_v3.1.R')
##################################################################################################################################################################
# Main function to analyze a series of libraries individually in parallel.
# output:
#   1. export top aptamer sequences from each library
#   2. plot speces vs count distribution for each library at a chosen cutoff
#   3. plot lorenz curve for each libray at various cutoffs
# input_filelist: two columns data.frame, column V1, contains the csv file names from Sting's output, 
#                 column V2, contains the preferred sample names
##################################################################################################################################################################
Main_individual_library <- function(input_filelist, outfilename, 
                                    # export top sequences into .csv files
                                    ToDo_topSequence = T, 
                                    numberOfTopSpecies = 10000,
                                    # plot cumulative count vs cutoffs
                                    ToDo_cumCountAtCutoffs = T,
                                    # plot species vs count histograms, number of rows: grid_panel[1], number of columns: grid_panel[2]
                                    ToDo_speciesVsCount = T,
                                    spvscount_cutoff = 5, scales = list(c(1000, 1000), c(10000, 5000), c(50000, 5000)), 
                                    # plot lorenz curves for individual library
                                    ToDo_lorenzCurve = T,
                                    cutoffs_lorenz = c(5, 10, 20, 50, 100)
                                    ){
  library(doParallel)
  library(foreach)
  library(gridExtra)
  library(methods)
  library(ggplot2)

  nLib <- nrow(input_filelist)
  cl <- makeCluster(nLib)
  registerDoParallel(cl)
  cat('start processing libraries in parallel...\n')
  allLibraryPlots <- foreach(
    i = 1:nLib, .combine = 'c', 
    .export = c('topSeq', 'lorenzCurve', 'cum_count_cutoffs', 'speciesVsCountPlot'), 
    .packages = c('data.table', 'ggplot2')
  ) %dopar% {
    log_file <- paste(Sys.Date(), '_log_', input_filelist$V2[i], '.txt', sep = '')
    sink(log_file, append = T, type = 'output')
    
    z <- list()
    # read individual library
    cat('\nreading files...\n')
    libraryDF <- fread(input_filelist$V1[i], data.table = F, sep = ' ', head = F, select = c(1,2), stringsAsFactors = F)
    names(libraryDF) <- c('count', 'seq')
    libraryDF <- data.frame(libraryDF, stringsAsFactors = F, row.names = NULL)
    cat('finished reading files!\n')

    # export top sequences into .csv files
    if(ToDo_topSequence){
      cat('\nwriting top sequences into .csv files...\n')
      topSequencesPerLibrary <- topSeq(libraryDF, numberOfTopSpecies = numberOfTopSpecies)
      filename = paste(Sys.Date(), '_', input_filelist$V2[i], '_top', numberOfTopSpecies, '.csv', sep = '')
      write.csv(topSequencesPerLibrary, file = filename)
      cat('finished writing top sequences!\n')
    }
    
    # plot cumulative count vs cutoffs
    if(ToDo_cumCountAtCutoffs){
      cat('\nplotting cumulative count vs cutoffs...\n')
      cum_count_plots <- cum_count_cutoffs(libraryDF, main = input_filelist$V2[i], 
                                           cutoffs = c(0, 1, 5, 10, 15, 20, 30, 40, 50, 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000, 10000),
                                           plot_only = T)
      # bar_cum_count <- cum_count_plots$cumCountPlot
      # bar_cum_percent <- cum_count_plots$cumPercentPlot
      z[[length(z)+1]] <- cum_count_plots
      names(z)[length(z)] <- 'cum_count_plots'
      cat('finished plotting cutmulative count vs cutoffs!\n')
    }
    
    # plot species vs count histograms
    if(ToDo_speciesVsCount){
      cat('\nplotting species vs count histograms with various scales...\n')
      libs_cutoff <- libraryDF[libraryDF$count > spvscount_cutoff, ]
      # rm(libraryDF)
      main_1 <- paste(input_filelist$V2[i], 'cutoff', spvscount_cutoff, sep = '_')
      histPlot_species_count <- lapply(1:length(scales), function(y){
        cat('plotting x, y scales', scales[[y]], '...\n')
        x_scale = c(0, scales[[y]][1])
        y_scale = c(0, scales[[y]][2])
        speciesVsCountPlot(libs_cutoff, x_scale = x_scale, y_scale = y_scale, main = main_1)
      })
      names(histPlot_species_count) <- paste('scales_', scales, sep = '')
      
      z[[length(z)+1]] <- histPlot_species_count
      names(z)[length(z)] <- 'histPlot_species_count'
      cat('finished plotting histograms!\n')
    }
    
    # plot lorenz curves for individual library
    if(ToDo_lorenzCurve){
      cat('\nplotting lorenz curves for each library at each cutoff...\n')
      
      lib_cutoffs_lorenz <- lapply(1:length(cutoffs_lorenz), function(y){
        featureIDX <- which(libraryDF$count > cutoffs_lorenz[y])
        list(lib = libraryDF[featureIDX, ], cutoff = cutoffs_lorenz[y])
      })
      names(lib_cutoffs_lorenz) <- paste('cutoff_', cutoffs_lorenz, sep = '')
      rm(libraryDF)
      
      lorenzPlot <- lapply(1:length(lib_cutoffs_lorenz), function(y){
        cat('plotting cutoff', cutoffs_lorenz[y], '...\n')
        main_2 <- paste(input_filelist$V2[i], names(lib_cutoffs_lorenz)[y], sep = '_')
        lorenzCurve(lib_cutoffs_lorenz[[y]]$lib, main = main_2, cutoff = lib_cutoffs_lorenz[[y]]$cutoff, plot_only = T)
        # ggplotGrob(lorenzCurve(lib_cutoffs_lorenz[[y]]$lib, plot = T, main = names(lib_cutoffs_lorenz)[y], cutoff = lib_cutoffs_lorenz[[y]]$cutoff, plot_only = T))
      })
      names(lorenzPlot) <- paste('cutoff_', cutoffs_lorenz, sep = '')
      
      z[[length(z)+1]] <- lorenzPlot
      names(z)[length(z)] <- 'lorenzPlot'
      cat('finished plotting lorenz curves!\n')
    }
    # cat('\nsaving plotting data...\n')
    # save(z, file = paste(Sys.Date(), '_', input_filelist$V2[i], '.RData', sep = ''))
    z <- list(z)
  }
  names(allLibraryPlots) <- input_filelist$V2
  stopCluster(cl)
  # save(allLibraryPlots, file = paste(Sys.Date(), outfilename, 'ggplots.RData', sep = '_'))
  
  if(ToDo_cumCountAtCutoffs){
    cat('\nplotting cumulative count vs cutoffs...\n')
    cumCountPlotList <- lapply(allLibraryPlots, function(x)x$cum_count_plots)
    cumCountPlotList <- unlist(cumCountPlotList, recursive = F)
    png(paste(Sys.Date(), outfilename, 'cumCount_vs_cutoffs.png', sep = '_'), 1500*4, 1000*nLib, res = 200)
    grid_plot_ggplot2(cumCountPlotList, nrow = nLib, ncol = 4)
    dev.off()
    cat('finished!\n')
  }
  
  if(ToDo_speciesVsCount){
    cat('\nplotting species vs count histograms with various scales...\n')
    n_scales = length(scales)
    histPlotList <- lapply(allLibraryPlots, function(x)x$histPlot_species_count)
    histPlotList <- unlist(histPlotList, recursive = F)
    png(paste(Sys.Date(), outfilename, 'species-vs-count.png', sep = '_'), 1500*n_scales, 1000*nLib, res = 200)
    grid_plot_ggplot2(histPlotList, nrow = nLib, ncol = n_scales)
    dev.off()
    cat('finished!\n')
  }
  
  if(ToDo_lorenzCurve){
    cat('\nplotting lorenz curves for each library at each cutoff...\n')
    n_cutoffs_lorenz <- length(cutoffs_lorenz)
    lorenzPlotList <- lapply(allLibraryPlots, function(x)x$lorenzPlot)
    lorenzPlotList <- unlist(lorenzPlotList, recursive = F)
    png(paste(Sys.Date(), outfilename, 'lorenz.png', sep = '_'), 1500*n_cutoffs_lorenz, 1500*nLib, res = 200)
    grid_plot_ggplot2(lorenzPlotList, nrow = nLib, ncol = n_cutoffs_lorenz)
    dev.off()
    cat('finished!\n')
  }
}


##################################################################################################################################################################
# Main function to analyze the kmers from one given library, with the consideration of count as weight
# output: 
#    1. barplots of kmer count for each k at each chosen cutoff
#    2. histograms of kmer count for each k at each chosen cutoff
#    3. top kmers for each k at each chosen cutoff
# input_filelist: one row two columns data.frame, column V1, contains the csv file name from Sting's output, 
#                 column V2, contains the preferred sample name
##################################################################################################################################################################
Main_kmer_analysis <- function(input_filelist, outfilename,  
                               cutoffs_kmer = c(20, 50), k = c(1, 2, 3)){
  library(doParallel)
  library(foreach)
  library(WriteXLS)
  
  cat('\nreading file...\n')
  libraryDF <- read.table(input_filelist$V1, header = F, stringsAsFactors = F)
  names(libraryDF) <- c('count', 'seq')
  libraryDF <- data.frame(libraryDF, stringsAsFactors = F)
  cat('finished reading files!\n')
  
  system(paste('mkdir', input_filelist$V2))
  setwd(input_filelist$V2)
  
  n_k <- length(k)
  
  # preparing libraries with various cutoffs
  lib_cutoffs_kmer <- lapply(1:length(cutoffs_kmer), function(y){
    featureIDX <- which(libraryDF$count > cutoffs_kmer[y])
    list(lib = libraryDF[featureIDX, ], cutoff = cutoffs_kmer[y])
  })
  names(lib_cutoffs_kmer) <- paste('cutoff_', cutoffs_kmer, sep = '')
  # examine and removing empty libraries after cutoff filtering
  empty_lib_idx <- which(sapply(lib_cutoffs_kmer, function(x)nrow(x$lib) == 0))
  lib_cutoffs_kmer[empty_lib_idx] <- NULL
  cat('The following cutoffs yielded empty library: \n', cutoffs_kmer[empty_lib_idx], '\n')
  # update the number of loops
  n_cutoffs_kmer <- length(lib_cutoffs_kmer)
  
  # remove initial library object from memory
  rm(libraryDF)
  
  # compute the plots in parallel by k at each cutoff, store kmer_distribution objects in kmerPlots 
  cat('\ncomputing kmer counts distribution for each library at each cutoff...\n')
  kmerPlots <- list()
  for(cutoff_idx in 1:n_cutoffs_kmer){
    cat('computing cutoff', cutoffs_kmer[cutoff_idx], '...\n')
    cl <- makeCluster(n_k)
    registerDoParallel(cl)
    kmerPlot <- foreach(i = 1:n_k, .combine = 'c', .export = c('kmer_distribution')) %dopar% {
      log_file <- paste(Sys.Date(), '_log_', input_filelist$V2, '_', k[i], 'kmer.txt', sep = '')
      sink(log_file, append = T, type = 'output')
      
      main_title <- paste(input_filelist$V2, '_', names(lib_cutoffs_kmer)[cutoff_idx], '_', k[i], 'mer', sep = '')
      eachPlot <- kmer_distribution(lib_cutoffs_kmer[[cutoff_idx]]$lib, k = k[i], weighted = T, plot = T, main = main_title)
      eachPlot <- list(eachPlot)
    }
    names(kmerPlot) <- paste('cutoff', cutoffs_kmer[cutoff_idx], '_', k, 'mer', sep = '')
    stopCluster(cl)
    kmerPlots[[cutoff_idx]] <- kmerPlot
  }
  names(kmerPlots) <- paste('cutoff', cutoffs_kmer, sep = '_')
  cat('finished computing...\n')
  
  # start plotting and exporting files
  cat('\nplotting kmer counts distribution for each library at each cutoff...\n')
  # extract histgram objects from kmerPlots
  kmerPlotList_hist <- lapply(1:n_cutoffs_kmer, function(x){
    eachCutoffPlots <- kmerPlots[[x]]
    lapply(1:n_k, function(y){
      eachCutoffPlots[[y]]$kmerHistogram
    })
  })
  kmerPlotList_hist <- unlist(kmerPlotList_hist, recursive = F)
  png(paste(Sys.Date(), outfilename, input_filelist$V2, 'kmer_histgram.png', sep = '_'), 1500*n_k, 1000*n_cutoffs_kmer, res = 200)
  grid_plot_ggplot2(kmerPlotList_hist, nrow = n_cutoffs_kmer, ncol = n_k)
  dev.off()
  cat('\nfinished plotting kmer counts distribution!\n')
  
  cat('\nplotting kmer counts barplots for each library at each cutoff...\n')
  # extract barplot objects from kmerPlots
  kmerPlotList_bar <- lapply(1:n_cutoffs_kmer, function(x){
    eachCutoffPlots <- kmerPlots[[x]]
    lapply(1:n_k, function(y){
      eachCutoffPlots[[y]]$kmerBarplot
    })
  })
  kmerPlotList_bar <- unlist(kmerPlotList_bar, recursive = F)
  png(paste(Sys.Date(), outfilename, input_filelist$V2, 'kmer_barplot.png', sep = '_'), 2000*n_k, 1000*n_cutoffs_kmer, res = 200)
  grid_plot_ggplot2(kmerPlotList_bar, nrow = n_cutoffs_kmer, ncol = n_k)
  dev.off()
  cat('\nfinished plotting kmer counts barplots!\n')
  
  cat('\nsaving top kmers from each condition...\n')
  # extract top sequence objects from kmerPlots
  kmerTop <- lapply(1:n_cutoffs_kmer, function(x){
    eachCutoffPlots <- kmerPlots[[x]]
    lapply(1:n_k, function(y){
      eachCutoffPlots[[y]]$kmer_count_top
    })
  })
  kmerTop <- unlist(kmerTop, recursive = F)
  names(kmerTop) <- unlist(lapply(kmerPlots, function(x)names(x)), recursive = F)
  WriteXLS(kmerTop, ExcelFileName = paste(Sys.Date(), outfilename, input_filelist$V2, 'top_kmers.xls', sep = '_'), SheetNames = names(kmerTop))
  
  save(kmerPlots, file = paste(Sys.Date(), outfilename, input_filelist$V2, 'kmerAnalysis.RData', sep = '_'))
}


##################################################################################################################################################################
# Main function to plot pairwised scatterplots:
# output: 
#    1. pairwised scatter plots
# input_data_table_csv: the file name for the data table, output table from AptamerGroupAnalysis_v4.1
##################################################################################################################################################################
Main_scatter_plot <- function(input_data_table_csv, outfilename, cutoff = 5){
  library(GGally)
  library(ggplot2)

  cat('\nreading and processing files...\n')
  dataFrame <- read.csv(input_data_table_csv, header = F, stringsAsFactors = F)
  dataFrame <- DF_prep(dataFrame)
  dataFrame$AllDF <- DF_top(dataFrame$AllDF, cutoff = cutoff)$newDF
  nLib <- ncol(dataFrame$AllDF)
  
  cat('\nplotting...\n')
  pair_scatter <- ggpairs(dataFrame$AllDF, lower = list(continuous = "points"), title = outfilename) + 
    theme(title = element_text(size = 20), axis.title = element_text(size = 20), axis.text = element_text(size = 20),
          panel.grid = element_line(size = 1),
          panel.background = element_rect(fill = "white"))
  png(paste(Sys.Date(), outfilename, 'scatter_plots.png'), 1000*nLib, 1000*nLib, res = 200)
  pair_scatter
  dev.off()
  pair_scatter
}


##################################################################################################################################################################
# Main process
# Usage
#   To plot species vs count, lorenz curve. 
#       Rscript Main_enrich_v1.R lib_dist
#   To perform kmer analysis. 
#       Rscript Main_enrich_v1.R kmer [lib_idx] 
#   To plot pairwise scatter plots. 
#       Rscript Main_enrich_v1.R scatter XX.csv 
##################################################################################################################################################################
inputArguments <- commandArgs(trailingOnly = T)
analysisType <- inputArguments[1]
outfilename <- 'ExpID5315_LibSeq'

outputDIR <- getwd()
if(analysisType == 'lib_dist'){
  setwd(outputDIR)
  input_filelist <- inputArguments[2]
  filelist <- read.csv(input_filelist, header = F, stringsAsFactors = F)
  Process1 <- Main_individual_library(filelist, outfilename = outfilename, 
                                      # export top sequences into .csv files
                                      ToDo_topSequence = T, numberOfTopSpecies = 100000,
                                      # plot cumulative count vs cutoffs
                                      ToDo_cumCountAtCutoffs = T,
                                      # plot species vs count histograms
                                      ToDo_speciesVsCount = T, spvscount_cutoff = 1, scales = list(c(1000, 1000), c(10000, 500), c(50000, 500), c(10000, 5000), c(50000, 5000)), 
                                      # plot lorenz curves for individual library
                                      ToDo_lorenzCurve = T, cutoffs_lorenz = c(5, 10, 20, 50, 100))
}

if(analysisType == 'scatter'){
  setwd(outputDIR)
  filename <- inputArguments[2]
  Process2 <- Main_scatter_plot(filename, outfilename = outfilename, cutoff = 5)
}

if(analysisType == 'kmer'){
  setwd(outputDIR)
  lib_idx <- inputArguments[2]
  filelist <- read.csv('filelist.csv', header = F, stringsAsFactors = F)
  filelist <- filelist[lib_idx, ]
  # Process3 <- Main_kmer_analysis(filelist, outfilename = outfilename, cutoffs_kmer = c(100), k = c(1, 2, 3, 4, 5, 6))
  Process4 <- Main_kmer_analysis(filelist, outfilename = outfilename, cutoffs_kmer = c(5), k = c(1, 2))
}


















