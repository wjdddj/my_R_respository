##################################################################################################################################################################
## Description: Perform library analysis for each individual library
##################################################################################################################################################################
individual_library_analysis <- function(
  file_name, sample_name,
  # export top sequences into .csv files
  ToDo_topSequence = T, 
  numberOfTopSpecies = 10000,
  # plot cumulative count vs cutoffs
  ToDo_cumCountAtCutoffs = T, ToPlot_cumCountAtCutoffs = T,
  # plot species vs count histograms, number of rows: grid_panel[1], number of columns: grid_panel[2]
  ToDo_speciesVsCount = T, ToPlot_speciesVsCount = T,
  spvscount_cutoff = 5, scales = list(c(1000, 1000), c(10000, 500), c(50000, 500), c(10000, 5000), c(50000, 5000)), 
  # plot lorenz curves for individual library
  ToDo_lorenzCurve = T, ToPlot_lorenzCurve = T, 
  cutoffs_lorenz = c(5, 10, 20, 50, 100)
){
  #library(gridExtra)
  #library(methods)
  library(ggplot2)
  library(reshape2)
  
  now <- Sys.Date()
  # read file
  libraryDF <- read_lib(file_name, select = c(1,2))
  
  z <- list()
  # prepare top seq
  if(ToDo_topSequence){
    cat('\nwriting top sequences into .csv files...\n')
    topSequencesPerLibrary <- topSeq(libraryDF, numberOfTopSpecies = numberOfTopSpecies)
    topseq_filename = paste0(now, '_', sample_name, '_top', numberOfTopSpecies, '.csv')
    write.csv(topSequencesPerLibrary, file = topseq_filename)
    cat('finished writing top sequences!\n')
  }
  
  # cum-count cum-species analysis, output table
  if(ToDo_cumCountAtCutoffs){
    cat('\npreparing cumulative count vs cutoffs...\n')
    cum_count_plots <- cum_count_cutoffs(
      libraryDF, main = sample_name, 
      cutoffs = c(0, 1, 5, 10, 15, 20, 30, 40, 50, 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000, 10000),
      plot = T, plot_only = F
    )
    z[[length(z)+1]] <- cum_count_plots$plots
    names(z)[length(z)] <- 'cum_count_plots'
    head(cum_count_plots$cumDF)
    cumDF <- melt(cum_count_plots$cumDF, id.vars = 'cutoffs')
    cumDF$cutoffs <- paste0('>', cumDF$cutoffs)
    cumDF <- data.frame(variable = cumDF$variable, cutoffs = cumDF$cutoffs, value = cumDF$value)
    gs_filename = paste0(now, '_', sample_name, '_generalstats.csv')
    write.csv(cumDF, file = gs_filename, row.names = F)
    z[[length(z)+1]] <- cumDF
    names(z)[length(z)] <- 'cumDF'
    cat('finished preparing cutmulative count vs cutoffs!\n')
    
    if(ToPlot_cumCountAtCutoffs){
      cat('\nplotting cumulative count vs cutoffs...\n')
      gsplot_filename <- paste0(now, '_', sample_name, '_generalstats.png')
      png(gsplot_filename, 1500*4, 1000, res = 200)
      grid_plot_ggplot2(cum_count_plots$plots, nrow = 1, ncol = 4)
      dev.off()
      cat('finished plotting cutmulative count vs cutoffs!\n')
    }
  }
  
  # species vs count
  if(ToDo_speciesVsCount){
    cat('\npreparing species vs count histograms with various scales...\n')
    libs_cutoff <- libraryDF[libraryDF$count > spvscount_cutoff, ]
    # rm(libraryDF)
    n_scales <- length(scales)
    main_1 <- paste(sample_name, 'cutoff', spvscount_cutoff, sep = '_')
    histPlot_species_count <- lapply(1:n_scales, function(y){
      cat('plotting x, y scales', scales[[y]], '...\n')
      x_scale = c(0, scales[[y]][1])
      y_scale = c(0, scales[[y]][2])
      speciesVsCountPlot(libs_cutoff, x_scale = x_scale, y_scale = y_scale, main = main_1)
    })
    names(histPlot_species_count) <- paste('scales_', scales, sep = '')
    z[[length(z)+1]] <- histPlot_species_count
    names(z)[length(z)] <- 'histPlot_species_count'
    cat('finished preparing histograms!\n')
    
    if(ToPlot_speciesVsCount){
      cat('\nplotting species vs count histograms with various scales...\n')
      scplot_filename <- paste0(now, '_', sample_name, '_species_vs_count.png')
      png(scplot_filename, 1500*n_scales, 1000, res = 200)
      grid_plot_ggplot2(histPlot_species_count, nrow = 1, ncol = n_scales)
      dev.off()
      cat('finished plotting histograms!\n')
    }
  }
  
  # lorenz curve
  if(ToDo_lorenzCurve){
    cat('\npreparing lorenz curves at each cutoff...\n')
    n_cutoff_lorenz <- length(cutoffs_lorenz)
    lib_cutoffs_lorenz <- lapply(1:n_cutoff_lorenz, function(y){
      featureIDX <- which(libraryDF$count > cutoffs_lorenz[y])
      list(lib = libraryDF[featureIDX, ], cutoff = cutoffs_lorenz[y])
    })
    names(lib_cutoffs_lorenz) <- paste('cutoff_', cutoffs_lorenz, sep = '')
    rm(libraryDF)
    lib_cutoffs_lorenz
    
    lorenzPlot <- lapply(1:n_cutoff_lorenz, function(y){
      cat('plotting cutoff', cutoffs_lorenz[y], '...\n')
      main_2 <- paste(sample_name, names(lib_cutoffs_lorenz)[y], sep = '_')
      lorenzCurve(lib_cutoffs_lorenz[[y]]$lib, main = main_2, cutoff = lib_cutoffs_lorenz[[y]]$cutoff, plot_only = T)
      # ggplotGrob(lorenzCurve(lib_cutoffs_lorenz[[y]]$lib, plot = T, main = names(lib_cutoffs_lorenz)[y], cutoff = lib_cutoffs_lorenz[[y]]$cutoff, plot_only = T))
    })
    names(lorenzPlot) <- paste('cutoff_', cutoffs_lorenz, sep = '')
    z[[length(z)+1]] <- lorenzPlot
    names(z)[length(z)] <- 'lorenzPlot'
    cat('finished preparing lorenz curves!\n')
    
    if(ToPlot_lorenzCurve){
      cat('\nplotting lorenz curves at each cutoff...\n')
      lcplot_filename <- paste0(now, '_', sample_name, '_lorenz_curve.png')
      png(lcplot_filename, 1500*n_scales, 1500, res = 200)
      grid_plot_ggplot2(lorenzPlot, nrow = 1, ncol = n_scales)
      dev.off()
      cat('finished plotting lorenz curves!\n')
    }
  }
  
  cat('\nsaving Rda file...\n')
  rda_filename <- paste0(now, '_', sample_name, '_robj.Rda' )
  save(z, file = rda_filename)
  cat('finished saving Rda file!\n')
  
  z
}


##################################################################################################################################################################
##################################################################################################################################################################
individual_kmer_analysis <- function(
  file_name, sample_name,
  cutoff_kmer = 5, k = c(1, 2, 3, 4, 5, 6, 7)
){
  library(doParallel)
  library(foreach)
  library(WriteXLS)
  
  now <- Sys.Date()
  n_k <- length(k)
  
  # read file
  libraryDF <- read_lib(file_name, select = c(1,2))
  
  # preparing libraries with various cutoffs
  featureIDX <- which(libraryDF$count > cutoff_kmer)
  if(length(featureIDX) == 0)
    stop('The chosen cutoff yielded empty library! Please choose a lower cutoff!\n')
  libraryDF <- libraryDF[featureIDX, ]

  cl <- makeCluster(n_k)
  registerDoParallel(cl)
  kmerPlot <- foreach(i = 1:n_k, .combine = 'c', .export = c('kmer_distribution')) %dopar% {
    log_file <- paste(now, '_log_', sample_name, '_', k[i], 'kmer.txt', sep = '')
    sink(log_file, append = T, type = 'output')
    
    main_title <- paste0(sample_name, '_cutoff', cutoff_kmer, '_', k[i], 'mer')
    eachPlot <- kmer_distribution(libraryDF, k = k[i], weighted = T, plot = T, main = main_title)
    eachPlot <- list(eachPlot)
  }
  names(kmerPlot) <- paste('cutoff', cutoff_kmer, '_', k, 'mer', sep = '')
  stopCluster(cl)
  cat('finished computing...\n')
  
  # start plotting and exporting files
  cat('\nplotting kmer counts distribution for each library at cutoff', cutoff_kmer, '...\n')
  # extract histgram objects from kmerPlots
  kmerPlotList_hist <- lapply(kmerPlot, function(x) x$kmerHistogram)
  png(paste(now, sample_name, 'cutoff', cutoff_kmer, 'kmer_histgram.png', sep = '_'), 1500*n_k, 1000*1, res = 200)
  grid_plot_ggplot2(kmerPlotList_hist, nrow = 1, ncol = n_k)
  dev.off()
  cat('\nfinished plotting kmer counts distribution!\n')
  
  cat('\nplotting kmer counts barplots for each library at cutoff', cutoff_kmer, '...\n')
  # extract barplot objects from kmerPlots
  kmerPlotList_bar <- lapply(kmerPlot, function(x) x$kmerBarplot)
  png(paste(now, sample_name, 'cutoff', cutoff_kmer, 'kmer_barplot.png', sep = '_'), 2000*n_k, 1000*1, res = 200)
  grid_plot_ggplot2(kmerPlotList_bar, nrow = 1, ncol = n_k)
  dev.off()
  cat('\nfinished plotting kmer counts barplots!\n')
  
  cat('\nsaving top kmers...\n')
  # extract top sequence objects from kmerPlots
  kmerTop <- lapply(kmerPlot, function(x) x$kmer_count_top)
  names(kmerTop) <- names(kmerPlot)
  WriteXLS(kmerTop, ExcelFileName = paste(now, sample_name, 'cutoff', cutoff_kmer, 'top_kmers.xls', sep = '_'), SheetNames = names(kmerTop))
  
  save(kmerPlot, file = paste(now, sample_name, 'cutoff', cutoff_kmer, 'kmerAnalysis.RData', sep = '_'))
}





