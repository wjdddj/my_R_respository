##################################################################################################################################################################
## read demultiplexed 'enSeq' or 'top seq' output files using data.table
##################################################################################################################################################################
read_lib <- function(file_name, ...){
  library(data.table)
  cat('reading file', file_name, '...\n')
  df_lib <- fread(file_name, stringsAsFactors = F, data.table = F, ...)
  if(class(df_lib[,1]) %in% c('numeric', 'integer')){
    colnames(df_lib) <- c('count', 'seq')
  }else if(class(df_lib[,2]) %in% c('numeric', 'integer')){
    colnames(df_lib) <- c('seq', 'count')
    df_lib <- data.frame(
      count = df_lib$count,
      seq = df_lib$seq,
      stringsAsFactors = F
    )
  }else{
    stop('Wrong file format, please check! File should have two columns (count and seq)\n')
  }
  cat('finished reading file!\n')
  df_lib
}

##################################################################################################################################################################
##################################################################################################################################################################
meme_norm_to_min_count <- function(df_seq){
  min_count <- min(df_seq$count)
  df_seq$count <- round(df_seq$count/min_count, 0)
  df_seq
}

##################################################################################################################################################################
## input: df_seq a data.frame with two columns, seq and count.
## output: save fasta format
##################################################################################################################################################################
to_fasta <- function(
  df_seq, 
  out_file_name
){
  library(seqinr)
  # df_seq <- df_lib
  # now <- format(Sys.time(), '%Y%m%d-%H%M%S')
  # out_file_name <- paste0(now, '.fa')
  df_seq$seq <- as.character(df_seq$seq)
  
  n <- nrow(df_seq)
  seq <- lapply(1:n, function(i){
    rep(df_seq$seq[i], df_seq$count[i])
  })
  seq <- do.call(c, seq)
  
  write.fasta(as.list(seq), names = 1:length(seq), file = out_file_name)
}

##################################################################################################################################################################
# summary the top numberOfTopSpecies of a give library
# input is a data.frame: first column is count, second column is species id
# output is a vector
##################################################################################################################################################################
topSeq <- function(libDataFrame, numberOfTopSpecies = 10000){
  topIDX <- order(libDataFrame[,1], decreasing = T)[1:numberOfTopSpecies]
  S00Out <- libDataFrame[topIDX, 1]
  names(S00Out) <- libDataFrame[topIDX, 2]
  S00Out
}

##################################################################################################################################################################
# generate a panel of plots from a list of ggplot objects
##################################################################################################################################################################
grid_plot_ggplot2 <- function(ggplotList, ncol, nrow){
  library(gridExtra)
  library(methods)
  plotList <- c(ggplotList, nrow = nrow, ncol = ncol)
  do.call(grid.arrange, plotList)
}

##################################################################################################################################################################
# For a given library, compute cumulative counts, cumulative count percent, number of species, and percent of total species above a series cutoffs.
# input: libDataFrame, a data.frame of an Aptamer library, first column is the count, second column is the sequence
# output: four ggplots object and a data.frame of the computed counts and percentages.
##################################################################################################################################################################
cum_count_cutoffs <- function(
  libDataFrame, main, 
  cutoffs = c(0, 1, 2, 5, 10, 15, 20, 30, 40, 50, 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000, 10000),
  plot = T, plot_only = T
){
  
  library(parallel)
  library(ggplot2)
  NumberOfValSpecies <- mclapply(1:length(cutoffs), function(x){
    sum(libDataFrame[,1] > cutoffs[x])
  }, mc.cores = 4)
  NumberOfValSpecies <- unlist(NumberOfValSpecies)
  total_species <- nrow(libDataFrame)
  
  cum_count <- mclapply(1:length(cutoffs), function(x){
    sum(libDataFrame[,1][libDataFrame[,1] > cutoffs[x]])
  }, mc.cores = 4)
  cum_count <- unlist(cum_count)
  total_count <- sum(libDataFrame[,1])
  
  cumDF <- data.frame(
    cutoffs = cutoffs,
    NumberOfValSpecies = NumberOfValSpecies, 
    Pct_NumberOfValSpecies = NumberOfValSpecies/total_species, 
    NumberOfReadsFromSpecies = cum_count, 
    Pct_NumberOfReadsFromSpecies = cum_count/total_count
  )
  
  if(!plot){
    z <- list(cumDF = cumDF)
  }else{
    cumSpeciesCountPlot <- ggplot(cumDF, aes(y = NumberOfValSpecies, x = factor(cutoffs))) + 
      geom_bar(stat = 'identity') +
      labs(x = 'cutoffs', y = 'number of species above cutoff', title = main) + 
      scale_y_log10() +
      theme(title = element_text(size = 20), axis.title = element_text(size = 20), axis.text = element_text(size = 15),
            axis.text.x = element_text(angle = 90, hjust = 1),
            panel.grid = element_line(size = 1),
            panel.background = element_rect(fill = "white"),
            legend.position='none')
    
    cumSpeciesPercentPlot <- ggplot(cumDF, aes(y = Pct_NumberOfValSpecies, x = factor(cutoffs))) + 
      geom_bar(stat = 'identity') +
      labs(x = 'cutoffs', y = 'percent of species above cutoff', title = main) + 
      theme(title = element_text(size = 20), axis.title = element_text(size = 20), axis.text = element_text(size = 15),
            axis.text.x = element_text(angle = 90, hjust = 1),
            panel.grid = element_line(size = 1),
            panel.background = element_rect(fill = "white"))
    
    cumCountPlot <- ggplot(cumDF, aes(y = NumberOfReadsFromSpecies, x = factor(cutoffs))) + 
      geom_bar(stat = 'identity') +
      labs(x = 'cutoffs', y = 'cumulative count above cutoff', title = main) + 
      # scale_y_log10() +
      theme(title = element_text(size = 20), axis.title = element_text(size = 20), axis.text = element_text(size = 15),
            axis.text.x = element_text(angle = 90, hjust = 1),
            panel.grid = element_line(size = 1),
            panel.background = element_rect(fill = "white"))
    
    cumPercentPlot <- ggplot(cumDF, aes(y = Pct_NumberOfReadsFromSpecies, x = factor(cutoffs))) + 
      geom_bar(stat = 'identity') +
      labs(x = 'cutoffs', y = 'cumulative count percent above cutoff', title = main) + 
      theme(title = element_text(size = 20), axis.title = element_text(size = 15), axis.text = element_text(size = 15),
            axis.text.x = element_text(angle = 90, hjust = 1),
            panel.grid = element_line(size = 1),
            panel.background = element_rect(fill = "white"))
    
    z <- list(plots = list(cumSpeciesCountPlot = cumSpeciesCountPlot, 
                           cumSpeciesPercentPlot = cumSpeciesPercentPlot,
                           cumCountPlot = cumCountPlot, 
                           cumPercentPlot = cumPercentPlot),
              cumDF = cumDF)
    
    if(plot_only){
      z <- list(cumSpeciesCountPlot = cumSpeciesCountPlot, 
                cumSpeciesPercentPlot = cumSpeciesPercentPlot,
                cumCountPlot = cumCountPlot, 
                cumPercentPlot = cumPercentPlot)
    }
  }
  z
}

##################################################################################################################################################################
# For a given library, plot the species ver count histgrams.
# input: 
#       libDataFrame, a data.frame, first column is count, second column is species id
#       x_scale, y_scale, the axis limits. e.g. x_scale = c(0, 10000)  
# output: 
#       a ggplot object of histgram
##################################################################################################################################################################
speciesVsCountPlot <- function(LibDataFrame, x_scale, y_scale, main){
  library(ggplot2)
  names(LibDataFrame) <- c('count', 'seq')
  histPlot_species_count <- ggplot(LibDataFrame, aes(count)) + 
    coord_cartesian(xlim = x_scale, ylim = y_scale) +
    geom_histogram(binwidth = x_scale[2]/100) +
    labs(x = 'count', y = 'number of species', title = main) + 
    theme(title = element_text(size = 20), axis.title = element_text(size = 20), axis.text = element_text(size = 20),
          line = element_line(linetype = 'solid', size = 4),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_line(size = 1),
          panel.background = element_rect(fill = "white"))
  histPlot_species_count
}

##################################################################################################################################################################
# Analysis of the species and count distribution of a given library
# Plot lorenz curve which shows the proportion of overall counts cumulated by the bottom % of species
# compute the cumulative count percentage of all sequences from low count to high count
# compute the cumulative number of species percentage of all sequences from low count to high count
# input arguments: 
#         libDataFrame, a data.frame, first column is count, second column is species id
# output:
#         libDataFrame, the same data.frame with addition of cumulative count percentage and cumulative species percentage
#         GiniCoef: Gini coefficient of the library, the smaller the more skewed distribution of the species over counts, meaning top species acounts for 
#                   more propotion of the total counts.
##################################################################################################################################################################
lorenzCurve <- function(libDataFrame, main = NULL, cutoff = 0, plot_only = T){
  library(ineq)
  library(ggplot2)
  
  libDataFrame <- libDataFrame[libDataFrame[,1] > cutoff,]
  
  if(nrow(libDataFrame) == 0){
    cat('No Data Above this Cutoff!\n')
    cum_speciesPercent <- c(0, 1)
    cum_countPercent <- c(0, 1)
    GiniCoef <- NULL
  }else{
    dataFrameSort <- libDataFrame[order(libDataFrame[,1], decreasing = F), ]
    countVec <- dataFrameSort[,1]
    GiniCoef <- Gini(countVec)
    totalCount <- sum(countVec)
    totalSpecies <- nrow(dataFrameSort)
    
    cum_count <- cumsum(countVec)
    cum_countPercent <- cum_count/totalCount
    cum_speciesPercent <- c(1:totalSpecies)/totalSpecies
  }
  
  if(nrow(libDataFrame) == 0){
    cat('plotting blank...\n')
    p <- qplot(x = cum_speciesPercent, y = cum_countPercent, geom = 'blank', main = paste(main, ' No Data Above Cutoff!', sep = ''),
               xlab = 'cumulative species (percent)', ylab = 'cumulative count (percent)')
    lorenz_curve <- p + 
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = 'longdash') + 
      theme(title = element_text(size = 20), axis.title = element_text(size = 20), axis.text = element_text(size = 20),
            line = element_line(linetype = 'solid', size = 4),
            panel.grid = element_line(size = 1))
  }else{
    p <- qplot(x = cum_speciesPercent, y = cum_countPercent, geom = 'path', main = paste(main, ' Gini =', round(GiniCoef, 3), sep = ''),
               xlab = 'cumulative species (percent)', ylab = 'cumulative count (percent)')
    lorenz_curve <- p + 
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = 'longdash') + 
      theme(title = element_text(size = 20), axis.title = element_text(size = 20), axis.text = element_text(size = 20),
            line = element_line(linetype = 'solid', size = 4),
            panel.grid = element_line(size = 1),
            panel.background = element_rect(fill = "white"))
    libDataFrame <- data.frame(libDataFrame, cum_countPercent = cum_countPercent, cum_speciesPercent = cum_speciesPercent)
  }
  
  if(plot_only){
    z <- lorenz_curve
  }else{
    z <- list(libDataFrame = libDataFrame, GiniCoef = GiniCoef, lorenz_curve = lorenz_curve)
  }
  z
}

##################################################################################################################################################################
# analysis of kmer count and plot its distribution in the library
# input is a vector of the sequences from a library, and k of kmer
# output is a table of all possible kmers and their appearance frequency in the library
##################################################################################################################################################################
kmer_distribution <- function(libDataFrame, k = 4, weighted = T, plot = T, main){
  library(Biostrings)
  library(ggplot2)
  # library(doParallel)
  # library(foreach)
  
  parameters <- list(k = k, weighted = weighted)
  if(!is.vector(libDataFrame)){
    cat('sequence and count information are available...\n')
    DNA_vector <- libDataFrame[,2]
    DNA_count <- as.numeric(libDataFrame[,1])
  }else{
    cat('only sequence information is available...\n')
    DNA_vector <- libDataFrame
  }
  seqList <- sapply(DNA_vector, DNAString)
  
  if(weighted){
    kmer_count <- lapply(1:length(seqList), function(x){
      kmer_count_per_seq <- oligonucleotideFrequency(seqList[[x]], k) * DNA_count[x]
    })
  }else{
    kmer_count <- lapply(1:length(seqList), function(x){
      kmer_count_per_seq <- oligonucleotideFrequency(seqList[[x]], k) * 1
    })
  }
  
  kmer_count <- do.call(cbind, kmer_count)
  kmer_count_all <- data.frame(count = apply(kmer_count, 1, sum),
                               seq = rownames(kmer_count))
  rownames(kmer_count_all) <- NULL
  kmer_count_all_sort <- kmer_count_all[order(kmer_count_all$count, decreasing = T), ]
  kmer_count_top <- kmer_count_all_sort[1:min(100, 4^k), ]
  
  if(plot){
    expected_frequency <- sum(kmer_count_all$count) / 4^k
    
    # histogram
    cat('computing histogram of kmer count distribution...\n')
    binwidth = max(max(kmer_count_all$count)/10, max(kmer_count_all$count)/length(kmer_count_all$count))
    kmerHistogram <- ggplot(kmer_count_all, aes(count)) + 
      geom_histogram(binwidth = binwidth) + 
      geom_vline(xintercept = expected_frequency, color = "red") + 
      labs(x = 'count of each kmer', y = 'number of kmers', title = main) + 
      theme(title = element_text(size = 20), axis.title = element_text(size = 20), axis.text = element_text(size = 20),
            line = element_line(linetype = 'solid', size = 4),
            panel.grid = element_line(size = 1),
            panel.background = element_rect(fill = "white"))
    
    # barplot
    cat('computing barplot of kmer count...\n')
    kmerBarplot <- ggplot(kmer_count_all_sort, aes(x = seq, y = count)) + 
      geom_bar(stat = 'identity') +
      scale_x_discrete(limits = kmer_count_all_sort$seq) +
      geom_hline(yintercept = expected_frequency, color = "red") + 
      labs(x = 'kmers', y = 'counts', title = main)
    
    if(k < 3){
      theme_plot <- theme(title = element_text(size = 20), axis.title = element_text(size = 20), axis.text = element_text(size = 20),
                          line = element_line(linetype = 'solid', size = 4),
                          panel.grid = element_line(size = 1),
                          panel.background = element_rect(fill = "white"))
    }else{
      theme_plot <- theme(title = element_text(size = 20), axis.title = element_text(size = 20), axis.text = element_text(size = 20),
                          axis.text.x = element_blank(), 
                          line = element_line(linetype = 'solid', size = 4),
                          panel.grid = element_line(size = 1),
                          panel.background = element_rect(fill = "white"))
    }
    kmerBarplot <- kmerBarplot + theme_plot
    
  }
  
  z <- list(kmer_count_all = kmer_count_all, 
            kmer_count_top = kmer_count_top, 
            parameters = parameters, 
            kmerHistogram = kmerHistogram, 
            kmerBarplot = kmerBarplot)
}


##################################################################################################################################################################
# Analysis of fractions of 'A', 'T', 'C', 'G' of each sequence from a given library.
# plot the histogram of each nucleotide fraction ditribution of all sequences in the library
# input is a vector of the sequences from a library
# out put is a table of nucleotide fractions of each sequence
##################################################################################################################################################################
nucFractionAnalysis <- function(DNA_vector, plot = T){
  library(parallel)
  nuc <- c('A', 'T', 'C', 'G')
  fraction <- mclapply(nuc, mc.cores = 4, function(x){
    nucleotide <- x
    sapply(DNA_vector, function(x){
      length(gregexpr(nucleotide, x)[[1]]) / nchar(x)
    })
  })
  
  fractionBySeq <- c()
  for(i in 1:4){
    fractionBySeq <- cbind(fractionBySeq, fraction[[i]])
  }
  colnames(fractionBySeq) <- paste(nuc, 'fraction', sep = '_')
  
  if(plot){
    par(mfrow = c(4,1))
    sapply(1:4, function(x){
      hist(fractionBySeq[,x], breaks = seq(0, 1, 0.01), xlim = c(0, 1), main = colnames(fractionBySeq)[x],
           xlab = 'fraction', ylab = 'number of species')
    })
  }
  
  z <- fractionBySeq
}


##################################################################################################################################################################
# compute the alignment score of a given sequence to a series of libraries
# input arguments: 
#        seqToAlign, a list of libraries, each contain a vector of sequence strings
#          e.g. seqToAlign <- list(S001_top = names(S001_top), S002_top = names(S002_top), S003_top = names(S003_top))
#
#        posString, a sequence to align to
#          e.g. posString <- 'TGCACTTGTCATTTTGTATATGTATTTGGTTTTTGGCTCT'
#
#        type = c('global', 'local', 'overlap', 'global-local', 'local-global')
# output: a list of libraries with alignment scores
##################################################################################################################################################################
alignToPos <- function(seqToAlign, posString, type = 'global', scoreOnly = TRUE){
  library(doParallel)
  library(foreach)
  nthread = length(seqToAlign) * length(posString)
  cl <- makeCluster(nthread)
  registerDoParallel(cl)
  
  alignMatrix <- expand.grid(names(seqToAlign), posString)
  colnames(alignMatrix) <- c('seqToAlign', 'posString')
  alignScore <- foreach(i = 1:nthread, .combine = 'c') %dopar% {
    library(Biostrings)
    idx_lib <- which(names(seqToAlign) == alignMatrix$seqToAlign[i])
    align <- sapply(seqToAlign[[idx_lib]], function(x){
      pairwiseAlignment(alignMatrix$posString[i], x, type = type, scoreOnly = scoreOnly)
    })
    z <- list(align)
  }
  stopCluster(cl)
  names(alignScore) <- paste(alignMatrix$seqToAlign, alignMatrix$posString, sep = '_')
  z <- list(alignScore = alignScore, alignMatrix = alignMatrix)
}






##################################################################################################################################################################
# species vs count plots with various scales
# input_filelist: could be either a data.frame, first column is a list of filenames in the current directory, second column is the names of the library; 
#                 or a list of libraries, each library is already structured so that count as the first column, seq as the second.
# output a list of ggplot objects, one for each library
# to plot, use grid_plot_ggplot2()
# only for interactive jobs
##################################################################################################################################################################
speciesVsCountPlot_arkiv <- function(input_filelist, cutoff = 50, x_scale = c(0, 500), y_scale = c(0, 500)){
  library(doParallel)
  library(foreach)
  
  if(is.data.frame(input_filelist)){
    # if input is a data.frame of filenames and preferred library names
    nLib <- nrow(input_filelist)
  }else{
    # if input is a list of libraries
    nLib <- length(input_filelist)
  }
  cl <- makeCluster(nLib)
  registerDoParallel(cl)
  
  histPlotList <- foreach(i = 1:nLib, .combine = 'c') %dopar% {
    library(ggplot2)
    
    if(is.data.frame(input_filelist)){
      libs <- read.table(input_filelist$V1[i], head = F, stringsAsFactors = F)
      main_1 <- paste(input_filelist$V2[i], 'cutoff', cutoff, sep = '_')
    }else{
      libs <- input_filelist[[i]]
      main_1 <- paste(names(input_filelist)[i], 'cutoff', cutoff, sep = '_')
    }
    
    names(libs) <- c('count', 'seq')
    libs <- data.frame(libs, stringsAsFactors = F)
    libs_cutoff <- libs[libs$count > cutoff, ]
    rm(libs)
    
    histPlot_lib <- ggplot(libs_cutoff, aes(count)) + 
      coord_cartesian(xlim = x_scale, ylim = y_scale) +
      geom_histogram(binwidth = 1) +
      labs(x = 'count', y = 'number of species', title = main_1) + 
      theme(title = element_text(size = 20), axis.title = element_text(size = 20), axis.text = element_text(size = 20),
            line = element_line(linetype = 'solid', size = 4),
            panel.grid = element_line(size = 1),
            panel.background = element_rect(fill = "white"))
    z <- list(histPlot_lib)
    names(z) <- main_1
    z
  }
  stopCluster(cl)
  histPlotList
}

##################################################################################################################################################################
# species vs count plots with various scales
# input is a list of libraries to plot
# not suitable for large files
##################################################################################################################################################################
# speciesVsCountPlot_arkiv <- function(libListToPlot, xlimits = c(1000, 5000, 10000), ylimits = c(50, 100, 200, 400, 800, 1600), filename){
#  library(ggplot2)
#  # library(parallel)
#  # libListToPlot <- libList_cutoff20

#  nlib <- length(libListToPlot)
#  scales <- expand.grid(xlimits, ylimits)
#  numberOfScale <- nrow(scales)
#  plotsList <- lapply(1:nlib, function(x){
#    libToPlot <- libListToPlot[[x]]$lib
#    main_1 <- names(libListToPlot)[x]
#    plotsByScales <- lapply(1:numberOfScale, function(y){
#      histPlot <- ggplot(libToPlot, aes(count)) + 
#        coord_cartesian(xlim = c(0, scales[y, 1]), ylim = c(0, scales[y, 2])) +
#        geom_histogram(binwidth = 1) +
#        labs(x = 'count', y = 'number of species', title = main_1) + 
#        theme(title = element_text(size = 20), axis.title = element_text(size = 20), axis.text = element_text(size = 20),
#              line = element_line(linetype = 'solid', size = 4),
#              panel.grid = element_line(size = 1))
#    })
#  })
#  plotsList <- unlist(plotsList, recursive = F)

#  png(paste(Sys.Date(), filename, 'specesVsCount.png', sep = '_'), 1000*numberOfScale, 1000*nlib, res = 200)
#  grid_plot_ggplot2(plotsList, nrow = nlib, ncol = numberOfScale)
#  dev.off()
# }


##################################################################################################################################################################
# kmer distribution by library by cutoff
##################################################################################################################################################################
kmer_panel_plots <- function(libList_cutoffs, cutoff = 1, kmers = c(2, 3, 4, 5, 6, 7, 8), filename = 'kmer_hist_cutoff50'){
  # library(parallel)
  nkmers <- length(kmers)
  nlib = length(libList_cutoffs)
  
  kmer_plots <- lapply(1:nlib, function(x){
    cutoffs <- unlist(lapply(libList_cutoffs[[x]], function(lib)lib$cutoff))
    libToPlot <- libList_cutoffs[[x]][[which(cutoffs == cutoff)]]
    main_1 <- names(libList_cutoffs[[x]])[which(cutoffs == cutoff)]
    lapply(1:nkmers, function(y){
      kmer_distribution(libToPlot$lib, k = kmers[y], weighted = F, plot = T, main = paste(main_1, ' k=', kmers[y], sep = ''))
    })
  })
  kmer_plots <- unlist(kmer_plots, recursive = F)
  
  kmer_hist <- lapply(kmer_plots, function(x)x$kmerHistogram)
  png(paste(Sys.Date(), filename, 'hist.png', sep = '_'), 1000*nkmers, 1000*nlib, res = 200)
  grid_plot_ggplot2(kmer_hist, nrow = nlib, ncol = nkmers)
  dev.off()
  
  kmer_bar <- lapply(kmer_plots, function(x)x$kmerBarplot)
  png(paste(Sys.Date(), filename, 'bar.png', sep = '_'), 1000*nkmers, 1000*nlib, res = 200)
  grid_plot_ggplot2(kmer_bar, nrow = nlib, ncol = nkmers)
  dev.off()
  
  kmer_plots
}

##################################################################################################################################################################
# calculated counts of overlapping sequences (minimum count of each overlapping sequence)
##################################################################################################################################################################
overlap_weighted <- function(libList){
  #libList <- ls_lib
  n_lib <- length(libList)
  ## find and compute common sequences from libraries
  common_seq <- Reduce(function(x, y) merge(x, y, by = 'seq', all = F), libList)
  common_count <- common_seq[,2:ncol(common_seq)]
  common_count[is.na(common_count)] <- 0
  colnames(common_count) <- names(libList)
  common_seq <- common_seq$seq
  common_count$overlap_count <- sapply(1:length(common_seq), function(i) min(common_count[i,], na.rm = T))
  overlap_sum <- sum(common_count$overlap_count, na.rm = T)
}

##################################################################################################################################################################
# draw venn diagram considering count information, depend on overlap_weighted()
##################################################################################################################################################################
vennD_weighted <- function(
  libList, 
  print.mode = 'raw', colors = colors, 
  plot = T
){
  # libList <- ls_lib[1:3]
  
  library(VennDiagram)
  n_lib <- length(libList)
  
  if(n_lib == 2){
    area12 <- overlap_weighted(libList)
    areas12 <- sapply(libList, function(x) sum(x$count))
    vennD <- draw.pairwise.venn(
      area1 = areas12[1], area2 = areas12[2], cross.area = area12,
      category = names(libList),
      print.mode = print.mode, colors = colors
    )
    
  }else if(n_lib == 3){
    areas123 <- sapply(libList, function(x) sum(x$count))
    area123 <- overlap_weighted(libList)
    area12 <- overlap_weighted(libList[c(1,2)])
    area13 <- overlap_weighted(libList[c(1,3)])
    area23 <- overlap_weighted(libList[c(2,3)])
    
    vennD <- draw.triple.venn(
      area1 = areas123[1], area2 = areas123[2], area3 = areas123[3],
      n12 = area12, n13 = area13, n23 = area23,
      n123 = area123,
      category = names(libList),
      print.mode = print.mode, colors = colors
    )
  }else if(n_lib == 4){
    areas1234 <- sapply(libList, function(x) sum(x$count))
    area1234 <- overlap_weighted(libList)
    ol3 <- combn(4, 3)
    areas3 <- sapply(1:ncol(ol3), function(i) overlap_weighted(libList[ol3[,i]]))
    ol2 <- combn(4, 2)
    areas2 <- sapply(1:ncol(ol2), function(i) overlap_weighted(libList[ol2[,i]]))
    
    vennD <- draw.quad.venn(
      area1 = areas1234[1], area2 = areas1234[2], area3 = areas1234[3], area4 = areas1234[4],
      n12 = areas2[1], n13 = areas2[2], n14 = areas2[3], n23 = areas2[4], n24 = areas2[5], n34 = areas2[6],
      n123 = areas3[1], n124 = areas3[2], n134 = areas3[3], n234 = areas3[4], 
      n1234 = area1234,
      category = names(libList),
      print.mode = print.mode, colors = colors
    )
  }
  
  if(plot) grid.draw(vennD)
  
  z <- vennD
}



##################################################################################################################################################################
# obetain weighted overlapping/unique items (similar to venn diagram weighted by count)
##################################################################################################################################################################
get_weighted <- function(libList){
  n_lib <- length(libList)
  totalcount <- sum(sapply(libList, function(x) sum(x$count)))
  
  ## common items
  seqList <- lapply(libList, function(x) x$seq)
  common_item <- Reduce(intersect, seqList)
  common_count <- sum(sapply(libList, function(x){
    sum(x$count[which(x$seq %in% common_item)])
  }))
  common_percent <- common_count/totalcount*100
  
  ## unique items to each
  unique_item <- lapply(1:n_lib, function(i) setdiff(seqList[[i]], Reduce(union, seqList[-i])))
  unique_count <- lapply(1:n_lib, function(i){
    count <- libList[[i]]$count
    seq <- libList[[i]]$seq
    sum(count[which(seq %in% unique_item[[i]])])
  })
  names(unique_count) <- names(libList)
  unique_percent <- lapply(unique_count, function(x) x/totalcount*100)
  
  unique_count[[length(unique_count)+1]] <- common_count
  names(unique_count)[length(unique_count)] <- 'common'
  unique_percent[[length(unique_percent)+1]] <- common_percent
  z <- data.frame(
    name = c(names(unique_count)),
    count = unlist(unique_count),
    percent = unlist(unique_percent),
    row.names = NULL
  )
}

##################################################################################################################################################################
# Generate random sequence library
##################################################################################################################################################################
randomLibrary <- function(librarySize = 10000, seqLength = 35){
  nuc <- c('A', 'T', 'C', 'G')
  randomSeq <- sapply(1:librarySize, function(x){
    a <- nuc[sample(1:4, seqLength, replace = T)]
    paste(a, collapse = '')
  })
  randomSeq <- data.frame(count = 1, seq = randomSeq)
}

##################################################################################################################################################################
# Generate simulated sequence library with spike in motifs
##################################################################################################################################################################
spikeInMotif <- function(BackgroundSeq, motif, freq = 0.1){
  librarySize <- length(BackgroundSeq)
  spikeInSize <- round(freq * librarySize, 0)
  spikeInLength <- trunc(mean(sapply(BackgroundSeq, nchar)))
  
  motif_vec <- strsplit(motif, split = '')[[1]]
  motifLength <- length(motif_vec)
  
  spikeInSeq <- sapply(1:spikeInSize, function(x){
    a <- nuc[sample(1:4, spikeInLength, replace = T)]
    randomPosition <- sample(1:(spikeInLength - motifLength + 1), 1)
    a[randomPosition:(randomPosition + motifLength - 1)] <- motif_vec
    paste(a, collapse = '')
  })
  
  replaceIDX <- sample(librarySize, spikeInSize)
  PositiveSeq <- BackgroundSeq
  PositiveSeq[replaceIDX] <- spikeInSeq
  
  parameters <- list(motif = motif, freq = freq, motifLength = motifLength)
  
  z <- list(PositiveSeq = PositiveSeq, spikeInSeq = spikeInSeq, BackgroundSeq = BackgroundSeq, parameters = parameters)
}


##################################################################################################################################################################
# read files according to a filelist
# input: files, a dataframe with first column as filenames, second column as the prefered names of this sample
# output: A list of tables from the directory with filnames matching input first column
##################################################################################################################################################################
readFiles <- function(files){
  filenames <- files$V1
  fileList <- lapply(filenames, function(x){
    z <- read.table(x, head = F, stringsAsFactors = F)
    names(z) <- c('count', 'seq')
    z <- data.frame(z, stringsAsFactors = F)
    z
  })
  names(fileList) <- files$V2
  fileList
}

readFiles_para <- function(files){
  library(doParallel)
  library(foreach)
  nLib <- nrow(files)
  filenames <- files$V1
  
  cl <- makeCluster(nLib)
  registerDoParallel(cl)
  fileList <- foreach(i = 1:nLib, .combine = c('c')) %dopar% {
    libraryDF <- read.table(filenames[i], head = F, stringsAsFactors = F)
    names(libraryDF) <- c('count', 'seq')
    libraryDF <- data.frame(libraryDF, stringsAsFactors = F, row.names = NULL)
  }
  stopCluster(cl)
  names(fileList) <- files$V2
  fileList
}

##################################################################################################################################################################
# process library with only sequences above a certain cutoff
##################################################################################################################################################################
libList_cut <- function(libList, cutoff = 0){
  libList_out <- lapply(libList, function(x){
    featureIDX <- which(x$count > cutoff)
    x[featureIDX, ]
  })
}

##################################################################################################################################################################
# prepare list of libraries from a data.frame, output table from AptamerGroupAnalysis_v4.1
##################################################################################################################################################################
DFtoList_libAnalysis <- function(dataFrame, cutoff = 0){
  listNames <- colnames(dataFrame)
  nLib <- length(listNames)
  libList <- c()
  for(i in 1:nLib){
    featureIDX <- which(dataFrame[,i] > cutoff)
    libList[[i]] <- data.frame(count = dataFrame[,i][featureIDX],
                               seq = rownames(dataFrame)[featureIDX],
                               stringsAsFactors = F)
  }
  names(libList) <- listNames
  z <- libList
}








