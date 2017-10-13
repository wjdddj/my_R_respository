suppressWarnings(library(data.table))
suppressWarnings(library(gsubfn))

print_help <- function(){
  cat('This is a script to plot histograms of aptamer sequence length.\n')
  cat('Usage:\n')
  cat('    Rscript getSeqLengthHist.R [inputFile] [outputDir]\n')
  cat('inputFile is the standard \'enSeq\' file output from the aptamer pipeline.\n')
}

getSeqFromEnSeqFile <- function(
  seqFile
){
  seqDF <- fread(seqFile, header = F, stringsAsFactors = F, data.table = F, select = c(2))
  seq <- seqDF$V2
}

getSeqLength <- function(
  seq
){
  seqLength <- nchar(seq)
}

getSampleNameFromEnSeqFile <- function(
  seqFile
){
  fileName <- basename(seqFile)
  sampleName <- strapply(fileName, '^enSeq_([^_]*)_.*$', simplify = T)
}

makeSeqLengthHist <- function(
  seqFile, outDir
){
  
  now <- format(Sys.Date(), '%Y%m%d')
  
  seq <- getSeqFromEnSeqFile(seqFile)
  seqLength <- nchar(seq)
  sampleName <- getSampleNameFromEnSeqFile(seqFile)
  outName <- paste0(outDir, now, '_', sampleName, '.pdf')
  
  pdf(outName)
  hist(
    seqLength, breaks = seq(20, 50, 1),
    main = sampleName
  )
  dev.off()
}


inputArguments <- commandArgs(trailingOnly = TRUE)
if(length(inputArguments) == 0){
  print_help()
  stop()
}

seqFile <- inputArguments[1]
outDir <- inputArguments[2]
# seqFile <- '/gpfs1/dev/bioserv/home/jwang/Aptamer_out/061417-SeqofEnrichedlibHiseq1A-FTRin35nB-7pM-0pctPhiX-S001S005-ExpID6692/enSeq_S001-HJC7KBCXY_S00_L00_R1_none.V.44046.SUM.csv'
makeSeqLengthHist(seqFile, outDir)



