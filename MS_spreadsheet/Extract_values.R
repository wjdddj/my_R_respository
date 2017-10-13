######################################################################################################################################################
# Date: 2016-4-26
# Author: Ryan Wang
# Description: 
# This module aims to provide tools for extracting values from Proteome Discovery software outputs (either csv or xlsx files).
# currently the implemented functions are 
#     1. check the file type; 
#     2. examine number of files and sheets to analyze;
#     3. read in csv or xlsx files;
#     4. extract values
#
# The general work flow is:
#     1. organize Proteome Discovery output files in a single directory under the same file format (either csv or xlsx).
#     2. run this script's main_value_extract() function.
#     3. output csv files can be found in the working directory.
######################################################################################################################################################


######################################################################################################################################################
## check filetype
######################################################################################################################################################
get_file_type <- function(filedir){
  filenames <- dir(filedir)
  file_exten <- unique(sapply(strsplit(filenames, split = '\\.'), function(x)x[length(x)]))
  # file_exten <- unique(strsplit(filenames, split = '\\.')[[1]][2])
  if(length(file_exten) > 1){
    stop('more than one type of files in the directory, please maintain only one file type!\n\'xlsx\' or \'csv\'')
  }else{
    cat('detected filetype:', file_exten, '\n')
  }
  file_exten
}

######################################################################################################################################################
## check number of sheets of each sheet
######################################################################################################################################################
check_file_sheets <- function(filedir, filetype){
  filenames <- dir(filedir)
  cat('Detected', length(filenames), 'files.\n')
  if(length(filenames) <= 30){
    cat('detected files are:\n')
    print(filenames)
    cat('\n')
  }
  z <- list(filenames = filenames)
  # setwd(filedir)
  ## check number of sheets for xlsx filetype
  if(filetype == 'xlsx'){
    library(openxlsx)
    library(data.table)
    sheetNames <- lapply(paste(filedir, filenames, sep = ''), getSheetNames)
    names(sheetNames) <- filenames
    n_sheet <- sapply(sheetNames, length)
    output <- data.frame(filenames = filenames, n_sheet)
    rownames(output) <- NULL
    cat('sheets per table:\n')
    print(output)
    cat('\n')
    z <- list(filenames = filenames, sheetNames = sheetNames)
  }
  z
}

######################################################################################################################################################
## read protein discovery output spreadsheets
######################################################################################################################################################
read_files <- function(filedir, file_sheet, filetype){
  
  cat('reading files...\n')
  filenames <- file_sheet$filenames
  if(filetype == 'xlsx'){
    library(openxlsx)
    library(data.table)
    fileList <- lapply(1:length(filenames), function(x){
      sheetNames <- file_sheet$sheetNames[[x]]
      file_sheet_list <- lapply(1:length(sheetNames), function(y){
        read.xlsx(paste(filedir, filenames[x], sep = ''), sheet = y)
      })
      names(file_sheet_list) <- sheetNames
      file_sheet_list
    })
    fileList <- unlist(fileList, recursive = F)
  }
  if(filetype == 'csv'){
    fileList <- lapply(paste(filedir, filenames, sep = ''), function(x){
      z <- read.csv(x, header = T, stringsAsFactors = F)
    })
    names(fileList) <- sapply(strsplit(filenames, split = '\\.'), function(x)x[1])
  }
  fileList
}


######################################################################################################################################################
## function to extract values of a certain attribute from a series of output spreadsheets from Proteome Discovery, 
## by matching all(Accession, Description, and Gene_ID).
## input: file_list, a list of spreadsheets with preferrably similar structures
##        attr_name, name of the attributes to extract from the file_list; 
##        Exact, indicate whether the input attr_name is the exact colname or not, for attr_name = '#.PSM', usually use Exact = T; 
##               for attr_name = 'Area', usually use Exact = F
##        cutoff_coverage, protein coverage cutoff filter
## output: the extracted data.frame
######################################################################################################################################################
extract_value <- function(file_list, attr_name, Exact = T, cutoff_coverage = 0, empty_cell = 'n.d.', write_table = T, filename = ''){
  library(data.table)
  # cat('extracting values...\n')
  cat('extracting values...\n')
  # obtain unique ids of all spreadsheets after filtering
  uniqComb <- unique(rbindlist(lapply(file_list, function(x){
    # filter by coverage
    if(cutoff_coverage > 0){
      x <- x[x$Coverage > cutoff_coverage, ]
    }
    data.frame(x$Accession, x$Description, x$Gene.ID)
  })))
  colnames(uniqComb) <- c('Accession', 'Description', 'Gene_ID')
  
  # extracting information based on matching between uniqComb and uniqComb_file
  attr_value_list <- lapply(1:length(file_list), function(x){
    attr_value <- rep(empty_cell, length(uniqComb$Description))
    
    # obtain unique ids of each spreadsheets
    uniqComb_file <- data.frame(file_list[[x]]$Accession, file_list[[x]]$Description, file_list[[x]]$Gene.ID)
    colnames(uniqComb_file) <- c('Accession', 'Description', 'Gene_ID')
    # matching each entry based on combined ids (Accession, Description, Gene_ID)
    uniqIDs_file <- apply(uniqComb_file, 1, function(y){paste(y, collapse = '_')})
    uniqIDs_to_match <- apply(uniqComb, 1, function(y){paste(y, collapse = '_')})
    interAccess <- intersect(uniqIDs_to_match, uniqIDs_file)
    idx1 <- match(interAccess, uniqIDs_to_match)
    idx2 <- match(interAccess, uniqIDs_file)
    if(Exact){
      if(sum(colnames(file_list[[x]]) == attr_name) != 1){
        stop('Cannot identify unique column. Please input new attr_name!\n')
      }
      attr_value[idx1] <- file_list[[x]][, colnames(file_list[[x]]) == attr_name][idx2]
    }else{
      if(sum(grepl(attr_name, colnames(file_list[[x]]))) != 1){
        stop('Cannot identify unique column. Please input new attr_name!\n')
      }
      attr_value[idx1] <- file_list[[x]][, grepl(attr_name, colnames(file_list[[x]]))][idx2]
    }
    # attr_value <- as.numeric(attr_value)
    attr_value
  })
  names(attr_value_list) <- names(file_list)
  attr_value_table <- data.frame(do.call(cbind, attr_value_list))
  # attr_value_table <- attr_value_table[, order(colnames(attr_value_table))]
  attr_value_table <- data.frame(Accession = uniqComb$Accession, Description = uniqComb$Description, Gene_ID = uniqComb$Gene_ID, attr_value_table)
  
  if(write_table){
    write.csv(attr_value_table, file = paste(Sys.Date(), filename, '_coverageAbove', cutoff_coverage, '_', attr_name, '.csv', sep = ''))
  }
  
  attr_value_table
}



######################################################################################################################################################
## Main function for extracting PSM and Area from Proteome Discovery output files
######################################################################################################################################################
main_value_extract <- function(filedir){
  workingDIR <- getwd()
  file_type = get_file_type(filedir)
  file_sheet = check_file_sheets(filedir, filetype = file_type)
  fileList <- read_files(filedir, file_sheet, filetype = file_type)
  setwd(workingDIR)
  PSMTab <- extract_value(fileList, attr_name = 'PSMs', Exact = F, cutoff_coverage = 0, write_table = T, empty_cell = 'n.d.')
  AreaTab <- extract_value(fileList, attr_name = 'Area', Exact = F, cutoff_coverage = 0, write_table = T, empty_cell = 'n.d.')
  cat('finished!\n')
}




######################################################################################################################################################
## to combine .csv files into one .xls file
## input is the directory that only contain .csv files to combine into an .xls file.
######################################################################################################################################################
combine_csv <- function(file_dir){
  library(WriteXLS)
  library(plyr)
  library(gsubfn)
  files <- dir(file_dir)
  ls_files <- llply(1:length(files), function(i){
    read.csv(paste0(file_dir, files[i]), head = T, stringsAsFactors = F)
  })
  names(ls_files) <- strapplyc(files, '^(.+)\\.csv$', simplify = T)
  WriteXLS(ls_files, ExcelFileName = paste0(Sys.Date(), '_byPhysician.xls'), SheetNames = names(ls_files))
}








