# source('~/R_modules/clinical_data/clinical_data_pipeline.R')

readMA <- function(){
  dfMA <- read.csv(
    '~/Documents/Project8_Evofosfamide/MAESTRO_trial/MAESTRO speadsheet - Caris - 10-16.csv', 
    header = T, stringsAsFactors = F
  )
  colnames(dfMA) <- gsub('\\.+', '_', colnames(dfMA))
  dfMA$Row <- NULL
  
  dfMA_v2 <- read.csv(
    '~/Documents/Project8_Evofosfamide/MAESTRO_trial/MAESTRO speadsheet - Caris - ver2 - 10-31-16.csv', 
    header = T, stringsAsFactors = F
  )
  colnames(dfMA_v2) <- gsub('\\.+', '_', colnames(dfMA_v2))
  dfMA_v2$Row <- NULL
  
  dfMA_v2$is_treated <- dfMA$Group != ''
  dfMA_v2
}

dfMA <- readMA()
