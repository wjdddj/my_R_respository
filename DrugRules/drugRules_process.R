source('~/R_modules/DrugRules/drugRules.R')

## IHC status: c('Positive', 'Negative', 'Equivocal', 'Other', 'Test not performed')
## SEQ status: c('Mutated', 'Mutated - Other', 'Mutation Not Detected', 'No Result', 'QNS', 'Test not performed')
convert_SEQ_result <- function(SEQ_results){
  #SEQ_results <- c('No Result', 'QNS', 'Test not performed')
  SEQ_results <- gsub('Mutated - Other|No Result|QNS|Test not performed', 'Nodata', SEQ_results)
}

convert_IHC_result <- function(IHC_result){
  #IHC_result <- c('Equivocal', 'Other', 'Test not performed')
  IHC_result <- gsub('Equivocal|Other|Test not performed', 'Nodata', IHC_result)
}

match_markers <- function(marker_status, markers_in_rule){
  #markers_in_rule <- c('BRCA1_SEQ', 'BRCA2_SEQ', 'ERCC1_IHC')
  #marker_status <- data.frame(biomarker_tech = c('BRCA1_SEQ', 'BRCA2_SEQ'),
  #                            status = c('Mutated - Other', 'No Result'))
  marker_status_out <- data.frame(biomarker_tech = markers_in_rule, stringsAsFactors = F)
  match_idx <- match(markers_in_rule, marker_status$biomarker_tech)
  marker_status_out$status <- as.character(marker_status$status[match_idx])
  marker_status_out$status[which(is.na(match_idx))] <- 'Nodata'
  marker_status_out
}


read_cmi_mstatus <- function(dir_cmi_file){
  #dir_cmi_file <- '~/Documents/Project3_Drug_Biomarker_Survival/RProjects/20160613_drug_rule_explore/drugrule_files/20160616_cmi_ova.csv'
  df_mstatus <- read.csv(dir_cmi_file, head = T, stringsAsFactors = F, row.names = 1)
  df_mstatus
}

marker_format <- function(df_mstatus){
  library(dplyr)
  df_mstatus <- df_mstatus %>%
    mutate(biomarker_tech = paste(biomarkerName, Tech, sep = '_'),
           status = conclusionresult) %>%
    select(MasterDeID, biomarker_tech, status) %>%
    group_by(MasterDeID)
  df_mstatus
}

drug_recommend <- function(df_mstatus_by_case, drugRules){
  library(dplyr)
  library(plyr)
  #library(multidplyr)
  list_drugrules <- as.list(as.character(df_drugrules$functions))
  n_rules <- length(list_drugrules)
  ls_recom <- llply(1:n_rules, function(x){
    dr_fun <- get(list_drugrules[[x]])
    df_mstatus_by_case %>%
      do(dr_fun(.))
  })
  recom <- do.call(rbind, ls_recom)
  #recom <- data.frame(recom)
}

write_recommend <- function(recom){
  now <- format(Sys.time(), '%Y-%m-%d-%H-%M-%S')
  out_file_name <- paste0(now, '_', 'recommend', '_', drug_rule_version, '.csv')
  write.csv(recom, out_file_name)
}








