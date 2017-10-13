################################################################################################################
## Obtain query from mySQL databases using Ryan's account
################################################################################################################
get_ryan_query <- function(database, query){
  library(RMySQL)
  drv <- dbDriver("MySQL") 
  dbcon <- dbConnect(MySQL(), user="jwang", password="Pass1234", host = 'bioserv.clsnet.intranet', dbname = database, port = 3307) 
  rs1 = dbSendQuery(dbcon, query)
  df_out <- dbFetch(rs1, n = -1)
  # dbListConnections(dbDriver( drv = "MySQL"))
  suppressWarnings(dbDisconnect(dbcon))
  # lapply(dbListConnections(dbDriver( drv = "MySQL")), dbDisconnect)
  df_out
}


########################################################################################################################################################################
## connect to database
## obtain case_reports_normal table
########################################################################################################################################################################
get_case_reports_normal <- function(){
  library(RMySQL)
  drv <- dbDriver("MySQL") 
  brds <- dbConnect(MySQL(), user="jwang", password="Pass1234", host = 'bioserv.clsnet.intranet', dbname = 'brds', port = 3307) 
  rs1 = dbSendQuery(brds, 'select * from case_reports_normal')
  case_reports_normal <- fetch(rs1, n = -1)
  dbListConnections( dbDriver( drv = "MySQL"))
  suppressWarnings(dbDisconnect(brds))
  case_reports_normal
}


########################################################################################################################################################################
## connect to database
## obtain cln_brca_case table
########################################################################################################################################################################
get_cln_brca_case <- function(){
  library(RMySQL)
  drv <- dbDriver("MySQL") 
  brds <- dbConnect(MySQL(), user="jwang", password="Pass1234", host = 'bioserv.clsnet.intranet', dbname = 'brds', port = 3307) 
  rs1 = dbSendQuery(brds, 'select * from cln_brca_case')
  cln_brca_case <- fetch(rs1, n = -1)
  dbListConnections( dbDriver( drv = "MySQL"))
  suppressWarnings(dbDisconnect(brds))
  cln_brca_case
}

########################################################################################################################################################################
## connect to database
## obtain case_specimen_status table
########################################################################################################################################################################
get_case_specimen_status <- function(){
  library(RMySQL)
  drv <- dbDriver("MySQL") 
  brds <- dbConnect(MySQL(), user="jwang", password="Pass1234", host = 'bioserv.clsnet.intranet', dbname = 'brds', port = 3307) 
  rs1 = dbSendQuery(brds, 'select * from case_specimen_status')
  case_specimen_status <- fetch(rs1, n = -1)
  dbListConnections( dbDriver( drv = "MySQL"))
  suppressWarnings(dbDisconnect(brds))
  case_specimen_status
}

########################################################################################################################################################################
## connect to database
## obtain brds_specimen_location table
########################################################################################################################################################################
get_brds_specimen_location <- function(){
  library(RMySQL)
  drv <- dbDriver("MySQL") 
  brds <- dbConnect(MySQL(), user="jwang", password="Pass1234", host = 'bioserv.clsnet.intranet', dbname = 'brds', port = 3307) 
  rs2 = dbSendQuery(brds, 'select * from brds_specimen_location')
  brds_specimen_location <- fetch(rs2, n = -1)
  dbListConnections( dbDriver( drv = "MySQL"))
  suppressWarnings(dbDisconnect(brds))
  brds_specimen_location
}


########################################################################################################################################################################
## function to examine current vial status in the database, especially useful to check whether the recently selected vials are updated or not; if updated, we can 
##  directly select samples from the database. There is no need to check again.
## input are CASE_Barcodes, Vial number, and the case_specimen_status table from brds
## output are currentStatusVials
########################################################################################################################################################################
get_vialStatus <- function(specimen_barcodes, case_specimen_status = NULL){
  if(is.null(case_specimen_status))
    case_specimen_status <- get_case_specimen_status()
  
  CASE_Barcodes <- substr(specimen_barcodes, 1, 9)
  Vials <- substr(specimen_barcodes, 11, 12)
  currentStatusSamples <- case_specimen_status[match(CASE_Barcodes, case_specimen_status$CASE_Barcode), ]
  VialList <- c('1A','1B','1C','2A','2B','2C','3A','3B','3C','4A','4B','4C')
  vialIDX <- match(Vials, VialList) + 5
  
  currentStatusVials <- unlist(lapply(1:nrow(currentStatusSamples), function(x){
    currentStatusSamples[x, vialIDX[x]]
  }))
  
  #currentStatusVials <- c()
  #for( i in 1: length(vialIDX) ){
  #  currentStatusVials <- c(currentStatusVials, currentStatusSamples[i, vialIDX[i]])
  #}
  
  z <- data.frame(
    CASE_Barcode = currentStatusSamples$CASE_Barcode, 
    Vials = Vials, 
    specimen_barcode = paste(currentStatusSamples$CASE_Barcode, Vials, sep = '-'), 
    current_vial_status_in_brds = currentStatusVials
  )
}

########################################################################################################################################################################
## For the case_specimen_status table, manually change vials status to 'RLV' for the selected, but not yet updated in the database. 
## The new table can then be used for new selection.
########################################################################################################################################################################
alter_case_specimen_status <- function(case_specimen_status, selected_specimen_barcode){
  #selected_specimen_barcode <- as.character(selected_vials$specimen_barcode)
  if(!is.vector(selected_specimen_barcode))
    stop('Please input a vector of specimen_barcode!\n')
  n_vials <- length(selected_specimen_barcode)
  selected_vials <- data.frame(
    CASE_Barcodes = substr(selected_specimen_barcode, 1, 9),
    vials = substr(selected_specimen_barcode, 11, 12)
  )
  VialList <- c('1A','1B','1C','2A','2B','2C','3A','3B','3C','4A','4B','4C')
  for(i in 1:n_vials){
    row_idx <- which(case_specimen_status$CASE_Barcode == selected_vials$CASE_Barcode[i])
    col_idx <- which(VialList == selected_vials$vials[i]) + 5
    case_specimen_status[row_idx, col_idx] <- 'RLV_selected'
  }
  case_specimen_status$rem_vials <- apply(case_specimen_status, 1, function(x) sum(x == 'Yes', na.rm = T))
  case_specimen_status
}


########################################################################################################################################################################
## Obtain consent level and HIPAA level from case_tbl table, simplify them, and merge to case_specimen_status table. 
## The new table can be then used for new selection.
########################################################################################################################################################################
add_consent_case_specimen_status <- function(case_specimen_status){
  case_tbl <- get_ryan_query('brds', 'select Barcode, ConsentLevel, LevelOfHIPAAAuthorization from case_tbl;')
  consent <- sapply(case_tbl$ConsentLevel, function(x){
    if(grepl('Tier 2', x)){
      'Tier 2'
    }else if(grepl('Tier 3', x)){
      'Tier 3'
    }else if(grepl('Tier 1', x)){
      'Tier 1'
    }else{
      x
    }
  })
  
  HIPAA <- sapply(case_tbl$LevelOfHIPAAAuthorization, function(x){
    if(grepl('Tier 2', x)){
      'Tier 2'
    }else if(grepl('Tier 1', x)){
      'Tier 1'
    }else{
      x
    }
  })
  case_simpl <- data.frame(CASE_Barcode = case_tbl$Barcode,
                           consent = consent, 
                           HIPAA = HIPAA)
  case_specimen_status <- unique(merge(case_specimen_status, case_simpl, by = 'CASE_Barcode', all.x = T))
}


########################################################################################################################################################################
## update case_specimen_status table by resitriction on hemolysis, only used in select_samples()
##    realizing this function by merging brds_specimen_location and df_selected_samples, and update the column 6:17 with hemo_fail, as well as add
## a rem_vial_hemo_restricted column
########################################################################################################################################################################
update_by_hemo <- function(
  # a dataframe with the same structure as case_specimen_status table
  df_selected_samples,
  # hemolysis option in the brds_specimen_location table
  hemolysis = c(1, 2)
){
  brds_specimen_location <- get_brds_specimen_location()
  brds_specimen_location$CASE_Barcode = substr(brds_specimen_location$specimen_barcode, 1, 9)
  df_selected_samples_vials <- unique(merge(df_selected_samples, brds_specimen_location, by = 'CASE_Barcode'))
  df_selected_samples_vials$vial <- substr(df_selected_samples_vials$specimen_barcode, 11, 12)

  VialList <- c('1A','1B','1C','2A','2B','2C','3A','3B','3C','4A','4B','4C')
  df_selected_samples_new <- lapply(1:nrow(df_selected_samples), function(x){
    #x <- 11
    #cat(x, '\n')
    each_sample <- df_selected_samples_vials[df_selected_samples_vials$CASE_Barcode == df_selected_samples$CASE_Barcode[x], ]
    
    ## update df_selected_samples column 6:17 with hemo_fail labels when hemolysis criteria are not met
    vial <- each_sample$vial[which(!each_sample$hemolysis %in% hemolysis)]
    vial_idx <- match(vial, VialList) + 5
    vial_idx <- vial_idx[!is.na(vial_idx)]
    each_selected_sample <- df_selected_samples[x, ]
    each_selected_sample[, vial_idx] <- 'hemo_fail'
    
    ## compute remaining number of vials after restricted on hemolysis status
    each_selected_sample$rem_vial_hemo_restricted <- length(which(each_sample$hemolysis %in% hemolysis))
    each_selected_sample
  })
  df_selected_samples_new <- do.call(rbind, df_selected_samples_new)
  df_selected_samples_new
}

########################################################################################################################################################################
## select samples according to case_specimen_status
## two modes:
##   1. select samples based on predefined input CASE_Barcode;
##      usage: select_samples(case_specimen_status, CASE_Barcode = input_case_barcode, n_vials = 1, min_rem_vials = 2)
##   2. select samples based on restrictions (e.g. n_samples, PROTOCOL_Name, CaseStatus, hemolysis, consent)
##      usage: select_samples(case_specimen_status, n_samples = 10, n_vials = 1, PROTOCOL_Name = 'Breast', CaseStatus = 'Non-Cancer', min_rem_vials = 2)
########################################################################################################################################################################
select_samples <- function(
  case_specimen_status,
  CASE_Barcode = NULL, 
  n_samples_to_select = NULL, n_vials_to_select = 1, 
  PROTOCOL_Name = NULL, 
  CaseStatus = NULL, 
  hemolysis = NULL, 
  PROTOCOL_Type = NULL,
  # need to add consent level and HIPAA level
  consent = NULL, HIPAA = NULL,
  # minimum remaining vials after current selection
  min_rem_vials = 2, 
  VAL = F
){
  
  # n_samples = 5; n_vials_to_select = 8; PROTOCOL_Name = NULL; CaseStatus = c('Normal_Male'); min_rem_vials = 4; hemolysis = c(1,2)
  input_criteria <- data.frame(
    PROTOCOL_Name = ifelse(is.null(PROTOCOL_Name), NA, PROTOCOL_Name), 
    CaseStatus = ifelse(is.null(CaseStatus), NA, CaseStatus), 
    min_rem_vials = ifelse(is.null(min_rem_vials), NA, min_rem_vials),
    hemolysis = ifelse(is.null(hemolysis), NA, paste(hemolysis, collapse = ',')),
    PROTOCOL_Type = ifelse(is.null(PROTOCOL_Type), NA, paste(PROTOCOL_Type, collapse = ',')),
    consent = ifelse(is.null(consent), NA, paste(consent, collapse = ',')),
    HIPAA = ifelse(is.null(HIPAA), NA, paste(HIPAA, collapse = ','))
  )
  
  if(VAL){
    case_specimen_status$rem_vials <- apply(
      case_specimen_status[, 6:17], 1, function(x){
        sum(x %in% c('Yes', 'Yes-VAL'))
      }
    )
  }
  
  if(all(is.null(CASE_Barcode), is.null(n_samples_to_select)))
    stop('Arguments missing! Please input at least CASE_Barcode or n_samples_to_select!\n')
  # case_specimen_status <- get_case_specimen_status()
  
  # selecting samples by input case_barcode
  if(!is.null(CASE_Barcode)){
    cat('selecting samples according to input case_barcodes...\nignoring n_samples_to_select...\n')
    df_selected_samples <- case_specimen_status[match(CASE_Barcode, case_specimen_status$CASE_Barcode), ]
    df_low_samples <- df_selected_samples[df_selected_samples$rem_vials - n_vials_to_select < min_rem_vials, ]
    if(nrow(df_low_samples) > 0){
      cat('not enough vials for sample:\n')
      print(df_low_samples$CASE_Barcode)
    }
  # select samples based on restrictions
  }else{
    cat('Applying input restrictions to select samples\n...\n')
    df_selected_samples <- case_specimen_status
    ## manditory: filtered by min_rem_vials
    df_selected_samples <- df_selected_samples[which(df_selected_samples$rem_vials - n_vials_to_select >= min_rem_vials), ]
    ## optional: filtered by PROTOCOL_Name
    if(!is.null(PROTOCOL_Name)){
      df_selected_samples <- df_selected_samples[which(df_selected_samples$PROTOCOL_Name %in% PROTOCOL_Name), ]
    }
    ## optional: fitered by CaseStatus
    if(!is.null(CaseStatus)){
      df_selected_samples <- df_selected_samples[which(df_selected_samples$CaseStatus %in% CaseStatus), ]
    }
    ## optional: fitered by hemolysis
    if(!is.null(hemolysis)){
      df_selected_samples <- update_by_hemo(df_selected_samples, hemolysis = hemolysis)
      # df_selected_samples$rem_vial_hemo_restricted - apply(df_selected_samples, 1, function(x)sum(x == 'Yes', na.rm = T))
      df_selected_samples <- df_selected_samples[which(df_selected_samples$rem_vial_hemo_restricted - n_vials_to_select >= min_rem_vials), ]
    }
    ## optional: filtered by PROTOCOL_Type
    if(!is.null(PROTOCOL_Type)){
      if(!'PROTOCOL_Type' %in% colnames(case_specimen_status)){
        stop('Please join case_specimen_status with local PROTOCOL_INFO table on COLLECTION_STRATEGY_Name!')
      }
      df_selected_samples <- df_selected_samples[which(df_selected_samples$PROTOCOL_Type %in% PROTOCOL_Type), ]
    }
    ## optional: filtered by consent level
    if(!is.null(consent)){
      df_selected_samples <- df_selected_samples[which(df_selected_samples$consent %in% consent), ]
    }
    ## optional: filtered by HIPAA level
    if(!is.null(HIPAA)){
      df_selected_samples <- df_selected_samples[which(df_selected_samples$HIPAA %in% HIPAA), ]
    }
    
    if(nrow(df_selected_samples) < n_samples_to_select){
      cat('Not enough samples to select, Please loose the input criteria!\n')
    }else{
      cat('Found', nrow(df_selected_samples),'samples matching the input criteria! ^_^ \n')
      cat('\nAfter selection,', nrow(df_selected_samples) - n_samples_to_select,'will be left under the input criteria:\n')
      print(input_criteria)
      ## randomly select samples from candidate pools, with rem_vials as weighted probability
      rem_vials <- df_selected_samples$rem_vials
      df_selected_samples <- df_selected_samples[sample(1:nrow(df_selected_samples), n_samples_to_select, replace = F, prob = rem_vials), ]
    }
  }
  
  df_selected_samples <- df_selected_samples[order(df_selected_samples$rem_vials, decreasing = T), ]
  df_selected_samples
}



########################################################################################################################################################################
## select vials from databases
########################################################################################################################################################################
select_vials <- function(df_selected_samples, n_vials_to_select = 1, VAL = F){
  # n_vials_to_select = 5
  # df_selected_samples = df_select_sample1
  
  VialList <- c('1A','1B','1C','2A','2B','2C','3A','3B','3C','4A','4B','4C')
  brds_specimen_location <- get_brds_specimen_location()
  # case_specimen_status <- get_case_specimen_status()
  
  if(VAL){
    status = c('Yes', 'Yes-VAL')
  }else{
    status = 'Yes'
  }
  
  list_vialIDX <- lapply(1 : nrow(df_selected_samples), function(x){
    idx_to_choose <- which(df_selected_samples[x,6:17] %in% status)
    if(length(idx_to_choose) == 1){
      idx_to_choose
    }else{
      sample(idx_to_choose, n_vials_to_select, replace = F)
    }
  })
  selected_vials <- lapply(list_vialIDX, function(x)VialList[x])
  selected_specimen_barcode <- lapply(1:nrow(df_selected_samples), function(x){
    data.frame(specimen_barcode = paste(df_selected_samples$CASE_Barcode[x], selected_vials[[x]], sep = '-'),
               CASE_Barcode = df_selected_samples$CASE_Barcode[x],
               vial = selected_vials[[x]],
               CaseStatus = df_selected_samples$CaseStatus[x])
  })
  selected_specimen_barcode <- do.call(rbind, selected_specimen_barcode)
  selectedSpecimen_location <- merge(selected_specimen_barcode, brds_specimen_location, by.x = 'specimen_barcode', by.y = 'specimen_barcode')
  
  # internal check vial status in the database
  vialStatus <- get_vialStatus(selectedSpecimen_location$specimen_barcode)
  if(any(!vialStatus$current_vial_status_in_brds %in% status)){
    cat('errors in selection!\nfollowing vials not available in database!\n')
    print(vialStatus[vialStatus$current_vial_status_in_brds!='Yes',])
    stop()
  }else{
    vialStatus <- vialStatus[, c('specimen_barcode', 'current_vial_status_in_brds')]
    selectedSpecimen_location <- merge(selectedSpecimen_location, vialStatus, by.x = 'specimen_barcode', by.y = 'specimen_barcode')
  }
  selectedSpecimen_location
}


########################################################################################################################################################################
## examine vial status of recently selected samples
########################################################################################################################################################################
get_vialStatus_selected <- function(file_URL){
  recently_selected <- read.csv(file_URL, head = T, stringsAsFactors = F)
  vialStatus <- get_vialStatus(recently_selected$specimen_barcode)
  vialStatus_notRLV <- vialStatus[vialStatus$current_vial_status_in_brds != 'RLV', ]
  z <- list(vialStatus = vialStatus, vialStatus_notRLV = vialStatus_notRLV)
}


########################################################################################################################################################################
## update local recent_selected folder with a new .csv file
## append selected vials into the file
########################################################################################################################################################################
update_local_record <- function(select_vials, last_file_URL){
  recently_selected <- read.csv(last_file_URL, head = T, stringsAsFactors = F, row.names = 1)
  current_vialStatus <- get_vialStatus(recently_selected$specimen_barcode)
  # update with current status in brds
  recently_selected$current_vial_status_in_brds <- current_vialStatus$current_vial_status_in_brds
  # add Date_selected column
  select_vials$Date_selected <- format(Sys.Date(), '%m/%d/%Y')
  # append the newly selected vials
  recently_selected <- rbind(recently_selected, select_vials)
  write.csv(recently_selected, file = paste('~/Documents/Project5_Biorepository/recently_selected_samples/', Sys.Date(), '_updated.csv', sep = ''))
}

########################################################################################################################################################################
# demographic summarize module
########################################################################################################################################################################

## trim off heading/trailing spaces and/or punctuations from string
trim <- function (x) gsub("^(\\s|[[:punct:]])+|([[:punct:]]|\\s)+$", "", x)

# ER
get_ER_status <- function(df_clinical_info){
  ER <- trim(df_clinical_info$EstrogenReceptorStatus)
  SF_ER <- trim(df_clinical_info$SFEstrogenReceptorStatus)
  ER_out <- sapply(1:length(ER), function(x){
    ifelse(ER[x]=='Unknown' | ER[x]=='Not Performed' | is.na(ER[x]), SF_ER[x], ER[x])
  })
  ER_out[ER_out == 'Not Performed' | ER_out == 'NULL' | is.na(ER_out)] <- 'Unknown'
  ER_out
}

# PR
get_PR_status <- function(df_clinical_info){
  PR <- trim(df_clinical_info$ProgesteroneReceptorStatus)
  SF_PR <- trim(df_clinical_info$SFProgesteroneReceptorStatus)
  PR_out <- sapply(1:length(PR), function(x){
    ifelse(PR[x]=='Unknown' | PR[x]=='Not Performed' | is.na(PR[x]), SF_PR[x], PR[x])
  })
  PR_out[PR_out == 'Not Performed' | PR_out == 'NULL' | is.na(PR_out)] <- 'Unknown'
  PR_out
}

# HER2
get_HER2_status <- function(df_clinical_info){
  HER2 <- trim(df_clinical_info$HER2GeneAmplificationStatusByFISH)
  SF_HER2 <- trim(df_clinical_info$SFHER2GeneAmplificationStatusByFISH)
  HER2_out <- sapply(1:length(HER2), function(x){
    ifelse(HER2[x]=='Unknown' | HER2[x]=='Not Performed' | is.na(HER2[x]), SF_HER2[x], HER2[x])
  })
  HER2_out[HER2_out == 'Not Performed' | HER2_out == 'NULL' | is.na(HER2_out)] <- 'Unknown'
  HER2_out
}

# Stage
get_stage <- function(df_clinical_info, collapse = T){
  StagePathol <- trim(df_clinical_info$AJCCGroupStagePathology)
  StageClinic <- trim(df_clinical_info$AJCCGroupStageClinical)
  stage_out <- sapply(1:length(StagePathol), function(x){
    ifelse(StagePathol[x] == 'NULL', StageClinic[x], StagePathol[x])
  })
  if(collapse){
    stage_out[stage_out %in% c('IA', 'IB', 'IC')] <- 'I'
    stage_out[stage_out %in% c('IIA', 'IIB', 'IIC')] <- 'II'
    stage_out[stage_out %in% c('IIIA', 'IIIB', 'IIIC')] <- 'III'
    stage_out[stage_out %in% c('IVA', 'IVB', 'IVC')] <- 'IV'
  }
  stage_out[stage_out == 'Not Performed' | stage_out == 'NULL'] <- 'Unknown'
  stage_out
}

# Stage T
get_stage_T <- function(df_clinical_info){
  Stage_T_Pathol <- trim(df_clinical_info$AJCCStageTPathology)
  Stage_T_Clinic <- trim(df_clinical_info$AJCCStageTClinical)
  stage_T_out <- sapply(1:length(Stage_T_Pathol), function(x){
    ifelse(Stage_T_Pathol[x] == 'NULL', Stage_T_Clinic[x], Stage_T_Pathol[x])
  })
  stage_T_out[stage_T_out == 'Not Performed' | stage_T_out == 'NULL'] <- 'Unknown'
  stage_T_out
}

# Stage N
get_stage_N <- function(df_clinical_info){
  Stage_N_Pathol <- trim(df_clinical_info$AJCCStageNPathology)
  Stage_N_Clinic <- trim(df_clinical_info$AJCCStageNClinical)
  stage_N_out <- sapply(1:length(Stage_N_Pathol), function(x){
    ifelse(Stage_N_Pathol[x] == 'NULL', Stage_N_Clinic[x], Stage_N_Pathol[x])
  })
  stage_N_out[stage_N_out == 'Not Performed' | stage_N_out == 'NULL'] <- 'Unknown'
  stage_N_out
}

# Stage M
get_stage_M <- function(df_clinical_info){
  Stage_M_Pathol <- trim(df_clinical_info$AJCCStageMPathology)
  Stage_M_Clinic <- trim(df_clinical_info$AJCCStageMClinical)
  stage_M_out <- sapply(1:length(Stage_M_Pathol), function(x){
    ifelse(Stage_M_Pathol[x] == 'NULL', Stage_M_Clinic[x], Stage_M_Pathol[x])
  })
  stage_M_out[stage_M_out == 'Not Performed' | stage_M_out == 'NULL'] <- 'Unknown'
  stage_M_out
}

# age
get_cat_age <- function(age){
  cut(age, c(0, 20, 30, 40, 50, 60, 70, 80, 100))
}




