################################################################################################################
## detach specific loaded package
################################################################################################################
detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}


################################################################################################################
## detach all loaded packages
################################################################################################################
detachAllPackages <- function() {
  basic_packages <- c("package:stats",
                      "package:graphics",
                      "package:grDevices",
                      "package:utils",
                      "package:datasets",
                      "package:methods",
                      "package:base")
  package_list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package_list <- setdiff(package_list, basic_packages)
  if (length(package_list)>0)  
    lapply(package_list, detach, character.only=TRUE)
}

################################################################################################################
## Obtain table from mySQL databases using Sting's account
################################################################################################################
get_sting_table <- function(database, table_name){
  library(RMySQL)
  drv <- dbDriver("MySQL") 
  dbcon <- dbConnect(MySQL(), user="sting", password="Pass1234", host = 'bioserv.clsnet.intranet', dbname = database, port = 3307) 
  rs1 = dbSendQuery(dbcon, paste0('select * from ', table_name))
  df_out <- dbFetch(rs1, n = -1)
  # dbListConnections(dbDriver( drv = "MySQL"))
  suppressWarnings(dbDisconnect(dbcon))
  # lapply(dbListConnections(dbDriver( drv = "MySQL")), dbDisconnect)
  df_out
}

################################################################################################################
## Obtain query from mySQL databases using Sting's account
################################################################################################################
get_sting_query <- function(database, query){
  library(RMySQL)
  drv <- dbDriver("MySQL") 
  dbcon <- dbConnect(MySQL(), user="sting", password="Pass1234", host = 'bioserv.clsnet.intranet', dbname = database, port = 3307) 
  rs1 = dbSendQuery(dbcon, query)
  df_out <- dbFetch(rs1, n = -1)
  # dbListConnections(dbDriver( drv = "MySQL"))
  suppressWarnings(dbDisconnect(dbcon))
  # lapply(dbListConnections(dbDriver( drv = "MySQL")), dbDisconnect)
  df_out
}

################################################################################################################
## Obtain table from mySQL databases using Ryan's account
################################################################################################################
get_ryan_table <- function(database, table_name){
  library(RMySQL)
  drv <- dbDriver("MySQL") 
  dbcon <- dbConnect(MySQL(), user="jwang", password="Pass1234", host = 'bioserv.clsnet.intranet', dbname = database, port = 3307) 
  rs1 = dbSendQuery(dbcon, paste0('select * from ', table_name))
  df_out <- dbFetch(rs1, n = -1)
  # dbListConnections(dbDriver( drv = "MySQL"))
  suppressWarnings(dbDisconnect(dbcon))
  # lapply(dbListConnections(dbDriver( drv = "MySQL")), dbDisconnect)
  df_out
}

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

################################################################################################################
## Compute tnt from therapy data of a SINGLE subject, using input of line numbers of start_line and end_line. 
## input:
##    1. df_therapy_by_subject, a data.frame of therapy of a single subject. It MUST contain these columns: 
##  'masterdeid', 'line_num', 'regimen_start', 'deathdate', 'lastcontactdate'.
##    2. start_line, line number of the start line.
##    3. end_line, line number of the end line, it could also be 'max', equivalent to overall survival from start
##  of the start line.
##
## output:
##    1. tnt_event, a data.frame contains tnt and event columns with only one row of data.
##    2. df_therapy_by_subject, the updated input data.frame with tnt and event columns added. The tnt and event
##  will have the same value among rows.
##
## rule: 
##    1. event is considered as either start of the end_line, or death when end_line is not available; 
##  censored when no deathdate available.
##    2. if event, tnt is the time difference between start of start_line and start of end_line, or time difference 
##  between start of start_line and death; if censored, tnt is the time difference between start of start_line and 
##  last contact date.
################################################################################################################
get_tnt <- function(df_therapy_by_subject, start_line = 1, end_line = 'max'){
  library(plyr)
  #start_line = 2
  #end_line = 3
  #df_therapy_by_subject <- df_therapy_ken[df_therapy_ken$masterdeid == '9632', ]
  #df_therapy_by_subject <- dlply(df_therapy, .variables = 'masterdeid')[[1]]
  regimen_start <- as.Date(df_therapy_by_subject$regimen_start, format = '%Y-%m-%d')
  deathdate <- as.Date(df_therapy_by_subject$deathdate, format = '%Y-%m-%d')
  lastcontactdate <- as.Date(df_therapy_by_subject$lastcontactdate, format = '%Y-%m-%d')
  line_num <- df_therapy_by_subject$line_num
  
  ## obtain the start time point 
  start_idx <- which(line_num == start_line)
  if(length(start_idx) != 1 | is.na(regimen_start[start_idx])){
    #cat('start point: regimen_start of line', start_line, 'is not available for case', df_therapy_by_subject$masterdeid[1], 
    #    '. NA will be returned\n')
    startpoint <- NA
  }else{
    startpoint <- regimen_start[start_idx]
  }
  
  ## obtain the end time point
  if(end_line == 'max'){
    end_idx <- which.max(line_num)
    if(!is.na(deathdate[end_idx])){
      ## there is an event characterized by death
      endpoint <- deathdate[end_idx]
      event <- 1
    }else{
      ## censored at last contact date
      endpoint <- lastcontactdate[end_idx]
      event <- 0
    }
  }else{
    end_idx <- which(line_num == end_line)
    if(length(end_idx) != 1){
      #cat('end point: regimen_start of line', end_line, 'is not available for case', df_therapy_by_subject$masterdeid[1], 
      #    '. lastcontact date or death date will be returned\n')
      end_idx <- which.max(line_num)
      if(!is.na(deathdate[end_idx])){
        ## there is an event characterized by death
        endpoint <- deathdate[end_idx]
        event <- 1
      }else{
        ## censored at last contact date
        endpoint <- lastcontactdate[end_idx]
        event <- 0
      }
    }else if(is.na(regimen_start[end_idx])){
      #cat('end point: regimen_start of line', end_line, 'is not available for case', df_therapy_by_subject$masterdeid[1], 
      #    '. lastcontact date or death date will be returned\n')
      end_idx <- which.max(line_num)
      if(!is.na(deathdate[end_idx])){
        ## there is an event characterized by death
        endpoint <- deathdate[end_idx]
        event <- 1
      }else{
        ## censored at last contact date
        endpoint <- lastcontactdate[end_idx]
        event <- 0
      }
    }else{
      endpoint <- regimen_start[end_idx]
      event <- 1
    }
  }
  
  ## any error will be returned as NA, and can be filtered out after
  tnt <- tryCatch(endpoint - startpoint + 1, error = function(e)return(NA))
  tnt_event <- data.frame(tnt = as.numeric(tnt), event = event)
  df_therapy_by_subject <- data.frame(df_therapy_by_subject, 
                                      tnt_event)
  z <- list(tnt_event = tnt_event, df_therapy_by_subject = df_therapy_by_subject)
}


################################################################################################################
## Obtain therapy combinations from start_line to end_line.
## input:
##    1. df_therapy_by_subject, a data.frame of therapy of a single subject. It MUST contain these columns: 
##  'masterdeid', 'line_num'; optional for tnt calculation: 'regimen_start', 'deathdate', 'lastcontactdate'.
##    2. start_line, line number of the start line.
##    3. end_line, line number of the end line, it could also be 'max', the last line will be included.
##
## output: a character string of all drugs used between start_line and end_line, collapsed with ','.
################################################################################################################
get_therapy <- function(df_therapy_by_subject, start_line, end_line){
  therapy <- df_therapy_by_subject$therapies
  line_num <- df_therapy_by_subject$line_num
  if(end_line == 'max'){
    end_line <- which.max(line_num)
    lines_idx <- which(line_num >= start_line & line_num <= end_line)
  }else{
    lines_idx <- which(line_num >= start_line & line_num < end_line)
  }
  if(length(lines_idx) == 0){
    cat('no treatment is found between', start_line, 'and', end_line, 'inclusively. NA will be returned as therapy\n')
    out_therapy <- NA
  }else{
    out_therapy <- paste(therapy[lines_idx], collapse = ',')
    out_therapy <- paste(sort(unique(strsplit(out_therapy, ',')[[1]])), collapse = ',')
  }
  out_therapy
}

################################################################################################################
## A wrapper function to obtain therapies and tnts for multiple subjects using get_tnt() and get_therapy().
## input:
##    1. df_therapy, a data.frame of therapy of a single or multiple subjects. It MUST contain these columns: 
##  'masterdeid', 'line_num', 'regimen_start', 'deathdate', 'lastcontactdate'.
##    2. start_line, line number of the start line.
##    3. end_line, line number of the end line, it could also be 'max', the last line will be included.
##
## output:
##    1. df_therapy_out, a data.frame with three columns, masterdeid, therapies, tnt, event.
################################################################################################################
get_therapy_tnt_df <- function(df_therapy, start_line = 1, end_line = 'max'){
  # df_therapy <- df_therapy_2nd
  # start_line = 2
  # end_line = 3
  library(plyr)
  ls_therapy <- dlply(df_therapy, .variables = 'masterdeid', function(x){
    # cat(x$masterdeid, '\n')
    # x <- dlply(df_therapy, .variables = 'masterdeid')[[1]]
    tnt_event <- get_tnt(x, start_line = start_line, end_line = end_line)$tnt_event
    therapies <- get_therapy(x, start_line = start_line, end_line = end_line)
    out <- data.frame(masterdeid = x$masterdeid[1], therapies = therapies, tnt_event)
  })
  df_therapy_out <- do.call(rbind, ls_therapy)
  df_therapy_out <- data.frame(df_therapy_out)
  rownames(df_therapy_out) <- NULL
  df_therapy_out
}


################################################################################################################
## Compute tnt of each line from therapy data of each subject. 
## input: 
##    1. df_therapy_by_subject, a data.frame of therapy by subject. It MUST contain these columns with name 
##  'line_num', 'regimen_start', 'deathdate', 'lastcontactdate'.
##
## output:
##    1. df_therapy_by_subject, the updated input data.frame with tnt and event columns added. 
##
## rule: 
##    1. event is considered as either start of the nextline, or death when nextline is not available; 
##  censored when no deathdate available.
##    2. if event, tnt is the time difference between start of current and start of next line or time difference 
##  between start of current line and death; if censored, tnt is the time difference between start of current and 
##  last contact date.
##    3. if nextline is not available while not max line, NA will be returned.
################################################################################################################
get_tnt_by_line <- function(df_therapy_by_subject){
  library(plyr)
  regimen_start <- as.Date(df_therapy_by_subject$regimen_start, format = '%Y-%m-%d')
  deathdate <- as.Date(df_therapy_by_subject$deathdate, format = '%Y-%m-%d')
  lastcontactdate <- as.Date(df_therapy_by_subject$lastcontactdate, format = '%Y-%m-%d')
  line_num <- df_therapy_by_subject$line_num
  maxline <- max(line_num)
  
  ls_tnt_event <- llply(df_therapy_by_subject$line_num, function(line){
    line_idx <- which(line_num == line)
    if(line < maxline){
      ## there is an event characterized by next line of therapy
      nextline_idx <- which(line_num == (line + 1))
      if(length(nextline_idx) != 1){
        #cat('nextline: regimen_start of line', (line + 1), 'is not available for case', df_therapy_by_subject$masterdeid[1], 
        #    'NA will be returned\n')
        tnt <- NA
        event <- NA
      }else{
        tnt <- regimen_start[nextline_idx] - regimen_start[line_idx] + 1
        event <- 1
      }
    }else if(!is.na(deathdate[line_idx])){
      ## there is an event characterized by death
      tnt <- deathdate[line_idx] - regimen_start[line_idx] + 1
      event <- 1
    }else{
      ## censored at last contact date
      tnt <- lastcontactdate[line_idx] - regimen_start[line_idx] + 1
      event <- 0
    }
    z <- data.frame(tnt = as.numeric(tnt), event = event)
  })
  
  tnt_event <- do.call(rbind, ls_tnt_event)
  df_therapy_by_subject <- data.frame(df_therapy_by_subject, tnt_event)
}


################################################################################################################
## Compute tnt and event by single line in batch
## A wrapper function to compute multiple subjects for get_tnt_by_line()
################################################################################################################
get_tnt_df_by_line <- function(df_therapy){
  library(plyr)
  ls_therapy <- dlply(df_therapy, .variables = 'masterdeid', function(x){
    # cat(x$masterdeid[1], '\n')
    x_with_tnt <- get_tnt_by_line(x)
  })
  df_therapy <- do.call(rbind, ls_therapy)
  rownames(df_therapy) <- NULL
  df_therapy
}


################################################################################################################
## Convert drug rule output file format such that columns are drugs in registry_name, rows are masterdeid
## input: 
##      1. recommend_table. The long format of drug recommendations
##      column1, masterdeid; column2, registry_name; column3, recommendation.
##
## output: 
##      1. df_recommend. A matrix, columns are individual drugs in registry_name (n = 100), rows are subjects' 
##      masterdeid.
################################################################################################################
format_recommend <- function(recommend_table){
  # recommend_table <- legacy_table
  library(plyr)
  
  query <- 'select * from cmi_to_registry_agent_lookup'
  reg_cmi_mapping <- get_ryan_query('registry_freeze_18_aug_2015', query)
  # drug_cmi_name <- data.frame(cmi_name = sort(as.character(reg_cmi_mapping$cmi_name)))
  drug_registry_name <- data.frame(registry_name = sort(unique(as.character(reg_cmi_mapping$registry_name))))
  
  recommend_table$prediction_num[recommend_table$recommendation == 'Benefit'] <- 1
  recommend_table$prediction_num[recommend_table$recommendation == 'Lack Of Benefit'] <- -1
  recommend_table$prediction_num[recommend_table$recommendation == 'Indeterminate'
                               |recommend_table$recommendation == 'DoNotReport'] <- 0  

  ls_byPatient_recommend <- dlply(recommend_table, .variables = 'masterdeid', function(masterdeid){
    out <- data.frame(registry_name = masterdeid$registry_name,
                      prediction = masterdeid$prediction_num)
    out <- merge(out, drug_registry_name, by = 'registry_name', all.y = T)
    out$registry_name <- as.character(out$registry_name)
    out <- out[order(out$registry_name), ]
    dim(out)
    ## @@ some drugs are in different drug-rule groups, remove if duplicated for now, assuming there is no conflict
    if(any(duplicated(out$registry_name))){
      idx_dups <- which(duplicated(out$registry_name))
      out <- out[-idx_dups, ]
    }
    out
  })

  ls_recommend <- llply(ls_byPatient_recommend, function(x){
    as.integer(as.character(x$prediction))
  })
  
  df_recommend <- do.call(rbind, ls_recommend)
  df_recommend[is.na(df_recommend)] <- 0
  colnames(df_recommend) <- drug_registry_name$registry_name
  rownames(df_recommend) <- names(ls_byPatient_recommend)
  df_recommend
}



################################################################################################################
## convert AccessionNumber to masterdeid
################################################################################################################
access_to_masterdeid <- function(AccessionNumber){
  query <- 'select ms.accessionnumber, cp.masterpatientid as masterdeid
            from ms1.casepatient cp 
            join ms1.ms1case ms on ms.casepatientid=cp.casepatientid;'
  lookup_table <- get_ryan_query('ms1', query)
  masterdeid <- lookup_table$masterdeid[match(AccessionNumber, lookup_table$accessionnumber)]
}

################################################################################################################
## convert cmi_names to registry_names
## cmi to registry: n to 1 mapping. This conversion will NOT yield conflict.
################################################################################################################
drugname_cmi_to_registry <- function(cmi_names){
  query <- 'select * from cmi_to_registry_agent_lookup'
  reg_cmi_mapping <- get_ryan_query('registry_freeze_18_aug_2015', query)
  registry_name <- reg_cmi_mapping$registry_name[match(cmi_names, reg_cmi_mapping$cmi_name)]
}

################################################################################################################
## convert registry_names to cmi_names
## cmi to registry: n to 1 mapping. This conversion will potentially yield conflict.
################################################################################################################
drugname_registry_to_cmi <- function(registry_names){
  query <- 'select * from cmi_to_registry_agent_lookup'
  reg_cmi_mapping <- get_ryan_query('registry_freeze_18_aug_2015', query)
  cmi_names <- reg_cmi_mapping$cmi_name[match(registry_names, reg_cmi_mapping$registry_name)]
}

################################################################################################################
## Obtain matched and unmatched information for one patient. Need to make sure the naming for drugs are 
## consistent among cmi and registry and ion.
## input: 
##    1. therapy_line, a vector of drugs (registry_name) used in the line of therapy
##    2. recommend_line, a vector of recommendation output from drug rules, names are drugs (registry_name). 
##  -1, lack of benefit; 0, neither; 1, benefit.
##
## rule: if any drug the subject used is recommended as Lack Of Benefit (-1), return unmatched
##       otherwise if any drug the subject used is recommended as Benefit (1), return matched
##       otherwise if all drugs the subject used are recommended as either Indeterminant or DoNotReport (0), 
##  return neither.
################################################################################################################
is_matched <- function(therapy_line, recommend_line){
  # therapy_line <- strsplit(df_therapy[120, ]$therapy, split = ',')[[1]]
  # recommend_line <- df_recommend1[120, ]
  idx <- match(therapy_line, names(recommend_line))
  match_drugs <- recommend_line[idx]
  if(any(match_drugs == -1, na.rm = T)){
    is_matched <- 'unmatched'
  }else if(sum(match_drugs, na.rm = T) == 0){
    is_matched <- 'neither'
  }else{
    is_matched <- 'matched'
  }
  is_matched
}

################################################################################################################
## compute matched unmatched info and add a column to the therapy table. 
## input: 
##    1. df_therapy: rows are individual cases. It must have a column of 'therapies' that contains drugs used in 
##  this line, collapsed by ','.
##    2. df_recommend: It MUST be a matrix!(colnames of a dataframe will be corrupted). Rows are individual cases 
##  ordered in the same way as df_therapy. rownames is masterdeid, colnames are registry_name of the drugs.
################################################################################################################
get_matched_df <- function(df_therapy, df_recommend){
  library(plyr)
  n <- nrow(df_therapy)
  followup <- unlist(llply(1:n, function(i){
    # cat(i, ' ')
    therapy <- strsplit(as.character(df_therapy$therapies[i]), ',')[[1]]
    recommend <- df_recommend[i,]
    names(recommend) <- colnames(df_recommend)
    is_matched(therapy, recommend)
  }))
  
  df_therapy <- data.frame(df_therapy, 
                           followed = followup)
}

################################################################################################################
## To compute matched/unmatched/neither for a single case. An improved version of is_matched(). Match on 
## druggroup instead of registry_name. 
## input:
##      1. received <- c('druggroup1', 'druggoup2', ...) 
##      2. recommend <- data.frame(druggroup, recommendation)
## rule: if any drug the subject used is recommended as Lack Of Benefit (-1), return unmatched
##       otherwise if any drug the subject used is recommended as Benefit (1), return matched
##       otherwise if all drugs the subject used are recommended as either Indeterminant or DoNotReport (0), 
##  return neither.
################################################################################################################
get_match_binary <- function(received, recommend){
  #received <- c('F', 'G')
  #recommend <- data.frame(druggroup = c('A', 'B', 'C', 'D', 'E'),
  #                        recommendation = c(1,1,0,-1,0))
  
  match_idx1 <- match(received, recommend$druggroup)
  recommendation <- recommend$recommendation[match_idx1]
  if(any(recommendation < 0, na.rm = T)){
    is_matched <- 'unmatched'
  }else if(sum(recommendation, na.rm = T) == 0){
    is_matched <- 'neither'
  }else{
    is_matched <- 'matched'
  }
  is_matched
}

################################################################################################################
## To compute match_score for a single case by matching weights and recommendations to received drugs. weights 
## were normalized within received drugs.
## input:
##      1. received <- c('druggroup1', 'druggoup2', ...)
##      2. recommend <- data.frame(druggroup, recommendation)
##      3. wt_star <- data.frame(druggroup, weight)
##      4. adjust_recommend = c(T, F) # whether to adjust the weight by benefit driven or lackOfBenefit driven.
##      5. if adjust_recommend == T, recom_weight = c(lackofbenefit, benefit) #lackofbenefit + benefit = 1
## output:
##      match_score. A numeric variable range [-1, 1]
################################################################################################################
get_match_score <- function(received, recommend, wt_star, adjust_recommend, recom_weight){
  #received <- df_therapy_2nd$druggroups[1]
  match_idx1 <- match(received, recommend$druggroup)
  match_idx2 <- match(received, wt_star$druggroup)
  if(length(match_idx1) == 0 |length(match_idx2) == 0){
    recommendation <- NA
    weight <- NA
  }else{
    recommendation <- recommend$recommendation[match_idx1]
    weight <- wt_star$weight[match_idx2]
  }
  
  if(adjust_recommend){
    weight <- sapply(1:length(weight), function(i){
      if(is.na(recommendation[i])){
        z <- NA
      }else if(recommendation[i] < 0){
        z <- weight[i] * recom_weight[1] 
      }else{
        z <- weight[i] * recom_weight[2]
      }
      z
    })
  }
  weight <- convert_weight(weight)
  match_score = sum(recommendation * weight, na.rm = T)
  match_score
}


################################################################################################################
## A wrapper function to compute match_score in batches using get_match_score(). Only suitable for a single 
## version of drugrule at a single version of weight combinations. A bigger wrapper is needed to run different 
## combinations. It also scale the match_scores from [-1, 1] to [0, 1].
## IMPORTANT!! match the order between df_therapy and ls_recommend for input.
## input:
##      1. df_therapy. contains a column named 'druggroups'. 
##      2. ls_recommend. contains drug rule recommendations organized by masterdeid. It MUST be in the same order
##      as df_therapy. ls_recommend should contain columns, named 'drugRule', 'recom'.
##      3. method = c('binary', 'score')
##      4. if method == 'score', adjust_recommend and recom_weight are inherited from get_match_score()
## output:
##      df_therapy. the same data.frame as input with a new column named 'match_scores' appended.
################################################################################################################
get_match_score_df <- function(df_therapy, ls_recommend, method, wt_star, adjust_recommend, recom_weight){
  # ls_recommend <- ls_recommend[match(inter_masterdeid, names(ls_recommend))]
  # df_therapy <- df_therapy_2nd[match(inter_masterdeid, df_therapy_2nd$masterdeid), ]
  # wt_star <- data.frame(druggroup = colnames(weight_matrix),
  #                       weight = convert_weight(c(0.3,0.3,0.6,0.3,0.3)), stringsAsFactors = F)
  
  n_sample <- length(ls_recommend)
  match_scores <- llply(1:n_sample, function(i){
    # cat(i,'\n')
    received <- strsplit(df_therapy$druggroups[i], ',')[[1]]
    recommend <- data.frame(druggroup = ls_recommend[[i]]$drugRule,
                            recommendation = ls_recommend[[i]]$recom, 
                            stringsAsFactors = F)
    if(method == 'score'){
      match_score <- get_match_score(received, recommend, wt_star, adjust_recommend = adjust_recommend, recom_weight = recom_weight)
      match_score <- (match_score - -1)/2 ## recommendation score is between -1 and 1, scale to 0 and 1.
    }else if(method == 'binary'){
      match_score <- get_match_binary(received, recommend)
    }
  })
  df_therapy$match_score <- unlist(match_scores)
  df_therapy
}

################################################################################################################
## text processing functions to clean up clinical information
################################################################################################################
## trim off heading/trailing spaces and/or punctuations from string
trim <- function (x) gsub("^(\\s|[[:punct:]])+|([[:punct:]]|\\s)+$", "", x)

## clean up unknowns
convert_unkowns <- function(x) gsub('NA|NULL|(^[Uu]nknown.*$)|(^[Nn]ot.*$)', 'Unknown', x)

## clean up grade
convert_grade <- function(x) {
  x <- gsub('\\s*/\\s*', '/', trim(x))
  x <- convert_unkowns(x)
  x
}

## clean up stage
## consolidate I, II, and IV; only remove trailing numbers for III
convert_stage <- function(x, collapseIII = T) {
  x <- gsub('^I[A-D1-9]+$', 'I', x)
  x <- gsub('^II[A-D1-9]+$', 'II', x)
  if(collapseIII){
    x <- gsub('^III[A-D1-9]+$', 'III', x)
  }else if(any(grepl('^III$', x))){
    cat('detected \'III\', set collapseIII to TRUE\n')
    x <- gsub('^III[A-D1-9]+$', 'III', x)
  }else{
    x <- ifelse(grepl('^III[A-D1-9]+$', x), gsub('[1-9]+$', '', x), x)
  }
  x <- gsub('^IV[A-D1-9]+$', 'IV', x)
  x <- convert_unkowns(x)
  x
}


################################################################################################################
## generate the survfit and coxph objects
################################################################################################################
surv_fit <- function(df_clinical_info, plot_covariate){
  library(survival)
  # df_clinical_info <- df_plot1
  # plot_covariate <- 'followed'
  is_unknown <- any(sapply(df_clinical_info, function(x) grepl('Unknown', x)))
  if(is_unknown){
    stop('unknowns detected, please remove the corresponding samples and rerun!\n')
  }
  ## with covariates intended for adjustment
  fitCox <- coxph(Surv(time, event) ~ ., data = df_clinical_info)
  ## KM plot just for 'followed'
  surv_formula <- as.formula(paste('Surv(time, event) ~ ', plot_covariate, sep = ''))
  fitSurv <- survfit(surv_formula, data = df_clinical_info)
  ## need to pass colid_covariate to the kmplots function to locate HR and p-values from the coef table
  colid_covariate <- which(colnames(df_clinical_info) == plot_covariate)
  z <- list(fitCox = fitCox, 
            fitSurv = fitSurv, 
            colid_covariate = colid_covariate)
}

################################################################################################################
## plot KM curves using the surv_fit output
## input: 
##      1. surv_fit, output from surv_fit(), contains coxph(), survfit(), as well as colid_covariate the column
##      index of the target covariate
##      2. display_coef_idx, when levels of target covariates are more than two, pick the desired contrast to
##      display by inputing the index of them from the summary(fitCox)$coef table
################################################################################################################
kmplots <- function(
  surv_fit, main, display_coef_idx = NULL, simplify = F, drawSegment = T, xlim = NULL
){
  library(survival)
  fitCox <- surv_fit$fixCox
  fitSurv <- surv_fit$fitSurv
  colid_covariate <- surv_fit$colid_covariate
  fitCox_sum <- summary(surv_fit$fitCox)
  fitSurv_sum <- summary(surv_fit$fitSurv)
  
  ## prepare for KM plot
  line_colors <- c('red', 'blue', 'green', 'pink', 'orange', 'purple', 'cyan')
  surv_strata <- toupper(sapply(names(fitSurv$strata), function(x)strsplit(x, '=')[[1]][2]))
  n_levels <- length(fitSurv$strata)
  if(n_levels > 7) stop('plot_covaraite has more than 7 levels, plot will be ugly, please check covariates!\n')
  if(n_levels > 2 & is.null(display_coef_idx)){
    stop('more than two levels in the selected covariate, please manually input the row numbers of the desired coeffients to display by setting display_coef_idx!\n')
  }else if(!is.null(display_coef_idx)){
    rowid <- display_coef_idx
  }else{
    rowid <- colid_covariate - 2
  }
  
  ## start plotting
  par(mar = c(4,4,5,2))
  plot(
    fitSurv, col = line_colors[1:n_levels], lwd = 4, 
    mark.time = T,
    xlab = 'Time', ylab = 'Proportion Event Free', 
    main = main, 
    cex.main = 1, cex = 1, cex.lab = 1,
    xaxs = 'i', yaxs = 'i', # remove extra space between plot and axis
    xlim = xlim
  )
  if(simplify == F){
    legend('topright', 
           legend = c(paste0(surv_strata, ' (n = ', fitSurv$n, ', event = ', fitSurv_sum$table[, 4], ')'),
                      paste0(
                        'HR = ', round(fitCox_sum$coef[rowid,2], 3),
                        ' (', round(fitCox_sum$conf.int[3], 3), '-',  round(fitCox_sum$conf.int[4], 3),')'
                        #', p-value =', round(fitCox_sum$coef[rowid,5], 3)
                      ),
                      paste0('log-rank test p-value =', round(fitCox_sum$sctest[3], 3))),
           col = c(line_colors[1:n_levels], rep('white', length(rowid) + 1)), lty = 1, lwd = 2, bty='n', cex = 0.8)
    ## add segments and text of median survival time
    x_del_position <- (-1)^(1:7) * 0.1
    # y_del_position <- floor(((1:7)-1)/2) * 0.1
    y_del_position <- (1:7) * 0.1 - 0.1
    idx = which(colnames(fitSurv_sum$table) == 'median')
    lapply(1:n_levels, function(x){
      segments(fitSurv_sum$table[x,idx], -0.1, fitSurv_sum$table[x,idx], 0.5, col = line_colors[x], lty = 2, lwd = 2)
      text(fitSurv_sum$table[x,idx] + x_del_position[x], y_del_position[x], round(fitSurv_sum$table[x,idx], 2), col = line_colors[x])
    })
    segments(0, 0.5, max(fitSurv_sum$table[,idx]), 0.5, col = 'black', lty = 2, lwd = 2)
  }else{
    if(drawSegment){
      ## add segments and text of median survival time
      x_del_position <- (-1)^(1:7) * 0.1
      # y_del_position <- floor(((1:7)-1)/2) * 0.1
      y_del_position <- (1:7) * 0.1 - 0.1
      idx = which(colnames(fitSurv_sum$table) == 'median')
      lapply(1:n_levels, function(x){
        segments(fitSurv_sum$table[x,idx], -0.1, fitSurv_sum$table[x,idx], 0.5, col = line_colors[x], lty = 2, lwd = 2)
      })
      segments(0, 0.5, max(fitSurv_sum$table[,idx]), 0.5, col = 'black', lty = 2, lwd = 2)
    }
  }
}


################################################################################################################
## Demographic table generator
## input:
##      1. dat, a dataframe with all demographic covariates to summarize, column names will be
##    displayed in the final table, pick your favourite as input
##      2. y_covariate, the target covariate to split as columns, must be two levels
##      3. x_covariate, a vector of demographic covariates of interest
##      4. category_test, whether to carry out proportion test on each category
## output:
##      a demographic table
##example function usage:
## get_demogr_table(regova,"firstline_aftercollection",c("Age","Grade","Stage","Race"),TRUE)
################################################################################################################
get_demogr_tina <- function(dat, y_covariate, x_covariate, category_test=TRUE){
  
  # dat=regova
  # y_covariate="firstline_aftercollection"
  # x_covariate=c("agecat","grade","stage","racecat")
  # category_test=TRUE
  
  #library(xtable)
  library(plyr)
  #y_idx <- which(colnames(df_demographics) == y_covariate)
  y_cov <- dat[, y_covariate]
  x_cov <- dat[,x_covariate]
  n_y_level <- nlevels(y_cov)
  nColumns <- n_y_level + 2
  
  ##
  n_table <- table(y_cov)
  ## process numeric variables
  x_num_idx <- which(sapply(x_cov, is.numeric))
  if (length(x_num_idx)!=0) {
    x_num <- data.frame(x_cov[, x_num_idx])
    x_num_names <- colnames(x_cov)[x_num_idx]
    colnames(x_num) <- x_num_names
    n_x_num <- ncol(x_num)
    
    if(n_y_level < 2){
      stop('y_covariate must have at least two levels to compute the p values.\n')
    }else if(n_y_level == 2){
      num_p <- round(sapply(1:n_x_num, function(x_idx){
        #wilcox.test(x_num[, x_idx] ~ y_cov)$p.value
        t.test(x_num[, x_idx] ~ y_cov)$p.value
      }), 4)
    }else{
      num_p <- round(sapply(1:n_x_num, function(x_idx){
        aov_sum <- summary(aov(x_num[, x_idx] ~ y_cov))
        aov_sum[[1]]$`Pr(>F)`[1]
      }))
    }
    
    ls_num_tables <- lapply(1:n_x_num, function(x_idx){
      #tab <- tapply(x_num[,x_idx], y_cov, median)
      #iqr_low <- tapply(x_num[,x_idx], y_cov, function(x)quantile(x, 0.25))
      #iqr_high <- tapply(x_num[,x_idx], y_cov, function(x)quantile(x, 0.75))
      #tab <- paste0(tab, '(', iqr_low, '-', iqr_high, ')')
      
      num_mean <- tapply(x_num[,x_idx], y_cov, mean, na.rm = T)
      num_sd <- tapply(x_num[,x_idx], y_cov, sd, na.rm = T)
      num_mean <- paste(num_mean, '(',num_sd ,')')
    })
    names(ls_num_tables) <- x_num_names
    x_num_table <- do.call(rbind, ls_num_tables)
    x_num_table <- cbind(x_num_table, num_p)
  }
  
  if (ncol(x_cov)>length(x_num_idx)) {
    ## process categorical variables
    if (length(x_num_idx)==0) {
      x_cat <- data.frame(x_cov)
      x_cat_names <- colnames(x_cov)
      colnames(x_cat) <- x_cat_names
    } else {
      x_cat <- data.frame(x_cov[, -x_num_idx])
      x_cat_names <- colnames(x_cov[,-x_num_idx])
      colnames(x_cat) <- x_cat_names
    }
    
    n_x_cat <- ncol(x_cat)
    ls_cat_tables <- lapply(1:n_x_cat, function(x_idx){
      tab1 <- table(x_cat[,x_idx], y_cov)
      tab2 <- rbind(NA, tab1)
      rownames(tab2)[1] <- x_cat_names[x_idx]
      tab2
    })
    fisher_p <- sapply(ls_cat_tables, function(x) {
      p_value <- tryCatch(fisher.test(x[-1,], workspace = 2e8)$p.value, error = function(e)return(NA))
      p_value <- round(p_value, 4)
    })
    
    ls_cat_tables <- lapply(1:n_x_cat, function(x){
      nx <- nrow(ls_cat_tables[[x]])
      if (category_test) {
        prop_p <- sapply(2:nx, function(x_row) {
          p_value <- tryCatch(prop.test(ls_cat_tables[[x]][x_row,],n_table)$p.value, error = function(e)return(NA))
          p_value <- round(p_value, 4)
        })
        out <- cbind(ls_cat_tables[[x]], c(fisher_p[x], prop_p))
        out <- rbind(NA, out)
        out
      } else {
        out <- cbind(ls_cat_tables[[x]], c(fisher_p[x], rep(NA, nx-1)))
        out <- rbind(NA, out)
        out
      }
    })
    
    x_cat_table_tmp <- do.call(rbind, ls_cat_tables)
    
    x_cat_table <- sapply(1:(ncol(x_cat_table_tmp)-1), function(x){
      paste0(x_cat_table_tmp[, x], '(', round(x_cat_table_tmp[, x]/n_table[x]*100, 1), '%)')
    })
    x_cat_table <- cbind(x_cat_table, x_cat_table_tmp[,ncol(x_cat_table_tmp)])
    x_cat_table[grep('NA', x_cat_table)] <- NA
    
  }
  
  n_table <- c('n', paste0(n_table), NA)
  x_num_table <- cbind(rownames(x_num_table), x_num_table)
  x_cat_table <- cbind(rownames(x_cat_table_tmp), x_cat_table)
  x_table <- rbind(n_table, x_num_table, x_cat_table)
  x_table <- rbind(c(NA, levels(y_cov), 'P-value'), x_table)
  rownames(x_table) <- NULL
  colnames(x_table) <- NULL
  x_table[grep('^$', x_table)] <- NA
  
  return(x_table)
  #align <- paste0(c('ll', rep('c', ncol(x_table)-1)), collapse = '')
  #print(xtable(x_table, align=align, caption = caption), include.rownames=FALSE, include.colnames=FALSE,hline.after=c(0, 2, nrow(x_table)), file = outfile, caption.placement = 'top', size = 'scriptsize')
  
}



################################################################################################################
## Format individual subject for drug time line plotting
################################################################################################################
format_drug_line <- function(df_therapy_by_subject, unit = 'month'){
  # df_therapy_by_subject <- ls_drug_by_subject[[13]]
  df_therapy_by_subject$drugstart <- as.Date(df_therapy_by_subject$drugstart, format = '%Y-%m-%d')
  df_therapy_by_subject$drugend <- as.Date(df_therapy_by_subject$drugend, format = '%Y-%m-%d')
  df_therapy_by_subject$collectiondate <- as.Date(df_therapy_by_subject$collectiondate, format = '%Y-%m-%d')
  df_therapy_by_subject$determinationdate <- as.Date(df_therapy_by_subject$determinationdate, format = '%Y-%m-%d')
  
  if(unit == 'month'){
    unitToDivide = 30
  }else if(unit == 'day'){
    unitToDivide = 1
  }
  
  df_therapy_by_subject <- data.frame(
    masterdeid = df_therapy_by_subject$masterdeid,
    registry_name = df_therapy_by_subject$registry_name,
    collectionPoint = 0,
    drugStartPoint = as.numeric(df_therapy_by_subject$drugstart - df_therapy_by_subject$collectiondate)/unitToDivide,
    drugEndPoint = as.numeric(df_therapy_by_subject$drugend - df_therapy_by_subject$collectiondate)/unitToDivide,
    drugDuration = as.numeric(df_therapy_by_subject$drugend - df_therapy_by_subject$drugstart)/unitToDivide,
    determinationPoint = as.numeric(df_therapy_by_subject$determinationdate - df_therapy_by_subject$collectiondate)/unitToDivide
  )
  df_therapy_by_subject <- arrange(df_therapy_by_subject, drugStartPoint, drugEndPoint)
  df_therapy_by_subject$unknownDuration <- ifelse(is.na(df_therapy_by_subject$drugDuration) | df_therapy_by_subject$drugDuration == 0,
                                                  T, F)
  
  df_therapy_by_subject$plotEndPoint <- ifelse(df_therapy_by_subject$unknownDuration == T,
                                               df_therapy_by_subject$drugStartPoint + 15/unitToDivide,
                                               df_therapy_by_subject$drugEndPoint)
  levels_order = df_therapy_by_subject$registry_name[order(df_therapy_by_subject$drugStartPoint, 
                                                           df_therapy_by_subject$drugEndPoint,
                                                           decreasing = T)]
  levels_order <- levels_order[which(!duplicated(levels_order, fromLast = T))]
  df_therapy_by_subject$registry_name <- factor(
    df_therapy_by_subject$registry_name,
    levels = levels_order)
  df_format <- df_therapy_by_subject
}

################################################################################################################
## plot the drug time line using the output from format_drug_line() as input
################################################################################################################
plot_drug_line <- function(df_format){
  library(ggplot2)
  plot_drugs <- ggplot(df_format, aes(colour = unknownDuration)) + 
    geom_segment(aes(x = drugStartPoint, xend = plotEndPoint, y = registry_name, yend = registry_name), 
                 size = 5) +
    labs(x = 'Time (months)', y = NULL, 
         title = paste0('masterdeid: ', df_format$masterdeid[1])) + 
    geom_vline(xintercept = 0, color = 'purple') + 
    geom_vline(xintercept = df_format$determinationPoint[1], color = 'black') + 
    theme(title = element_text(size = 20), axis.title = element_text(size = 20), axis.text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid = element_line(size = 1),
          legend.position='none')  
  z <- list(plot_drugs = plot_drugs)
}



################################################################################################################
## ovarian cancer registry ion data process appended functions
################################################################################################################


################################################################################################################
## To generate all combinations of weight matrix for given options and druggroups
## input:
##      1. weight_options <- c(1, 0.6, 0.3)
##      2. druggroup <- c('dg1', 'dg2', 'dg3', 'dg4', 'dg5')
## output: 
##      a matrix. columns are druggroup names, rows are weight combinations.
################################################################################################################
generate_weight_matrix <- function(weight_options, druggroup){
  library(gtools)
  num_values <- length(weight_options)
  num_drugs <- length(druggroup)
  wt_matrix <- permutations(n = num_values, r = num_drugs, weight_options, repeats.allowed = T)
  colnames(wt_matrix) <- druggroup
  cat(num_drugs, 'drug groups\n')
  cat(num_values, 'values for each drug group\n')
  cat(nrow(wt_matrix), 'combinations\n')
  wt_matrix
}


################################################################################################################
## weight conversion to sum to 1
################################################################################################################
convert_weight <- function(wt){
  wt_star <- wt/sum(wt, na.rm = T)
}


################################################################################################################
## convert masterpatientid to masterdeid
################################################################################################################
convert_masterdeid <- function(masterpatientid){
  query <- 'select * from ion_reg_ova_case_lst;'
  masterdeid_mapping <- get_sting_query('bioinfo_projects', query)
  masterdeid <- masterdeid_mapping$masterdeid[match(masterpatientid, masterdeid_mapping$masterpatientid)]
}

################################################################################################################
## convert registry_names to druggroups
## mapping_table <- data.frame(registry_name, druggroup_name)
################################################################################################################
convert_registry_to_druggroup <- function(registry_names, mapping_table){
  #registry_names <- strsplit(df_therapy_ken$therapies[1], ',')[[1]]
  #mapping_table <- reg_druggroup_mapping
  druggroup <- as.character(mapping_table$druggroup_name[match(registry_names, mapping_table$registry_name)])
  #if(any(is.na(druggroup))){
  #  druggroup <- druggroup[-which(is.na(druggroup))]
  #}
}


################################################################################################################
## Followings are archived functions. DO NOT use
################################################################################################################




################################################################################################################
## Convert drug rule output file format such that columns are drugs in cmi_name, rows are AccessionNumber
## @@@ Important to note that there are cases where a single drug recommendation might be derived from two or more 
## @@@ independent drug rules, which might lead to conflict recommendation. Need to resolve this later.
################################################################################################################
format_recommend_akiv <- function(drug_rule_out){
  library(plyr)
  
  query <- 'select * from cmi_to_registry_agent_lookup'
  reg_cmi_mapping <- get_ryan_query('registry_freeze_18_aug_2015', query)
  drug_cmi_name <- data.frame(cmi_name = sort(as.character(reg_cmi_mapping$cmi_name)))
  #drug_registry_name <- data.frame(registry_name = sort(as.character(reg_cmi_mapping$registry_name)))
  
  # colnames(drug_rule_out) <- c('AccessionNumber', 'drug_group', 'xx', 'xxx', 'prediction', 'marker', 'cmi_name', 'xxxx', 'xxxxx', 'xxxxxx')
  drug_rule_out$prediction_num <- drug_rule_out$prediction
  drug_rule_out$prediction_num[drug_rule_out$prediction == 'Benefit'] <- 1
  drug_rule_out$prediction_num[drug_rule_out$prediction == 'Lack Of Benefit'] <- -1
  drug_rule_out$prediction_num[drug_rule_out$prediction == 'Indeterminate'
                               |drug_rule_out$prediction == 'DoNotReport'] <- 0
  
  ls_byPatient_recommend <- dlply(drug_rule_out, .variables = 'AccessionNumber', function(AccessionNumber){
    
    out <- data.frame(cmi_name = AccessionNumber$cmi_name,
                      prediction = AccessionNumber$prediction_num)
    out <- merge(out, drug_cmi_name, by = 'cmi_name', all.y = T)
    out$cmi_name <- as.character(out$cmi_name)
    out <- out[order(out$cmi_name), ]
    ## @@ some drugs are in different drug-rule groups, remove duplicated for now, assuming there is no conflict
    if(any(duplicated(out$cmi_name))){
      idx_dups <- which(duplicated(out$cmi_name))
      out <- out[-idx_dups, ]
    }
    out
  })
  
  ls_recommend <- llply(ls_byPatient_recommend, function(x){
    as.integer(as.character(x$prediction))
  })
  
  df_recommend <- do.call(rbind, ls_recommend)
  df_recommend[is.na(df_recommend)] <- 0
  colnames(df_recommend) <- drug_cmi_name$cmi_name
  rownames(df_recommend) <- names(ls_byPatient_recommend)
  df_recommend
}


################################################################################################################
## Demographic table generator using xtable
## input:
##      1. df_demographics, a dataframe with all demographic covariates to summarize, column names will be
##    displayed in the final table, pick your favourite as input
##      2. y_covariate, the target covariate to split as columns, must be at least two levels
##      3. x_covariate, a vector of demographic covariates of interest
##      4. category_test, whether to carry out proportion test on each category
##      5. caption, table name
##      6. outfile, output directory and filename 
################################################################################################################
get_demogr_ryan <- function(
  df_demographics, y_covariate, x_covariate, category_test=TRUE,
  caption = caption,
  outfile = ''
){
  # caption <- 'my table'
  # df_demographics <- df_plot4[, -c(1:2)]
  # df_demographics <- data.frame(age1 = sample(df_plot4$age),
  #                              df_demographics)
  # df_demographics <- df_demographics[, -(4:5)]
  # y_covariate <- 'followed'
  library(xtable)
  library(plyr)
  y_idx <- which(colnames(df_demographics) == y_covariate)
  y_cov <- df_demographics[, y_idx]
  x_cov <- df_demographics[, x_covariate]
  n_y_level <- nlevels(y_cov)
  nColumns <- n_y_level + 2
  
  ##
  n_table <- table(y_cov)
  ## process numeric variables
  x_num_idx <- which(sapply(x_cov, is.numeric))
  x_num <- data.frame(x_cov[, x_num_idx])
  x_num_names <- colnames(x_cov)[x_num_idx]
  colnames(x_num) <- x_num_names
  n_x_num <- ncol(x_num)
  
  if(n_y_level < 2){
    stop('y_covariate must have at least two levels to compute the p values.\n')
  }else if(n_y_level == 2){
    num_p <- round(sapply(1:n_x_num, function(x_idx){
      wilcox.test(x_num[, x_idx] ~ y_cov)$p.value
    }), 4)
  }else{
    num_p <- round(sapply(1:n_x_num, function(x_idx){
      aov_sum <- summary(aov(x_num[, x_idx] ~ y_cov))
      aov_sum[[1]]$`Pr(>F)`[1]
    }))
  }
  
  ls_num_tables <- lapply(1:n_x_num, function(x_idx){
    tab <- tapply(x_num[,x_idx], y_cov, median, na.rm = T)
    iqr_low <- tapply(x_num[,x_idx], y_cov, function(x)quantile(x, 0.25, na.rm = T))
    iqr_high <- tapply(x_num[,x_idx], y_cov, function(x)quantile(x, 0.75, na.rm = T))
    tab <- paste0(tab, '(', iqr_low, '-', iqr_high, ')')
  })
  names(ls_num_tables) <- x_num_names
  x_num_table <- do.call(rbind, ls_num_tables)
  x_num_table <- cbind(x_num_table, num_p)
  
  ## process categorical variables
  x_cat <- data.frame(x_cov[, -x_num_idx])
  x_cat_names <- colnames(x_cov[-x_num_idx])
  colnames(x_cat) <- x_cat_names
  n_x_cat <- ncol(x_cat)
  ls_cat_tables <- lapply(1:n_x_cat, function(x_idx){
    tab1 <- table(x_cat[,x_idx], y_cov)
    tab2 <- rbind(NA, tab1)
    rownames(tab2)[1] <- x_cat_names[x_idx]
    tab2
  })
  fisher_p <- sapply(ls_cat_tables, function(x) {
    p_value <- tryCatch(fisher.test(x[-1,], workspace = 2e8)$p.value, error = function(e)return(NA))
    p_value <- round(p_value, 4)
  })
  ls_cat_tables <- lapply(1:n_x_cat, function(x){
    nx <- nrow(ls_cat_tables[[x]])
    if (category_test) {
      prop_p <- sapply(2:nx, function(x_row) {
        p_value <- tryCatch(prop.test(ls_cat_tables[[x]][x_row,],n_table)$p.value, error = function(e)return(NA))
        p_value <- round(p_value, 4)
      })
      out <- cbind(ls_cat_tables[[x]], c(fisher_p[x], prop_p))
      out <- rbind(NA, out)
      out
    } else {
      out <- cbind(ls_cat_tables[[x]], c(fisher_p[x], rep(NA, nx-1)))
      out <- rbind(NA, out)
      out
    }
  })
  x_cat_table_tmp <- do.call(rbind, ls_cat_tables)
  x_cat_table <- sapply(1:(ncol(x_cat_table_tmp)-1), function(x){
    paste0(x_cat_table_tmp[, x], '(', round(x_cat_table_tmp[, x]/n_table[x]*100, 1), '%)')
  })
  x_cat_table <- cbind(x_cat_table, x_cat_table_tmp[,ncol(x_cat_table_tmp)])
  x_cat_table[grep('NA', x_cat_table)] <- NA
  
  n_table <- c('n', paste0('(', n_table, ')'), NA)
  x_num_table <- cbind(rownames(x_num_table), x_num_table)
  x_cat_table <- cbind(rownames(x_cat_table_tmp), x_cat_table)
  x_table <- rbind(n_table, x_num_table, x_cat_table)
  x_table <- rbind(c(NA, levels(y_cov), 'P-value'), x_table)
  rownames(x_table) <- NULL
  colnames(x_table) <- NULL
  x_table[grep('^$', x_table)] <- NA
  
  align <- paste0(c('ll', rep('c', ncol(x_table)-1)), collapse = '')
  print(xtable(x_table, align=align, caption = caption), 
        include.rownames=FALSE, include.colnames=FALSE,
        hline.after=c(0, 2, nrow(x_table)), 
        file = outfile, 
        #tabular.environment="longtable",
        caption.placement = 'top', size = 'scriptsize')
  
}






