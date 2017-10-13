rm(list = ls())

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
## obtain all COE_NAME
################################################################################################################
get_coe_names <- function(){
  query <- "
  SELECT DISTINCT
    COE_NAME
  FROM
  bioinfo_projects.coe_summary_case_lst"
  coe_names <- get_sting_query('bioinfo_projects', query)
  coe_names <- coe_names$COE_NAME
}

################################################################################################################
## obtain number of patients from coe site and physician
################################################################################################################
get_site_physician_summary <- function(coe_name, summary_cutoff = 20){
  query <- sprintf("
  SELECT DISTINCT
    OrderingPhysician, COUNT(*) AS 'n_sample'
  FROM
  bioinfo_projects.coe_summary_case_lst cscl
  WHERE COE_NAME = '%s'
  GROUP BY OrderingPhysician;
  ", coe_name)
  coe_physician_tab <- get_sting_query('bioinfo_projects', query)
  coe_physician_tab$Physician <- ifelse(
    coe_physician_tab$n_sample > summary_cutoff, 
    coe_physician_tab$OrderingPhysician,
    'Others'
  )
  coe_physician_tab <- data.frame(
    OrderingPhysician = coe_physician_tab$OrderingPhysician,
    Physician = coe_physician_tab$Physician,
    n_sample = coe_physician_tab$n_sample
  )
}

################################################################################################################
## obtain test results for each COE site
################################################################################################################
get_test_result <- function(coe_name, summary_cutoff){
  sqlcmd <- sprintf(
    "SELECT 
    *
    FROM
    bioinfo_projects.coe_summary_test cst
    JOIN
    (SELECT 
    AccessionNumber, Lineage, OrderingPhysician
    FROM
    bioinfo_projects.coe_summary_case_lst
    WHERE
    COE_NAME = '%s') x using(AccessionNumber);",
    coe_name
  )
  df_test <- get_sting_query('bioinfo_projects', sqlcmd)
  coe_physician_tab <- get_site_physician_summary(coe_name, summary_cutoff = summary_cutoff)
  df_test <- merge(df_test, coe_physician_tab, by = 'OrderingPhysician', all.x = T)
}

################################################################################################################
## summarize the frequency of Positive result by each lineage
## type = c('SEQ', 'CNV')
################################################################################################################
get_ngs_summary_by_case <- function(df_cmi_by_site, type){
  library(dplyr)
  library(reshape2)
  #type <- 'CNV'
  #df_cmi_by_site <- df_test
  df_cmi_by_site$Lineage[is.na(df_cmi_by_site$Lineage)] <- 'NA'
  df_cmi_by_site <- df_cmi_by_site[df_cmi_by_site$Tech == type, ]
  Test_summary <- df_cmi_by_site %>% 
    group_by(AccessionNumber, BiomarkerName, Tech, Lineage) %>% 
    summarize(
      #len = length(AccessionNumber), 
      unique_normResult = ifelse(any(NormResult == 'Positive'), 'Positive', NormResult)
    ) %>%
    group_by(BiomarkerName, Tech, Lineage) %>%
    summarize(
      Positive_count = sum(unique_normResult == 'Positive')
    ) %>%
    group_by(Lineage) %>%
    arrange(desc(Positive_count)) %>%
    slice(1:10) %>%
    group_by(Lineage) %>%
    mutate(len = length(Lineage)) %>%
    mutate(
      Test = if(len[1] < 10){
        paste0('TOP', paste0('0', 1:length(Lineage)))
      }else{
        paste0('TOP', c(paste0('0', 1:9), 10))
      },
      Value = paste0(BiomarkerName, '[', Positive_count, ']')
    )
  Test_summary <- Test_summary[, -which(colnames(Test_summary) %in% c('Positive_count', 'BiomarkerName', 'len'))]
  colnames(Test_summary)[3] <- 'BiomarkerName'
  Test_summary
  #out_summary <- dcast(Test_summary, Tech + BiomarkerName ~ Lineage, value.var = 'Top10')
}


################################################################################################################
## summarize the frequency of Positive result by each physician and each lineage
################################################################################################################
get_other_summary_by_case <- function(df_cmi_by_physician){
  library(dplyr)
  
  #df_cmi_by_physician <- df_test[df_test$Physician == 'Philip Philip', ]
  
  ## total number of patient by lineage is calculated before filtering by technology
  df_cmi_by_physician$Lineage[is.na(df_cmi_by_physician$Lineage)] <- 'NA'
  
  total_by_lineage <- df_cmi_by_physician %>%
    group_by(Lineage) %>%
    summarize(N = length(unique(AccessionNumber)))
  total_by_lineage <- rbind(total_by_lineage, c('N_test', sum(total_by_lineage$N)))
  
  exclusion_tech_list <- c('SEQ', 'CNV', 'Microarray', 'RT-PCR', 'qPCR')
  df_cmi_by_physician <- df_cmi_by_physician[!df_cmi_by_physician$Tech %in% exclusion_tech_list, ]
  
  total_by_test <- df_cmi_by_physician %>%
    group_by(Tech, BiomarkerName) %>%
    summarize(N = length(unique(AccessionNumber)))
  
  Test_summary <- df_cmi_by_physician %>% 
    group_by(AccessionNumber, BiomarkerName, Tech, Lineage) %>% 
    summarize(
      #len = length(AccessionNumber), 
      unique_normResult = ifelse(any(NormResult == 'Positive'), 'Positive', NormResult)
    ) %>%
    group_by(Tech, Lineage, BiomarkerName) %>%
    summarize(
      Value = sum(unique_normResult == 'Positive')
    )
  
  z <- list(Test_summary = Test_summary, 
            total_by_test = total_by_test,
            total_by_lineage = total_by_lineage)
}

################################################################################################################
## format the summary results tables
################################################################################################################
format_out_summary <- function(out_summary, total_by_test, total_by_lineage){
  tech_order <- c('NGS', 'CNV', 'Fusion', 'IHC', 'ISH', 'Methylation', 'Rearrangement', 'FA')
  format_summary <- c()
  
  ## format rows
  for(tech in tech_order){
    if(tech %in% out_summary$Tech){
      tech_summary <- out_summary[out_summary$Tech == tech, ]
      format_summary <- rbind(format_summary, tech_summary)
    }
  }
  
  ## format columns, put 'None Of These Apply' and 'NA' to the end of the columns
  if('None Of These Apply' %in% colnames(format_summary)){
    format_summary <- cbind(
      format_summary[, -which(colnames(format_summary) == 'None Of These Apply')],
      format_summary[, c('None Of These Apply')]
    )
    colnames(format_summary)[ncol(format_summary)] <- 'None Of These Apply'
    format_summary$`None Of These Apply` <- as.character(format_summary$`None Of These Apply`)
  }
  
  if('NA' %in% colnames(format_summary)){
    format_summary <- cbind(
      format_summary[, -which(colnames(format_summary) == 'NA')],
      format_summary[, c('NA')]
    )
    colnames(format_summary)[ncol(format_summary)] <- 'NA'
    format_summary$`NA` <- as.character(format_summary$`NA`)
  }
  
  ## add total number of test
  format_summary <- cbind(
    format_summary,
    N_test = c(
      total_by_test$N[match(paste0(format_summary$Tech, format_summary$BiomarkerName),
                            paste0(total_by_test$Tech, total_by_test$BiomarkerName))]
    )
  )
  ## add total number of patients
  format_summary <- rbind(
    format_summary,
    c('Total', '#Patient', 
      total_by_lineage$N[match(colnames(format_summary)[3:ncol(format_summary)], total_by_lineage$Lineage)])
  )
  colnames(format_summary)[length(format_summary)] <- 'Count of Tests'
  
  format_summary
}



################################################################################################################
## Main wrapper of functions, call for each COE site
################################################################################################################
coe_combine_summary <- function(coe_name, summary_cutoff){
  library(reshape2)
  library(WriteXLS)

  now <- format(Sys.Date(), '%Y%m%d')
  
  df_test <- get_test_result(coe_name, summary_cutoff = summary_cutoff)
  
  ## obtain overall summary of SEQ, CNV, and all other tests
  cnv_summary <- get_ngs_summary_by_case(df_test, type = 'CNV')
  seq_summary <- get_ngs_summary_by_case(df_test, type = 'SEQ')
  seq_summary$Tech <- 'NGS'
  test_all_summary <- get_other_summary_by_case(df_test)
  
  ## obtain by physician summary
  ls_test_summary <- dlply(
    df_test, 
    .variables = 'Physician',
    function(x){
      get_other_summary_by_case(x)
    }
  )
  ls_test_summary[[length(ls_test_summary)+1]] <- test_all_summary
  names(ls_test_summary)[length(ls_test_summary)] <- 'All'
  
  ## combined CNV, NGS, and others into one sheet for each physician
  ls_combined_summary <- llply(
    ls_test_summary, 
    function(x){
      #x <- ls_test_summary$`Ma Patrick`
      total_by_test <- x$total_by_test
      total_by_lineage <- x$total_by_lineage
      lineages <- unique(x$Test_summary$Lineage)
      Test_summary <- rbind(x$Test_summary, cnv_summary, seq_summary)
      Test_summary <- Test_summary[Test_summary$Lineage %in% lineages, ] 
      
      # format the output tabs
      out_summary <- dcast(Test_summary, Tech + BiomarkerName ~ Lineage, value.var = 'Value')
      out_summary <- format_out_summary(out_summary, total_by_test, total_by_lineage)
    }
  )
  
  # prepare the summary tab
  Summary <- get_site_physician_summary(coe_name = coe_name, summary_cutoff = summary_cutoff)
  Summary$OrderingPhysician <- as.character(Summary$OrderingPhysician)
  Summary <- rbind(Summary, c('Total', NA, sum(Summary$n_sample)))
  colnames(Summary) <- c('Physician', 'Category', '#Patients')
  ls_combined_summary[[length(ls_combined_summary) + 1]] <- Summary
  names(ls_combined_summary)[length(ls_combined_summary)] <- 'Summary'
  
  # reorder the tabs
  n_tabs <- length(ls_combined_summary)
  ls_combined_summary_out <- list()
  ls_combined_summary_out[[1]] <- ls_combined_summary[[n_tabs]]
  names(ls_combined_summary_out)[1] <- names(ls_combined_summary)[n_tabs]
  ls_combined_summary_out[[2]] <- ls_combined_summary[[n_tabs - 1]]
  names(ls_combined_summary_out)[2] <- names(ls_combined_summary)[n_tabs - 1]
  for(i in 1:(n_tabs - 2)){
    ls_combined_summary_out[[2+i]] <- ls_combined_summary[[i]]
    names(ls_combined_summary_out)[2+i] <- names(ls_combined_summary)[i]
  }
  
  # write excel file
  WriteXLS(ls_combined_summary_out, ExcelFileName = paste0(now, '_', gsub('\\/', '_', coe_name), '_byPhysician.xls'), SheetNames = names(ls_combined_summary_out))
  
  ls_combined_summary_out
}

################################################################################################################################################################################################################################
################################################################################################################################################################################################################################
library(plyr)

## configurations
now <- format(Sys.Date(), '%Y%m%d')
out_dir <- paste0('/Users/jwang/Documents/Project4_CMI_Request/COE_detailed_report/', now, '_COE_detailed_report')
dir.create(out_dir)
setwd(out_dir)
summary_cutoff <- 20


coe_names <- get_coe_names()

#working_dir <- getwd()
#out_dir <- paste0(working_dir, '/COE_summary')
#dir.create(path = out_dir)
#setwd(out_dir)

ls_summary <- llply(
  coe_names, 
  function(x){
    cat('summarizing', x, '...\n')
    coe_combine_summary(x, summary_cutoff)
  }
)


# coe_name <- 'Valley Medical Oncology Consortium'
# summary_cutoff <- 20
# test <- coe_combine_summary('Valley Medical Oncology Consortium', 20)

#setwd(working_dir)



