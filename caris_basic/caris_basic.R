################################################################################################################
## detach specific loaded packages
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
## Obtain table from mySQL bioserv02 databases using Ryan's account
################################################################################################################
get_ryan_query_bioserv02 <- function(database, query){
  library(RMySQL)
  drv <- dbDriver("MySQL") 
  dbcon <- dbConnect(MySQL(), user="jwang", password="Pass1234", host = 'bioserv02.clsnet.intranet', dbname = database, port = 3306) 
  rs1 = dbSendQuery(dbcon, query)
  df_out <- dbFetch(rs1, n = -1)
  # dbListConnections(dbDriver( drv = "MySQL"))
  suppressWarnings(dbDisconnect(dbcon))
  # lapply(dbListConnections(dbDriver( drv = "MySQL")), dbDisconnect)
  df_out
}

################################################################################################################
## Close all DB connection
################################################################################################################
closeAllConnections <- function(){
  lapply(dbListConnections(MySQL()), dbDisconnect)
}

################################################################################################################
## trim off heading/trailing spaces and/or punctuations from string
################################################################################################################
trim <- function (x) gsub("^(\\s|[[:punct:]])+|([[:punct:]]|\\s)+$", "", x)



################################################################################################################
## convert different format of NAs into NA
################################################################################################################
convertNA <- function(df, na = c('NA', 'na', 'n.a.')){
  dfNew <- df
  dfNew[is.na(dfNew)] <- NA
  for(i in 1:ncol(dfNew)){
    dfNew[dfNew[,i] %in% na, i] <- NA
  }
#   dfNew <- sapply(
#     dfNew,
#     function(x){
#       x[x %in% na] <- NA
#       x
#     }
#   )
  dfNew
}

################################################################################################################
## collapse a vector to accommodate mysql query where in '()'.
################################################################################################################
formQueryItemArray <- function(
  queryItems
){
  out <- paste0(
    "('",
    paste(
      queryItems,
      collapse = "','"
    ),
    "')"
  )
}

################################################################################################################
## to convert IDs.
################################################################################################################
retrieveIDs <- function(
  inputID, 
  inidType = 'AccessionNumber'
#  outidType = 'masterdeid'
){
#   inputID <- clinDatasetDF$AccessionNumber
  mySQLIDs <- formQueryItemArray(inputID)
  IDtable <- get_sting_query(
    database = 'bioinfo',
    query = paste0(
      'select * from deidlookup\n',
      'where ', inidType, ' in ', mySQLIDs, ';'
    )
  )
#  IDtable[, c(inidType, outidType)]
  IDtable
}

