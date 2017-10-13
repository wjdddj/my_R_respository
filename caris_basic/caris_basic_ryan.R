source('~/R_modules/caris_basic/caris_basic_pub.R')
library(RMySQL)

################################################################################################################
## Obtain query from bioserv:3307 using Sting's account
## depends on caris_basic_pub.R
################################################################################################################
get_sting_bioserv <- function(database, query){
  df_out <- get_query(
    server = 'bioserv.clsnet.intranet', 
    port = 3307,
    username = "sting",
    password = "Pass1234",
    database = database,
    query = query
  )
  df_out
}

################################################################################################################
## Obtain query from bioserv:3307 using Ryan's account
## depends on caris_basic_pub.R
################################################################################################################
get_ryan_bioserv <- function(database, query){
  df_out <- get_query(
    server = 'bioserv.clsnet.intranet', 
    port = 3307,
    username = "jwang",
    password = "Pass1234",
    database = database,
    query = query
  )
  df_out
}

################################################################################################################
## Obtain query from bioserv02:3306 using Ryan's account
## depends on caris_basic_pub.R
################################################################################################################
get_ryan_bioserv02_3306 <- function(database, query){
  df_out <- get_query(
    server = 'bioserv02.clsnet.intranet', 
    port = 3306,
    username = "jwang",
    password = "Pass1234",
    database = database,
    query = query
  )
  df_out
}

################################################################################################################
## Obtain query from bioserv02:3307 using Ryan's account
## depends on caris_basic_pub.R
################################################################################################################
get_ryan_bioserv02_3307 <- function(database, query){
  df_out <- get_query(
    server = 'bioserv02.clsnet.intranet', 
    port = 3307,
    username = "jwang",
    password = "wryan",
    database = database,
    query = query
  )
  df_out
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
  IDtable <- get_sting_bioserv(
    database = 'bioinfo',
    query = paste0(
      'select * from deidlookup\n',
      'where ', inidType, ' in ', mySQLIDs, ';'
    )
  )
  #  IDtable[, c(inidType, outidType)]
  IDtable
}

