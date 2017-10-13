########################################################################################################################################################################
## connect to CMI (bioinfo) database
## obtain all_test_cases table
########################################################################################################################################################################
get_all_test_cases <- function(SQL_query = 'select * from all_test_cases'){
  library(RMySQL)
  drv <- dbDriver("MySQL") 
  bioinfo <- dbConnect(MySQL(), user="jwang", password="Pass1234", host = 'bioserv.clsnet.intranet', dbname = 'bioinfo', port = 3307) 
  rs1 = dbSendQuery(bioinfo, SQL_query)
  all_test_cases <- fetch(rs1, n = -1)
  dbListConnections( dbDriver( drv = "MySQL"))
  lapply(dbListConnections( dbDriver( drv = "MySQL")), dbDisconnect)
  all_test_cases
}



convert_IHC_result <- function(IHC_result){
  IHC_result <- gsub('Above Threshold', 'Positive', IHC_result)
  IHC_result <- gsub('Test not performed', 'Other', IHC_result)
}

