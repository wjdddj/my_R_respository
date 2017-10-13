library(RMySQL)

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
## Obtain query from mySQL databases on server
################################################################################################################
get_query <- function(
  server,
  port,
  username,
  password,
  database,
  query
){
  drv <- dbDriver("MySQL") 
  dbcon <- dbConnect(
    MySQL(), 
    user = username, 
    password = password, 
    host = server, 
    dbname = database, 
    port = port
  ) 
  rs1 = dbSendQuery(dbcon, query)
  df_out <- dbFetch(rs1, n = -1)
  suppressWarnings(dbDisconnect(dbcon))
  df_out
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


