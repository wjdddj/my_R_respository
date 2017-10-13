rm(list=ls())
library(reshape2)
library(plyr)
library(RMySQL)

## file configurations
today <- format(Sys.time(), '%Y%m%d')

row_name_file <- '/Users/jwang/Documents/Project4_CMI_Request/CMI_high_level_summary/configuration/20160725_match_rownames.csv'
col_name_file <- '/Users/jwang/Documents/Project4_CMI_Request/CMI_high_level_summary/configuration/20160725_match_colnames.csv'
dirpath = paste0('/Users/jwang/Documents/Project4_CMI_Request/CMI_high_level_summary/', today, '_cmi_high_level')
dir.create(dirpath)
setwd(dirpath)

#options(width=240)


## db configurations
user="sting"
password="Pass1234"
host="bioserv.clsnet.intranet"
dbname="bioinfo_projects"
port=3307

#now <- format(Sys.time(), '%Y%m%d-%H%M%S')

con <- dbConnect(MySQL(),user = user, password = password, dbname = dbname, host = host, port = port)

### Overall Summary
sqlcmd = "
select 'Total N' as 'Header', Lineage, count(*) as 'value' from cmi_highlevel_case_summary group by Lineage UNION
select 'IHC', Lineage, count(*) as 'value' from cmi_highlevel_case_summary where IHC!=0 group by Lineage UNION 
select 'ISH', Lineage, count(*) as 'value' from cmi_highlevel_case_summary where ISH!=0 group by Lineage UNION 
select 'Microarray', Lineage, count(*) as 'value' from cmi_highlevel_case_summary where Microarray!=0 group by Lineage UNION 
select '45Gene', Lineage, count(*) as 'value' from cmi_highlevel_case_summary where 45Gene!=0 group by Lineage UNION 
select 'BRCA', Lineage, count(*) as 'value' from cmi_highlevel_case_summary where BRCA!=0 group by Lineage UNION 
select '600Gene', Lineage, count(*) as 'value' from cmi_highlevel_case_summary where 600Gene!=0 group by Lineage UNION 
select 'Fusion', Lineage, count(*) as 'value' from cmi_highlevel_case_summary where Fusion!=0 group by Lineage UNION
select 'CNV', Lineage, count(*) as 'value' from cmi_highlevel_case_summary where CNV!=0 group by Lineage UNION
select 'RegION', Lineage, count(*) as 'value' from cmi_highlevel_case_summary where RegION in ('Reg','ION') and isexcluded = 0 group by Lineage UNION
select RegION, Lineage, count(*) as 'value' from cmi_highlevel_case_summary where RegION in ('Reg','ION') and isexcluded = 0 group by Lineage, RegION UNION
select concat('Outcome_', Outcome) as 'Outcome', Lineage, count(*) as 'value' from cmi_highlevel_case_summary where Outcome is not null and isexcluded = 0 group by Lineage, Outcome UNION
select concat('Race_',Race) as 'Race', Lineage, count(*) as 'value' from cmi_highlevel_case_summary where Race is not null and isexcluded = 0 group by Lineage, Race UNION
select if(Sex='f', 'Female', if(Sex='m', 'Male', 'NA')) as 'Sex', Lineage, count(*) as 'value' from cmi_highlevel_case_summary where Sex is not null and Sex in ('f','m') group by Lineage, Sex UNION
select concat('SlidesRemaining_', SlidesRemaining) as SlidesRemaining, Lineage, count(*) as 'value' from cmi_highlevel_case_summary where isexcluded = 0 group by Lineage,SlidesRemaining UNION
select concat('BlockRemaining_', BlockRemaining) as BlockRemaining, Lineage, count(*) as 'value' from cmi_highlevel_case_summary where isexcluded = 0 group by Lineage, BlockRemaining
;"

dbout <- dbGetQuery(con, sqlcmd)
dbout$Lineage[is.na(dbout$Lineage)] <- 'NA'

pivot.all <- dcast(dbout, Lineage ~ Header, sum)

#### ION/Reg
sqlcmd = "
select RegION, MatchedCat, Lineage, count(*) as 'value' from cmi_highlevel_case_summary where RegION in ('Reg','ION') and MatchedCat is not null and isexcluded = 0 group by Lineage, RegION, MatchedCat
;"

dbout <- dbGetQuery(con, sqlcmd)
dbout$Lineage[is.na(dbout$Lineage)] <- 'NA'

pivot <- dcast(dbout, Lineage ~ RegION + MatchedCat, sum)

pivot.all <- merge(pivot.all, pivot, by = "Lineage", all = TRUE, incomparables = 0)

pivot.all[is.na(pivot.all)] <- 0

#### ION/Reg Stage
sqlcmd = "
select concat('Stage_',Stage) as 'Stage', Lineage, count(*) as 'value' from cmi_highlevel_case_summary 
where Stage is not null and isexcluded = 0 
group by Lineage, Stage;"

dbout <- dbGetQuery(con, sqlcmd)
dbout$Lineage[is.na(dbout$Lineage)] <- 'NA'

pivot <- dcast(dbout, Lineage ~ Stage, sum)

Stage_0 <- apply(pivot[,2:3], 1, sum)
Stage_I <- apply(pivot[,4:9], 1, sum)
Stage_II <- apply(pivot[,10:13], 1, sum)
Stage_III <- pivot[,14:16]
Stage_IIIC <- apply(pivot[,17:18], 1, sum)
Stage_IV <- apply(pivot[,19:22], 1, sum)
Stage_Unknown <- pivot[,23]
pivot <- data.frame(Lineage = pivot$Lineage,
                    Stage_0,
                    Stage_I, Stage_II, Stage_III, Stage_IIIC, Stage_IV, Stage_Unknown)

pivot.all <- merge(pivot.all, pivot, by = "Lineage", all = TRUE, incomparables = 0)

pivot.all[is.na(pivot.all)] <- 0

#### Age
sqlcmd="
select
Lineage,
if(age<40,'Age_39-',
if(age>=70,'Age_70+',
concat('Age_',floor(age/10),'0-',floor(age/10),'9'))) as 'AgeRange', count(*) as 'value'
from cmi_highlevel_case_summary
group by
Lineage,
if(age<40,'Age_39-',
if(age>=70,'Age_70+',
concat('Age_',floor(age/10),'0-',floor(age/10),'9')));"

dbout <- dbGetQuery(con, sqlcmd)
dbout$Lineage[is.na(dbout$Lineage)] <- 'NA'

pivot <- dcast(dbout, Lineage ~ AgeRange, sum)

pivot.all <- merge(pivot.all, pivot, by = "Lineage", all = TRUE, incomparables = 0)

pivot.all[is.na(pivot.all)] <- 0

#### SampleTested
sqlcmd = "
select concat('SampleTested_',SampleTested) as 'SampleTested', Lineage, count(*) as 'value' 
from cmi_highlevel_case_summary 
where SampleTested is not null and SampleTested!='NA' and isexcluded = 0 
group by Lineage, SampleTested;"

dbout <- dbGetQuery(con, sqlcmd)
dbout$Lineage[is.na(dbout$Lineage)] <- 'NA'
dbout_consolidate <- dlply(dbout, .variables = 'Lineage')
priority <- c('Primary', 'Metastatic', 'Local Recurrence', 'Unknown')
dbout_clean <- llply(1:length(dbout_consolidate), function(x){
  n <- x
  x <- dbout_consolidate[[x]]
  count.all <- c()
  for(each_priority in priority){
    rowid <- grep(each_priority, x$SampleTested)
    count <- sum(x$value[rowid])
    x <- x[-rowid, ]
    count.all <- c(count.all, count)
  }
  count.all <- data.frame(SampleTested = priority, 
                          Lineage = names(dbout_consolidate)[n], 
                          value = count.all)
})
dbout_clean <- do.call(rbind, dbout_clean)
dbout_clean$Lineage <- as.character(dbout_clean$Lineage)
pivot <- dcast(dbout_clean, Lineage ~ SampleTested, sum)
pivot$Lineage[is.na(pivot$Lineage)] <- 'NA'

pivot.all<-merge(pivot.all, pivot, by = "Lineage", all = TRUE, incomparables = 0)

pivot.all[is.na(pivot.all)] <- 0
pivot.all <- rbind(pivot.all, c('Total', sapply(pivot.all[, -1], sum, na.rm = T)))

##### Last process and output
pivot.all <- t(pivot.all)

## match to the row orders
row_mat <- read.csv(row_name_file, head = T, stringsAsFactors = F)
pivot.all <- pivot.all[match(row_mat$db_name, rownames(pivot.all)), ]
col_mat <- read.csv(col_name_file, head = T, stringsAsFactors = F)
col_mat[is.na(col_mat)] <- 'NA'
pivot.all <- pivot.all[, match(as.character(col_mat$db_name), pivot.all[1,])]
pivot.all <- cbind(row_mat$report_name, pivot.all)
colnames(pivot.all) <- NULL
ofilename = paste(dirpath, paste(today, '_', 'CMI_HighLevel_Summary', '.csv',sep=''), sep="/")

write.table(pivot.all, file = ofilename, sep = ',', row.names = F, col.names = F)

dbDisconnect(con)

