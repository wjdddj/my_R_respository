rm(list=ls())
library(reshape2)
library(plyr)
library(RMySQL)

## file configurations
today <- format(Sys.time(), '%Y%m%d')

row_name_file <- '/Users/jwang/Documents/Project4_CMI_Request/CMI_high_level_summary/configuration/20170303_byhistology_match_rownames.csv'
col_name_file <- '/Users/jwang/Documents/Project4_CMI_Request/CMI_high_level_summary/configuration/20170303_byhistology_match_colnames.csv'
dirpath = paste0('/Users/jwang/Documents/Project4_CMI_Request/CMI_high_level_summary/', today, '_cmi_high_level')
dir.create(dirpath)
setwd(dirpath)

## db configurations
user="sting"
password="Pass1234"
host="bioserv.clsnet.intranet"
dbname="bioinfo_projects"
port=3307
con <- dbConnect(MySQL(),user = user, password = password, dbname = dbname, host = host, port = port)

### Overall Summary

sqlcmd = "
select 'Total N' as 'Header', Lineage, Histology, count(*) as 'value' from cmi_highlevel_case_summary group by Lineage, Histology UNION
select 'IHC', Lineage, Histology, count(*) as 'value' from cmi_highlevel_case_summary where IHC!=0 group by Lineage, Histology UNION 
select 'ISH', Lineage, Histology, count(*) as 'cases' from cmi_highlevel_case_summary where ISH!=0 group by Lineage, Histology UNION 
select 'Microarray', Lineage, Histology, count(*) as 'cases' from cmi_highlevel_case_summary where Microarray!=0 group by Lineage, Histology UNION 
select '45Gene', Lineage, Histology, count(*) as 'cases' from cmi_highlevel_case_summary where 45Gene!=0 group by Lineage, Histology UNION 
select 'BRCA', Lineage, Histology, count(*) as 'cases' from cmi_highlevel_case_summary where BRCA!=0 group by Lineage, Histology UNION 
select '600Gene', Lineage, Histology, count(*) as 'cases' from cmi_highlevel_case_summary where 600Gene!=0 group by Lineage, Histology UNION 
select 'Fusion', Lineage, Histology, count(*) as 'cases' from cmi_highlevel_case_summary where Fusion!=0 group by Lineage, Histology UNION
select 'CNV', Lineage, Histology, count(*) as 'value' from cmi_highlevel_case_summary where CNV!=0 group by Lineage, Histology UNION
select 'RegION', Lineage, Histology, count(*) as 'value' from cmi_highlevel_case_summary where RegION in ('Reg','ION') and isexcluded = 0 group by Lineage, Histology UNION
select RegION, Lineage, Histology, count(*) as 'value' from cmi_highlevel_case_summary where RegION in ('Reg','ION') group by Lineage, Histology, RegION UNION
select concat('Outcome_',Outcome) as 'Outcome', Lineage, Histology, count(*) as 'value' from cmi_highlevel_case_summary where Outcome is not null group by Lineage, Histology, Outcome UNION
select concat('Race_',Race) as 'Race', Lineage, Histology, count(*) as 'value' from cmi_highlevel_case_summary where Race is not null group by Lineage, Histology, Race UNION
select if(Sex='f','Female',if(Sex='m','Male','NA')) as 'Sex', Lineage, Histology, count(*) as 'value' from cmi_highlevel_case_summary where Sex is not null and Sex in ('f','m') group by Lineage, Histology, Sex UNION
select concat('SlidesRemaining_', SlidesRemaining) as SlidesRemaining, Lineage, Histology, count(*) as 'value' from bioinfo_projects.cmi_highlevel_case_summary group by Lineage, Histology, SlidesRemaining UNION
select concat('BlockRemaining_', BlockRemaining) as BlockRemaining, Lineage, Histology, count(*) as 'value' from bioinfo_projects.cmi_highlevel_case_summary group by Lineage, Histology, BlockRemaining;
"

dbout <- dbGetQuery(con, sqlcmd)
dbout$Lineage[is.na(dbout$Lineage)] <- 'NA'
dbout$Histology[is.na(dbout$Histology)] <- 'NA'

pivot.all <- dcast(dbout, Lineage + Histology ~ Header, sum)

#### ION/Reg MatchedCat
sqlcmd = "
select RegION, MatchedCat, Lineage, Histology, count(*) as 'value' from cmi_highlevel_case_summary 
where RegION in ('Reg','ION') and MatchedCat is not null
group by Lineage, Histology, RegION, MatchedCat;"

dbout <- dbGetQuery(con, sqlcmd)
dbout$Lineage[is.na(dbout$Lineage)] <- 'NA'
dbout$Histology[is.na(dbout$Histology)] <- 'NA'

pivot <- dcast(dbout, Lineage + Histology ~ RegION + MatchedCat, sum)
pivot.all <- merge(
  pivot.all,
  pivot,
  by = c('Lineage', 'Histology'),
  all = TRUE
)
pivot.all[is.na(pivot.all)] <- 0


#### ION/Reg Stage
sqlcmd = "
select concat('Stage_',Stage) as 'Stage', Lineage, Histology, count(*) as 'value' from cmi_highlevel_case_summary 
where Stage is not null
group by Lineage, Histology, Stage;"

dbout <- dbGetQuery(con, sqlcmd)
dbout$Lineage[is.na(dbout$Lineage)] <- 'NA'
dbout$Histology[is.na(dbout$Histology)] <- 'NA'

pivot <- dcast(dbout, Lineage + Histology ~ Stage, sum)
# colnames(pivot)
# Stage_0 <- apply(pivot[,3:4], 1, sum)
# Stage_I <- apply(pivot[,5:10], 1, sum)
# Stage_II <- apply(pivot[,11:14], 1, sum)
# Stage_III <- pivot[,15:17]
# Stage_IIIC <- apply(pivot[,18:19], 1, sum)
# Stage_IV <- apply(pivot[,20:23], 1, sum)
# Stage_Unknown <- pivot[,24]

Stage_0 <- apply(pivot[, grep('Stage_0', colnames(pivot))], 1, sum)
Stage_I <- apply(pivot[, grep('Stage_I[A-C0-9]*$', colnames(pivot))], 1, sum)
Stage_II <- apply(pivot[, grep('Stage_II[A-C0-9]*$', colnames(pivot))], 1, sum)
Stage_III <- pivot[, grep('Stage_III[A-B0-9]*$', colnames(pivot))]
Stage_IIIC <- pivot[, grep('Stage_IIIC[0-9]*$', colnames(pivot))]
Stage_IV <- apply(pivot[, grep('Stage_IV[A-C0-9]*$', colnames(pivot))], 1, sum)
Stage_Unknown <- pivot[, grep('Stage_Unknown$', colnames(pivot))]

pivot <- data.frame(
  pivot[,c('Lineage', 'Histology')],
  Stage_0, Stage_I, Stage_II, Stage_III, Stage_IIIC, Stage_IV, Stage_Unknown
)

pivot.all <- merge(
  pivot.all,
  pivot,
  by = c('Lineage', 'Histology'),
  all = TRUE
)
pivot.all[is.na(pivot.all)] <- 0

#### Age
sqlcmd = "
select 
Lineage, Histology,
if(age<40,'Age_39-',
if(age>=70,'Age_70+',
concat('Age_',floor(age/10),'0-',floor(age/10),'9'))) as 'AgeRange',count(*) as 'value'
from cmi_highlevel_case_summary
group by
Lineage, Histology,
if(age<40,'Age_39-',
if(age>=70,'Age_70+',
concat('Age_',floor(age/10),'0-',floor(age/10),'9')))
;"

dbout <- dbGetQuery(con, sqlcmd)
dbout$Lineage[is.na(dbout$Lineage)] <- 'NA'
dbout$Histology[is.na(dbout$Histology)] <- 'NA'

pivot <- dcast(dbout, Lineage + Histology ~ AgeRange, sum)

pivot.all <- merge(
  pivot.all,
  pivot,
  by = c('Lineage', 'Histology'),
  all = TRUE
)
pivot.all[is.na(pivot.all)] <- 0


#### SampleTested
sqlcmd = "
select concat('SampleTested_',SampleTested) as 'SampleTested',
Lineage, Histology,
count(*) as 'value' 
from cmi_highlevel_case_summary 
where SampleTested is not null and SampleTested!='NA'
group by Lineage, Histology, SampleTested;"

dbout <- dbGetQuery(con, sqlcmd)
dbout$Lineage[is.na(dbout$Lineage)] <- 'NA'
dbout$Histology[is.na(dbout$Histology)] <- 'NA'
dbout_consolidate <- dlply(dbout, .variables = c('Lineage','Histology'))
priority <- c('Primary', 'Metastatic', 'Local Recurrence', 'Unknown')
dbout_clean <- llply(1:length(dbout_consolidate), function(x){
  n <- x
  x <- dbout_consolidate[[x]]
  Lineage <- x$Lineage[1]
  Histology <- x$Histology[1]
  count.all <- c()
  for(each_priority in priority){
    rowid <- grep(each_priority, x$SampleTested)
    count <- sum(x$value[rowid])
    x <- x[-rowid, ]
    count.all <- c(count.all, count)
  }
  count.all <- data.frame(
    SampleTested = paste('SampleTested_', priority, sep = ''), 
    Lineage = Lineage,
    Histology = Histology,
    value = count.all
  )
})
dbout_clean <- do.call(rbind, dbout_clean)

pivot <- dcast(dbout_clean, Lineage + Histology ~ SampleTested, sum)
pivot.all <- merge(
  pivot.all,
  pivot,
  by = c('Lineage', 'Histology'),
  all = TRUE
)
pivot.all[is.na(pivot.all)] <- 0

##### Last process and output
pivot.all <- t(pivot.all)
## match to the row orders
row_mat <- read.csv(row_name_file, head = T, stringsAsFactors = F)
pivot.all <- pivot.all[match(row_mat$db_name, rownames(pivot.all)), ]
## match to the col orders
col_mat <- read.csv(col_name_file, head = T, stringsAsFactors = F)
col_mat[is.na(col_mat)] <- 'NA'
pivot.all <- pivot.all[, match(
  paste0(col_mat$db_lineage[2:nrow(col_mat)], col_mat$db_histology[2:nrow(col_mat)]), 
  paste0(pivot.all[1,], pivot.all[2,])
)]
pivot.all <- cbind(row_mat$report_name, pivot.all)

colnames(pivot.all) <- NULL

ofilename <- paste(dirpath, paste(today,'_','CMI_HighLevel_Summary_with_histology','.csv',sep=''), sep="/")
write.table(pivot.all, file = ofilename, sep = ',', row.names = F, col.names = F)

## only 'Female Genital Tract Malignancy' and 'Ovarian Surface Epithelial Carcinomas'
ofilename2 <- paste(dirpath, paste(today,'_','CMI_HighLevel_Summary_with_histology_ovarian&femaleGenital','.csv',sep=''), sep="/")
write.table(
  pivot.all[, c(1, which(pivot.all[1,]%in%c(
    'Female Genital Tract Malignancy', 
    'Ovarian Surface Epithelial Carcinomas'
  )))], 
  file = ofilename2, 
  sep = ',', row.names = F, col.names = F
)


dbDisconnect(con)
