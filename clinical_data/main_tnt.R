################################################################################################################
## main function to calculate match/unmatch, tnt, and perform survival analysis and kmplot
## input:
##    1. drug_rule_recommend, a data.frame with three columns 'masterdeid', 'registry_name'(drug name), 
##  'recommendation'('Benefit', 'Lack Of Benefit', 'Indeterminate', 'DoNotReport')
##    2. df_therapy, a data.frame of therapy of a single or multiple subjects. It MUST contain these columns: 
##  'masterdeid', 'line_num', 'regimen_start', 'deathdate', 'lastcontactdate'.
##    3. start_line, line number of the start line.
##    4. end_line, line number of the end line, it could also be 'max', the last line will be included.
##
## example work flow:
##    1. prepare/format the drug_rule_recommend and df_therapy objects.
##    2. call main_tnt().
##       secondline_restricted_legacy <- main_tnt(legacy_table, df_therapy, start_line = 2, end_line = 3)
##    3. plot.
##       pdf(paste0(Sys.Date(), 'filename.pdf'), 7, 7)
##       kmplots(secondline_restricted_legacy$fit_surv, main = 'xxx')
##       dev.off()
################################################################################################################
main_tnt <- function(drug_rule_recommend, df_therapy, start_line = 2, end_line = 3){
  # drug_rule_recommend <- newrule_table
  mat_recommend <- format_recommend(drug_rule_recommend)
  # compute the tnt and drug combinations of each line period of interest
  df_therapy <- get_therapy_tnt_df(df_therapy_ken, start_line = start_line, end_line = end_line)
  # convert therapies to lower case
  df_therapy$therapies <- tolower(df_therapy$therapies)
  # examine matching between recommendations and therapies
  # whether recommendations covered all the query cases
  if(sum(!df_therapy$masterdeid%in%rownames(mat_recommend)) > 0){
    missing_idx <- which(!df_therapy$masterdeid%in%rownames(mat_recommend))
    cat('\nThe following cases (in masterdeid) do not have recommendations. They will be removed from the analysis. 
        Please check your drug_rule_recommend table!\n')
    print(as.character(df_therapy$masterdeid[missing_idx]))
  }
  
  inter_masterdeid <- intersect(df_therapy$masterdeid, rownames(mat_recommend))
  mat_recommend_matched <- mat_recommend[match(inter_masterdeid, rownames(mat_recommend)), ]
  df_therapy_matched <- df_therapy[match(inter_masterdeid, df_therapy$masterdeid), ]
  df_therapy_matched <- get_matched_df(df_therapy_matched, mat_recommend_matched)
  
  # prepare for survival analysis
  df_surv <- data.frame(time = as.numeric(df_therapy_matched$tnt)/365,
                        censored = df_therapy_matched$event,
                        followed = df_therapy_matched$followed)
  
  df_surv <- df_surv[df_surv$followed != 'neither', ]
  df_surv$followed <- factor(df_surv$followed, levels = c('unmatched', 'matched'))
  
  fit_surv <- surv_fit(df_surv, plot_covariate = 'followed')
  
  z <- list(fit_surv = fit_surv, 
            df_therapy_matched = df_therapy_matched,
            mat_recommend_matched = mat_recommend_matched)
}











