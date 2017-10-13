################################################################################################################
## main function to calculate match/unmatch,and perform survival analysis and kmplot
## Target: first line after collection, tnt is already calculated by Ken, therefore, no need to call 
##        get_therapy_tnt_df().
## input:
##      1. drug_rule_recommend. MUST follow the input format of format_recommend(). 
##        column1, masterdeid; column2, registry_name; column3, recommendation
##      2. df_therapy. MUST contain columns: masterdeid, therapies, event, tnt
################################################################################################################
main_tnt_1st_line_after_coll <- function(drug_rule_recommend, df_therapy){
  # drug_rule_recommend <- newrule_6804_rerunBRCA
  mat_recommend <- format_recommend(drug_rule_recommend)
  # compute the tnt and drug combinations of each line period of interest
  # df_therapy <- get_therapy_tnt_df(df_therapy_ken, start_line = 2, end_line = 3)
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











