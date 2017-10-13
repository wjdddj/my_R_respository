################################################################################################################
## fit coxph model: surv ~ match_score + registry_case
## and extract HR and 95%CI and p value
################################################################################################################
get_cox_hr <- function(df_clinical_info){
  #df_clinical_info <- df_surv
  
  library(survival)
  # df_clinical_info <- df_plot1
  # plot_covariate <- 'followed'
  is_unknown <- any(sapply(df_clinical_info, function(x) grepl('Unknown', x)))
  if(is_unknown){
    stop('unknowns detected, please remove the corresponding samples and rerun!\n')
  }
  ## with covariates intended for adjustment
  fitCox_sum <- summary(coxph(Surv(time, censored) ~ ., data = df_clinical_info))
  z <- list(HR_match = round(fitCox_sum$coef[1, 2], 3), 
            CI_95_match = paste0('(',round(fitCox_sum$conf.int[1, 3], 3), ',', round(fitCox_sum$conf.int[1, 4], 3), ')'),
            p_value_match = round(fitCox_sum$coef[1, 5], 3),
            HR_reg = round(fitCox_sum$coef[2, 2], 3), 
            CI_95_reg = paste0('(',round(fitCox_sum$conf.int[2, 3], 3), ',', round(fitCox_sum$conf.int[2, 4], 3), ')'),
            p_value_reg = round(fitCox_sum$coef[2, 5], 3),
            p_value_lr = round(fitCox_sum$sctest[3], 3))
}


################################################################################################################
## main function for first line after collection
## fit coxph model: surv ~ match_score + registry_case
################################################################################################################
main_per_weight_1stline <- function(df_therapy, ls_recommend, wt_star, reg_druggroup_mapping){
  #wt_star <- data.frame(druggroup = colnames(wt_star_matrix),
  #                      weight = wt_star_matrix[1,])
  df_therapy$therapies <- tolower(df_therapy$therapies)
  
  ## convert to druggroups
  df_therapy$druggroups <- sapply(df_therapy$therapies, function(x){
    paste(convert_registry_to_druggroup(strsplit(x, ',')[[1]], reg_druggroup_mapping), collapse = ',')
  })
  
  ## compute matchness
  inter_masterdeid <- intersect(df_therapy$masterdeid, names(ls_recommend))
  ls_recommend <- ls_recommend[match(inter_masterdeid, names(ls_recommend))]
  df_therapy <- df_therapy[match(inter_masterdeid, df_therapy$masterdeid), ]
  # df_therapy$masterdeid == names(ls_recommend)
  df_therapy <- get_match_score_df(df_therapy, ls_recommend, wt_star, 
                                   adjust_recommend = adjust_recommend, recom_weight = recom_weight)
  #head(df_therapy)
  ## prepare for survival analysis
  df_surv <- data.frame(time = as.numeric(df_therapy$tnt)/365,
                        censored = df_therapy$event,
                        match_score = df_therapy$match_score,
                        registry_case = df_therapy$registry_case)
  df_surv$registry_case <- factor(df_surv$registry_case, levels = c('N', 'Y'))
  HR_p <- get_cox_hr(df_surv)
  
  ## output
  match_scores <- df_surv$match_score
  df_summary <- data.frame(drugRule_version = ls_recommend[[1]]$Version[1],
                           weight_combo = paste(apply(wt_star, 1, function(x)paste(x, collapse = '=')), collapse = '|'),
                           HR_match = HR_p$HR_match,
                           CI_95_match = HR_p$CI_95_match,
                           p_value_match = HR_p$p_value_match,
                           HR_reg = HR_p$HR_reg,
                           CI_95_reg = HR_p$CI_95_reg,
                           p_value_reg = HR_p$p_value_reg,
                           p_value_lr = HR_p$p_value_lr)
  z <- list(match_scores = match_scores,
            df_summary = df_summary)
}
