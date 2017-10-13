rm(list = ls())
source('~/R_modules/clinical_data/clinical_data_pipeline.R')
source('~/R_modules/DrugRules/drugRules_process.R')

main <- function(){
  df_mstatus <- read_cmi_mstatus(dir_cmi_file)
  df_mstatus <- marker_format(df_mstatus)
  df_recommend <- drug_recommend(df_mstatus)
  write_recommend(df_recommend)
}

####
dir_cmi_file <- '~/Documents/Project3_Drug_Biomarker_Survival/RProjects/20160613_drug_rule_explore/drugrule_files/20160616_cmi_ova.csv'
#dir_cmi_file <- '~/Documents/Project3_Drug_Biomarker_Survival/RProjects/20160613_drug_rule_explore/drugrule_files/20160616_cmi_ova_brca_rerun.csv'

system.time(main())
