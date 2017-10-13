## Antimetabolites (gemcitabine)
dr_gemcit <- function(marker_status){
  #marker_status <- df_mstatus_test1
  dr_gemcit_version <- 'ver0'
  
  markers_in_rule <- c('RRM1_IHC')
  marker_status <- match_markers(marker_status, markers_in_rule)
  RRM1_IHC <- convert_IHC_result(as.character(marker_status$status[marker_status$biomarker_tech == 'RRM1_IHC']))
  if(RRM1_IHC == 'Positive'){
    recom <- -1
  }else if(RRM1_IHC == 'Negative'){
    recom <- 1
  }else{
    recom <- 0
  }
  z <- data.frame(DrugLineage = 'Ovarian Surface Epithelial Carcinoma (EOC)',
                  drugRule = 'Antimetabolites (gemcitabine)',
                  recom = recom,
                  version = dr_gemcit_version)
}