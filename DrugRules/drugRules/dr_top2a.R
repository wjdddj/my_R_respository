## Anthracyclines and related substances
dr_top2a <- function(marker_status){
  dr_top2a_version <- 'ver0'
  
  markers_in_rule <- c('TOP2A_IHC')
  marker_status <- match_markers(marker_status, markers_in_rule)
  TOP2A_IHC <- convert_IHC_result(as.character(marker_status$status[marker_status$biomarker_tech == 'TOP2A_IHC']))
  if(TOP2A_IHC == 'Positive'){
    recom <- 1
  }else if(TOP2A_IHC == 'Negative'){
    recom <- -1
  }else{
    recom <- 0
  }
  z <- data.frame(DrugLineage = 'Ovarian Surface Epithelial Carcinoma (EOC)',
                  drugRule = 'Anthracyclines and related substances',
                  recom = recom,
                  version = dr_top2a_version)
}