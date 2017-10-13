## Topo1 inhibitors
dr_topo1 <- function(marker_status){
  dr_topo1_version <- 'ver0'
  
  markers_in_rule <- c('TOPO1_IHC')
  marker_status <- match_markers(marker_status, markers_in_rule)
  TOPO1_IHC <- convert_IHC_result(as.character(marker_status$status[marker_status$biomarker_tech == 'TOPO1_IHC']))
  if(TOPO1_IHC == 'Positive'){
    recom <- 1
  }else if(TOPO1_IHC == 'Negative'){
    recom <- -1
  }else{
    recom <- 0
  }
  z <- data.frame(DrugLineage = 'Ovarian Surface Epithelial Carcinoma (EOC)',
                  drugRule = 'Topo1 inhibitors',
                  recom = recom,
                  version = dr_topo1_version)
}