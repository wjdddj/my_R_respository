## Taxanes
dr_taxanes <- function(marker_status){
  #marker_status <- data.frame(biomarker_tech = c('BRCA1_SEQ', 'BRCA2_SEQ'),
  #                            status = c('Mutated - Other', 'No Result'))  
  dr_taxanes_version <- 'ver0'
  
  markers_in_rule <- c('TLE3_IHC', 'TUBB3_IHC')
  marker_status <- match_markers(marker_status, markers_in_rule)
  TLE3_IHC <- convert_IHC_result(as.character(marker_status$status[marker_status$biomarker_tech == 'TLE3_IHC']))
  TUBB3_IHC <- convert_IHC_result(as.character(marker_status$status[marker_status$biomarker_tech == 'TUBB3_IHC']))
  
  if(TLE3_IHC == 'Positive'){
    recom <- 1
  }else if(TUBB3_IHC == 'Negative'){
    recom <- -1
  }else{
    recom <- 0
  }
  z <- data.frame(DrugLineage = 'Ovarian Surface Epithelial Carcinoma (EOC)',
                  drugRule = 'Taxanes',
                  recom = recom,
                  version = dr_taxanes_version)
}