## Platinum compounds
dr_platinum <- function(marker_status){
  #marker_status <- data.frame(biomarker_tech = c('BRCA1_SEQ', 'BRCA2_SEQ', 'ERCC1_IHC'),
  #                            status = c('Mutated - Other', 'No Result', 'Negative'))
  dr_platinum_version <- 'ver0'
  
  markers_in_rule <- c('BRCA1_SEQ', 'BRCA2_SEQ', 'ERCC1_IHC')
  marker_status <- match_markers(marker_status, markers_in_rule)
  BRCA1_SEQ <- convert_SEQ_result(as.character(marker_status$status[marker_status$biomarker_tech == 'BRCA1_SEQ']))
  BRCA2_SEQ <- convert_SEQ_result(as.character(marker_status$status[marker_status$biomarker_tech == 'BRCA2_SEQ']))
  ERCC1_IHC <- convert_IHC_result(as.character(marker_status$status[marker_status$biomarker_tech == 'ERCC1_IHC']))
  if((BRCA1_SEQ == 'Mutation Not Detected' & BRCA2_SEQ == 'Mutation Not Detected') | (ERCC1_IHC == 'Positive')){
    recom <- -1
  }else if(BRCA1_SEQ == 'Nodata' & BRCA2_SEQ == 'Nodata' & ERCC1_IHC == 'Nodata'){
    recom <- 0
  }else{
    recom <- 1
  }
  z <- data.frame(DrugLineage = 'Ovarian Surface Epithelial Carcinoma (EOC)',
                  drugRule = 'Platinum compounds',
                  recom = recom,
                  version = dr_platinum_version)
}