## DB info, Column Names Follows Bioserv Formatting
## AccessionNumber as the primary Key

source('~/R_modules/caris_basic/caris_basic.R')
source('~/R_modules/ADAPT_tissue/ADAPT_tissue.R')

################################################################################################################
################################################################################################################
getBiomarkerSResultFromDB <- function(
  AccessionNumbers, Biomarker, Technology, outBiomarkerName, fuzzyMatch = F # very dangerous to set fuzzyMatch to T
){
  ## query from bioserv:bioinfo.all_test_cases
  mySQLAccessionNumber <- formQueryItemArray(AccessionNumbers)
  if(!fuzzyMatch){
    Q = paste0(
      "SELECT AccessionNumber, SResult as ", outBiomarkerName, " \n ",
      "FROM bioinfo.all_test_cases \n ",
      "WHERE Test = '", Biomarker, "' and technology = '", Technology, "' and AccessionNumber in ", mySQLAccessionNumber, ";"
    )
  }else{
    warning('DANGEROUS! set fuzzyMatch = T; Please make sure db contains the right Test names.')
    Q = paste0(
      "SELECT AccessionNumber, SResult as ", outBiomarkerName, " \n ",
      "FROM bioinfo.all_test_cases \n ",
      "WHERE Test like '%", Biomarker, "%' and technology = '", Technology, "' and AccessionNumber in ", mySQLAccessionNumber, ";"
    )
  }
  z <- get_sting_query('bioinfo', Q)
}

################################################################################################################
# IHC and IHC_IA should be mutually exclusive; When conflicted, IHC take priority
################################################################################################################
getCombinedIHCfromIA <- function(
  AccessionNumbers, Biomarker, outBiomarkerName
){

  IAcolName <- paste0(outBiomarkerName, '_IA')
  resultIHCIA <- getBiomarkerSResultFromDB(
    AccessionNumbers, 
    Biomarker = paste0(Biomarker, ' IHC-IA'), Technology = 'IHC-IA', outBiomarkerName = IAcolName,
    fuzzyMatch = F
  )
  
  nonIAcolName <- paste0(outBiomarkerName, '_nonIA')
  resultIHC <- getBiomarkerSResultFromDB(
    AccessionNumbers, 
    Biomarker = Biomarker, Technology = 'IHC', outBiomarkerName = nonIAcolName,
    fuzzyMatch = F
  )
  
  z <- merge(
    resultIHCIA,
    resultIHC,
    by = 'AccessionNumber',
    all = T
  )
  # IHC and IHC_IA should be mutually exclusive; When conflicted, IHC take priority
  z <- data.frame(
    z,
    ifelse(
      is.na(z[, nonIAcolName]), z[, IAcolName], z[, nonIAcolName]
    )
  )
  colnames(z)[4] <- outBiomarkerName
  z
}


################################################################################################################
################################################################################################################
getHER2ISHFromDB <- function(
  AccessionNumbers
){
  # AccessionNumbers = clinDatasetDF_20170406$AccessionNumber
  HER2_FISH <- getBiomarkerSResultFromDB(
    AccessionNumbers, 
    Biomarker = 'Her2/Neu FISH', Technology = 'FISH', outBiomarkerName = 'HER2_FISH',
    fuzzyMatch = F
  )
  HER2_CISH <- getBiomarkerSResultFromDB(
    AccessionNumbers, 
    Biomarker = 'Her2 CISH', Technology = 'CISH', outBiomarkerName = 'HER2_CISH',
    fuzzyMatch = F
  )
  z <- merge(
    HER2_FISH,
    HER2_CISH,
    by = 'AccessionNumber',
    all = T
  )
  # FISH and CISH should be mutually exclusive; When conflicted, FISH take priority
  z$HER2_ISH <- ifelse(
    is.na(z$HER2_FISH), z$HER2_CISH, z$HER2_FISH
  )
  z
}

################################################################################################################
################################################################################################################
getIHCIntensityFromDB <- function(
  AccessionNumbers, Biomarker, outBiomarkerName, fuzzyMatch = F
){
  ## query from bioserv:bioinfo.ihc_analysis_reports
  mySQLAccessionNumber <- formQueryItemArray(AccessionNumbers)
  if(!fuzzyMatch){
    Q = paste0(
      "SELECT AccessionNumber as AccessionNumber, \n",
      "StainingIntensity as ", outBiomarkerName, "_Intenstiy, \n",
      "PercentStaining as ", outBiomarkerName, "_PercentStaining \n",
      "FROM bioinfo.ihc_analysis_reports
      WHERE Biomarker = '", Biomarker, "' and AccessionNumber in ", mySQLAccessionNumber, ";"
    )
  }else{
    warning('DANGEROUS! set fuzzyMatch = T; Please make sure db contains the right Test names.')
    Q = paste0(
      "SELECT AccessionNumber as AccessionNumber, \n",
      "StainingIntensity as ", outBiomarkerName, "_Intenstiy, \n",
      "PercentStaining as ", outBiomarkerName, "_PercentStaining \n",
      "FROM bioinfo.ihc_analysis_reports
      WHERE Biomarker like '%", Biomarker, "%' and AccessionNumber in ", mySQLAccessionNumber, ";"
    )
  }
  z <- get_sting_query('bioinfo', Q)
}

################################################################################################################
################################################################################################################
getFISHIntensityFromDB <- function(
  AccessionNumbers, Biomarker, outBiomarkerName, fuzzyMatch = F
){
  ## query from bioserv:bioinfo.fish_analysis_reports
  mySQLAccessionNumber <- formQueryItemArray(AccessionNumbers)
  if(!fuzzyMatch){
    Q = paste0(
      "SELECT AccessionNumber as AccessionNumber, \n",
      "`Total/Avg Gene Copy Number` as ", outBiomarkerName, "_Gene_Copy, \n",
      "`Total/Avg Control Copy Number` as ", outBiomarkerName, "_Control_Copy \n",
      "FROM bioinfo.fish_analysis_reports
      WHERE Gene = '", Biomarker, "' and AccessionNumber in ", mySQLAccessionNumber, ";"
    )
  }else{
    warning('DANGEROUS! set fuzzyMatch = T; Please make sure db contains the right Test names.')
    Q = paste0(
      "SELECT AccessionNumber as AccessionNumber, \n",
      "`Total/Avg Gene Copy Number` as ", outBiomarkerName, "_Gene_Copy, \n",
      "`Total/Avg Control Copy Number` as ", outBiomarkerName, "_Control_Copy \n",
      "FROM bioinfo.fish_analysis_reports
      WHERE Gene like '%", Biomarker, "%' and AccessionNumber in ", mySQLAccessionNumber, ";"
    )
  }
  z <- get_sting_query('bioinfo', Q)
}

################################################################################################################
################################################################################################################
getCISHIntensityFromDB <- function(
  AccessionNumbers, Biomarker, outBiomarkerName, fuzzyMatch = F
){
  ## query from bioserv:bioinfo.fish_analysis_reports
  mySQLAccessionNumber <- formQueryItemArray(AccessionNumbers)
  if(!fuzzyMatch){
    Q = paste0(
      "SELECT AccessionNumber as AccessionNumber, \n",
      "`Total/Avg Gene Copy Number` as ", outBiomarkerName, "_Gene_Copy, \n",
      "`Total/Avg Control Copy Number` as ", outBiomarkerName, "_Control_Copy \n",
      "FROM bioinfo.cish_analysis_reports
      WHERE Gene = '", Biomarker, "' and AccessionNumber in ", mySQLAccessionNumber, ";"
    )
  }else{
    warning('DANGEROUS! set fuzzyMatch = T; Please make sure db contains the right Test names.')
    Q = paste0(
      "SELECT AccessionNumber as AccessionNumber, \n",
      "`Total/Avg Gene Copy Number` as ", outBiomarkerName, "_Gene_Copy, \n",
      "`Total/Avg Control Copy Number` as ", outBiomarkerName, "_Control_Copy \n",
      "FROM bioinfo.cish_analysis_reports
      WHERE Gene like '%", Biomarker, "%' and AccessionNumber in ", mySQLAccessionNumber, ";"
    )
  }
  z <- get_sting_query('bioinfo', Q)
}

#########################################################################################################
## calculate HER2 status from IHC and ISH results
## input df must contain column "HER2_ISH", "HER2_IHC"
#########################################################################################################
calculateHER2StatusFromISHIHC <- function(df){
  #df <- df_train3
  HER2_Status <- ifelse(
    df$HER2_IHC == 'Positive' | df$HER2_ISH == 'Positive',
    'Positive',
    ifelse(
      df$HER2_IHC == 'Negative' | df$HER2_ISH == 'Negative',
      'Negative',
      ifelse(
        df$HER2_IHC == 'Equivocal' | df$HER2_ISH == 'Equivocal',
        'Equivocal',
        NA
      )
    )
  )
  HER2_Status
}


#########################################################################################################
## calculate HER2 IHC from IHC intensity and percent
## input df must contain column "HER2_IHC_Intensity", "HER2_IHC_Percent"
#########################################################################################################
calculateHER2IHCfromScore <- function(df){
  #df <- df_train3
  HER2_IHC = ifelse(
    df$HER2_IHC_Intensity == 3 & df$HER2_IHC_Percent >= 10,
    'Positive', 
    ifelse(
      df$HER2_IHC_Intensity == 2 & df$HER2_IHC_Percent >= 10,
      'Equivocal',
      'Negative'
    )
  )
  
  HER2_IHC
}






