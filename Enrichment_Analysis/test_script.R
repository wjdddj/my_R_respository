source('~/R_modules/Enrichment_Analysis/enrichment_analysis_v2.0.R')

## cumulative count percentile as cutoff
df_lib <- read_lib('enSeq_R7lib-HVHNHBCXX_S00_L00_R1_none.V.44046.SUM.csv', select = c(1,2))
head(df_lib)
df_lib[min(which(cumsum(df_lib$count)/sum(df_lib$count) > 0.8)), ]
cum_p <- 0.9

get_cutoff_by_cum_p <- function(df_lib, cum_p, sorted = T){
  if(sorted == F){
    df_lib <- df_lib[order(df_lib$count, decreasing = T), ]
  }
  cutoff <- df_lib$count[min(which(cumsum(df_lib$count)/sum(df_lib$count) > cum_p))]
}

get_cum_p_by_cutoff <- function(df_lib, cutoff){
  
}
