## split only experimental arm into responder and non-responder
rm(list = ls())
source('read_MAESTRO.R')
source('simulate_functions.R')
source('summary_functions.R')
library(lattice)
library(survival)
library(dplyr)

now <- format(Sys.time(), '%Y%m%d-%H%M%S')

df_OS <- data.frame(
  time = as.numeric(dfMA$OS),
  event = ifelse(dfMA$OSC_1_censored == 1, 0, 1),
  Group = factor(dfMA$Group, levels = c('PLACEBO', 'TH-302'))
)

paramSets <- generateParameters(
  surv_cutoff = c(100, 150, 200), # cutoffs for determining true label
  arms = c('single', 'double', 'combine'), # arms = c('single', 'double', 'combine')
  sesp = list(
    c(0.95, 0.1),  
    c(0.95, 0.2), c(0.9, 0.2), 
    c(0.95, 0.3), c(0.9, 0.3), c(0.8, 0.3),
    c(0.95, 0.4), c(0.9, 0.4), c(0.8, 0.4), c(0.7, 0.4), 
    c(0.95, 0.5), c(0.9, 0.5), c(0.8, 0.5), c(0.7, 0.5), c(0.6, 0.5), 
    c(0.95, 0.6), c(0.9, 0.6), c(0.8, 0.6), c(0.7, 0.6), c(0.6, 0.6), c(0.5, 0.6), 
    c(0.95, 0.7), c(0.9, 0.7), c(0.8, 0.7), c(0.7, 0.7), c(0.6, 0.7), c(0.5, 0.7), c(0.4, 0.7),
    c(0.95, 0.8), c(0.9, 0.8), c(0.8, 0.8), c(0.7, 0.8), c(0.6, 0.8), c(0.5, 0.8), c(0.4, 0.8), c(0.3, 0.8), 
    c(0.95, 0.9), c(0.9, 0.9), c(0.8, 0.9), c(0.7, 0.9), c(0.6, 0.9), c(0.5, 0.9), c(0.4, 0.9), c(0.3, 0.9), c(0.2, 0.9),
    c(0.95, 0.95), c(0.9, 0.95), c(0.8, 0.95), c(0.7, 0.95), c(0.6, 0.95), c(0.5, 0.95), c(0.4, 0.95), c(0.3, 0.95), c(0.2, 0.95), c(0.1, 0.95)
  ), # sensitivity/specificity for experimental arm, regardless of arms parameter.
  sp_plb = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2) # sensitivity/specificity for placebo arm (rate between placebo and experimental arm), when arms == 'double'.
  #pos_rate = c(0.5) # positive rate for placebo arm, when arms == 'single'
)

df_param <- do.call(rbind, paramSets)

system.time(mySimulation <- performSimulation(
  df_OS, paramSets, r = 100
))
format(object.size(mySimulation), units = 'Mb')

simSum <- summarizeSimulation(mySimulation)

currentDir <- getwd()
newDir <- paste0(currentDir, '/', now, '_plots')
dir.create(newDir)
setwd(newDir)
write.csv(simSum, paste0(now, '_simulation1.csv'), row.names = F)
plotSampleSurv(df_OS, mySimulation, plot_type = 'R')
plotSampleSurv(df_OS, mySimulation, plot_type = 'both')
save(mySimulation, file = paste0(now, '_mySimulation.Rda'))
setwd(currentDir)

