# 01/11/2014
# Team PREDICTS-PD
# How does PD respond to human impacts?

# START
cat (paste0 ('\nStage 6 started at [', Sys.time (), ']'))

# LIBS
library (dplyr)

# DIRS
input.dir <- '5_metrics'
output.dir <- '6_analysis'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\nReading in data ....')
# read in RDS
predicts.data <- readRDS (file.path(input.dir, 'predictsdata_wpd.rds'))
# convert to dplyr format
predicts.data <- tbl_df (predicts.data)

# PROCESS
# how many species made it into the trees?
cat ('Mean [', 1 - mean (predicts.data$PD_pdropped), '] species in tree distritbution', sep = '')
# how does PD respond to use intensity?
plot (predicts.data$Est_mean_PD~predicts.data$Use_intensity)  # VERY crude plot

# OUTPUT

# FINISH
cat (paste0 ('\nStage 6 finished at [', Sys.time (), ']'))