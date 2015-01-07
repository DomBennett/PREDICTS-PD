# 01/11/2014
# Team PREDICTS-PD
# Calculate phylo metrics per site

# START
cat (paste0 ('\nStage 5 started at [', Sys.time (), ']'))

# LIBS
# always source before dplyr
source (file.path ('functions', 'tools.R'))
library (dplyr)
library (ape)

# DIRS
input.dirs <- c ('0_data', '3_parse')
output.dir <- '5_metrics'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\nReading in data ....')
# read in RDS
predicts.data <- readRDS (file.path (input.dirs[1], 'diversity-2014-10-29-03-40-20.rds'))
# convert to dplyr format
predicts.data <- tbl_df (predicts.data)
# find appropriate sources, stick Source_ID and Study_number together
predicts.data$SSID <- paste0(predicts.data$Source_ID, '_', predicts.data$Study_number)
# read in trees
all.trees <- readInTrees (folder = input.dirs[2])
cat ('\nDone.')

# PROCESS
cat ('\nCalculating PD estimates ....')
study.counter <- 0
res <- tbl_df (data.frame ())
for (i in 1:length (studies)) {
  # progress
  study <- studies[i]
  cat ('\n.... study [', study, '], [', i, '/', length (studies), ']',
       sep = '')
  
  # get study data
  study.data  <- filter (predicts.data, SSID == study)
  
  # get tree distribution
  treedist <- all.trees[study][[1]]
  
  # extract community matrix, necessary for commPD
  cmatrix <- getCommunityMatrix(study.data)
  
  # # estimate spp dropped
  tip.labels <- getNames (treedist)
  pdropped <- sum (!colnames (cmatrix) %in% tip.labels)/
    ncol(cmatrix)
  
  # calc PD by site for all trees in dist
  multi.comm.pds <- multiCommPD (phylos=treedist, comm.data=cmatrix,
                                 type=2, min.spp=0)
  
  # get mean and sd PD per site
  comm.pds <- data.frame (mean = apply (multi.comm.pds, 2, mean, na.rm =TRUE),
                     sd = apply (multi.comm.pds, 2, sd, na.rm =TRUE))
  
  # add results to study.data (this is pretty ugly... oh well)
  site.counts <- table (select (study.data, Site_number)[ ,1])
  study.data$Est_mean_PD <- rep (comm.pds[ ,'mean'], site.counts)
  study.data$Est_sd_PD <- rep (comm.pds[ ,'sd'], site.counts)
  study.data$PD_pdropped <- pdropped
  
  # bind study.data to res
  res <- rbind (res, study.data)
}
cat ('\nDone. Calculated PD estimates for [', study.counter, '] studies.',
     sep='')

# OUTPUT
saveRDS (res, file=file.path(output.dir, 'predictsdata_wpd.rds'))

# FINISH
cat (paste0 ('\nStage 5 finished at [', Sys.time (), ']'))