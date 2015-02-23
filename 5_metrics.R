# 01/11/2014
# Team PREDICTS-PD
# Calculate phylo metrics per site

# START
cat (paste0 ('\nStage 5 started at [', Sys.time (), ']'))

# LIBS
source (file.path ('tools', 'tree_tools.R'))
source (file.path ('tools', 'community_tools.R'))

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
# find appropriate sources, stick Source_ID and Study_number together
predicts.data$SSID <- paste0(predicts.data$Source_ID, '_', predicts.data$Study_number)
# read in trees
all.trees <- readInTrees (folder = input.dirs[2])
cat ('\nDone.')

# PROCESS
cat ('\nCalculating PD estimates ....')
study.counter <- 0
res <- data.frame ()
for (i in 1:length (all.trees)) {
  # progress
  study <- names (all.trees)[i]
  cat ('\n.... study [', study, '], [', i, '/', length (all.trees), ']',
       sep = '')
  # get study data
  study.data  <- predicts.data[predicts.data$SSID == study, ]
  # get tree distribution
  treedist <- all.trees[[i]]
  # extract community matrix, necessary for commPD
  cmatrix <- getCommunityMatrix(study.data)
  # # estimate spp dropped
  tip.labels <- getNames (treedist)
  pdropped <- sum (!colnames (cmatrix) %in% tip.labels)/
    ncol(cmatrix)
  # calc phylogenetic metrics per site for all trees in dist
  phymets.res <- data.frame (PD1.mean=rep(NA, nrow(cmatrix)),
                             PD1.sd=rep(NA, nrow(cmatrix)),
                             PSV.mean=rep(NA, nrow(cmatrix)),
                             PSV.sd=rep(NA, nrow(cmatrix)),
                             PSE.mean=rep(NA, nrow(cmatrix)),
                             PSE.sd=rep(NA, nrow(cmatrix)))
  for (metric in c( 'PD1', 'PSV', 'PSE')) {
    multi <- multiCommPhyMets (trees=treedist, cmatrix=cmatrix, metric=metric)
    # get mean and sd per site
    phymets.res[ ,paste0 (metric, '.mean')] <-  apply (multi, 2, mean, na.rm =TRUE)
    phymets.res[ ,paste0 (metric, '.sd')] <-  apply (multi, 2, sd, na.rm =TRUE)
  }
  # add results to study.data
  site.counts <- table (study.data$Site_number)
  study.data$Est_mean_PD <- rep (phymets.res$PD1.mean, site.counts)
  study.data$Est_sd_PD <- rep (phymets.res$PD1.sd, site.counts)
  study.data$Est_mean_PSV <- rep (phymets.res$PSV.mean, site.counts)
  study.data$Est_sd_PSV <- rep (phymets.res$PSV.sd, site.counts)
  study.data$Est_mean_PSE <- rep (phymets.res$PSE.mean, site.counts)
  study.data$Est_sd_PSE <- rep (phymets.res$PSE.sd, site.counts)
  # add pdropped to study.data
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