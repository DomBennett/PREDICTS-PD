# 01/11/2014
# Team PREDICTS-PD
# Calculate phylo metrics per site

# START
cat (paste0 ('\nMetrics started at [', Sys.time (), ']\n'))

# LIBS
source (file.path ('tools', 'tree_tools.R'))
source (file.path ('tools', 'community_tools.R'))

# PARAMETERS
#normalise <- 'age'  # how to normalise the branch lengths (age, total, none)

# DIRS
parse.dir <- 'B1_parse.R'
predicts.dir <- file.path ('0_data', 'PREDICTS-DATA')
output.dir <- 'B3_metrics'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\nReading in data ....')
# read in RDS
predicts.data <- readRDS (file.path (predicts.dir,
                                     'diversity-2014-10-29-03-40-20.rds'))
# get SSID
predicts.data$SSID <- paste0 (predicts.data$Source_ID, '_',
                              predicts.data$Study_number, '_',
                              predicts.data$Class)
# read in trees
trees <- readInTrees (folder=parse.dir, recursive=TRUE)
cat ('\nDone.')

# PROCESS
cat ('\nCalculating PD estimates ....')
study.counter <- 0
ntaxa <- list ('predicts'=rep (NA, length (trees)),
               'mapped'=rep (NA, length (trees)),
               'pglt'=rep (NA, length (trees)))
res <- data.frame ()
for (i in 1:length (trees)) {
  # progress
  study <- names (trees)[i]
  cat ('\n.... study [', study, '], [', i, '/', length (trees), ']',
       sep = '')
  # get study data
  study.data  <- predicts.data[predicts.data$SSID == study, ]
  # extract community matrix, necessary for commPD
  cmatrix <- getCommunityMatrix(study.data)
  ntaxa[['predicts']][i] <- ncol (cmatrix)
  # calculate stats for each set of trees
  phymets.res <- list ()
  for (j in 1:length (trees[[i]])) {
    treedist <- trees[[i]][[j]]
    prefix <- names (trees[[i]])[j]
    phymets.res[[prefix]] <- list ()
    # estimate spp dropped
    tip.labels <- getNames (treedist)
    ntaxa[[prefix]][i] <- length (tip.labels)
    pdropped <- sum (!colnames (cmatrix) %in% tip.labels)/
      ncol(cmatrix)
    phymets.res[[prefix]][['pdropped']] <- pdropped
    # calc phylogenetic metrics per site for all trees in dist
    for (metric in c( 'PD1', 'PSV', 'PSE')) {
      multi <- multiCommPhyMets (trees=treedist, cmatrix=cmatrix,
                                 metric=metric, normalise=normalise)
      # get mean and sd per site
      phymets.res[[prefix]][[paste0 (metric, '_mean')]] <-
        apply (multi, 2, mean, na.rm =TRUE)
      phymets.res[[prefix]][[paste0 (metric, '_sd')]] <-
        apply (multi, 2, sd, na.rm =TRUE)
    }
  }
  # add results to study.data
  # TODO: is there a better way of doing this?
  #  This is something I don't like about R
  site.counts <- table (study.data$Site_number)
  if ('pglt' %in% names (phymets.res)) {
    study.data$pglt_mean_PD <-
      rep (phymets.res$pglt$PD1_mean, site.counts)
    study.data$pglt_sd_PD <-
      rep (phymets.res$pglt$PD1_sd, site.counts)
    study.data$pglt_mean_PSV <-
      rep (phymets.res$pglt$PSV_mean, site.counts)
    study.data$pglt_sd_PSV <-
      rep (phymets.res$pglt$PSV_sd, site.counts)
    study.data$pglt_mean_PSE <-
      rep (phymets.res$pglt$PSE_mean, site.counts)
    study.data$pglt_sd_PSE <-
      rep (phymets.res$pglt$PSE_sd, site.counts)
    study.data$pglt_pdropped <-
      phymets.res$pglt$pdropped
  } else {
    study.data$pglt_mean_PD <- rep (NA, nrow (study.data))
    study.data$pglt_sd_PD <- rep (NA, nrow (study.data))
    study.data$pglt_mean_PSV <- rep (NA, nrow (study.data))
    study.data$pglt_sd_PSV <- rep (NA, nrow (study.data))
    study.data$pglt_mean_PSE <- rep (NA, nrow (study.data))
    study.data$pglt_sd_PSE <- rep (NA, nrow (study.data))
    study.data$pglt_pdropped <- rep (NA, nrow (study.data))
  }
  if ('mapped' %in% names (phymets.res)) {
    study.data$mapped_mean_PD <-
      rep (phymets.res$mapped$PD1_mean, site.counts)
    study.data$mapped_sd_PD <-
      rep (phymets.res$mapped$PD1_sd, site.counts)
    study.data$mapped_mean_PSV <-
      rep (phymets.res$mapped$PSV_mean, site.counts)
    study.data$mapped_sd_PSV <-
      rep (phymets.res$mapped$PSV_sd, site.counts)
    study.data$mapped_mean_PSE <-
      rep (phymets.res$mapped$PSE_mean, site.counts)
    study.data$mapped_sd_PSE <-
      rep (phymets.res$mapped$PSE_sd, site.counts)
    study.data$mapped_pdropped <-
      phymets.res$mapped$pdropped
  } else {
    study.data$mapped_mean_PD <- rep (NA, nrow (study.data))
    study.data$mapped_sd_PD <- rep (NA, nrow (study.data))
    study.data$mapped_mean_PSV <- rep (NA, nrow (study.data))
    study.data$mapped_sd_PSV <- rep (NA, nrow (study.data))
    study.data$mapped_mean_PSE <- rep (NA, nrow (study.data))
    study.data$mapped_sd_PSE <- rep (NA, nrow (study.data))
    study.data$mapped_pdropped <- rep (NA, nrow (study.data))
  }
  # bind study.data to res
  res <- rbind (res, study.data)
  study.counter <- study.counter + 1
}
ntaxa.predicts <- sum (ntaxa[['predicts']], na.rm=TRUE)
ntaxa.mapped <- sum (ntaxa[['mapped']], na.rm=TRUE)
ntaxa.pglt <- sum (ntaxa[['pglt']], na.rm=TRUE)
pp.pglt <-mean (ntaxa[['pglt']]*100/ntaxa[['predicts']], na.rm=TRUE)
pp.mapped <-mean (ntaxa[['mapped']]*100/ntaxa[['predicts']], na.rm=TRUE)
pglt.more <- sum (ntaxa$pglt - ntaxa$mapped, na.rm=TRUE)
pglt.more.p <- mean (ntaxa$pglt/ntaxa$mapped, na.rm=TRUE)
# work out total names between both
pglt.pull <- ntaxa$pglt > ntaxa$mapped
pglt.pull[is.na (pglt.pull)] <- FALSE
pglt.pull[is.na(ntaxa$mapped)] <- TRUE
tot.taxa <- sum (c (ntaxa$pglt[pglt.pull], ntaxa$mapped[!pglt.pull]))
cat ('\nDone. Calculated PD estimates for [',
     study.counter, '] studies. Found [', ntaxa.predicts,
     '] PREDICTS names of which [', ntaxa.pglt,
     '] were represented in pglt trees and [', ntaxa.mapped,
     '] were represented in mapped trees.
     For studies that shared both mapped and pglt trees there were [', pglt.more,
     '] more in pglt trees, represeneting [', pglt.more.p, 'x] more names on average.
     Mean [', pp.pglt, '%] of study names in studies for pglt trees.
     Mean [', pp.mapped, '%] of study names in studies for mapped trees.
     In total, a max [', tot.taxa, '] names is repsented by both.',  sep='')

# OUTPUT
outfile <- paste0 ('predictsdata_wpd_', normalise, '.rds')
saveRDS (res, file=file.path(output.dir, outfile))

# FINISH
cat (paste0 ('\nStage 6 finished at [', Sys.time (), ']'))