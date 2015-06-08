# 01/11/2014
# Team PREDICTS-PD
# Calculate phylo metrics per site

# START
cat (paste0 ('\nStage 6 started at [', Sys.time (), ']\n'))

# LIBS
source (file.path ('tools', 'tree_tools.R'))
source (file.path ('tools', 'community_tools.R'))

# DIRS
input.dirs <- c ('0_data', '4_parse')
predicts.dir <- file.path ('0_data', 'PREDICTS-DATA')
output.dir <- '6_metrics'
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
trees <- readInTrees (folder=input.dirs[2], recursive=TRUE)
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
                                 metric=metric)
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
mean.pglt.ntaxa <- mean (ntaxa[['pglt']], na.rm=TRUE)
mean.mapped.ntaxa <- mean (ntaxa[['mapped']], na.rm=TRUE)
ntaxa.predicts <- sum (ntaxa[['predicts']], na.rm=TRUE)
ntaxa.mapped <- sum (ntaxa[['mapped']], na.rm=TRUE)
ntaxa.pglt <- sum (ntaxa[['pglt']], na.rm=TRUE)
pdf (file.path(output.dir, 'names.pdf'))
barplot (c (ntaxa.pglt, ntaxa.mapped, ntaxa.predicts),
         names.arg=c ('pG-lt', 'Mapped', 'PREDICTS'),
         col=rainbow(3), ylab='N names')
barplot (c (mean.pglt.ntaxa, mean.mapped.ntaxa),
         names.arg=c ('pG-lt', 'Mapped'),
         col=rainbow(2), ylab='Mean N names per study')
dev.off ()
cat ('\nDone. Calculated PD estimates for [',
     study.counter, '] studies. Found [', ntaxa.predicts,
     '] names of which [', ntaxa.pglt,
     '] were represented in pglt trees and [', ntaxa.mapped,
     '] were represented in mapped trees with mean [', mean.pglt.ntaxa,
     '] and [', mean.mapped.ntaxa, '] number of taxa per study, respectively.',  sep='')

# OUTPUT
saveRDS (res, file=file.path(output.dir, 'predictsdata_wpd.rds'))

# FINISH
cat (paste0 ('\nStage 6 finished at [', Sys.time (), ']'))