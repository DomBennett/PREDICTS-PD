# 27/04/2015
# Team PREDICTS-PD
# Plot communities by site

# START
cat (paste0 ('\nPlots started at [', Sys.time (), ']\n'))

# LIBS
source (file.path ('tools', 'tree_tools.R'))
source (file.path ('tools', 'community_tools.R'))
source (file.path ('tools', 'plotting_tools.R'))

# DIRS
predicts.dir <- file.path ('0_data', 'PREDICTS-DATA')
pglt.dir <- '0_pglt'
output.dir <- 'B4_plots'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\nReading in data ....')
# read in RDS
predicts.data <- readRDS (file.path (predicts.dir,
                                     'diversity-2014-10-29-03-40-20.rds'))
predics.data <- as.data.frame (predicts.data, stringsAsFactors=FALSE)
# get SSID
predicts.data$SSID <- paste0 (predicts.data$Source_ID, '_',
                              predicts.data$Study_number, '_',
                              predicts.data$Class)
studies <- unique (predicts.data$SSID)
cat ('\nDone.')

# PROCESS
cat ('\nReading in trees and plotting ....')
study.counter <- 0
for (i in 1:length (studies)) {
  # progress
  study <- studies[i]
  cat ('\n.... study [', study, '], [', i, '/', length (studies), ']',
       sep = '')
  fpath <- file.path (pglt.dir, study, '4_phylogeny', 'consensus.tre')
  tree <- suppressWarnings (try (read.tree (fpath), silent = TRUE))
  if (class (tree) == 'try-error') {
    next
  }
  tree$tip.label <- gsub ('_', ' ', tree$tip.label)
  # get study data
  study.data  <- predicts.data[predicts.data$SSID == study, ]
  # extract community matrix
  cmatrix <- getCommunityMatrix (study.data)
  # open pdf
  pdf (file.path (output.dir, paste0 (study, '.pdf')))
  # set margin area, 5 space at bottom for legend
  par(mar= c (5, 0, 0, 0) + 0.1, xpd=FALSE)
  # find site indexes
  site.index <- cumsum (table (study.data$Site_number))
  # get uses, make sure in order and complot abundance and incidence
  uses <- as.vector (study.data$Use_intensity[site.index])
  cmatrix.uses <- cmatrix[order (uses), ]
  uses <- uses[order (uses)]
  commplot (cmatrix=cmatrix, tree=tree, groups=uses, no.margin=FALSE)
  plotLegend (groups=uses)
  commplot (cmatrix=cmatrix>0, tree=tree, groups=uses, no.margin=FALSE)
  plotLegend (groups=uses)
  # repeat for prim habs
  habs <- as.vector (study.data$Predominant_habitat[site.index])
  cmatrix.habs <- cmatrix[order (habs), ]
  habs <- habs[order (habs)]
  commplot (cmatrix=cmatrix, tree=tree, groups=habs, no.margin=FALSE)
  plotLegend (groups=habs)
  commplot (cmatrix=cmatrix>0, tree=tree, groups=habs, no.margin=FALSE)
  plotLegend (groups=habs)
  dev.off ()
  study.counter <- study.counter + 1
}
cat ('\nDone. Plotted for [', study.counter, '] studies.', sep='')

# FINISH
cat (paste0 ('\nStage 7 finished at [', Sys.time (), ']'))