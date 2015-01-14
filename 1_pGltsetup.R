# 01/11/2014
# Team PREDICTS-PD
# Create setup files for pG-lt

# START
cat (paste0 ('\nStage 1 started at [', Sys.time (), ']'))

# LIBS
# source before dplyr! Else, a plyr and dplyr conflict.
source (file.path ('functions', 'tools.R'))
library (dplyr)

# PARAMETERS
easy.names <- TRUE  # setup for 'easy' names?
non.sp.tol <- 0.5  # the maximum proportion of names resolved above the species level
non.sci.tol <- 0.5  # the maximum proportion of non-scientific names

# DIRS
input.dir <- '0_data'
output.dir <- '1_pGltsetup'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\nReading in data ....')
# read in RDS
x <- readRDS (file.path (input.dir, 'diversity-2014-10-29-03-40-20.rds'))  # TODO: what's the difference between this and sites-*.rds?
# x <- readRDS (file.path (input.dir, 'djb208-2014-08-04-08-40-03.rds'))  # urban data
# convert to dplyr format
x <- tbl_df (x)
cat ('\nDone.')

# MANIPULATE
cat ('\nManipulating data ....')
# find appropriate sources, stick Source_ID and Study_number together
x$SSID <- paste0(x$Source_ID, '_', x$Study_number)
sources <- group_by (x, SSID)
sources <- summarise (sources, N.names = n_distinct (Parsed_name))
# take only studies with more than 5 names
sources <- filter (sources, N.names >= 5)
cat ('\nDone.')

# EXTRACT NAMES AND SAVE
cat ('\nExtracting names ....')
#TODO: deal with studies holding the same names
study.counter <- 0
names.counter <- 0
for (i in 1:nrow (sources)) {
  # study ID
  study <- sources$SSID[i]
  
  # progress
  cat ('\n.... study [', study, '], [', i, '/', nrow(sources), ']', sep = '')
  
  # filter out those studies with too few species level names
  species <- as.character (select (filter (x, SSID == study), Species)[ ,1])
  # calc proportion of absent species names
  if (sum ('' == species)/length (species) > non.sp.tol) {
    next
  }
  
  # filter out non-scientific names (at the moment pglt can't handle common names)
  resolutions <- as.character (select (filter (x, SSID == study), Resolution_entered)[ ,1])
  # calc proportion of non-scientific names
  if (sum ('Scientific' != resolutions)/length (resolutions) > non.sci.tol) {
    next
  }
  
  
  # create taxonomic tbl_df
  taxon <- select (filter (x, SSID == study), Genus, Family, Order, Class)
  
  # EASY NAMES filter
  if (easy.names) {
    # only easy groups: Aves (class), Mammalia (class), Hymenoptera (class)
    continue <- FALSE
    if (all (taxon$Order == 'Hymenoptera')) {
      continue <- TRUE
    }
    if (all (taxon$Class == 'Mammalia')) {
      continue <- TRUE
    }
    if (all (taxon$Class == 'Aves')) {
      continue <- TRUE
    }
    if (!continue) {
      next
    }
  }
  
  # find lowest shared taxonomic group -- the parent
  for (j in 1:ncol (taxon)) {
    if (any (taxon[ ,j] == '') ||
          length (unique (taxon[ ,j])) > 1) {
      next
    }
    parent <- unique (taxon[ ,j])
    break
  }
  if (length (parent) > 1) {
    cat (paste0 ('\nMore than 1 higher taxon for [',
                 study,'] -- skipping'))
    next
  }
  
  # get parentid -- this is necessary for names resolution in pglt
  resolved <- taxaResolve (names = as.character (parent))
  parentid <- resolved$taxid
  if (length (parentid) > 1) {
    cat (paste0 ('\nMore than 1 higher taxon txid for [',
                 study,'] -- skipping'))
    next
  }
  
  # write out parameters.csv and names.txt for pglt
  parameters <- data.frame (Parameter = 'parentid',
                            Value = parentid, Note = parent)
  temp.dir <- file.path (output.dir, study)
  if (!file.exists (temp.dir)) {
    dir.create (temp.dir)
  }
  names <- select (filter (x, SSID == study), Best_guess_binomial)[ ,1]
  names <- unique (as.character (names))
  write.table (names, file = file.path (temp.dir, 'names.txt'),
               col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table (parameters, file = file.path (temp.dir, 'parameters.csv'),
               col.names = TRUE, row.names = FALSE, quote = FALSE,
               sep = ',')
  study.counter <- study.counter + 1
  names.counter <- names.counter + length (names)
}
cat ('Done. Setup for [', names.counter, '] names and [', study.counter, '] studies.',
     sep='')

# FINISH
cat (paste0 ('\nStage 1 finished at [', Sys.time (), ']'))