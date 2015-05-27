# 01/11/2014
# Team PREDICTS-PD
# Parse published and pG-lt trees

# START
cat (paste0 ('\nStage 4 started at [', Sys.time (), ']'))

# PARAMETERS
use.unconstrained <- FALSE

# LIBS
source (file.path ('tools', 'tree_tools.R'))

# DIRS
data.dir <- '0_data'
pglt.dir <- '2_pGltrun'
mapped.dir <- '3_map'
output.dir <- '4_parse'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\n Reading in pglt and mapped trees for each study ....')
# counters
pglt.counter <- map.counter <- 0
studies <- list.files ('1_pGltsetup')
trees <- list ()
for (i in 1:length (studies)) {
  study <- studies[i]
  cat ('\n.... study [', study, '], [', i, '/', length (studies), ']',
       sep = '')
  # init container
  study.tree <- list ()
  # pglt-trees
  if (use.unconstrained) {
    fpath <- file.path (pglt.dir, study, '4_phylogeny', 'distribution_unconstrained.tre')
  } else {
    fpath <- file.path (pglt.dir, study, '4_phylogeny', 'distribution.tre')
  }
  pglt.trees <- suppressWarnings (try (read.tree (fpath), silent = TRUE))
  if (class (pglt.trees) != 'try-error') {
    cat ('\n.... attempting to rate-smooth')
    pglt.trees <- runChronos (pglt.trees)
    study.tree['pglt.trees'] <- list (pglt.trees)
    pglt.counter <- pglt.counter + 1
  }
  # mapped trees
  fpath <- file.path (mapped.dir, paste0 (study, '.tre'))
  mapped.trees <- suppressWarnings (try (read.tree (fpath), silent = TRUE))
  if (class (mapped.trees) != 'try-error') {
    study.tree['mapped.trees'] <- list (mapped.trees)
    map.counter <- map.counter + 1
  }
  if (length (study.tree) > 0) {
    # add study container to all container
    trees[study] <- list (study.tree)
  }
}
cat ('\nDone.')

# OUTPUT
cat ('\nWriting out ....')
for (i in 1:length (trees)) {
  pglt.trees <- trees[[i]][['pglt.trees']]
  mapped.trees <- trees[[i]][['mapped.trees']]
  folder.path <- file.path (output.dir, paste0 (names (trees)[i]))
  if (!file.exists (folder.path)) {
    dir.create (folder.path)
  }
  if (!is.null (mapped.trees)) {
    write.tree (file = file.path (folder.path, 'mapped.tre'), phy = mapped.trees)
  }
  if (!is.null (pglt.trees)) {
    write.tree (file = file.path (folder.path, 'pglt.tre'), phy = pglt.trees)
  }
}
cat ('Done. Discovered [', length (studies),
     '] studies of which [', pglt.counter,
     '] had pglt distirbutions and [', map.counter ,
     '] had mapped distributions', sep='')

# FINISH
cat (paste0 ('\nStage 4 finished at [', Sys.time (), ']'))