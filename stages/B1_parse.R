# 01/11/2014
# Team PREDICTS-PD
# Parse mapped and pG-lt trees

# START
if (use.unconstrained) {
  cat (paste0 ('\nParse (unconstrained) started at [', Sys.time (), ']'))
} else {
  cat (paste0 ('\nParse started at [', Sys.time (), ']'))
}

# PARAMETERS
#use.unconstrained <- FALSE

# LIBS (UNIX ONLY)
library (foreach)
library (doMC)
source (file.path ('tools', 'tree_tools.R'))

writeTree <- function (tree, filename, folder.path) {
  # Create folder and writeout
  if (!file.exists (folder.path)) {
    dir.create (folder.path)
  }
  write.tree (file=file.path (folder.path, filename), phy=tree)
}

# PARAMETERS
registerDoMC (ncpus)

# DIRS
setup.dir <- 'A3_pgltsetup'
pglt.dir <- '0_pglt'
mapped.dir <- 'A4_map'
if (use.unconstrained) {
  output.dir <- 'B1_parse_unconstrained'
} else {
  output.dir <- 'B1_parse'
}
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\nReading in pglt and mapped trees for each study ....')
studies <- list.files (setup.dir)
counter <- foreach (i=1:length (studies)) %dopar% {
  counter <- NULL  # record p or m
  study <- studies[i]
  cat ('\n.... study [', study, '], [', i, '/', length (studies), ']',
       sep = '')
  # output path
  folder.path <- file.path (output.dir, studies[i])
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
    pglt.trees <- runRateSmoother (pglt.trees, i)
    writeTree (pglt.trees, 'pglt.tre', folder.path)
    counter <- c (counter, 'p')
  }
  # mapped trees
  fpath <- file.path (mapped.dir, paste0 (study, '.tre'))
  mapped.trees <- suppressWarnings (try (read.tree (fpath), silent = TRUE))
  if (class (mapped.trees) != 'try-error') {
    writeTree (mapped.trees, 'mapped.tre', folder.path)
    counter <- c (counter, 'm')
  }
  counter
}
counter <- table (unlist (counter))
pglt.counter <- ifelse (is.na (counter['p']), 0, counter['p'])
map.counter <- ifelse (is.na (counter['m']), 0, counter['m'])
cat ('Done. Discovered [', length (studies),
     '] studies of which [', pglt.counter,
     '] had pglt distirbutions and [', map.counter ,
     '] had mapped distributions', sep='')

# FINISH
cat (paste0 ('\nB1 finished at [', Sys.time (), ']'))