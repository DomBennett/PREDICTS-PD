# 16/01/2015
# D.J. Bennett
# Determining the RTT stat: use the root to tip distances
#  to determine the 'goodness' of a phylogenetic tree


# LIBS
library (treebase)

# FUNCTIONS
calcRTTStat <- function (tree) {
  # Return RTT stat -- the coefficient of variation of the
  #  normalised root to tip distances
  rtt.dists <- diag (vcv.phylo (tree))  # calc rtts
  norm.rtt.dists <- rtt.dists/max (rtt.dists)  # normalise
  sd (norm.rtt.dists) / mean (norm.rtt.dists)  # return CoV
}
safeConnect <- function (expr, wait = 3600, trys = 10,
                         verbose = TRUE) {
  # Safely connect to online datasources using expr
  #  in large data pipelines without pipeline failing
  # I use an hour wait and ten trys because TreeBase can be
  #  temperamental
  for (i in 1:trys) {
    res <- try (expr, TRUE)
    if (any (class (res) == 'try-error')) {
      error.message <- geterrmessage ()
      if (verbose) {
        cat ('\n.... Connection error occurred:')
        cat (error.message)
        cat (paste0 ('\n.... waiting [', wait,'] seconds'))
      }
      Sys.sleep (wait)
    } else {
      return (res)
    }
  }
  stop ('Unable to connect')
}

# DEMONSTRATION WITH FAKE TREES
n <- 100
tree <- stree(n)  # fake tree
# 1. Best tree -- rtt stat = 0
tree$edge.length <- rep (1, n)
rtt.stat <- calcRTTStat(tree)
plot (tree, main = paste0 ('Best tree. RTT stat: ', rtt.stat))
# 2. Good tree -- rtt stat ~ 0
tree$edge.length <- rnorm(n, 100, 1)
rtt.stat <- calcRTTStat(tree)
plot (tree, main = paste0 ('Good tree. RTT stat: ', round (rtt.stat, 4)))
# 3. Bad tree -- rtt stat > 0.3
tree$edge.length <- rpois (n, lambda = 10)
rtt.stat <- calcRTTStat(tree)
plot (tree, main = paste0 ('Bad tree. RTT stat: ', round (rtt.stat, 4)))
# 3. Worst tree -- rtt stat > 1
tree$edge.length <- rpois (n, lambda = 0.5)
rtt.stat <- calcRTTStat(tree)
plot (tree, main = paste0 ('Worst tree. RTT stat: ', round (rtt.stat, 4)))

# CALCULATE RTT STAT REFERENCE (Note download takes a long time)
subsample <- 100  # how many published trees to download?
# search
cat ('\nSearching suitable trees ....')
meta <- safeConnect (expr = suppressWarnings (metadata ()))
# choose those with more than 5 taxa
meta <- meta[meta$ntaxa >= 5, ]
# choose only species trees
meta <- meta[meta$kind == 'Species Tree', ]
meta <- meta[sample (nrow (meta)), ]
# download
cat ('\nDownloading trees ....')
counter <- 0
i <- 0
trees <- list ()
while (counter < subsample) {
  i <- i + 1
  cat (paste0 ('\n.... tree [',counter,']\n'))
  tree.data <- meta[i, ]
  tree <- safeConnect (expr = {
    suppressWarnings (search_treebase (tree.data$Tree.id,
                                       by = 'id.tree',
                                       verbose = FALSE))})
  closeAllConnections ()
  if (length (tree) > 0) {
    tree <- tree[[1]]
    # filter out bad trees
    error <- try (expr = {
      plot (tree, plot = FALSE)
    }, silent = TRUE)
    # make sure tree has branch lengths and is rooted
    if (!is.null (tree$edge.length) & is.rooted (tree) &
          class (error) == 'list') {
      trees <- c (trees, list (tree))
      counter <- counter + 1
    }
  } else {
    cat ('\n........ no tree could be retrieved')
  }
}
# calculate RTT stat
pds <- ntips <- rtt.stats <- rep (NA, length (trees))
for (i in 1:length (trees)) {
  tree <- trees[[i]]
  print (i)
  ntips[i] <- length (tree$tip.label)
  pds[i] <- sum (tree$edge.length)
  rtt.stats[i] <- calcRTTStat (tree)
  plot (tree, main = paste0 ('Tree[', i, '] | [',
                             round (rtt.stats[i], 4), ']'),
        show.tip.label = FALSE)
}
ref.rtt.stat <- max (rtt.stats, na.rm = TRUE)
hist (rtt.stats)
plot (log (pds/ntips) ~ rtt.stats)  # there should be no relationship
plot (log (ntips) ~ rtt.stats)
# MAX RTT STAT FOR PGLT TREES: 0.5