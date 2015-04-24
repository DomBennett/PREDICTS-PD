# 01/11/2014
# Team PREDICTS-PD
# Compare pG-lt and published trees.

# PARAMETERS
min.tree <- 5  # minimum number of tips in a tree for reference
iterations <- 10  # how many samples to be randomly taken from the distributions?

# START
cat (paste0 ('\nStage 4 started at [', Sys.time (), ']'))

# LIBS
source (file.path ('tools', 'tree_tools.R'))

# DIRS
input.dirs <- c ('0_data', '4_parse')
output.dir <- '5_compare'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\nReading in trees ....')
# read in trees
trees <- readInTrees (folder=input.dirs[2], recursive=TRUE)
# read in pub trees again for reference
pub.trees <- readInTrees (folder=file.path (input.dirs[1], 'parsed_trees'))
cat ('\nDone.')

# FILTER
cat ('\nFiltering ....')
keep <- NULL
for (i in 1:length (trees)) {
  if (length (trees[[i]]) > 1) {
    keep <- append(keep, i)
  }
}
trees <- trees[keep]
cat ('\nDone')

# PROCESS
cat ('\nCalculating tree distance metrics ....')
# calculate distances between pglt and mapped trees
ph85 <- score <- dmat <- ntaxa <- etaxa <- ref.trees <-
  rep (NA, length (trees) * 100)
c <- 1
for (i in 1:length (trees)) {
  cat ('\n.... tree [', i, '/', length (trees), ']', sep = '')
  # get dists
  pglt.trees <- trees[[i]][['pglt']]
  mapped.trees <- trees[[i]][['mapped']]
  # find from whence the mapped trees came
  tip.labels <- getNames (pglt.trees)
  # find best reference tree
  ref.tree <- names (findBestRef (tip.labels, pub.trees))[1]
  # quick check that ref.tree is returned
  if (is.null (ref.tree)) {
    next
  }
  ref.trees[c:(c + 99)] <- rep (ref.tree, 100)
  for (j in 1:iterations) {
    pglt.tree <- pglt.trees[[sample (1:length (pglt.trees), 1)]]
    mapped.tree <- mapped.trees[[sample (1:length (mapped.trees), 1)]]
    shared.ntaxa <- sum (duplicated (c (pglt.tree$tip.label,
                                        mapped.tree$tip.label)))
    if (shared.ntaxa >= min.tree) {
      res <- calcDist (pglt.tree, mapped.tree)
      ph85[c] <- res[['PH85']]
      score[c] <- res[['score']]
      dmat[c] <- res[['dmat']]
      ntaxa[c] <- shared.ntaxa
      etaxa[c] <- getSize (pglt.tree) - shared.ntaxa
    }
    c <- c + 1
  }
}
p.ph85.mean <- mean (ph85, na.rm = TRUE)
p.ph85.sd <- sd (ph85, na.rm = TRUE)
p.score.mean <- mean (score, na.rm = TRUE)
p.score.sd <- sd (score, na.rm = TRUE)
p.dmat.mean <- mean (dmat, na.rm = TRUE)
p.dmat.sd <- sd (dmat, na.rm = TRUE)
p.etaxa.mean <- mean (etaxa / (ntaxa + etaxa), na.rm = TRUE)
p.etaxa.sd <- mean (etaxa / (ntaxa + etaxa), na.rm = TRUE)
cat ('\nDone....')
cat ('[', p.ph85.mean, '±', p.ph85.sd, '] PH85', sep='')
cat ('[', p.score.mean, '±', p.score.sd, '] score', sep='')
cat ('[', p.dmat.mean, '±', p.dmat.sd, '] dmat', sep='')
cat ('[', p.etaxa.mean, '±', p.etaxa.sd, '] extra taxa in pG-lt trees', sep='')

# OUTPUT
cat ('\nPlotting and summarising ....')
# quick ggplot boxplot function
ggBoxplot <- function (dist.metric, ylab) {
  pull <- !is.na (get (dist.metric))
  p <- ggplot (plot.data[pull, ], aes_string ('Comparison_Tree', dist.metric))
  p <- p + geom_boxplot (aes (fill = Comparison_Tree)) +
    ylab (ylab) + theme_bw () +
    theme (axis.text.x = element_blank(), axis.title.x = element_blank())
  print (p)
}
# quick ggplot scatter for showing relation to ntaxa
ggScatter <- function (dist.metric, ylab) {
  i <- which (colnames(plot.data) == dist.metric)
  colnames(plot.data)[i] <- 'dist.metric'
  temp <- ddply (plot.data, .variables=.(study, Comparison_Tree, ntaxa),
                 .fun=summarize, mean = mean (dist.metric, na.rm = TRUE),
                 se = sd (dist.metric, na.rm = TRUE) / sqrt (length (dist.metric)))
  limits <- aes (colour = Comparison_Tree, ymax = mean + se,
                 ymin = mean - se)
  pull <- !is.na (temp$mean)
  p <- ggplot (temp[pull,], aes (x = ntaxa, y = mean))
  p <- p + geom_point (aes (colour = Comparison_Tree)) +
    geom_errorbar(limits, width=0.2) +
    stat_smooth (method = 'lm') + xlab ('N. taxa') +
    ylab (ylab) +
    theme_bw ()
  print (p)
  bytaxa <- ddply (temp, .(Comparison_Tree), summarize,
                   mean = mean (mean, na.rm = TRUE))
  cat ('\n.... results for [', dist.metric, ']', sep = '')
  for (i in 1:nrow (bytaxa)) {
    cat ('\n........ ', as.character (bytaxa[[i, 'Comparison_Tree']]),
         ' = [', bytaxa[[i, 'mean']], ']',
         sep = '')
  }
}
# open filehandle
pdf (file.path (output.dir, 'comparisons.pdf'))
# make gg dataset and print boxplots
study <- rep (1:length (trees), each = 100)
plot.data <- data.frame (ph85, score, dmat, ntaxa, study,
                         Comparison_Tree = factor (ref.trees))
ggBoxplot ('ph85', 'PH85 (P. Internal Branches)')
ggBoxplot ('score', 'Score (P. Internal Branches w/ Length)')
ggBoxplot ('dmat', 'Distance of cophenetic matrices')
ggScatter ('ph85', 'PH85 (P. Internal Branches) (±SE)')
ggScatter ('score', 'Score (P. Internal Branches w/ Length) (±SE)')
ggScatter ('dmat', 'Distance of cophenetic matrices (±SE)')
dev.off ()
save (ph85, score, dmat, ref.trees, ntaxa,
      file = file.path (output.dir, 'results.RData'))
cat ('\nDone.')

# FINISH
cat (paste0 ('\nStage 4 finished at [', Sys.time (), ']'))