# 01/11/2014
# Team PREDICTS-PD
# Compare pG-lt and published trees.

# PARAMETERS
min.tree <- 5  # minimum number of tips in a tree for reference

# START
cat (paste0 ('\nStage 4 started at [', Sys.time (), ']'))

# LIBS
source (file.path ('tools', 'tree_tools.R'))

# DIRS
input.dirs <- c ('0_data', '3_parse')
output.dir <- '4_compare'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\nReading in trees ....')
# read in pub trees
pub.trees <- readInTrees (folder = file.path (input.dirs[1], 'pub_phylos'))
# read in pG-lt trees
pglt.trees <- readInTrees (folder = input.dirs[2])
cat ('\nDone.')

# PROCESS
cat ('\nCalculating shared nodes ....')
# calculate distances between pglt tree and the best pub tree
ph85 <- score <- dmat <- ntaxa <- etaxa <- ref.trees <-
  rep (NA, length (pglt.trees) * 100)
c <- 1
for (i in 1:length (pglt.trees)) {
  cat ('\n.... tree [', i, '/', length (pglt.trees), ']', sep = '')
  # get dist
  treedist <- pglt.trees[[i]]
  tip.labels <- getNames (treedist)
  # find best reference tree
  ref.tree <- findBestRef (tip.labels, pub.trees)
  ref.trees[c:(c + 99)] <- rep (names (ref.tree)[1], 100)
  ref.tree <- ref.tree[[1]]
  for (j in 1:length (treedist)) {
    tree <- treedist[[j]]
    shared.ntaxa <- sum (tree$tip.label %in% ref.tree$tip.label)
    if (shared.ntaxa >= min.tree) {
      res <- calcDist (tree, ref.tree)
      ph85[c] <- res[['PH85']]
      score[c] <- res[['score']]
      dmat[c] <- res[['dmat']]
      ntaxa[c] <- shared.ntaxa
      etaxa[c] <- getSize (tree) - shared.ntaxa
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
study <- rep (1:length (pglt.trees), each = 100)
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