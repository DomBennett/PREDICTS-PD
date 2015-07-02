# 03/07/2015
# D.J. Bennett
# Plotting the difference in RTTstats

# LIB
library (ape)

# FUNCTION
calcRTTStat <- function (tree) {
  # Return RTT stat -- the coefficient of variation of the
  #  normalised root to tip distances
  rtt.dists <- diag (vcv.phylo (tree))  # calc rtts
  norm.rtt.dists <- rtt.dists/max (rtt.dists)  # normalise
  sd (norm.rtt.dists) / mean (norm.rtt.dists)  # return CoV
}

# PROCESS
# get a good tree
tree <- compute.brlen (rtree (100))
deep.edges <- which (!tree$edge[,2] %in% 1:length (tree$tip.label))
is <- sample (deep.edges, 10)
tree$edge.length[is] <- tree$edge.length[is] + rnorm (length(is), sd=0.01)
plot (tree, show.tip.label=FALSE)
rtt.stat <- calcRTTStat(tree)
good.tree <- tree
# get a bad tree
tree <- rtree (100)
tip.edges <- which (tree$edge[,2] %in% 1:length (tree$tip.label))
is <- sample (tip.edges, 2)
tree$edge.length[is] <- tree$edge.length[is] + abs (rnorm (length(is), sd=50))
plot (tree, show.tip.label=FALSE)
rtt.stat <- calcRTTStat(tree)
bad.tree <- tree
# plot them together
par (mfrow = c(1,2))
plot (good.tree, show.tip.label=FALSE, main='RTT-stat < 0.1')
plot (bad.tree, show.tip.label=FALSE, main='RTT-stat > 1.0')
