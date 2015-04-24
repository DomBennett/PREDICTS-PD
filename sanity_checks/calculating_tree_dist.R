# 11/02/2015
# Exploring different tree distance metrics

# LIBS
library (ape)

################################################################################
# PH85, TOPOLOGICAL
# 'The topological distance is defined as twice the number of internal branches
# defining different bipartitions of the tips'
################################################################################
# Simulated trees
tree1 <- stree (10, 'left')
tree2 <- stree (10, 'left')
# 1 node change leads to a distance of 2
holder <- tree1$tip.label[10]
tree1$tip.label[10] <- tree1$tip.label[8]
tree1$tip.label[8] <- holder
plot (tree1, type = 'unrooted')  # ignores root
(dist.topo (tree1, tree2))  # should be 2
# normalise by dividing by the total number of internal branches in both trees
((tree1$Nnode + tree2$Nnode) - 4)  # should be 14
# random trees again
tree1 <- stree (10, 'left')
tree2 <- stree (10, 'left')
# switch nodes around to make them maximally different
tree1$tip.label <- tree1$tip.label[c (10,2,8,4,6,5,7,3,9,1)]
plot (tree1, type = 'unrooted')
# normalised measure
dist.topo (tree1, tree2) / ((tree1$Nnode + tree2$Nnode) - 4)  # should be 1
# random test, n. internal branches shoudl never be greater than distance
bools <- rep (NA, 100)
for (i in 1:100) {
  tree1 <- rtree (100)
  tree2 <- rtree (100)
  bools[i] <- dist.topo (tree1, tree2) > ((tree1$Nnode + tree2$Nnode) - 4)
}
(any (bools))  # should be FALSE

################################################################################
# SCORE, BRANCH LENGTH
# 'square root of the sum of the squared differences of the (internal) branch
# lengths'
################################################################################
# first let's get a function that extracts internal branch lengths
getInternalEdgeLengths <- function (tree) {
  internal.branches <- which (!tree$edge[ ,2] %in% 1:length (tree$tip.label))
  tree$edge.length[internal.branches]
}
# now let's generate some trees with edge lengths -- make sure they're unrooted
# if they're rooted internal branches can be confusing, dist.topo assumes unrooted
tree1 <- unroot(stree (10, 'left'))
tree2 <- unroot(stree (10, 'left'))
tree1$edge.length <- rep (1, nrow (tree1$edge))
tree2$edge.length <- rep (1, nrow (tree2$edge))
# switching t8 and t10 is a small change and only changes 1 internal edge length
# we should expect there to be a sinlge edge difference. Since internal
# branches are all equal to 1, this should be the sqrt(2) -- the difference 
# squared and square rooted.
holder <- tree1$tip.label[10]
tree1$tip.label[10] <- tree1$tip.label[8]
tree1$tip.label[8] <- holder
plot (tree1, type = 'unrooted')
(dist.topo (tree1, tree2, method = 'score'))
# now let's trying rearranging all the tips to get a different tree
tree1 <- compute.brlen (unroot(stree (10, 'left')))
tree2 <- compute.brlen (unroot(stree (10, 'left')))
tree1$edge.length <- rep (1, nrow (tree1$edge))
tree2$edge.length <- rep (1, nrow (tree2$edge))
tree1$tip.label <- tree1$tip.label[c (10,2,8,4,6,5,7,3,9,1)]
dist.topo (tree1, tree2, method = 'score')  # now all the internal edges are different
# this is equal to sqrt(14) -- twice the number of different internal edges
# so the max difference possible can be expressed as:
#  sqrt (sum (internal.edge.lengths.tree1^2, internal.edge.lengths.tree1^2))
# Let's test this with some different lengths of branch
tree1 <- compute.brlen (unroot(stree (10, 'left')))
tree2 <- compute.brlen (unroot(stree (10, 'left')))
tree1$edge.length <- tree1$edge.length/sum(tree1$edge.length)
tree2$edge.length <- tree2$edge.length/sum(tree2$edge.length)
tree1$edge.length <- tree1$edge.length*10
tree1$tip.label <- tree1$tip.label[c (10,2,8,4,6,5,7,3,9,1)]
max.d <- sqrt (sum (getInternalEdgeLengths (tree1)^2, getInternalEdgeLengths (tree2)^2))
(dist.topo (tree1, tree2, method = 'score'))/max.d  # 1, maximum distance
# random test, n. internal branches should never be greater than distance
bools <- rep (NA, 100)
for (i in 1:100) {
  tree1 <- unroot(compute.brlen(rtree (100)))
  tree2 <- unroot(compute.brlen(rtree (100)))
  max.d <- sqrt (sum (getInternalEdgeLengths (tree1)^2, getInternalEdgeLengths (tree2)^2))
  obs.d <- dist.topo (tree1, tree2, method = 'score')
  bools[i] <- round (obs.d, 5) > round (max.d, 5)  # 5 decimal tol
}
(any (bools))  # should be FALSE

################################################################################
# CORRELATION OF DISTANCE MATRICES
################################################################################
# first test with the same tree
tree1 <- compute.brlen (stree (10, 'left'))
tree2 <- compute.brlen (stree (10, 'left'))
dist1 <- cophenetic.phylo (tree1)
dist2 <- cophenetic.phylo (tree2)
dist1 <- dist1[order (colnames(dist1)), order (colnames(dist1))]  # ensure tips are in the same order
dist2 <- dist2[order (colnames(dist2)), order (colnames(dist2))]
model <- cor.test (x = dist1, y = dist2)
model$estimate  # Pearson's r == 1
# Now switch tips 8 and 10 around
holder <- tree1$tip.label[10]
tree1$tip.label[10] <- tree1$tip.label[8]
tree1$tip.label[8] <- holder
dist1 <- cophenetic.phylo (tree1)
dist2 <- cophenetic.phylo (tree2)
dist1 <- dist1[order (colnames(dist1)), order (colnames(dist1))]
dist2 <- dist2[order (colnames(dist2)), order (colnames(dist2))]
model <- cor.test (x = dist1, y = dist2)
model$estimate  # Still very high
# Now muddle the tree
tree1 <- compute.brlen (stree (10, 'left'))
tree2 <- compute.brlen (stree (10, 'left'))
tree1$edge.length <- tree1$edge.length/sum(tree1$edge.length)
tree2$edge.length <- tree2$edge.length/sum(tree2$edge.length)
tree1$tip.label <- tree1$tip.label[c (10,2,8,4,6,5,7,3,9,1)]
dist1 <- cophenetic.phylo (tree1)
dist2 <- cophenetic.phylo (tree2)
dist1 <- dist1[order (colnames(dist1)), order (colnames(dist1))]
dist2 <- dist2[order (colnames(dist2)), order (colnames(dist2))]
model <- cor.test (x = dist1, y = dist2)
model$estimate  # that's as low as the correlation can get for this tree
# Can random trees for anything bigger get lower?
random.rs <- rep (NA, 100)
random.ns <- round (runif (100, 10, 200))
for (i in 1:100) {
  tree1 <- compute.brlen(rtree (random.ns[i]))
  tree2 <- compute.brlen(rtree (random.ns[i]))
  dist1 <- cophenetic.phylo (tree1)
  dist2 <- cophenetic.phylo (tree2)
  dist1 <- dist1[order (colnames(dist1)), order (colnames(dist1))]
  dist2 <- dist2[order (colnames(dist2)), order (colnames(dist2))]
  model <- cor.test (x = dist1, y = dist2)
  random.rs[i] <- model$estimate
}
hist(random.rs)
plot (random.rs ~ random.ns)  # yes, obviously bigger trees are much less likely to be similar