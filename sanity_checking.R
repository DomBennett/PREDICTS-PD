

# What's an acceptable non-uniform rtt dist?
library (ape)  # library
n <- 10  # number of tips
mean.rttdist <- 100
sd.rttdist <- 50
tree <- stree (n)  # simple tree to visualise branch attraction
tree$edge.length <- rnorm (n = n, mean = mean.rttdist,
                           sd = sd.rttdist)
# calculate residuals
residuals <- abs (tree$edge.length - mean (tree$edge.length))
residual.sd <- sd (residuals)
test.stat <- 10
success <- residual.sd < test.stat
plot (tree, main = paste0 ('Test: ', success))




