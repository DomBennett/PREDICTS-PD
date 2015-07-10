# 04/11/2014
# Team PREDICTS-PD
# Community tools

# LIBS
library (MoreTreeTools)
library (reshape)

# FUNCTIONS
getCommunityMatrix <- function (study.data) {
  # tidyr has a function called spread() that would do this in a line -- but it's not out yet!
  # https://rpubs.com/m_dev/tidyr-intro-and-demos
  molten.data <- melt.data.frame (study.data,  measure.vars = 'Measurement', na.rm = TRUE)
  pull <- names (molten.data) %in% c ('Site_number', 'Parsed_name', 'variable', 'value')
  molten.data <- molten.data[ ,pull]
  cmatrix <- cast (molten.data, Site_number ~ Parsed_name, mean)
  # replace any missing values with 0
  cmatrix[is.na (cmatrix)] <- 0
  rownames (cmatrix) <- cmatrix[ ,1]
  cmatrix[ ,-1]
}

multiCommPhyMets <- function(trees, cmatrix, metric,
                             normalise=c ('age', 'total', 'none')) {
  # Wrapper for commPD that works with multiphylos
  .calc <- function(i) {
    tree <- trees[[i]]
    if (normalise == 'age') {
      # normlaise RTT
      tree$edge.length <- tree$edge.length / getSize (tree, 'rtt')
    } else if (normalise == 'total') {
      # normalise edge lengths to sum to 1
      tree$edge.length <- tree$edge.length/sum (tree$edge.length)
    } else if (normalise == 'none') {
      # do nothing
    } else {
      stop ('Unknown normalise arg')
    }
    # drop names not in community matrix
    cmatrix <- cmatrix[ ,colnames (cmatrix) %in% tree$tip.label]
    t (calcComPhyMets (cmatrix=cmatrix, tree=tree, metrics=c (metric))[ ,2])
  }
  normalise <- match.arg (normalise)
  mdply (.data=data.frame(i=1:length (trees)), .fun=.calc)[ ,-1]
}