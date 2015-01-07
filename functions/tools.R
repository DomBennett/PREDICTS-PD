# 04/11/2014
# Team PREDICTS-PD
# Pipeline tools

# LIBRARIES
library (ape)
library (stringr)  # for str_replace_all
library (plyr)  # compact
library (reshape) # melt and cast
library (RJSONIO)  # fromJSON

# TREE WRANGLING
readInTrees <- function (folder) {
  # point at a folder, it will return a list of trees in that folder
  filenames <- list.files (path = folder, pattern = '\\.tre')
  studies <- sub ('\\.tre', '', filenames)
  all.trees <- list ()
  for (i in 1:length (studies)) {
    treedist <- read.tree (file.path (folder, filenames[i]))
    # drop underscores from tip names
    treedist <- dropUnderscore(treedist)
    all.trees <- c (all.trees, list (treedist))
    names (all.trees)[i] <- studies[i]
  }
  all.trees
}

dropUnderscore <- function (phylos) {
  # drop _ in names of trees in a multiphylo
  if (class (phylos) == 'multiPhylo') {
    .drop <- function (i) {
      phylo <- phylos[[i]]
      phylo$tip.label <- gsub ('_', ' ', phylo$tip.label)
      res <<- c (res, list (phylo))
    }
    res <- list ()
    m_ply (.data=data.frame (i=1:length (phylos)), .fun=.drop)
    class (res) <- 'multiPhylo'
    return (res)
  } else {
    phylos$tip.label <- gsub ('_', ' ', phylos$tip.label)
    return (phylos)
  }
}

findBestRef <- function (tip.labels, ref.trees) {
  # find the best reference tree based on given tip.labels
  ptips <- rep (NA, length (ref.trees))
  for (i in 1:length (ref.trees)) {
    ptips[i] <- sum (ref.trees[[i]]$tip.label %in% tip.labels)/length (tip.labels)
  }
  besti <- which (ptips == max (ptips))[1]
  ref.trees[[besti]]
}

getNames <- function(phylos) {
  # get tip names from a multiphylo
  .get <- function (i) {
    res <<- c (res, phylos[[i]]$tip.label)
  }
  res <- NULL
  m_ply (.data=data.frame (i=1:length (phylos)), .fun=.get)
  unique (res)
}

# PREDICTS DATA WRANGLING TOOLS
getCommunityMatrix <- function (study.data) {
  # tidyr has a function called spread() that would do this in a line -- but it's not out yet!
  # https://rpubs.com/m_dev/tidyr-intro-and-demos
  molten.data <- melt.data.frame (study.data,  measure.vars = 'Measurement', na.rm = TRUE)
  pull <- names (molten.data) %in% c ('Site_number', 'Parsed_name', 'variable', 'value')
  molten.data <- molten.data[ ,pull]
  cmatrix <- cast (molten.data, Site_number ~ Parsed_name, mean)
  rownames (cmatrix) <- cmatrix[ ,1]
  cmatrix[ ,-1]
}


# SEACRH FOR TAXON IDS ONLINE
.safeFromJSON <- function (url, max.trys = 10) {
  # Safe wrapper for fromJSON
  trys <- 0
  while (trys < max.trys) {
    json.obj <- try (fromJSON (url), silent = TRUE)
    if (class (json.obj) == 'try-error') {
      cat ('---- Connection failed: trying again ----\n')
      trys <- trys + 1
      Sys.sleep (10)
    } else {
      return (json.obj)
    }
  }
  stop ("Failed to connect, server may be down.")
}

taxaResolve <- function (names, batch = 100, datasource = 4){
  # Resolve taxonomic names via the Global Names Resolver.
  #  Names that cannot be resolved are returned as NA.
  #
  # Args:
  #  names: vector of names
  #  batch: the max batch number of names to be searched
  #
  # Return:
  #  dataframe containing GNR metadata for each name
  batchResolve <- function (batch.names) {
    #create query from names
    url <- "http://resolver.globalnames.org/name_resolvers.json?"
    data_source_ids <- paste0 ("&data_source_ids=", datasource)
    names2 <- paste0 ("names=", paste0 (str_replace_all (
      batch.names, " ", "+"), collapse = "|"))
    query <- paste (compact (list (url, names2,
                                   data_source_ids)), collapse = "")
    #search via API
    data <- .safeFromJSON (query)$data
    return (data)
  }
  data <- list()
  # Split names into batch sized chunks
  #  http://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
  x <- seq_along (names)
  di <- split (names, ceiling (x/batch))
  for (d in di) {
    temp.data <- batchResolve (d)
    data <- c (data, temp.data)
  }
  #transform results into output
  search.name <- name.string <- canonical.form <-
    lineage <- lineage.ids <- rank <- taxid <-
    match.type <- prescore <- score <- rep (NA, length (names))
  for (i in 1:length (names)){
    #print(i)
    if (length (data[[i]]) == 1){
      search.name[i] <- data[[i]][[1]]
    } else {
      search.name[i] <- data[[i]][[1]]
      name.string[i] <-
        data[[i]]$results[[1]]$name_string
      canonical.form[i] <-
        data[[i]]$results[[1]]$canonical_form
      lineage[i] <-
        data[[i]]$results[[1]]$classification_path
      lineage.ids[i] <-
        data[[i]]$results[[1]]$classification_path_ids
      rank[i] <-
        data[[i]]$results[[1]]$classification_path_ranks
      taxid[i] <-
        data[[i]]$results[[1]]$taxon_id
      match.type[i] <-
        data[[i]]$results[[1]]$match_type
      prescore[i] <-
        data[[i]]$results[[1]]$prescore
      score[i] <- data[[i]]$results[[1]]$score
    }
  }
  res <- data.frame (search.name = search.name,
                     name.string = name.string,
                     canonical.form = canonical.form,
                     lineage = lineage, lineage.ids =
                       lineage.ids, rank = rank,
                     taxid = taxid, match.type =
                       match.type, prescore = prescore,
                     score = score)
  return (res)
}

# CALC PD
extractEdges <- function(phylo, taxa, type = 1) {
  # Extract edges from a phylo object using 1 of 3 methods
  #
  # Args:
  #  phylo: phylogeny (ape class)
  #  taxa: vector of taxon names
  #  type:
  #     1 -- phylogeny consisting solely of the taxa, default
  #     2 -- edges from taxon tips to terminal node
  #     3 -- edges unique to taxa
  #
  # Return:
  #  vector of edges
  # TODO(01/07/2013): this may be more achievable with a vegan matrix
  if (!type %in% c(1,2,3)) {
    stop("Type must be an integer: 1, 2 or 3.")
  }
  if (!is.vector(taxa) | !is.character(taxa)) {
    stop("Invalid or no taxa given.")
  }
  if (length(taxa) == length (phylo$tip.label)){
    return(phylo$edge)
  }
  if (type == 1 & length (taxa) == 1){
    stop("length(taxa) == 1 :
         Cannot return a single edge for type 1.")
  }
  # start at the tips and step back into the phylogeny ...
  # ...add all connecting edges to a vector...
  # stop when all paths have met at the same node (type = 1)
  # or when all paths have reached the root node (type = 2)
  # or when all the nodes are unique (type = 3)
  edges <- match (match (taxa, phylo$tip.label), phylo$edge[,2])
  end.nodes <- phylo$edge[edges, 1]
  term.node <- length (phylo$tip.label) + 1
  if (all(end.nodes %in% term.node)) {
    return(edges)
  } else {
    if (type == 3){
      while (any (duplicated (end.nodes))){
        start.node <- end.nodes[duplicated(end.nodes)][1]
        if (sum (phylo$edge[,1] %in% start.node) == sum (end.nodes %in% start.node)){
          edge <- match (start.node, phylo$edge[,2])
          end.node <- phylo$edge[edge,1]
          edges <- c(edges, edge)
          end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
        } else {
          end.nodes <- end.nodes[end.nodes != start.node]
        }
      }
    } else {
      while (TRUE){
        end.nodes <- sort (end.nodes, TRUE)
        start.node <- end.nodes[1]
        edge <- match (start.node, phylo$edge[,2])
        end.node <- phylo$edge[edge,1]
        edges <- c(edges, edge)
        end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
        if (type == 2){
          if (sum (term.node == end.nodes) == length (end.nodes)){
            break
          }
        } else {
          if (sum (end.nodes[1] == end.nodes) == length (end.nodes)){
            break
          }
        }
      }
    }
    return (edges)
  }
  }

multiCommPD <- function(phylos, comm.data, type = 1, min.spp = 2,
                   taxon.names = colnames(comm.data)) {
  # Wrapper for commPD that works with multiphylos
  .calc <- function(i) {
    phylo <- phylos[[i]]
    # normalise edge lengths to sum to 1
    phylo$edge.length <- phylo$edge.length/sum (phylo$edge.length)
    # drop names not in comm data
    comm.data <- comm.data[ ,colnames (comm.data) %in% phylo$tip.label]
    res <- commPD (phylo, comm.data, type, min.spp)
    t (res)
  }
  mdply (.data=data.frame(i=1:length(phylos)), .fun=.calc)[ ,-1]
}

commPD <- function(phylo, comm.data, type = 1, min.spp = 2,
                   taxon.names = colnames(comm.data)) {
  # Calculate Faith's PD for a community given a community phylogeny and a community
  #  matrix. Depends on extractEdges().
  #
  # Args:
  #  phylo: community phylogeny (ape class)
  #  comm.data: matrix of community data, species as cols sites as rows
  #  type: which branches to include for calculating PD, either 1, 2 or 3,
  #   default 2 -- see extractEdges()
  #  min.spp: minimum number of species to include for calculating PD, default 2
  #  taxon.names: if taxon.names are not the column names of comm.data specify
  #   here, else ignore.
  #
  # Return:
  # vector of PDs
  if (type == 1 & min.spp < 2) {
    stop("Cannot compute type 2 PD for fewer than 2 species")
  }
  calcPD <- function (row) {
    if (sum(row) > 1){
      taxa <- taxon.names[as.logical(row)]
      edges <- extractEdges(phylo, taxa, type)
      return (sum(phylo$edge.length[edges]))
    } else {
      return (0)
    }
  }
  # add site names if none
  if (is.null(rownames(comm.data))) {
    rownames(comm.data) <- 1:nrow(comm.data)
  }
  # convert to incidence
  presences <- comm.data[comm.data > 0]
  comm.data[comm.data > 0] <- rep(1, length(presences))
  # drop sites with too few species
  row.sums <- rowSums(comm.data) >= min.spp
  site.names <- rownames(comm.data)[row.sums]
  comm.data <- comm.data[row.sums, ]
  # calculate PD by site
  pds <- as.matrix(apply(comm.data, 1, calcPD))
  rownames(pds) <- site.names
  return (pds)
}

# PLOT A COMMUNITY
plotComm <- function(comm.data, phylo, groups = rep(1, nrow(comm.data)),
                     no.margin = TRUE, ...){
  # Plot community data on community phylogeny to visualise
  #  differences in community structure. Use colours to distinguish
  #  site groups and alpha to distinguish abundances (works like rainbow())
  #
  # Args:
  #  comm.data: community data matrix (cols taxa, rows sites)
  #  phylo: community phylogeny
  #  groups: site groups
  #
  # Return:
  #  a matrix of community data
  # ultrametricize tree FIRST
  phylo<- compute.brlen(phylo, method="Grafen")
  # plot phylogeny, allow space for points
  edges <- extractEdges(phylo, phylo$tip.label[1], type = 2)
  phylo$edge.length <- phylo$edge.length/sum(phylo$edge.length[edges])
  # make all phylos the same length before plotting i.e. all branches from terminal
  # node to tip equal 1
  # for some weird reason the rules of plotting are dyanmic!
  if (nrow(comm.data) < 20) {
    variable.max <- 1 + nrow(comm.data)/20
    spacing.start <- 0.55
    spacing.i <- 0.05
  } else {
    variable.max <- nrow(comm.data)/10
    spacing.i <- 0.1 - 1/nrow(comm.data)
    spacing.start <- 0.5 + spacing.i
  }
  plot(phylo, no.margin = no.margin, show.tip.label = FALSE,
       x.lim = c(0, variable.max), ...)
  
  # generate alphas based on abundances
  n <- length(unique(groups))
  hs <- seq.int(0, 1 + max(1, n - 1)/n, length.out = n)%%1
  alphas <- comm.data/max(comm.data)
  
  # loop init
  ntips <- length(phylo$tip.label)
  spacing <- spacing.start
  group <- groups[1]
  j <- 1
  
  # loop through sites and plot points for present species
  for(i in 1:nrow(comm.data)){
    j <- ifelse(group == groups[i], j, j + 1)
    pull <- as.logical(comm.data[i,])
    taxa <- phylo$tip.label[pull]
    abunds <- alphas[i, pull]
    tiplabels(tip = match(taxa, phylo$tip.label),
              pch = 19, adj = spacing, col = hsv(rep(hs[j], ntips), 1, 1, abunds))
    spacing <- spacing + spacing.i
    group <- groups[i]
  }
}