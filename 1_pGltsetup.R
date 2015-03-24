# 01/11/2014
# Team PREDICTS-PD
# Create setup files for pG-lt

# START
cat (paste0 ('\nStage 1 started at [', Sys.time (), ']'))

# LIBS
# source before dplyr! Else, a plyr and dplyr conflict.
library (MoreTreeTools)
library (dplyr)

# PARAMETERS
reference.batch <- FALSE  # setup for reference batch?
non.sp.tol <- 0.5  # the maximum proportion of names resolved above the species level
non.sci.tol <- 0.5  # the maximum proportion of non-scientific names

# DIRS
input.dir <- '0_data'
output.dir <- '1_pGltsetup'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\nReading in data ....')
# read in RDS
x <- readRDS (file.path (input.dir, 'diversity-2014-10-29-03-40-20.rds'))
# convert to dplyr format
x <- tbl_df (x)
cat ('\nDone.')

# MANIPULATE
cat ('\nManipulating data ....')
# make new study ID: stick Source_ID, Study_number and Class together
x$SSID <- paste0 (x$Source_ID, '_', x$Study_number, '_', x$Class)
# group
sources <- group_by (x, SSID)
# summarise
sources <- summarise (sources, N_names = n_distinct (Parsed_name),
                      Class_ex = first (Class), N_classes = n_distinct(Class),
                      Order_ex = first (Order), N_orders = n_distinct(Order),
                      Family_ex = first (Family), N_families = n_distinct(Family),
                      Genus_ex = first (Genus), N_genera = n_distinct(Genus),
                      P_sp_names = 1 - (sum ('' == Species)/n()),
                      P_sci_names = sum ('Scientific' == Resolution_entered)/n())
cat ('\nDone.')

# FILTER
cat ('\nFiltering data ....')
filtered <- filter (sources, P_sp_names >= non.sp.tol)
filtered <- filter (filtered, P_sci_names >= non.sci.tol)
filtered <- filter (filtered, N_names > 5)
if (reference.batch) {
  birds <- filter (filtered, Class_ex == 'Aves')
  bees <- filter (filtered, Order_ex == 'Hymenoptera')
  mammals <- filtered <- filter (filtered, Class_ex == 'Mammalia')
  filtered <- rbind (birds, mammals)
}
cat ('\nDone. [', nrow (filtered), '] suitable stuides identified.', sep = '')

# EXTRACT NAMES AND SAVE
cat ('\nExtracting names ....')
# use these to extract parent
taxa.counters <- c ('N_classes', 'N_orders', 'N_families', 'N_genera')
taxa.names <- c ('Class_ex', 'Order_ex', 'Family_ex', 'Genus_ex')
# counters for stats out
study.counter <- 0
names.counter <- 0
# convert dplyr to data.frame
metadata <- data.frame (filtered, stringsAsFactors = FALSE)
for (i in 1:nrow (metadata)) {
  # study ID
  study <- metadata$SSID[i]
  # progress
  cat ('\n.... study [', study, '], [', i, '/', nrow(metadata), ']', sep = '')
  # find lowest shared taxonomic group -- the parent
  pull <- which (metadata[i, taxa.counters] == 1)
  parent <- metadata[i, taxa.names][[pull[length (pull)]]]
  # get parentid -- this is necessary for names resolution in pglt
  resolved <- taxaResolve (names = as.character (parent))
  parentid <- resolved$taxid
  if (length (parentid) > 1) {
    cat ('\n.... unable to resolve parent taxonomic group')
    next
  }
  # write out parameters.csv and names.txt for pglt
  parameters <- data.frame (Parameter = 'parentid',
                            Value = parentid, Note = parent)
  temp.dir <- file.path (output.dir, study)
  if (!file.exists (temp.dir)) {
    dir.create (temp.dir)
  }
  names <- as.data.frame (filter (x, SSID == study))$Parsed_name
  names <- unique (as.character (names))
  write.table (names, file = file.path (temp.dir, 'names.txt'),
               col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table (parameters, file = file.path (temp.dir, 'parameters.csv'),
               col.names = TRUE, row.names = FALSE, quote = FALSE,
               sep = ',')
  study.counter <- study.counter + 1
  names.counter <- names.counter + length (names)
}
cat ('Done. Setup for [', names.counter, '] names and [', study.counter, '] studies.',
     sep='')

# FINISH
cat (paste0 ('\nStage 1 finished at [', Sys.time (), ']'))