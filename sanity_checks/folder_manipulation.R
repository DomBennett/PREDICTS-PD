# NOT FOR SOURCING

# 09/03/2015 -- MOVING FOLDERS GIVEN STUDY NAME CHANGE
easy.folders <- list.files ('FINAL-BATCH')
full.folders <- list.files ('FULL-BATCH')
# loop through each easy folder, find match in full and move contents
for (i in 1:length (easy.folders)) {
  easy.folder <- easy.folders[i]
  pattern <- paste0 ('^', easy.folder)
  pull <- grepl (pattern, full.folders)
  if (any (pull)) {
    full.folder <- full.folders[pull]
    system (paste ("cp -r", file.path ('FINAL-BATCH', easy.folder, ''),
                   file.path ('FULL-BATCH', full.folder, '')))
  } else {
    stop (easy.folder)
  }
}

# 23/02/2015 -- EXTRACTING EASY BATCH, GETTING FINAL BATCH
# all folders that can be run
setup.folders <- list.files ('1_pGltsetup')
# current folders that have been run (easy batch)
run.folders <- list.files ('2_pGltrun')
# folders that shouldn't have been run
not.picked <- run.folders[!run.folders %in% setup.folders]
# easy batch
easy.batch <- run.folders[run.folders %in% setup.folders]
# folders that still need to be run
final.batch <- setup.folders[!setup.folders %in% run.folders]
# move folders (UNIX only)
for (folder in easy.batch) {
  system (paste ("cp -r", file.path('2_pGltrun', folder), 'EASY-BATCH/'))
}
for (folder in final.batch) {
  system (paste ("cp -r", file.path('1_pGltsetup', folder), 'FINAL-BATCH/'))
}

# 23/02/2015 -- WORKING OUT RUN SUCCESSES
completed <- 0
for (folder in run.folders) {
  if (file.exists (file.path ('2_pGltrun', folder, '4_phylogeny',
                              'distribution.tre'))) {
    completed <- completed + 1
  }
}
pcomplete <- round (completed * 100 / length (run.folders), 2)
cat ('[', pcomplete, '%] completed', sep = '')