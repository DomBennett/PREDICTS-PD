# 23/02/2015
# Explore, move, and investigate pG-lt folders

# EXTRACTING EASY BATCH, GETTING FINAL BATCH
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

# WORKING OUT RUN SUCCESSES
completed <- 0
for (folder in run.folders) {
  if (file.exists (file.path ('2_pGltrun', folder, '4_phylogeny',
                              'distribution.tre'))) {
    completed <- completed + 1
  }
}
pcomplete <- round (completed * 100 / length (run.folders), 2)
cat ('[', pcomplete, '%] completed', sep = '')