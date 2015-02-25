# 25/02/2015
# Team PREDICTS-PD
# Ensure deps are installed

# START
user.input <- readline('WARNING: this script will automaticaly install packages.
                       Are you sure you want to continue? (y/n) ')
if (user.input != 'y') {
  stop ('Execution halted')
}

# GET INSTALLED PACKAGES
packages <- installed.packages ()[ ,1]

# CHECK AND INSTALL
counter <- 0
if (!'devtools' %in% packages) {
  install.packages ('devtools')
  counter <- counter + 1
}
if (!'MoreTreeTools' %in% packages) {
  install_github ('https://github.com/DomBennett/MoreTreeTools.git')
  counter <- counter + 1
}
if (!'ape' %in% packages) {
  install.packages('ape')
  counter <- counter + 1
}
if (!'plyr' %in% packages) {
  install.packages('plyr')
  counter <- counter + 1
}
if (!'reshape' %in% packages) {
  install.packages('reshape')
  counter <- counter + 1
}
if (!'dplyr' %in% packages) {
  install.packages('dplyr')
  counter <- counter + 1
}

# END
cat ('\nComplete! Installed [', counter, '] packages.', sep = '')