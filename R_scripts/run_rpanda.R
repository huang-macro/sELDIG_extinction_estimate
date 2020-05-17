### run RPANDA (by Hélène Morlon) 

# install.packages('RPANDA',repos='http://cran.uni-muenster.de/')

library(ape)
library(RPANDA)

setwd("~/Dropbox/sDiv_working_group/sELDIG_extinction_estimate")

# source("R_scripts/read_trees.R")

# Function to test multiple models (from Luke Harmon's tutorial)
source("R_scripts/multiRPANDA.R")

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
## The same tree as for the BAMM analysis
tree<-read.tree("bamm_build/good_trees/tree_396.tre")

trypar<-list(c(0.4,0),c(0.4,-0.05,0),c(0.4,0.1,0.05),c(0.4,-0.05,0.1,0.05))
