### Set-up for BAMM (Dan Rabosky, http://bamm-project.org)

# install.packages('BAMMtools',repos='http://cran.uni-muenster.de/')

library(ape)
library(BAMMtools)

setwd("~/Dropbox/sDiv_working_group/sELDIG_extinction_estimate")

# source("R_scripts/read_trees.R")

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
## BAMM on mamy trees
file_list<-list.files("bamm_build/good_trees")

# A template of the control file needed for each run of BAMM
temp<-readLines("bamm_build/ControlFile_temp.txt")

for(i in 1:length(file_list)){

	# point to the tree file in the control file
	i.control<-gsub(pattern="tree_396.tre", replace=file_list[i], x=temp)

	# read a tree file
	i.file<-paste("bamm_build/good_trees/", file_list[i], sep='')
	i.tree<-read.tree(i.file)

	# Find priors for the general prior block for BAMM analysis
	i.priors<-setBAMMpriors(i.tree, outfile=NULL)
	l<-paste("lambdaInitPrior = ", i.priors["lambdaInitPrior"], sep='')
	ls<-paste("lambdaShiftPrior = ", i.priors["lambdaShiftPrior"], sep='')
	m<-paste("muInitPrior = ", i.priors["muInitPrior"], sep='')

	# Fill the priors into the control file
	i.control<-gsub(pattern="lambdaInitPrior = 18.4", replace=l, x=i.control)
	i.control<-gsub(pattern="lambdaShiftPrior = 0.003", replace=ls, x=i.control)
	i.control<-gsub(pattern="muInitPrior = 18.4", replace=m, x=i.control)

	# change prefix for output files
	tag<-strsplit(file_list[i], split="\\.")[[1]][1]
	i.control<-gsub(pattern="outName = trial", 
					replace=paste("outName = ", tag, sep=''), 
					x=i.control)

	# Save the file
	i.con<-paste("bamm_build/control_files/ControlFile_", tag, ".txt", sep='')
	writeLines(i.control, con=i.con)
}

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
## BAMM trial on one tree
tree<-read.tree("bamm_build/good_trees/tree_396.tre")

# General prior block for BAMM analysis
setBAMMpriors(tree, outfile=NULL)
