### Prep data for BAMM by Dan Rabosky (http://bamm-project.org)

install.packages('BAMMtools',repos='http://cran.uni-muenster.de/')

library(BAMMtools)

setwd("~/Dropbox/sDiv_working_group/sELDIG_extinction_estimate")

# source("R_scripts/read_trees.R")

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
## General prior block for BAMM analysis
setBAMMpriors(p, outfile=NULL)

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
## Look at BAMM output
tree<-read.tree("bamm_build/good_trees/tree_396.tre")
events<-read.csv("bamm_build/event_data.txt")

# Make a bammdata object
ed<-getEventData(tree, events, burnin=0.25)

# Rates at the tips
tip.rates<-getTipRates(ed)

###########################################################################
## Plot a colorful tree
pdf("bamm_build/trial_plots.pdf", width=6, height=5)

# the tree
trial<-plot.bammdata(ed, lwd=2, labels=F, cex=0.5)
addBAMMshifts(ed, cex=2) #no shift for tree 396
addBAMMlegend(trial, nTicks=4, side=4, las=1)
legend("bottomleft", "Tree 396", bty="n")

# speciation rate
plotRateThroughTime(ed, intervalCol="red", avgCol="red")

# tip rates
hist(tip.rates$lambda.avg, breaks=9, 
	col="firebrick2", xlab="average lambda", main="")
hist(tip.rates$mu.avg, breaks=9,
	col="skyblue2", xlab="average mu", main="")

dev.off()