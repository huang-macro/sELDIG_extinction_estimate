### post-BAMM analyses (Dan Rabosky, http://bamm-project.org)

# install.packages('BAMMtools',repos='http://cran.uni-muenster.de/')

library(BAMMtools)
library(dplyr)

setwd("~/Dropbox/sDiv_working_group/sELDIG_extinction_estimate")

# source("R_scripts/read_trees.R")

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
## Look at BAMM output from many trees

# get some info of the original trees
ori.info<-read.csv("data/original_tree_info.csv")

###########################################################################
# extract some bamm results
prob.nos<-rep(NA, nrow(ori.info))

for(i in 1:nrow(ori.info)){
	# get output information
	i.eventfile<-paste("bamm_build/event_data/", 
					   ori.info$tree.id[i], 
					   "_event_data.txt", sep='')
	i.events<-read.csv(i.eventfile)

	# get corresponding tree
	i.treefile<-paste("bamm_build/good_trees/", 
					  ori.info$tree.id[i],
					  ".tre", sep='')
	i.tree<-read.tree(i.treefile)

	# make a bammdata object
	i.ed<-getEventData(i.tree, i.events, burnin=0.25)

	# prob. of 0 shift (because the simulationd didn't have any)
	i.prob<-summary(i.ed)
	prob.nos[i]<-i.prob[i.prob$shifts==0, "prob"]
}

###########################################################################
hist(prob.nos)
hist(ori.info$age)

pdf("results/bamm_prob_noshift_extinction.pdf", width=9, height=3)
par(mfrow=c(1,3), las=1, mar=c(5,5,1,1))
plot(ori.info$no.ex, prob.nos, 
	xlab="Number of extinct species",
	ylab="Prob. no rate shift")
plot(ori.info$prop.ex, prob.nos,
	xlab="Proportion of extinct species",
	ylab="Prob. no rate shift")
plot(ori.info$age, prob.nos,
	xlab="Max. age",
	ylab="Prob. no rate shift")
dev.off()


###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
## Look at BAMM output from one tree 
tree<-read.tree("bamm_build/good_trees/tree_396.tre")
events<-read.csv("bamm_build/event_data/tree_396_event_data.txt")

# Make a bammdata object
ed<-getEventData(tree, events, burnin=0.25)

# Rates at the tips
tip.rates<-getTipRates(ed)

###########################################################################
## Plot a colorful tree
pdf("results/bamm_trial_plots.pdf", width=6, height=5)

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