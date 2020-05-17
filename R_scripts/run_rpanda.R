### run RPANDA (by Hélène Morlon) 

# install.packages('RPANDA',repos='http://cran.uni-muenster.de/')

library(ape)
library(RPANDA)

setwd("~/Dropbox/sDiv_working_group/sELDIG_extinction_estimate")

# source("R_scripts/read_trees.R")

# Function to test multiple models (fit.multi.rpanda)
source("R_scripts/multiRPANDA.R")

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
## The same tree as for the BAMM analysis
tree<-read.tree("bamm_build/good_trees/tree_396.tre")

# parameters based on the BAMM output
trypar<-list(c(0.02,0), 
			 c(0.02,-0.001,0), 
			 c(0.02,0.001,0.01),
			 c(0.02,-0.001,0.01,0.001))

# test models
tryrun<-fit.multi.rpanda(tree, trypar)

# compare models 
for(i in 1:length(tryrun)){
	print(tryrun[[i]]$aicc)
}

# plot rates
tot_time<-max(node.age(tree)$ages)
t<-seq(0, tot_time, length.out = 100)

pdf("results/rpanda_trial_plot.pdf", width=9, height=9)
par(mfrow=c(2,2))
for(i in 1:length(tryrun)){
	mod<-tryrun[[i]]
	y1<-mod$f.lamb(t)
	y2<-mod$f.mu(t)
	y3<-y1-y2

	plot(-t, y1, type = "l", 
		ylim=range(c(y1, y2, y3)),
		xlab="Time", ylab="Speciation rate",
		main=paste("Model aicc: ", round(mod$aicc), sep=''))

	lines(-t, y2, col=2)
	lines(-t, y3, lty=2, col=4)
	if(i==1) legend("right", lty=c(1,1,2), col=c(1,2,4), bty='n',
					legend=c("speciation rate", 
								"extinction rate", 
								"net diversificationr rate"))
}
dev.off()



