### run RPANDA (by Hélène Morlon) 

# install.packages('RPANDA',repos='http://cran.uni-muenster.de/')

library(ape)
library(RPANDA)
library(BAMMtools)

setwd("~/Dropbox/sDiv_working_group/sELDIG_extinction_estimate")

# source("R_scripts/read_trees.R")

# Function to test multiple models (fit.multi.rpanda)
source("R_scripts/f_multiRPANDA.R")

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
## The same tree as for the BAMM trial
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

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
## Run RPANDA on all the trees that could be used in BAMM
file_list<-list.files("bamm_build/good_trees")

# Set up global initial models (output should not depend on these)
b<-0.05; bs<-0.001; d<-0.01; ds<-0.001
init.par<-list(c(b,d), 
		 c(b,bs,d), 
		 c(b,d,ds),
		 c(b,bs,d,ds))

# all the good trees

# Loop through all the trees to run RPANDA
out.list<-list()
for(i in 1:length(file_list)){
	# read a tree file
	i.file<-paste("bamm_build/good_trees/", file_list[i], sep='')
	i.tree<-read.tree(i.file)

	# run RPANDA 
	i.run<-fit.multi.rpanda(i.tree, init.par)
	out.list[[i]]<-i.run
}
saveRDS(out.list, "results/rpanda_output_list.rds")

###########################################################################
## compare models (constant rates in all cases)
best<-()
for(i in 1:length(out.list)){
	i.run<-out.list[[i]]
	i.aicc<-c()
	for(j in 1:length(i.run)) i.aicc[j]<-tryrun[[j]]$aicc
	best[i]<-which.min(i.aicc)
}

# compare rates from the constant-rate model
m1.b<-m1.d<-c()
for(i in 1:length(out.list)){
	i.run<-out.list[[i]][[1]]

	m1.b[i]<-i.run$lamb_par
	m1.d[i]<-i.run$mu_par
}

pdf("results/rpanda_constant_rates.pdf", height=7, width=5)
par(mfrow=c(2,1), mar=c(5,5,1,1), las=1)
hist(m1.b, breaks=15,
	xlab='Speciation rate', main='', col="firebrick3")
hist(m1.d, breaks=15,
	xlab='Extinction rate', main='', col="skyblue2")
dev.off()

