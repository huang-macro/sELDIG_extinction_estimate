### Read trees from Oskar

# install.packages('geiger',repos='http://cran.uni-muenster.de/')

library(ape)
library(geiger)

setwd("~/Dropbox/sDiv_working_group/sELDIG_extinction_estimate")

# trees in a list
trees<-readRDS('data/Hagen_phylogenies.rds')

# corresponding parameters in a dataframe
# paras<-readRDS('data/Hagen_parameters.rds')

# Not all items in the tree list are trees
is.tree<-c()
for (i in 1:length(trees)){
	if (class(trees[[i]]) == "phylo") 
}

phylos<-trees[is.tree]

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
## BAMM only takes ultrametric, fully bifurcating trees with all 
## branch lengths >0
## Here, I only select trees with >= 50 extant species

for (i in 1:length(phylos)){
	p<-drop.extinct(phylos[[i]])
	if (length(p$tip.label)>50 & min(p$edge.length)>0){
		# p$edge.length[which(p$edge.length==0)]<-0.01
		write.tree(p, paste("bamm_build/good_trees/tree_", i, ".tre", sep=''))
	}
}
