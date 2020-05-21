### Read trees from Oskar
### Time unit: 0.17 Myr --> need conversion for standard methods

# install.packages('geiger',repos='http://cran.uni-muenster.de/')

library(ape)
library(geiger)

setwd("~/Dropbox/sDiv_working_group/sELDIG_extinction_estimate")

# function for handling zero branch length
source("R_scripts/f_zero_branch.R")

###########################################################################
# First set of trees in a list
trees<-readRDS('data/Hagen_phylogenies.rds')

# corresponding parameters in a dataframe
# paras<-readRDS('data/Hagen_parameters.rds')

# Not all items in the tree list are trees
is.tree<-c()
for (i in 1:length(trees)){
	if (class(trees[[i]]) == "phylo") is.tree<-c(is.tree, i)
}

phylos<-trees[is.tree]

###########################################################################
## BAMM only takes ultrametric, fully bifurcating trees with all 
## branch lengths >0
## Here, I only select trees with >= 50 extant species

for (i in is.tree){
	p<-drop.extinct(trees[[i]])
	if (length(p$tip.label)>50){
		# rescale the branch length to Myr
		p$edge.length<-p$edge.length/0.17

		# get rid of zero branch length
		while(min(p$edge.length)==0) p<-noZ(p)

		# save the good trees
		write.tree(p, paste("bamm_build/good_trees/tree_", i, ".tre", sep=''))
	}
}

# The original tree 396
p<-trees[[396]]
p

pdf("results/tree_396.pdf", width=6, height=11)
plot(p)
dev.off()

###########################################################################
###########################################################################
###########################################################################
# Another set of trees
nex_list<-list.files("data/1d")

trees_1d<-list()
for(i in 1:length(nex_list)){
	i.file<-paste("data/1d/", nex_list[i], sep='')
	trees_1d[[i]]<-read.nexus(i.file)

	i.name<-paste("1d", strsplit(nex_list[i], '\\.')[[1]][1], sep='_')
	names(trees_1d)[[i]]<-i.name
}

###########################################################################
# Find trees for BAMM
for (i in 1:length(trees_1d)){
	p<-drop.extinct(trees_1d[[i]])
	if (length(p$tip.label)>50){
		# rescale the branch length to Myr
		p$edge.length<-p$edge.length/0.17 

		# get rid of zero branch length
		while(min(p$edge.length)==0) p<-noZ(p)

		# save the good trees
		i.file<-paste("bamm_build/good_trees/tree_", 
						names(trees_1d)[i], ".tre", sep='')
		write.tree(p, i.file)
	}
}

