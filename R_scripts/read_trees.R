### Read trees from Oskar

# install.packages('ape',repos='http://cran.uni-muenster.de/')

library(ape)

setwd("~/Dropbox/sDiv_working_group/sELDIG_extinction_estimate")

# trees in a list
trees<-readRDS('data/Hagen_phylogenies.rds')

# corresponding parameters in a dataframe
# paras<-readRDS('data/Hagen_parameters.rds')

# Not all items in the tree list are trees
is.tree<-c()
for (i in 1:length(trees)){
	if (class(trees[[i]]) == "phylo") is.tree<-c(is.tree, i)
}

phylos<-trees[is.tree]