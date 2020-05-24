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
age<-no.ex<-prop.ex<-tree.id<-c()

for (i in is.tree){
	p<-drop.extinct(trees[[i]])
	if (length(p$tip.label)>50){
		# rescale the branch length to Myr
		p$edge.length<-p$edge.length*0.17

		# get rid of zero branch length
		while(min(p$edge.length)==0) p<-noZ(p)

		# save the good trees
		write.tree(p, paste("bamm_build/good_trees/tree_", i, ".tre", sep=''))

		# save some info
		tree.id[i]<-paste("tree", i, sep='_')
		age[i]<-max(node.age(trees[[i]])$ages)*0.17
		no.ex[i]<-length(trees[[i]]$tip.label)-length(p$tip.label)
		prop.ex[i]<-no.ex[i]/length(trees[[i]]$tip.label)
	}
}

ori.info<-na.omit(data.frame(tree.id, age, no.ex, prop.ex))

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
age_1d<-no.ex_1d<-tree.id_1d<-prop.ex_1d<-c()

for (i in 1:length(trees_1d)){
	p<-drop.extinct(trees_1d[[i]])
	if (length(p$tip.label)>50){
		# rescale the branch length to Myr
		p$edge.length<-p$edge.length*0.17 

		# get rid of zero branch length
		while(min(p$edge.length)==0) p<-noZ(p)

		# save the good trees
		i.file<-paste("bamm_build/good_trees/tree_", 
						names(trees_1d)[i], ".tre", sep='')
		write.tree(p, i.file)

		# save some info
		tree.id_1d[i]<-paste("tree", names(trees_1d)[i], sep='_')
		age_1d[i]<-max(node.age(trees_1d[[i]])$ages)*0.17
		no.ex_1d[i]<-length(trees_1d[[i]]$tip.label)-length(p$tip.label)
		prop.ex_1d[i]<-no.ex_1d[i]/length(trees_1d[[i]]$tip.label)
	}
}

ori.info_1d<-na.omit(data.frame(tree.id=tree.id_1d, age=age_1d, 
								no.ex=no.ex_1d, prop.ex=prop.ex_1d))
ori.info<-rbind(ori.info, ori.info_1d)

write.csv(ori.info, "data/original_tree_info.csv", row.names=F)