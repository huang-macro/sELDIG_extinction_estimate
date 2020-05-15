### Read trees from Oskar

# install.packages('ape',repos='http://cran.uni-muenster.de/')

library(ape)

setwd("~/Dropbox/sDiv_working_group/sELDIG_extinction_estimate")

# trees in a list
trees<-readRDS('data/Hagen_phylogenies.rds')

# corresponding parameters in a dataframe
# paras<-readRDS('data/Hagen_parameters.rds')

i<-1
i.tree<-trees[[i]]