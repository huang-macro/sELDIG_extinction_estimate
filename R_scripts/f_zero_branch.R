### Function to handle zero branch length

noZ<-function(tree, add.length=0.001){
	
	require(ape)

	# the zero-length branches
	z<-which(tree$edge.length==0)
	tree$edge.length[z]<-add.length

	# their daughter branches
	d<-which(tree$edge[, 1] %in% tree$edge[z, 2])
	tree$edge.length[d]<-tree$edge.length[d]-add.length

	tree
}