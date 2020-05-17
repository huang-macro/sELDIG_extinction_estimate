## Function to test multiple models in RPANDA
## (based on the diversification tutorial here: 
## http://lukejharmon.github.io/ilhabela/) 

# different scenarios of b and d rates
lambda.cst <- function(x,y){y}
lambda.var <- function(x,y){y[1]*exp(y[2]*x)}
mu.cst <- function(x,y){y}
mu.var <- function(x,y){y[1]*exp(y[2]*x)}

# fit multiple models
fit.multi.rpanda <- function(tree,par){
    require(RPANDA)

    # constant b & d rates
    bcstdcst <- fit_bd(tree, max(branching.times(tree)), 
				    	f.lamb=lambda.cst, f.mu=mu.cst, 
				    	lamb_par=par[[1]][1], mu_par=par[[1]][2],
				    	cst.lamb=TRUE, cst.mu=TRUE, 
				    	cond="crown", f=87/89, dt=1e-3)

    # varying b but constant d
    bvardcst <- fit_bd(tree, max(branching.times(tree)), 
				    	f.lamb=lambda.var, f.mu=mu.cst, 
				    	lamb_par=par[[2]][c(1,2)], mu_par=par[[2]][3],
				    	expo.lamb=TRUE, cst.mu=TRUE, 
				    	cond="crown", f=87/89, dt=1e-3)

    # constant b but varying d
    bcstdvar <- fit_bd(tree, max(branching.times(tree)), 
    					f.lamb=lambda.cst, f.mu=mu.var, 
    					lamb_par=par[[3]][1], mu_par=par[[3]][c(2,3)],
    					cst.lamb=TRUE, expo.mu=TRUE,
    					cond="crown", f=87/89, dt=1e-3)

    # varying both b & d
    bvardvar <- fit_bd(tree, max(branching.times(tree)), 
    					f.lamb=lambda.var, f.mu=mu.var, 
    					lamb_par=par[[4]][c(1,2)], mu_par=par[[4]][c(3,4)],
    					expo.lamb=TRUE, expo.mu=TRUE, 
    					cond="crown", f=87/89, dt=1e-3)
    
    return(list("bcstdcst"=bcstdcst,"bvardcst"=bvardcst,"bcstdvar"=bcstdvar,"bvardvar"=bvardvar))
}

