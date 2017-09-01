
### fit a conditional model based on a previously fit model. Default does not re-estimate variance components, but it may. 
updateNullModCond <- function(nullmod, G, covMatList = NULL, group.idx = NULL, AIREML.tol = 1e-6, maxIter = 100, dropZeros = TRUE, verbose = TRUE){
    
    ## a few checks that may be transfered to wrapper function: 
    if (nullmod$hetResid & is.null(group.idx)) stop("group indices are required for updating the null model")

    if (updateVarComp){ ## this check may be pulled out for wrapper function.
    	if (is.null(covMatList)) stop("covMatList is needed for udpating variance components")
    	if (nullmod$hetResid & is.null(group.idx)) stop("group indices are required for updating variance components")
    }

    
    X = cbind(nullmod$model.matrix, G)
    
    if (!nullmod$family$mixedmodel){ ## if it is not a mixed model, re-fit the model. (This includes hetrogeneous residuals).
       
       return(fitNullModel(nullmod$outcome, X, covMatList = NULL, 
       						group.idx = group.idx, family = nullmod$family$family))  	
    }
    
    
    ### if the function reached this far, nullmod is a mixed model!
    ### Re-fit the model with the new design matrix and start point the varComp
	## from the provided nullmod object. 
    
	return(fitNullModel(nullmod$outcome, X, covMatList = covMatList, 
				group.idx = group.idx, family = nullmod$family$family, start = nullmod$varComp,
				AIREML.tol = AIREML.tol, maxIter= maxIter, dropZeros = dropZeros, verbose = verbose))     
	    
	     
    
}


