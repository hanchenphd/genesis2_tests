
### fit a conditional model based on a previously fit model. Default does not re-estimate variance components, but it may. 
updateNullModCond <- function(nullmod, G, updateVarComp = FALSE, covMatList = NULL, group.idx = NULL, AIREML.tol = 1e-6, maxIter = 100, dropZeros = TRUE, verbose = TRUE){
    
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
    ### if updateVarComp, re-fit the model with the new design matrix and start point the varComp
	## from the provided nullmod object. 
    
	if (updateVarComp){
		return(fitNullModel(nullmod$outcome, X, covMatList = covMatList, 
				group.idx = group.idx, family = nullmod$family$family, start = nullmod$varComp,
				AIREML.tol = AIREML.tol, maxIter= maxIter, dropZeros = dropZeros, verbose = verbose))
	}        
	    
	
	if (!updatedVarComp){        
		stop("model without updated variance component is not functional")
		
		#### to complete... 
	    ## if not updating variance components
	    ## all variance components are the same 
	    ## need to update all other parameters. 
	    m <- length(covMatList)
	    g <- length(group.idx) 
		
	    n <- nrow(X)
	    k <- ncol(X)
	    Y <- nullmod$outcome
	    
	    sq <- .computeSigmaQuantities(nullmod$varComp, covMatList = covMatList, family = nullmod$family$family,
	    									group.idx = group.idx, fitted.values = nullmod$fitted.values, 
	    									resid.marginal = nullmod$resid.marginal)
	    
	    
	    # update lq: 
	    lq <- .calcLikelihoodQuantities(Y = Y, X = X, n = n, k = k, 
	    								Sigma.inv = sq$Sigma.inv, diag(sq$cholSigma))
	     			
	    # update eta as well:
	 	VinvR <- crossprod(sq$Sigma.inv, nullmod$resid.marginal)
		eta <- fitted.values + crossprod(sq$Sigma, VinvR)   
	    
	    # update AI:
	    
	    
	    if (family$family == "gaussian"){
	    	vc.mod <- list(varComp = nullmod$varComp, AI = AI, converged = nullmod$converged, 
	    					zeroFLAG = nullmod$zeroFLAG, beta=lq$beta, eta=eta, 
	    					logLikR=lq$logLikR, logLik=lq$logLik, RSS=lq$RSS)	
	    }   
	    
	    if (family$family != "gaussian"){
	    	vc.mod <- list(allZero = FALSE, varComp = nullmod$varComp, AI = AI, 
	    					converged = nullmod$converged, zeroFLAG = nullmod$zeroFLAG, beta = lq$beta, 
	    					eta = eta, logLikR=lq$logLikR, logLik=lq$logLik, RSS=lq$RSS)
	    }   
		
	
		out <- .nullModOutMM(y = Y, workingY =  , vc.mod = vc.mod, family = nullmod$family, 
									covMatList = covMatList, group.idx = group.idx, vmu = , gmuinv = , 
									use.sparsity = FALSE, dropZeros = TRUE)
	   
	           
	    return(out)	 
    }  
     
    
}


