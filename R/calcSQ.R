

.computeSigmaQuantities <- function(varComp, covMatList, group.idx = NULL,  vmu = NULL, gmuinv = NULL){
	m <- length(covMatList)
	
	if (is.null(vmu)){ ## this means the family is "gaussian"
		if (is.null(group.idx)){
			group.idx <- list( 1:nrow(covMatList[[1]]))
			names(group.idx)[1]  <- names(varComp)[length(varComp)]
			g <- 1
		} else{
			g <- length(group.idx)
		}
		
	} else {g <- 0}
	
	Sigma <- Vre <- Reduce("+", mapply("*", covMatList, varComp[1:m], SIMPLIFY=FALSE))
    if(g > 0){
        diagV <- rep(0,nrow(covMatList[[1]]))
        for(i in 1:g){
            diagV[group.idx[[i]]] <- varComp[m+i]
        }
        diag(Sigma) <- diag(Sigma) + diagV
    }
    
    ### if non-gaussian family:
    if (!is.null(vmu)){
		Sigma <- Sigma + diag(as.vector(vmu)/as.vector(gmuinv)^2)
    }
    
    	# cholesky decomposition
	cholSigma <- chol(Sigma)
	# inverse
	Sigma.inv <- chol2inv(cholSigma)
	        


	return(list(Sigma = Sigma, cholSigma = cholSigma, Sigma.inv = Sigma.inv, Vre = Vre))
   
	
}

