
## takes a null model and prepre specific arguments to streamline the testing
## idx.exclude are indices of individuals that should be excluded (e.g. because of missing genotypes)

nullModelTestPrep <- function(nullmod, idx.exclude = NULL){
    
    if (is.null(idx.exclude)){
        Y <- nullmod$workingY
        W <- nullmod$model.matrix
        resid <- nullmod$resid.marginal
    } else{
        Y <- nullmod$workingY[-idx.exclude]
        W <- nullmod$model.matrix[-idx.exclude,]
        resid <- nullmod$resid.marginal[-idx.exclude]
    }
    
    
    n <- length(Y)
    
    if (nullmod$family$mixedmodel){  ## n by n cholSigmaInv
        if (is.null(idx.exclude)){
            C <- nullmod$cholSigmaInv
        } else{
            C <- subsetCholSigmaInv(nullmod$cholSigmaInv, idx.exclude)
        }
        CW <- crossprod(C, W)	
        Mt <- C - tcrossprod(tcrossprod(C, tcrossprod(chol2inv(chol(crossprod(CW))), CW)), CW)
        resid <- as.vector(Mt %*% crossprod(Mt, Y))
    }
 
	if (!nullmod$family$mixedmodel & (nullmod$family$family != "gaussian")){
		sigma <- sqrt(nullmod$varComp)
		C <- diag(sigma)
		CW <- W * sigma
		Mt <- C - tcrossprod(tcrossprod(C, tcrossprod(chol2inv(chol(crossprod(CW))),CW)), CW)
	}

    
    if (!nullmod$family$mixedmodel & (nullmod$family$family == "gaussian")){  ## a vector or scalar cholSigmaInv
        
		if (nullmod$hetResid)	{  ## cholSigmaInv is a vector
            if (is.null(idx.exclude)){
            	C <- nullmod$cholSigmaInv
            } else{
            	C <- nullmod$cholSigmaInv[-idx.exclude]
            	}
         }   
        
        if (!nullmod$hetResid) { ## family is "gaussian", cholSigmaInv is a scalar.
            C <- nullmod$cholSigmaInv
        }	  		
        CW <- W * C      ## this is equal to crossprod(diag(C), W) when C is a vector
        # Mt <- diag(C) - tcrossprod(t(tcrossprod(chol2inv(chol(crossprod(CW))), CW))*C, CW)
        Mt <- -tcrossprod(t(tcrossprod(chol2inv(chol(crossprod(CW))), CW))*C, CW)
        diag(Mt) <- diag(Mt) + C
    }	
	
    # phenotype adjusted for the covariates/correlation structure
    Ytilde <- crossprod(Mt, Y)
    sY2 <- sum(Ytilde^2)

    return(list(Mt = Mt, Ytilde = Ytilde, sY2 = sY2, k = ncol(W), resid = resid))
}



nullModelBRprep <- function(nullmod){
    if (nullmod$family$mixedmodel) stop("BinomiRare should be used for IID observations.")
    if (nullmod$family$family != "binomial") stop("BinomiRare should be used for disease (binomial) outcomes.")
    
    probs <- nullmod$fitted.values
    
    return(list(D =nullmod$outcome, probs = probs))
    
}


# this is a fancy way of getting the inverse of the subset without having to get the original matrix
# cholesky decomposition of sigma inverse (inverse phenotype covariance matrix)
subsetCholSigmaInv <- function(cholSigmaInv, chol.idx) {
    if(length(chol.idx) > 0){
        # subset cholSigmaInv
        SigmaInv <- tcrossprod(cholSigmaInv)
        for(i in sort(chol.idx, decreasing=TRUE)){
            SigmaInv <- SigmaInv[-i,-i] - tcrossprod(SigmaInv[-i,i])/SigmaInv[i,i]
        }
        cholSigmaInv <- t(chol(SigmaInv))
    }
    
    cholSigmaInv
}
