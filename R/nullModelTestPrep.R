
## takes a null model and prepre specific arguments to streamline the testing
## idx.exclude are indices of individuals that should be excluded (e.g. because of missing genotypes)

nullModelTestPrep <- function(nullmod, idx.exclude = NULL){
	
	if (is.null(idx.exclude)){
		Y <- nullmod$workingY
		W <- nullmod$model.matrix
	} else{
		Y <- nullmod$workingY[-idx.exclude]
		W <- nullmod$model.matrix[-idx.exclude,]
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
	}
	
	if (!nullmod$family$mixedmodel){  ## a vector or scalar cholSigmaInv
				
		if (nullmod$hetResid | nullmod$family$family != "gaussian")	{  ## cholSigmaInv is a vector
			C <- ifelse(is.null(idx.exclude), nullmod$cholSigmaInv , nullmod$cholSigmaInv[-idx.exclude])
		} else { ## not hetResid, family is "gaussian", cholSigmaInv is a scalar.
			C <- nullmod$cholSigmaInv
		}	  		
		CW <- W * C      ## this is equal to crossprod(diag(C), W) when C is a vector
		# Mt <- diag(C) - tcrossprod(t(tcrossprod(chol2inv(chol(crossprod(CW))), CW))*C, CW)
		Mt <- -tcrossprod(t(tcrossprod(chol2inv(chol(crossprod(CW))), CW))*C, CW)
		diag(Mt) <- diag(Mt) + C
	}	
	
	Ytilde <- crossprod(Mt, Y)
	sY2 <- sum(Ytilde^2)

	return(list(CW = CW, Mt = Mt, Ytilde = Ytilde, sY2 = sY2, k = ncol(W)))
}



nullModelBRprep <- function(nullmod){
	if (nullmod$family$mixedmodel) stop("BinomiRare should be used for IID observations.")
	if (nullmod$family$family != "binomial") stop("BinomiRare should be used for disease (binomial) outcomes.")
	
	probs <- nullmod$fitted.values
	
	return(probs)
	
}