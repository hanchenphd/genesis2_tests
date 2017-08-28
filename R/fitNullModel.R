### A function that fits a null model, assuming complete data - will have to be wrapped using a different function that performs checks. 
## or checks will be added... for a start - complete valid data. nrow(X) == length(y), dimensions of covMatList similarly appropriate. We
## do not deal with IDs, just indices, as everything is assumed to match. 
## X is assumed to have an intercept. 
## non-gaussian families can only be binomial and poisson. 

#' Fit null model
#'
#' @param y outcome vector
#' @param X model.matrix
#' @param covMatList A list of matrices specifying the covariance structures of the random effects terms
#' @param group.idx list of indices for each group level
#' @param family A description of the error distribution to be used in the model. The default "gaussian" fits a linear mixed model; see \code{\link{family}} for further options.
#' @param AIREML.tol The convergence threshold for the Average Information REML (AIREML) procedure used to estimate the variance components of the random effects.
#' @param maxIter The maximum number of iterations allowed in the AIREML procedure
#' @param dropZeros Logical indicator of whether variance component terms that converge to 0 should be removed from the model
#' @param verbose Logical indicator of whether updates from the function should be printed to the console
#' @return The null model (TODO fill in details)
#'
#' @importFrom stats glm lm
#' @export
fitNullModel <- function(y, X, covMatList = NULL, group.idx = NULL, family = "gaussian", 
                         AIREML.tol = 1e-6, maxIter = 100, dropZeros = TRUE, verbose = TRUE){
    
    if(!is.null(covMatList)){
        if (class(covMatList) == "matrix"){
            covMatList <- list(A = covMatList)
        }
    } 

    if (is.null(colnames(X))){
        colnames(X) <- paste0("X", 1:ncol(X))
    }
    
    ## may be transferred to the wrapper function. 
    if(is.character(family)){
        family <- get(family)
    }
    if(is.function(family)){
        family <- family()
    }
    if(is.null(family$family)){
        stop("'family' not recognized")
    }

    if (!is.element(family$family, c("gaussian", "binomial", "poisson"))){
        stop("family must be one of gaussian, binomial, or poisson")
    }
    
    if (family$family == "gaussian"){
        if (is.null(covMatList) & is.null(group.idx)) {
            mod <- lm(y ~ -1 + X)  ## prepare output based on that. 
            out <- .nullModOutReg(y, X, mod, family)
        }
        if (is.null(covMatList) & !is.null(group.idx)){
            vc.mod <- .runWLSgaussian(y, X, group.idx = group.idx, start = NULL, 
                                      AIREML.tol = AIREML.tol,   maxIter = maxIter,  verbose = verbose)
            out <- .nullModOutWLS(y, X, vc.mod = vc.mod, family =  family, group.idx = group.idx)
        }
        if (!is.null(covMatList)){
            vc.mod <- .runAIREMLgaussian(y, X, start = NULL, covMatList, group.idx, AIREML.tol, dropZeros,  maxIter,  verbose)
            out <- .nullModOutMM(y, X, vc.mod, family = family, covMatList = covMatList, 
                                 group.idx = group.idx, dropZeros = dropZeros)
        }
    } 
    if (family$family != "gaussian"){ # separate condition instead of "else" for readability. 
        mod <- glm(y ~ X, family = family)
	
        
        if (!is.null(covMatList)){ ## iterate between computing workingY and estimating VCs. 
            eta <- mod$linear.predictors
            working.y <- .calcWorkingYnonGaussian(y, eta, family)
            newstart <- NULL
            Yreps <- 0
            repeat({
                Yreps <- Yreps + 1
                if(verbose) message("Computing Variance Component Estimates...")
                if(verbose) message(paste(paste("Sigma^2_",c(names(covMatList)),sep="", collapse="     "), "log-lik", "RSS", sep="     "))
                
                # estimate variance components
                out <- .runAIREMLother(Y=working.y$Y, X=X, start=newstart, covMatList=covMatList, 
                                       AIREML.tol=AIREML.tol, dropZeros=dropZeros, maxIter=maxIter, 
                                       verbose=verbose, vmu=working.y$vmu, gmuinv=working.y$gmuinv)
                
                if (out$allZero == TRUE) {
                    message("All variance components estimated as zero, using glm...")
                    break()
                }
                # update parameters
                if(verbose) message("Updating WorkingY Vector...")
                working.y <- .calcWorkingYnonGaussian(y, out$eta, family)
                
                # current variance component estimate
                newstart <- out$sigma2.k
                newstart[out$zeroFLAG] <- AIREML.tol
                
                # test for convergence
                stat <- sqrt(sum((out$eta - eta)^2))
                if(verbose) message(paste("Checking for Convergence...", stat, sep = "\t"))
                eta <- out$eta
                if(stat < AIREML.tol){ break() }
                if(Yreps == maxIter){
                    out$converged <- FALSE
                    warning("Maximum number of iterations for workingY reached without convergence!")
                    break()
                }
            })
	    ## check whether all variance components were estimated as zero:
            if (out$allZero == TRUE){
                out <- .nullModOutReg(y, X, mod, family)
            } else{
                out <- .nullModOutMM(y, X, vc.mod, family, covMatList = covMatList, vmu=working.y$vmu, gmuinv=working.y$gmuinv, dropZeros = dropZeros)
            }	
        } else{
            out <- .nullModOutReg(y, X, mod, family)
        }
        
    }	
    
    
    ## prepare output arguments. 
    ## First put arguments that are outputted using all regression models. 
    ## Then add arguments to the list according to the type of model. 
    
    ## preparing outputs in dedicated functions, add match.call() (Anyting else?)
    out$call <- match.call()
    
    return(out)	
    
    
}



## fit here is an object return by the glm() function.
## eta is X %*% beta
.calcWorkingYnonGaussian <- function(y, eta, family){
    mu <- family$linkinv(eta) # exp(eta)/(1 + exp(eta)) for binomial
    # weights
    vmu <- family$variance(mu) # mu(1-mu) for binomial
    # inverse of g'(mu)
    gmuinv <- family$mu.eta(eta) # = vmu for canonical link
    # working vector
    Y <- eta + (y - mu)/gmuinv

    return(list(Y=Y, vmu=vmu, gmuinv=gmuinv))
}
