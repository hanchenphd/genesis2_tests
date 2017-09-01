
### add heterogeneous variants - covMatList == NULL
# n - number of people 
# g - the number of groups (1 to more). if g == 1 group.idx is ignored. 
#' @importFrom stats var
.runWLSgaussian <- function (Y, X, group.idx, start, AIREML.tol,
                             maxIter,  verbose){
    n <- length(Y)
    k <- ncol(X)
    g <- length(group.idx)
    if (g <= 1) stop("group.idx must have length > 1")
    ## initializing parameters
    sigma2.p <- var(Y)
    AIREML.tol <- AIREML.tol * sigma2.p
    val <- 2 * AIREML.tol
    if (is.null(start)) {
        sigma2.k <- rep(sigma2.p,  g)
    } else {
        sigma2.k <- as.vector(start)
    }
    reps <- 0
    
    ## starting iteration
    repeat ({
        reps <- reps + 1
        
        if (reps > maxIter) {
            converged <- FALSE
            warning("Maximum number of iterations reached without convergence!")
            (break)()
        }

        
        ### Sigma is the matrix of variances. All non-diagonals are zero. 
        ## initialize Sigma (we don't need Sigma = only its diagonal!)
        diagSigma <- rep(0, n)
        
        ## set the values of the diagonal to be the group-specific variances. 
        for (i in 1:g) {
            diagSigma[group.idx[[i]]] <- sigma2.k[i]
        }
        
        ## just the diagonal - squared root of diagonals, and inverse of diagonals of Sigma.
        cholSigma.diag <- sqrt(diagSigma)
        Sigma.inv.diag <- 1/diagSigma
    
        lq <- .calcLikelihoodQuantities(Y, X, n, k, diag(Sigma.inv.diag), cholSigma.diag)

        
        
        if (verbose)  print(c(sigma2.k, lq$logLikR, lq$RSS))
        

        ## Updating variances and calculating their covariance matrix
        if (reps > 1) {
            
            score.AI <- .calcAIhetvars(lq$P, lq$PY, group.idx)
            score    <- score.AI$score
            AI       <- score.AI$AI
            AIinvScore <- solve(AI, score)

            sigma2.kplus1 <- sigma2.k + AIinvScore
            
            tau <- 1
            while (!all(sigma2.kplus1 >= 0)) {
                tau <- 0.5 * tau
                sigma2.kplus1 <- sigma2.k + tau * AIinvScore

            }
            val <- sqrt(sum((sigma2.kplus1 - sigma2.k)^2))
            sigma2.k <- sigma2.kplus1

            if (val < AIREML.tol) {
            	converged <- TRUE
                (break)()
            }
        }
        else { # first rep
            sigma2.kplus1 <- rep(NA, g)

            for (i in 1:g) {
             	sigma2.kplus1[i] <- (1/n) * (sigma2.k[i]^2 * crossprod(lq$PY[group.idx[[i]]]) + 
                                                         n *sigma2.k[ i] - sigma2.k[i]^2 * sum(diag(lq$P)[group.idx[[i]]]))
            }

            sigma2.k <- sigma2.kplus1
        }
    })
    
    ## after convergence, updated sigma again
	for (i in 1:g) {
		diagSigma[group.idx[[i]]] <- sigma2.k[i]
	}
        
	## just the diagonal - squared root of diagonals, and inverse of diagonals of Sigma.
	cholSigma.diag <- sqrt(diagSigma)
	Sigma.inv.diag <- 1/diagSigma
        
	lq <- .calcLikelihoodQuantities(Y, X, n, k, diag(Sigma.inv.diag), diag(cholSigma.diag))
	score.AI <- .calcAIhetvars(lq$P, lq$PY, group.idx)
	AI       <- score.AI$AI


    eta <- lq$fits
    return(list(varComp = sigma2.k, AI = AI, converged = converged,
                Sigma.inv.diag = Sigma.inv.diag, beta = lq$beta, 
                residM = lq$residM, fits = lq$fits, eta = eta, logLikR = lq$logLikR,
                logLik = lq$logLik, RSS = lq$RSS))
}



