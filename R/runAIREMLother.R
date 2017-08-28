
.runAIREMLother <- function(Y, X, start, covMatList, AIREML.tol, dropZeros, maxIter, verbose, vmu, gmuinv){
    
    m <- length(covMatList)
    n <- length(Y)
    k <- ncol(X)
    
    # initial values for variance components
    if(is.null(start)){
        sigma2.k <- rep(10*sqrt(AIREML.tol), m)
    }else{
        sigma2.k <- as.vector(start)
    }
    
    reps <- 0
    repeat({
        reps <- reps+1
        
        zeroFLAG <- sigma2.k < AIREML.tol # which elements have converged to "0"
        sigma2.k[zeroFLAG] <- 0 # set these to 0
        
        if (sum(zeroFLAG) == m)  return(list(allZero = TRUE))
        
        # variance matrix
        Vre <- Reduce("+", mapply("*", covMatList, sigma2.k[1:m], SIMPLIFY=FALSE))
        V <- Vre + diag(as.vector(vmu)/as.vector(gmuinv)^2)

        # cholesky decomposition
        cholSigma <- chol(V)
        # inverse
        Sigma.inv <- chol2inv(cholSigma)
        
        lq <- .calcLikelihoodQuantities(Y, X, n, k, Sigma.inv, diag(cholSigma))

        
        # print current estimates
        if(verbose) print(c(sigma2.k, lq$logLikR, lq$RSS))

        if(reps > 1){
            # Average Information and Scores
            covMats.score.AI <- .calcAIcovMats(lq$P, lq$PY, covMatList)
            AI <- covMats.score.AI$AI
            score <- covMats.score.AI$score
            
            if(dropZeros){  ## here need to exit if all terms were zero!!
                # remove Zero terms
                AI <- AI[!zeroFLAG,!zeroFLAG]
                score <- score[!zeroFLAG]
            }
            
            # update
            AIinvScore <- solve(AI, score)
            
            if(dropZeros){
                sigma2.kplus1[!zeroFLAG] <- sigma2.k[!zeroFLAG] + AIinvScore
                sigma2.kplus1[zeroFLAG] <- 0
            }else{
                sigma2.kplus1 <- sigma2.k + AIinvScore
                sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
            }
            
            # step-halving if step too far
            tau <- 1
            while(!all(sigma2.kplus1 >= 0)){
                tau <- 0.5*tau
                if(dropZeros){
                    sigma2.kplus1[!zeroFLAG] <- sigma2.k[!zeroFLAG] + tau*AIinvScore
                    sigma2.kplus1[zeroFLAG] <- 0
                }else{
                    sigma2.kplus1 <- sigma2.k + tau*AIinvScore
                    sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
                }
            }
            
            # test for convergence
            stat <- sqrt(sum((sigma2.kplus1 - sigma2.k)^2))
            sigma2.k <- sigma2.kplus1
            if(stat < AIREML.tol){
                converged <- TRUE
                break()
            }
            if(reps == maxIter){
                converged <- FALSE
                warning("Maximum number of iterations reached without convergence!")
                break()
            }
            
        }else{
            # EM step
            sigma2.kplus1 <- rep(NA,m)
            for(i in 1:m){
                PAPY <- crossprod(lq$P,crossprod(covMatList[[i]],lq$PY))
                sigma2.kplus1[i] <- (1/n)*((sigma2.k[i])^2*crossprod(Y,lq$PAPY) + (n*sigma2.k[i] - (sigma2.k[i])^2*sum(lq$P*covMatList[[i]])))
            }
            sigma2.k <- sigma2.kplus1
        }
        
    })
    
    # linear predictor
    ## what is VinvR ???
    eta <- lq$fits + crossprod(Vre, VinvR) # X\beta + Zb
    
    return(list(allZero = FALSE, sigma2.k = sigma2.k, AI = AI, converged = converged, zeroFLAG = zeroFLAG, beta = lq$beta, eta = lq$eta, logLikR=lq$logLikR, logLik=lq$logLik, RSS=lq$RSS))

}


