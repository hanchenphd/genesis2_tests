

.calcLikelihoodQuantities <- function(Y, X, n, k, Sigma.inv, cholSigma.diag){
    
    ### Calulate the weighted least squares estimate
    Sigma.inv_X <- crossprod(Sigma.inv, X)
    chol.Xt_Sigma.inv_X <- chol(crossprod(X, Sigma.inv_X))
    Xt_Sigma.inv_X.inv <- chol2inv(chol.Xt_Sigma.inv_X)
    beta <- crossprod(Xt_Sigma.inv_X.inv, crossprod(Sigma.inv_X, Y))
    
    ## calculate the mean of the outcomes
    fits <- tcrossprod(X, t(beta))
    
    # obtain marginal residuals
    residM <- as.vector(Y - fits)
    Sigma.inv_R <- crossprod(Sigma.inv, residM)
    
    # calculate likelihood quantities
    Rt_Sigma.inv_R <- crossprod(residM, Sigma.inv_R)
    RSS <- as.numeric(Rt_Sigma.inv_R/(n - k)) 
    logLik <- as.numeric(-0.5 * n * log(2 * pi * RSS) - sum(log( cholSigma.diag)) -
                         0.5 * Rt_Sigma.inv_R/RSS)
    
    ## log likelihood- REML type, accounting for estimation of mean effects.    
    logLikR <- as.numeric(logLik + 0.5 * k * log(2 * pi *
                                                 RSS) - sum(log(diag(chol.Xt_Sigma.inv_X))))
    
    
    ## calculate projection matrix.
    P <- Sigma.inv - tcrossprod(tcrossprod(Sigma.inv_X, Xt_Sigma.inv_X.inv),
                                Sigma.inv_X)
    PY <- crossprod(P, Y)
    
    return(list(P= P, PY = PY, RSS = RSS, logLik = logLik, logLikR = logLikR, Sigma.inv_R  = Sigma.inv_R , Sigma.inv_X = Sigma.inv_X, Xt_Sigma.inv_X.inv = Xt_Sigma.inv_X.inv, beta = beta, fits = fits, residM = residM))

}
