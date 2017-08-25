
		
.calcLikelihoodQuantities <- function(Y, W, n, k, Sigma.inv, cholSigma.diag){
		
		### Calulate the weighted least squares estimate
        Sigma.inv_W <- crossprod(Sigma.inv, W)
        chol.Wt_Sigma.inv_W <- chol(crossprod(W, Sigma.inv_W))
        Wt_Sigma.inv_W.inv <- chol2inv(chol.Wt_Sigma.inv_W)
        beta <- crossprod(Wt_Sigma.inv_W.inv, crossprod(Sigma.inv_W, Y))
        
        ## calculate the mean of the outcomes
        fits <- tcrossprod(W, t(beta))
        
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
            RSS) - sum(log(diag(chol.Wt_Sigma.inv_W))))
            
         
        ## calculate projection matrix.
        P <- Sigma.inv - tcrossprod(tcrossprod(Sigma.inv_W, Wt_Sigma.inv_W.inv),
            Sigma.inv_W)
        PY <- crossprod(P, Y)
        
        return(list(P= P, PY = PY, RSS = RSS, logLik = logLik, logLikR = logLikR, Sigma.inv_R  = Sigma.inv_R , Sigma.inv_W = Sigma.inv_W, Wt_Sigma.inv_W.inv = Wt_Sigma.inv_W.inv, beta = beta, fits = fits, residM = residM))

}
