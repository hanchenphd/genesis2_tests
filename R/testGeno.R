
## function that gets an n\times p matrix of p genotypes of n individuals, and a null model, and tests the genotypes associations with the outcomes. 
## Genetic data are always assumed complete. 
## Types of tests: 
## Single variant: Wald, score, BinomiRare, interaction. 
## Variant set: SKAT, burden, SKAT-O. Multiple types of p-values. Default: Davis with Koenen if does not converge. 


# E an environmntal variable for optional GxE interaction analysis. 
testGenoSingleVar <- function(nullprep, G, E = NULL, test = c("Wald", "Score"), GxE.return.cov = FALSE){
    test <- match.arg(test)
    
    n <- length(nullprep$Ytilde)
    p <- ncol(G)
    
    if (test == "Wald" & is.null(E)){
        res <- .testGenoSingleVarWald(nullprep$Mt, G, nullprep$Ytilde, nullprep$sY2, n, nullprep$k)
    }
    
    if (test == "Wald" & !is.null(E)){
        res <- .testGenoSingleVarWaldGxE(nullprep$Mt, G, E, nullprep$Ytilde, nullprep$sY2, n, nullprep$ k)
    }
    
    if (test == "Score"){
        res <- .testGenoSingleVarScore(nullprep$Mt, G, nullprep$Ytilde)
    }

    return(res)
}



.testGenoSingleVarScore <- function(Mt, G, Ytilde){
    Xtilde <- crossprod(Mt, G) # adjust genotypes for correlation structure and fixed effects
    XtX <- colSums(Xtilde^2) # vector of X^T P X (for each SNP) b/c (M^T M) = P
    score <- as.vector(crossprod(Xtilde, Ytilde)) # X^T P Y
    Stat <- score/sqrt(XtX)
    
    res <- cbind(Score = score, Score.SE = sqrt(XtX), Score.Stat = Stat, 
                 Score.pval = pchisq(Stat^2, df = 1, lower.tail = FALSE) )
    
    return(res)	
}



.testGenoSingleVarWald <- function(Mt, G, Ytilde, sY2, n, k){
    Xtilde <- crossprod(Mt, G) # adjust genotypes for correlation structure and fixed effects
    XtX <- colSums(Xtilde^2) # vector of X^T SigmaInv X (for each SNP)
    XtY <- as.vector(crossprod(Xtilde, Ytilde))
    beta <- XtY/XtX
    RSS <- as.numeric((sY2 - XtY * beta)/(n - k - 1))
    Vbeta <- RSS/XtX
    Stat <- beta/sqrt(Vbeta)
    res <- cbind(Est = beta, Est.SE = sqrt(Vbeta), Wald.Stat = Stat, 
                 Wald.pval = pchisq(Stat^2, df = 1, lower.tail = FALSE))
    return(res)	
}



.testGenoSingleVarWaldGxE <- function(Mt, G, E, Ytilde, sY2, n, k, GxE.return.cov.mat = FALSE){

    E <- as.matrix(E)
    p <- ncol(G)
    v <- ncol(E) + 1
    
    if (GxE.return.cov.mat) {
        res.Vbetas <- vector(mode = "list", length = p)
    } else {
        res.Vbetas <- NULL
    }
    
    intE <- cbind(1, E) # add intercept the "Environmental" variable E.
    
    var.names <- c("G", paste0("G", colnames(E), sep = ":"))
    
    res <- matrix(NA, nrow = p, ncol = length(var.names)*2 + 2,
                  dimnames = list(NULL, 
                                  c(paste0("Est.", var.names), paste0("SE.", var.names), "GxE.Stat", "Joint.Stat" ) ))
    
    if (ncol(E) == 1) res[,"cov.G.E"] <- NA
    
    for (g in 1:p) {
        Xtilde <- crossprod(Mt, G[, g] * intE)
        XtX <- crossprod(Xtilde)
        XtXinv <- tryCatch(chol2inv(chol(XtX)), error = function(e) {TRUE}) # this is inverse A matrix of sandwich
        # check that the error function above hasn't been called (which returns TRUE instead of the inverse matrix)
        if (is.logical(XtXinv)) next
        
        XtY <- crossprod(Xtilde, Ytilde)
        betas <- crossprod(XtXinv, XtY)
        res[g, grep("Est", colnames(res))] <- betas
        
        RSS <- as.numeric((sY2 - crossprod(XtY, betas))/(n - k - v))
        Vbetas <- XtXinv * RSS
        
        if (GxE.return.cov.mat) {
            res.Vbetas[[g]] <- Vbetas
        }
        
        res[g, grep(colnames(res), "SE")] <- sqrt(diag(Vbetas))
        
        res[g, "GxE.stat"] <- tryCatch(crossprod(betas[-1],
                                                 crossprod(chol2inv(chol(Vbetas[-1, -1])),
                                                           betas[-1])), 
                                       error = function(e) { NA })
        
        res[g, "Joint.Stat"] <- tryCatch(crossprod(betas,
                                                   crossprod(XtX, betas))/RSS, 
                                         error = function(e) { NA })
    }
    
    res[,"GxE.pval"] <- pchisq(res[,"GxE.Stat"], df = (v - 1), lower.tail = FALSE)
    res[,"Joint.pval"] <- pchisq(res[,"Joint.Stat"], df = v, lower.tail = FALSE)

    return(list(res = res, GxEcovMatList = res.Vbetas))
}



## G is an n by v matrix of 2 or more columns, all representing allels of the same (multi-allelic) variant. 
.testSingleVarMultAlleles <- function(Mt, G, Ytilde, sY2, n, k){
    v <- ncol(G)
    
    var.names <- colnames(G)
    
    res <- matrix(NA, nrow = 1, ncol = length(var.names)*2 + 2,
                  dimnames = list(NULL, 
                                  c(paste0("Est.", var.names), paste0("SE.", var.names), "Joint.Stat", "Joint.Pval" ) ))
    
    
    Xtilde <- crossprod(Mt, G)
    XtX <- crossprod(Xtilde)
    XtXinv <- tryCatch(chol2inv(chol(XtX)), error = function(e) {TRUE})
    
    if (is.logical(XtXinv)) return(list(res = res, allelesCovMat = NA))
    
    XtY <- crossprod(Xtilde, Ytilde)
    betas <- crossprod(XtXinv, XtY) ## effect estimates of the various alleles
    res[1, grep("Est", colnames(res))] <- betas
    
    RSS <- as.numeric((sY2 - crossprod(XtY, betas))/(n - k - v))
    Vbetas <- XtXinv * RSS
    
    res[1, grep(colnames(res), "SE")] <- sqrt(diag(Vbetas))
    
    res[1, "Joint.Stat"] <- tryCatch(crossprod(betas,
                                               crossprod(XtX, betas))/RSS, 
                                     error = function(e) { NA })
    
    res[,"Joint.pval"] <- pchisq(res[,"Joint.Stat"], df = v, lower.tail = FALSE)
    return(list(res = res, allelesCovMat = Vbetas))
}





