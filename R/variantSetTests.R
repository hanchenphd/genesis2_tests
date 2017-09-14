
testVariantSet <- function(nullprep, G, weights, test = c("SKAT", "Burden"), 
							burden.test = c("Score", "Wald"),  rho = 0, pval.method = "davies", 
							return.scores = FALSE, return.scores.cov = FALSE){

	if(!is.element(test, c("SKAT", "Burden"))) stop("Only SKAT and Burden tests are implemented for testing variant sets")
	
	if (test == "SKAT") return(.testVariantSetSKAT(nullprep, G, weights, rho, pval.method, 
											return.scores, return.scores.cov))
											
	if (test == "Burden") return(.testVariantSetBurden(nullprep, G, weights, burden.test))
	
}



.testVariantSetSKAT <- function(nullprep, G, weights,  rho = 0,  pval.method = "davies", 
								return.scores = FALSE, return.scores.cov = FALSE){
	 
	 U <- as.vector(crossprod(G, nullprep$resid))
     G <- crossprod(nullprep$Mt, G)
     if (length(rho) == 1) {
		testout <- .runSKATTest(scores = U, geno.adj = geno,
					weights = weights, rho = rho, pval.method = pval.method,
                  	optimal = FALSE)
		} else {
         testout <- .runSKATTest(scores = U, geno.adj = geno,
                  weights = weights, rho = rho, pval.method = pval.method,
                  optimal = TRUE)
            }
	return(testout)
}


## create the burden score, than calls the appropriate single variant test function. 
## can easily implement GxE interaction with the burden score... later!
.testVariantSetBurden <- function(nullprep, G, weights, burden.test = c("Score", "Wald")){
		
	burden <- colSums(t(G) * weights)
	if (burden.test == "Score") {
		out <- .testGenoSingleVarScore(nullprep$Mt, G = matrix(burden), nullprep$Ytilde, maf = 1 ) ## just an arbitrary value for maf, just so results isn't set to NA. 	
	}
	if (burden.test == "Wald"){
		out <- .testGenoSingleVarWald(nullprep$Mt, G = matrix(burden), nullprep$Ytilde, sY2 = nullprep$sY2, 
										n = length(nullprep$Ytile),k = nullprep$k ,maf = 1 ) ## just an arbitrary value for maf, just so results isn't set to NA. 	
	}
	
	return(out)
}