

## a function that tests a variance component (e.g. genetic relatedness, for testing heritability) 
## using the likelihood ration test. 
testSingleVarianceComponent <- function(nullmod, varCompName, covMatList, group.idx = NULL, AIREML.tol = 1e-6, maxIter = 100, dropZeros = TRUE, verbose = TRUE){
	
	if (nullmod$hetResid & is.null(group.idx)) stop("missing group indices group.idx")
	
	ind.test <- match(varCompName, names(covMatList))
	
	
	if (length(covMatList) == 1){
		nullmod.noVarComp <- fitNullModel(nullmod$outcome, X = nullmod$model.matrix, covMatList = NULL, 
				group.idx = group.idx, family = nullmod$family$family, start = NULL,
				AIREML.tol = AIREML.tol, maxIter= maxIter, dropZeros = dropZeros, verbose = verbose)

	} else{
		nullmod.noVarComp <- fitNullModel(nullmod$outcome, X = nullmod$model.matrix, covMatList = covMatList[-ind.test], 
				group.idx = group.idx, family = nullmod$family$family, start = nullmod$varComp,
				AIREML.tol = AIREML.tol, maxIter= maxIter, dropZeros = dropZeros, verbose = verbose)	
		}
	

	
	
	test.stat <- -2*(nullmod.noVarComp$logLik - nullmod$logLik)			
	
	pval <- 0.5*pchisq(test.stat, df = length(nullmod$varComp), lower.tail = FALSE) + 0.5*pchisq(test.stat, df = length(nullmod$varComp) - 1, lower.tail = F)  
	
	return(pval)
}