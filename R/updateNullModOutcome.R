
### takes an existing linear model (family has to be "gaussian"!), rank-normalize residuals and scale, and re-fit null model. If group.idx are provided, optional scale by group.idx. 

## the names of items in the list group.idx have to match the names of the corresponding variance components!

updateNullModOutcome <- function(nullmod, covMatList = NULL, group.idx = NULL, rankNorm.option = c("by.group", "all"), rescale = TRUE, AIREML.tol = 1e-6, maxIter = 100, verbose = TRUE){
    
   resid <- nullmod$resid.marginal
   if (!is.null(group.idx)) g <- length(group.idx)
   if (!is.null(covMatList)) m <- length(covMatList)
   
   ## checks that may be put into wrapper:
   if (rescale & is.null(group.idx)) stop("Rescaling is only done by groups, and group indices are missing.") 
    
   if (rankNorm.option == "by.group" & is.null(group.idx)) stop("Cannot rank normalize by group, missing group indices.")
   if (rankNorm.option == "by.group"){
   		for (i in 1:g){
   			group.resids <- rankNorm(resid[group.idx[[i]]])
   			if (rescale){
   				group.var <- .averageGroupVar(nullmod$varComp, covMatList, group.idx[i])
   			} else{
   				group.var <- 1
   			}
   			resid[group.idx[[i]]] <- group.resids*sqrt(group.var)
   			
   		}
   }
   
   if (rankNorm.option == "all"){
   		resid <- rankNorm(resid)
   		if (rescale) { 
   			for (i in 1:g){
	   			group.var <- .averageGroupVar(nullmod$varComp, covMatList, group.idx[i])
   				resid[group.idx[[i]]] <- resid[group.idx[[i]]]*sqrt(group.var)
   			}
   		}    		
   } 
   
   ### not re-fit the null model:
   
   new.nullmod <- fitNullModel(y = resid, X = nullmod$model.matrix, covMatList = covMatList,
   								group.idx = group.idx, family = "gaussian", start = nullmod$varComp, 
   								AIREML.tol = AIREML.tol, maxIter = maxIter, dropZeros = TRUE, 
   								verbose = verbose)
   
    
    new.nullmod$call <- match.call()
    return(new.nullmod)	   
     
    
}


## here group.idx is a list with just one entry - a vector of the group indices.
## the name of this list entry is the name of the group, same as it appears in the varComp vector. 
.averageGroupVar <- function(varComp, covMatList = NULL, group.idx = NULL){

	if (is.null(group.idx)){
		stop("group indices are not provided, cannot calculate average variance in the group")
	} 
		
	if (is.null(covMatList)) {
		m <- 0	
	} else{
		m <- length(covMatList)
	}
	
	## initialize sum of variance components
	sum.var <- 0
	if (m > 0){
		for (i in 1:m){
			sum.var <- sum.var + varComp[i]*mean(diag(covMatList[[i]])[group.idx[[1]]])
		}
	}
	
	## now add residual variance:
	sum.var <- sum.var + varComp[paste0("V_", names(group.idx)[1])]
	
	return(sum.var)
}