



.calcAIcovMats <- function(Y, P, PY, m, covMatList){
	AI <- matrix(NA, nrow =  m, ncol = m)
	score <- rep(NA, m)
	
	for (i in 1:m){
		PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
		score[i] <- -0.5*(sum(P*covMatList[[i]]) - crossprod(Y, PAPY)) # tr(PA) - YPAPY
		AI[i,i] <- 0.5*crossprod(PY, crossprod(covMatList[[i]],PAPY)) # YPAPAPY
		if((i+1) <= m){
			for(j in (i+1):m){
				AI[i,j] <- 0.5*crossprod(PY, crossprod(covMatList[[j]],PAPY)) # YPDPAPY
				AI[j,i] <- AI[i,j]
			}
		}
	}
	return(list(AI = AI, score = score))
}	
	
	
.calcAIcovMatsResids <- function(P, PY, m, covMatList, g, group.idx){
	AI <- matrix(NA, nrow =  m, ncol = g)
	
	for(i in 1:m){
		PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
		if(g == 1){
			AI[i,1] <- 0.5*crossprod(PY, PAPY) # YPIPAPY
		}else{
			for(j in 1:g){
				AI[i,j] <- 0.5*crossprod(PY[group.idx[[j]]], PAPY[group.idx[[j]]]) # YP(I_group)PAPY
			}
		}
	}

	return(AI)
}
	





.calcAIhetvars <- function(P, PY, g, group.idx){
	
	if (g == 1){
		score <-  -0.5*(sum(diag(P)) - crossprod(PY)) 
		AI <- 0.5*crossprod(PY,crossprod(P,PY))
	} else{	
		AI <- matrix(NA, nrow =  g, ncol = g)
		score <- rep(NA, g)
	            
		PAPY <- crossprod(P, PY)
		for (i in 1:g) {
			PIPY <- crossprod(P[group.idx[[i]], ], PY[group.idx[[i]]]) 
			score[ i] <- -0.5 * (sum(diag(P)[group.idx[[i]]]) - crossprod(PY[group.idx[[i]]]))
			AI[ i,  i] <- 0.5 * crossprod(PY[group.idx[[i]]], PIPY[group.idx[[i]]])
			
			if ((i + 1) <= g) {
				for (j in (i + 1):g) {
					AI[i, j] <- 0.5 * crossprod(PY[group.idx[[j]]], PIPY[group.idx[[j]]]) 
					AI[j, i] <- AI[i, j]
				}
			}
		}
	}
	return(list(AI = AI, score = score))
}

