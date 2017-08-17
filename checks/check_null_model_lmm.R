### set working directory...
setwd()

source("calcAIcovMatsFunctions.R")
source("calcLQ.R")
source("fitNullModel.R")
source("prepareOutputForNullModel.R")
source("runAIREMLgaussian.R")
source("runAIREMLother.R")
source("runWLS.R")


require(GENESIS)
require(GWASTools)

### Checks for the linear regression case:
##### successful! 
n <- 100
X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))

group.idx <- list(G1 = c(1:(n/2)), G2 = c((n/2 + 1):n))



## compare to GENESIS. I need to create a jusk correlation matrix that will zero out as having variance component zero...
scanData <- ScanAnnotationDataFrame(data = data.frame(scanID = paste0("p", 1:n), y = y, X1 = X[,1], X2 = X[,2], X3 = X[,3], group = c(rep("G1", n/2), rep("G2", n/2))))

varCompJunk <- FALSE
while(!varCompJunk){

	cor.mat <- matrix(rnorm(n*n, sd = 0.05),n,n)
	cor.mat <- crossprod(cor.mat)
	dimnames(cor.mat) <- list(scanData$scanID, scanData$scanID)
	lmm.genesis <- GENESIS:::fitNullMM(scanData, "y", covars = c("X1", "X2", "X3"), covMatList = cor.mat,group.var = "group")
	if (lmm.genesis$varComp[1] != 0 ) varCompJunk <- TRUE
}

nullmod <- fitNullModel(y, X, group.idx = group.idx, cor.mat)


nullmod$family$family == "gaussian"
nullmod$family$mixedmodel == FALSE
nullmod$hetResid == TRUE
nullmod$converged == TRUE
is.null(nullmod$zeroFLAG) 
all(nullmod$workingY == y)
all(nullmod$outcome == y)
all(nullmod$model.matrix == X)





## checks - GENESIS:
all(abs(nullmod$fixef - lmm.genesis$fixef) < 1e-9)
all(abs(nullmod$betaCov - lmm.genesis$betaCov) < 1e-9)
all(abs(nullmod$resid.marginal - lmm.genesis$resid.response) < 1e-9)
all(abs(nullmod$logLik - lmm.genesis$logLik) < 1e-9)
all(abs(nullmod$logLikR - lmm.genesis$logLikR) < 1e-9)

## currently GENESIS has a mistake, in AUC calculation it uses the number of 
## matrices and groups used, but not the actual number of non-zero variance components. 
## so this "2" fixes for it. 
all(abs(nullmod$AIC - lmm.genesis$AIC) < 1e-9)
all(nullmod$workingY == lmm.genesis$workingY)
all(nullmod$model.matrix == lmm.genesis$model.matrix)
all(abs(nullmod$varComp - lmm.genesis$varComp) < 1e-9)
all(abs(nullmod$varCompCov - lmm.genesis$varCompCov) < 1e-9) 
nullmod$family$family == lmm.genesis$family$family
all(nullmod$zeroFLAG == lmm.genesis$zeroFLAG)
all(abs(nullmod$cholSigmaInv - lmm.genesis$cholSigmaInv) < 1e-9)
all(abs(nullmod$RSS - lmm.genesis$RSS) < 1e-9)

