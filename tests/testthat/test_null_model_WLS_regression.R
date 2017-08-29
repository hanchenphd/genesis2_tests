context("check null model WLS regression")
require(GENESIS)
require(GWASTools)

test_that("WLS", {
### Checks for the linear regression case:
##### successful! 
n <- 100
X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))

group.idx <- list(G1 = c(1:(n/2)), G2 = c((n/2 + 1):n))
nullmod <- fitNullModel(y, X, group.idx = group.idx, verbose=FALSE)


## compare to GENESIS. I need to create a jusk correlation matrix that will zero out as having variance component zero...
scanData <- ScanAnnotationDataFrame(data = data.frame(scanID = paste0("p", 1:n), y = y, X1 = X[,1], X2 = X[,2], X3 = X[,3], group = c(rep("G1", n/2), rep("G2", n/2))))

varCompJunk <- TRUE
while(varCompJunk){

	cor.mat <- diag(rnorm(n, sd = 0.01))
	dimnames(cor.mat) <- list(scanData$scanID, scanData$scanID)
	lmm.genesis <- fitNullMM(scanData, "y", covars = c("X1", "X2", "X3"), covMatList = cor.mat,group.var = "group", verbose=FALSE)
	if (lmm.genesis$varComp[1] == 0 ) varCompJunk <- FALSE
}

expect_equal(nullmod$family$family, "gaussian")
expect_false(nullmod$family$mixedmodel)
expect_true(nullmod$hetResid)
expect_true(nullmod$converged)
expect_null(nullmod$zeroFLAG) 
expect_equal(nullmod$workingY, y)
expect_equal(nullmod$outcome, y)
expect_true(all(nullmod$model.matrix == X))


## checks - GENESIS:
expect_true(all(abs(nullmod$fixef - lmm.genesis$fixef) < 1e-6))
expect_true(all(abs(nullmod$betaCov - lmm.genesis$betaCov) < 1e-6))
expect_true(all(abs(nullmod$resid.marginal - lmm.genesis$resid.response) < 1e-8))
expect_true(all(abs(nullmod$logLik - lmm.genesis$logLik) < 1e-6))
expect_true(all(abs(nullmod$logLikR == lmm.genesis$logLikR) < 1e-6))

## currently GENESIS has a mistake, in AUC calculation it uses the number of 
## matrices and groups used, but not the actual number of non-zero variance components. 
## so this "2" fixes for it. 
expect_true(all(abs(nullmod$AIC - (lmm.genesis$AIC - 2)) < 1e-6))
expect_true(all(nullmod$workingY == lmm.genesis$workingY))
expect_true(all(nullmod$model.matrix == lmm.genesis$model.matrix))
expect_true(all(abs(nullmod$varComp - lmm.genesis$varComp[2:3]) < 1e-6))
expect_true(all(abs(nullmod$varCompCov - lmm.genesis$varCompCov[2:3, 2:3]) < 1e-6))  
expect_equal(nullmod$family$family, lmm.genesis$family$family)
expect_null(nullmod$zeroFLAG)
expect_true(all(abs(nullmod$cholSigmaInv - lmm.genesis$cholSigmaInv) < 1e-6))
expect_true(all(abs(nullmod$RSS - lmm.genesis$RSS) < 1e-6))


})
