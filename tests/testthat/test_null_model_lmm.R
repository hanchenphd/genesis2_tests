context("check null model lmm")
require(GENESIS)
require(GWASTools)

test_that("lmm", {
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
	lmm.genesis <- fitNullMM(scanData, "y", covars = c("X1", "X2", "X3"), covMatList = cor.mat,group.var = "group", verbose=FALSE)
	if (lmm.genesis$varComp[1] != 0 ) varCompJunk <- TRUE
}

nullmod <- fitNullModel(y, X, group.idx = group.idx, cor.mat, verbose=FALSE)


expect_equal(nullmod$family$family, "gaussian")
expect_true(nullmod$family$mixedmodel)
expect_true(nullmod$hetResid)
expect_true(nullmod$converged)
#expect_null(nullmod$zeroFLAG)
expect_true(all(nullmod$workingY == y))
expect_true(all(nullmod$outcome == y))
expect_true(all(nullmod$model.matrix == X))





## checks - GENESIS:
expect_true(all(abs(nullmod$fixef - lmm.genesis$fixef) < 1e-9))
expect_true(all(abs(nullmod$betaCov - lmm.genesis$betaCov) < 1e-9))
expect_true(all(abs(nullmod$resid.marginal - lmm.genesis$resid.response) < 1e-9))
expect_true(all(abs(nullmod$logLik - lmm.genesis$logLik) < 1e-9))
expect_true(all(abs(nullmod$logLikR - lmm.genesis$logLikR) < 1e-9))

## currently GENESIS has a mistake, in AUC calculation it uses the number of 
## matrices and groups used, but not the actual number of non-zero variance components. 
## so this "2" fixes for it. 
expect_true(all(abs(nullmod$AIC - lmm.genesis$AIC) < 1e-9))
expect_true(all(nullmod$workingY == lmm.genesis$workingY))
expect_true(all(nullmod$model.matrix == lmm.genesis$model.matrix))
expect_true(all(abs(nullmod$varComp - lmm.genesis$varComp) < 1e-9))
expect_true(all(abs(nullmod$varCompCov - lmm.genesis$varCompCov) < 1e-9)) 
expect_equal(nullmod$family$family, lmm.genesis$family$family)
expect_true(all(nullmod$zeroFLAG == lmm.genesis$zeroFLAG))
expect_true(all(abs(nullmod$cholSigmaInv - lmm.genesis$cholSigmaInv) < 1e-9))
expect_true(all(abs(nullmod$RSS - lmm.genesis$RSS) < 1e-9))


### test updating a conditional model: 
G = matrix(rnorm(100, 100,1))
nullmod2 <- updateNullModCond(nullmod, G, covMatList = list(cor.mat), group.idx = group.idx, AIREML.tol = 1e-7, verbose=FALSE)
nullmod3 <- fitNullModel(y, cbind(X, G), group.idx = group.idx, cor.mat, AIREML.tol = 1e-7, verbose=FALSE)

expect_true(max(abs(nullmod2$varComp - nullmod3$varComp)) < 1e-6)
expect_true(max(abs(nullmod3$fixef - nullmod2$fixef)) < 1e-5)
expect_true(all(abs(nullmod2$cholSigmaInv - nullmod3$cholSigmaInv) < 1e-5))
expect_true(all(abs(nullmod2$varCompCov - nullmod3$varCompCov) < 1e-5))


})
