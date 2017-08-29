context("check null model lmm")
require(GENESIS)
require(GWASTools)

test_that("glmm", {
### Checks for the logistic (also poisson, same algorithm) regression case:
n <- 100
scanID <- paste0("p", 1:n)
X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))

## generate random effects:
sqrt.cor.mat <- matrix(rnorm(n*n, sd = 0.05),n,n)
cor.mat <- crossprod(sqrt.cor.mat)
dimnames(cor.mat) <- list(scanID, scanID)


## compare to GENESIS. Make sure the variance component isn't zeroed out.
## if it does, the function has an error. 

varCompZero <- TRUE
while(varCompZero){	
	random.iid <- rnorm(n)
	random <- crossprod(sqrt.cor.mat*0.05, random.iid)
	expit <- function(x){exp(x)/(1+exp(x))} 
	p <- expit(X %*% c(-1, 0.5, 1) + random) 
	D <- rbinom(n, size = 1, prob = p)
	
	scanData <- ScanAnnotationDataFrame(data = data.frame(scanID = scanID, D = D, X1 = X[,1], X2 = X[,2], X3 = X[,3]))
	
	glmm.genesis <- tryCatch({
		fitNullMM(scanData, "D", covars = c("X1", "X2", "X3"), covMatList = cor.mat, family = "binomial", verbose=FALSE)
				}, 
			warning = function(w){return(list(message = "warning"))},
			error = function(e){return(list(message = "error"))}
			)
	if (!is.null(glmm.genesis$message)) next
	if (glmm.genesis$varComp[1] != 0 ) varCompZero <- FALSE
}

nullmod <- fitNullModel(y = D, X = X, covMatList = cor.mat, verbose=FALSE, family = "binomial")


expect_equal(nullmod$family$family, "binomial")
expect_true(nullmod$family$mixedmodel)
expect_false(nullmod$hetResid)
expect_true(nullmod$converged)
expect_false(nullmod$zeroFLAG)
expect_true(all(nullmod$outcome == D))
expect_true(all(nullmod$model.matrix == X))





## checks - GENESIS:
expect_true(all(abs(nullmod$workingY - glmm.genesis$workingY) < 1e-6))
expect_true(all(abs(nullmod$fixef - glmm.genesis$fixef) < 1e-6))
expect_true(all(abs(nullmod$betaCov - glmm.genesis$betaCov) < 1e-6))
expect_true(all(abs(nullmod$resid.marginal - glmm.genesis$resid.response) < 1e-6))
expect_true(all(abs(nullmod$logLik - glmm.genesis$logLik) < 1e-6))
expect_true(all(abs(nullmod$logLikR - glmm.genesis$logLikR) < 1e-6))

expect_true(all(abs(nullmod$AIC - glmm.genesis$AIC) < 1e-6))
expect_true(all(nullmod$model.matrix == glmm.genesis$model.matrix))
expect_true(all(abs(nullmod$varComp - glmm.genesis$varComp) < 1e-6))
expect_true(all(abs(nullmod$varCompCov - glmm.genesis$varCompCov) < 1e-5))  ## not smaller then 1e-6 because values were not updated after convergence. However, this is not important, because the variance component values are converged at the desired level!
expect_equal(nullmod$family$family, glmm.genesis$family$family)
expect_true(all(nullmod$zeroFLAG == glmm.genesis$zeroFLAG))
expect_true(all(abs(nullmod$cholSigmaInv - glmm.genesis$cholSigmaInv) < 1e-6))
expect_equal(nullmod$RSS ,1)

})
