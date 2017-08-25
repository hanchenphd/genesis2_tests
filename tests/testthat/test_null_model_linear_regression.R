context("check null model linear regression")
require(GENESIS)
require(GWASTools)

test_that("linear", {
### Checks for the linear regression case:
##### successful! 
n <- 100
X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = 2)

nullmod <- fitNullModel(y, X, verbose=FALSE)

# compare to linear regression
lm.mod <- lm(y ~ -1 + X)

## compare to GENESIS:
scanData <- ScanAnnotationDataFrame(data = data.frame(scanID = paste0("p", 1:n), y = y, X1 = X[,1], X2 = X[,2], X3 = X[,3]))
lm.genesis <- fitNullReg(scanData, "y", covars = c("X1", "X2", "X3"), verbose=FALSE)



## checks - linear regression:
expect_equal(nullmod$fitted.values, fitted(lm.mod))
expect_equal(nullmod$family$family, "gaussian")
expect_false(nullmod$family$mixedmodel)
expect_false(nullmod$hetResid)
expect_true(all(nullmod$resid.marginal == lm.mod$resid))
expect_true(all(nullmod$fixef == summary(lm.mod)$coef))
expect_equal(nullmod$varComp, summary(lm.mod)$sigma^2)
expect_null(nullmod$varCompCov)
expect_true(all(nullmod$betaCov == vcov(lm.mod)))
expect_true(all(nullmod$fitted.values == fitted(lm.mod)))
expect_equal(nullmod$logLik, as.numeric(logLik(lm.mod)))
expect_equal(nullmod$AIC, AIC(lm.mod))
expect_true(all(nullmod$workingY == y))
expect_true(all(nullmod$outcome == y))
expect_true(all(nullmod$model.matrix == X))
expect_equal(nullmod$cholSigmaInv, 1/summary(lm.mod)$sigma)
expect_true(nullmod$converged)
expect_null(nullmod$zeroFLAG)
expect_equal(nullmod$RSS, sum(lm.mod$resid^2)/(summary(lm.mod)$sigma^2*(n - ncol(X))))



## checks - GENESIS:
expect_true(all(nullmod$fixef == lm.genesis$fixef))
expect_true(all(nullmod$betaCov == lm.genesis$betaCov))
expect_true(all(nullmod$resid.marginal == lm.genesis$resid.response))
expect_true(all(nullmod$logLik == lm.genesis$logLik))
expect_true(all(nullmod$AIC == lm.genesis$AIC))
expect_true(all(nullmod$workingY == lm.genesis$workingY))
expect_true(all(nullmod$model.matrix == lm.genesis$model.matrix))
expect_true(all(nullmod$varComp == lm.genesis$sigma^2))
expect_equal(nullmod$family$family, lm.genesis$family$family)
})

