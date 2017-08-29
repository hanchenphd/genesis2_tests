context("check null model logistic regression")
require(GENESIS)
require(GWASTools)

test_that("logistic", {
### Checks for the logistic regression case:
##### successful! 
n <- 100
X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))


expit <- function(x){exp(x)/(1+exp(x))}
p <- expit(X %*% c(-1, 0.5, 1))
D <- rbinom(n, size = 1, prob = p)

nullmod <- fitNullModel(D, X, family = "binomial", verbose=FALSE)

# compare to logistic regression
glm.mod <- glm(D ~ -1 + X, family = "binomial")

## compare to GENESIS:
scanData <- ScanAnnotationDataFrame(data = data.frame(scanID = paste0("p", 1:n), D = D, X1 = X[,1], X2 = X[,2], X3 = X[,3]))
glm.genesis <- fitNullReg(scanData, "D", covars = c("X1", "X2", "X3"), family = "binomial", verbose=FALSE)

## checks - logistic regression:
expect_equal(nullmod$family$family, "binomial")
expect_false(nullmod$hetResid)
expect_true(all(nullmod$fitted.values == fitted(glm.mod)))
expect_false(nullmod$family$mixedmodel)
expect_true(all(nullmod$resid.marginal == resid(glm.mod, type = "response")))
expect_true(all(nullmod$fixef == summary(glm.mod)$coef))
expect_true(all(nullmod$varComp == fitted(glm.mod)*(1-fitted(glm.mod))))
expect_null(nullmod$varCompCov)
expect_true(all(nullmod$betaCov == vcov(glm.mod)))
expect_true(all(nullmod$fitted.values == fitted(glm.mod)))
expect_equal(nullmod$logLik, as.numeric(logLik(glm.mod)))
expect_equal(nullmod$AIC, AIC(glm.mod))
expect_true(all(nullmod$workingY == D))
expect_true(all(nullmod$outcome == D))
expect_true(all(nullmod$model.matrix == X))
expect_true(all(abs(nullmod$cholSigmaInv -1/sqrt(fitted(glm.mod)*(1-fitted(glm.mod)))) < 1e-10))
expect_equal(nullmod$converged, glm.mod$converged)
expect_null(nullmod$zeroFLAG)
expect_equal(nullmod$RSS, 1)



## checks - GENESIS:
expect_true(all(nullmod$fixef ==  glm.genesis$fixef))
expect_true(all(nullmod$betaCov == glm.genesis$betaCov))
expect_true(all(nullmod$resid.marginal == glm.genesis$resid.response))
expect_true(all(nullmod$logLik == glm.genesis$logLik))
expect_true(all(nullmod$AIC == glm.genesis$AIC))
expect_true(all(nullmod$workingY == glm.genesis$workingY))
expect_true(all(nullmod$model.matrix == glm.genesis$model.matrix))
expect_equal(nullmod$family$family, glm.genesis$family$family)
expect_true(all(abs(nullmod$varComp - glm.genesis$sigma^2) < 1e-10))

})


