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
y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = 2)

nullmod <- fitNullModel(y, X)

# compare to linear regression
lm.mod <- lm(y ~ -1 + X)

## compare to GENESIS:
scanData <- ScanAnnotationDataFrame(data = data.frame(scanID = paste0("p", 1:n), y = y, X1 = X[,1], X2 = X[,2], X3 = X[,3]))
lm.genesis <- fitNullReg(scanData, "y", covars = c("X1", "X2", "X3"))



## checks - linear regression:
all(nullmod$fitted.values == fitted(lm.mod))
nullmod$family$family == "gaussian"
nullmod$family$mixedmodel == FALSE
nullmod$hetResid == FALSE
all(nullmod$resid.marginal == lm.mod$resid)
all(nullmod$fixef == summary(lm.mod)$coef)
nullmod$varComp == summary(lm.mod)$sigma^2
is.null(nullmod$varCompCov)
all(nullmod$betaCov == vcov(lm.mod))
all(nullmod$fitted.values == fitted(lm.mod))
nullmod$logLik == as.numeric(logLik(lm.mod))
nullmod$AIC == AIC(lm.mod)
all(nullmod$workingY == y)
all(nullmod$outcome == y)
all(nullmod$model.matrix == X)
nullmod$cholSigmaInv == 1/summary(lm.mod)$sigma
nullmod$converged == TRUE
is.null(nullmod$zeroFLAG) 
nullmod$RSS == sum(lm.mod$resid^2)/(summary(lm.mod)$sigma^2*(n - ncol(X)))



## checks - GENESIS:
all(nullmod$fixef ==  lm.genesis$fixef)
all(nullmod$betaCov == lm.genesis$betaCov)
all(nullmod$resid.marginal == lm.genesis$resid.response)
all(nullmod$logLik == lm.genesis$logLik)
all(nullmod$AIC == lm.genesis$AIC)
all(nullmod$workingY == lm.genesis$workingY)
all(nullmod$model.matrix == lm.genesis$model.matrix)
all(nullmod$varComp == lm.genesis$sigma^2)
nullmod$family$family == lm.genesis$family$family


