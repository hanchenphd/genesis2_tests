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


expit <- function(x){exp(x)/(1+exp(x))}
p <- expit(X %*% c(-1, 0.5, 1))
D <- rbinom(n, size = 1, prob = p)

nullmod <- fitNullModel(D, X, family = "binomial")

# compare to logistic regression
glm.mod <- glm(D ~ -1 + X, family = "binomial")

## compare to GENESIS:
scanData <- ScanAnnotationDataFrame(data = data.frame(scanID = paste0("p", 1:n), D = D, X1 = X[,1], X2 = X[,2], X3 = X[,3]))
glm.genesis <- fitNullReg(scanData, "D", covars = c("X1", "X2", "X3"), family = "binomial")

## checks - logistic regression:
nullmod$family$family == "binomial"
nullmod$hetResid == FALSE
all(nullmod$fitted.values == fitted(glm.mod))
nullmod$family$mixedmodel == FALSE
all(nullmod$resid.marginal == resid(glm.mod, type = "response"))
all(nullmod$fixef == summary(glm.mod)$coef)
all(nullmod$varComp == fitted(glm.mod)*(1-fitted(glm.mod)))
is.null(nullmod$varCompCov)
all(nullmod$betaCov == vcov(glm.mod))
all(nullmod$fitted.values == fitted(glm.mod))
nullmod$logLik == as.numeric(logLik(glm.mod))
nullmod$AIC == AIC(glm.mod)
all(nullmod$workingY == D)
all(nullmod$outcome == D)
all(nullmod$model.matrix == X)
all(abs(nullmod$cholSigmaInv -1/sqrt(fitted(glm.mod)*(1-fitted(glm.mod)))) < 1e-10)
nullmod$converged ==  glm.mod$converged
is.null(nullmod$zeroFLAG) 
nullmod$RSS == 1



## checks - GENESIS:
all(nullmod$fixef ==  glm.genesis$fixef)
all(nullmod$betaCov == glm.genesis$betaCov)
all(nullmod$resid.marginal == glm.genesis$resid.response)
all(nullmod$logLik == glm.genesis$logLik)
all(nullmod$AIC == glm.genesis$AIC)
all(nullmod$workingY == glm.genesis$workingY)
all(nullmod$model.matrix == glm.genesis$model.matrix)
nullmod$family$family == glm.genesis$family$family
all(abs(nullmod$varComp - glm.genesis$sigma^2) < 1e-10)




