context("check update of null model after ranknormalizing and re-scaling residuals")
require(GENESIS)
require(GWASTools)

test_that("updateOutcome", {
### Checks that updating the outcome and re-fitting the null model works okay.
n <- 100
X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))

scanID <- paste0("p", 1:n)
group.idx <- list(G1 = c(1:(n/2)), G2 = c((n/2 + 1):n))

cor.mat <- matrix(rnorm(n*n, sd = 0.05),n,n)
cor.mat <- crossprod(cor.mat)
dimnames(cor.mat) <- list(scanID, scanID)
covMatList <- list(A = cor.mat)

nullmod <- fitNullModel(y, X, group.idx = group.idx, covMatList, verbose=FALSE)

group.ind <- 1
expect_equal(.averageGroupVar(nullmod$varComp, covMatList = covMatList, group.idx = group.idx[group.ind]),
             nullmod$varComp[1]*mean(diag(cor.mat)[group.idx[[group.ind]]]) + nullmod$varComp[group.ind + 1])

group.ind <- 2
expect_equal(.averageGroupVar(nullmod$varComp, covMatList = covMatList, group.idx = group.idx[group.ind]),
	     nullmod$varComp[1]*mean(diag(cor.mat)[group.idx[[group.ind]]]) + nullmod$varComp[group.ind + 1])

expect_equal(.averageGroupVar(nullmod$varComp, covMatList = NULL, group.idx = group.idx[group.ind]),
             nullmod$varComp[group.ind + 1])
						
throws_error(.averageGroupVar(nullmod$varComp, covMatList = NULL, group.idx = group.idx))
throws_error(.averageGroupVar(nullmod$varComp, covMatList = NULL, group.idx = NULL))
throws_error(.averageGroupVar(nullmod$varComp, covMatList = NULL, group.idx = group.idx[[1]]))

nullmod2 <- updateNullModOutcome(nullmod, covMatList = covMatList, group.idx = group.idx,  rankNorm.option = c("by.group"), rescale = c("None"), verbose=FALSE)


expect_true(abs(.averageGroupVar(nullmod2$varComp, covMatList = covMatList, group.idx = group.idx[1]) - 1) < 0.1 )
expect_true(abs(.averageGroupVar(nullmod2$varComp, covMatList = covMatList, group.idx = group.idx[2]) -1) < 0.1 )

# why are these tests not < 0.1???
nullmod3 <- updateNullModOutcome(nullmod, covMatList = covMatList, group.idx = group.idx,  rankNorm.option = c("by.group"), rescale = c("residSD"), verbose=FALSE)
expect_true(abs(.averageGroupVar(nullmod3$varComp, covMatList = covMatList, group.idx = group.idx[1]) - var(nullmod$resid.marginal[group.idx[[1]]] ) ) < 0.1 )
expect_true(abs(.averageGroupVar(nullmod3$varComp, covMatList = covMatList, group.idx = group.idx[2]) - var(nullmod$resid.marginal[group.idx[[2]]] ) ) < 0.1 )

nullmod4 <- updateNullModOutcome(nullmod, covMatList = covMatList, group.idx = group.idx,  rankNorm.option = c("by.group"), rescale = c("model"), verbose=FALSE)
expect_true(abs(.averageGroupVar(nullmod4$varComp, covMatList = covMatList, group.idx = group.idx[1]) - .averageGroupVar(nullmod$varComp, covMatList = covMatList, group.idx = group.idx[1]) ) < 0.1 )
expect_true(abs(.averageGroupVar(nullmod4$varComp, covMatList = covMatList, group.idx = group.idx[2]) - .averageGroupVar(nullmod$varComp, covMatList = covMatList, group.idx = group.idx[2]) ) < 0.1 )



})
