context("check null model prep")

test_that("nullModelTestPrep vs calculateProjection", {
	n <- 100
	X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
	y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))
	
	group.idx <- list(G1 = c(1:(n/2)), G2 = c((n/2 + 1):n))
	
	cor.mat <- matrix(rnorm(n*n, sd = 0.05),n,n, dimnames=list(1:n, 1:n))
	cor.mat <- crossprod(cor.mat)
	covMatList <- list(A = cor.mat)

        dat <- as.data.frame(cbind(scanID=1:n, y, X, group=1))
        names(dat)[2:5] <- c("y", paste0("X",1:3))
        dat$group[group.idx[[2]]] <- 2
        
        # basic
	nullmod <- fitNullModel(y, X, verbose=FALSE)
	nullprep <- nullModelTestPrep(nullmod)
 
        nullmod.orig <- GENESIS::fitNullReg(dat, outcome="y", covars=paste0("X",1:3), verbose=FALSE)
        proj <- GENESIS:::.calculateProjection(nullmod.orig, test="", burden.test="")

        expect_true(all(abs(nullprep$Mt - proj$Mt) < 1e-9))
        #expect_true(all(abs(nullprep$resid - proj$resid) < 1e-9))
        
        # with covMatList
	nullmod <- fitNullModel(y, X, covMatList=covMatList, verbose=FALSE)
	nullprep <- nullModelTestPrep(nullmod)

        nullmod.orig <- GENESIS::fitNullMM(dat, outcome="y", covars=paste0("X",1:3), covMatList=covMatList, verbose=FALSE)
        proj <- GENESIS:::.calculateProjection(nullmod.orig, test="", burden.test="")

        expect_true(all(abs(nullprep$Mt - proj$Mt) < 1e-9))
        expect_true(all(abs(nullprep$resid - proj$resid) < 1e-9))
        
        # with group
	nullmod <- fitNullModel(y, X, group.idx = group.idx, covMatList=covMatList, verbose=FALSE)
	nullprep <- nullModelTestPrep(nullmod)

        nullmod.orig <- GENESIS::fitNullMM(dat, outcome="y", covars=paste0("X",1:3), covMatList=covMatList, group.var="group", verbose=FALSE)
        proj <- GENESIS:::.calculateProjection(nullmod.orig, test="", burden.test="")

        expect_true(all(abs(nullprep$Mt - proj$Mt) < 1e-9))
        expect_true(all(abs(nullprep$resid - proj$resid) < 1e-9))
        
}
