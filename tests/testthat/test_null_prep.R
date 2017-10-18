context("check null model prep")

test_that("nullModelTestPrep", {
	n <- 100
	X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
	y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))

	cor.mat <- matrix(rnorm(n*n, sd = 0.05),n,n, dimnames=list(1:n, 1:n))
	cor.mat <- crossprod(cor.mat)

        # basic
	nullmod <- fitNullModel(y, X, verbose=FALSE)
	nullprep <- nullModelTestPrep(nullmod)

	expect_equal(nullprep$k, ncol(X))
	expect_equal(dim(nullprep$Mt), c(n, n))
	expect_equal(dim(nullprep$Ytilde), dim(y))

        # with covMatList
	nullmod <- fitNullModel(y, X, covMatList=cor.mat, verbose=FALSE)
	nullprep <- nullModelTestPrep(nullmod)

	expect_equal(nullprep$k, ncol(X))
	expect_equal(dim(nullprep$Mt), c(n, n))
	expect_equal(dim(nullprep$Ytilde), dim(y))
})

test_that("nullModelTestPrep vs calculateProjection", {
	n <- 100
	X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
	y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))
	
	group.idx <- list(G1 = c(1:(n/2)), G2 = c((n/2 + 1):n))
	
	cor.mat <- matrix(rnorm(n*n, sd = 0.05),n,n, dimnames=list(1:n, 1:n))
	cor.mat <- crossprod(cor.mat)

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
	nullmod <- fitNullModel(y, X, covMatList=cor.mat, verbose=FALSE)
	nullprep <- nullModelTestPrep(nullmod)

        nullmod.orig <- GENESIS::fitNullMM(dat, outcome="y", covars=paste0("X",1:3), covMatList=cor.mat, verbose=FALSE)
        proj <- GENESIS:::.calculateProjection(nullmod.orig, test="", burden.test="")

        expect_true(all(abs(nullprep$Mt - proj$Mt) < 1e-9))
        expect_true(all(abs(nullprep$resid - proj$resid) < 1e-9))
        
        # with group
	nullmod <- fitNullModel(y, X, group.idx = group.idx, covMatList=cor.mat, verbose=FALSE)
	nullprep <- nullModelTestPrep(nullmod)

        nullmod.orig <- GENESIS::fitNullMM(dat, outcome="y", covars=paste0("X",1:3), covMatList=cor.mat, group.var="group", verbose=FALSE)
        proj <- GENESIS:::.calculateProjection(nullmod.orig, test="", burden.test="")

        expect_true(all(abs(nullprep$Mt - proj$Mt) < 1e-9))
        expect_true(all(abs(nullprep$resid - proj$resid) < 1e-9))
        
})

test_that("nullModelTestPrep vs calculateProjection - binary", {
	n <- 100
	X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
	
	sqrt.cor.mat <- matrix(rnorm(n*n, sd = 0.05),n,n, dimnames=list(1:n, 1:n))
	cor.mat <- crossprod(sqrt.cor.mat)

varCompZero <- TRUE
while(varCompZero){	
	random.iid <- rnorm(n)
	random <- crossprod(sqrt.cor.mat*0.05, random.iid)
	expit <- function(x){exp(x)/(1+exp(x))} 
	p <- expit(X %*% c(-1, 0.5, 1) + random) 
	y <- rbinom(n, size = 1, prob = p)
	
        dat <- as.data.frame(cbind(scanID=1:n, y, X))
        names(dat)[2:5] <- c("y", paste0("X",1:3))
        
	
	glmm.genesis <- tryCatch({
		fitNullMM(dat, "y", covars = c("X1", "X2", "X3"), covMatList = cor.mat, family = "binomial", verbose=FALSE)
				}, 
			warning = function(w){return(list(message = "warning"))},
			error = function(e){return(list(message = "error"))}
			)
	if (!is.null(glmm.genesis$message)) next
	if (glmm.genesis$varComp[1] != 0 ) varCompZero <- FALSE
}

        # basic
	nullmod <- fitNullModel(y, X, family="binomial", verbose=FALSE)
	nullprep <- nullModelTestPrep(nullmod)
 
        nullmod.orig <- GENESIS::fitNullReg(dat, outcome="y", covars=paste0("X",1:3), family="binomial", verbose=FALSE)
        proj <- GENESIS:::.calculateProjection(nullmod.orig, test="", burden.test="")

        expect_true(all(abs(nullprep$Mt - proj$Mt) < 1e-9))
        expect_true(all(abs(nullprep$resid - proj$resid) < 1e-9))
        
        # with covMatList
	nullmod <- fitNullModel(y, X, covMatList=cor.mat, family="binomial", verbose=FALSE)
	nullprep <- nullModelTestPrep(nullmod)

        nullmod.orig <- GENESIS::fitNullMM(dat, outcome="y", covars=paste0("X",1:3), covMatList=cor.mat, family="binomial", verbose=FALSE)
        proj <- GENESIS:::.calculateProjection(nullmod.orig, test="", burden.test="")

        expect_true(all(abs(nullprep$Mt - proj$Mt) < 1e-8))
        expect_true(all(abs(nullprep$resid - proj$resid) < 1e-7))
})
