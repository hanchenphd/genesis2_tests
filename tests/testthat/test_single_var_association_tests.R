context("check single variant association tests")

test_that("singleVarTest", {
	n <- 100
	X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
	y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))
	
	group.idx <- list(G1 = c(1:(n/2)), G2 = c((n/2 + 1):n))
	
	nullmod <- fitNullModel(y, X, group.idx = group.idx, verbose=FALSE)
	
	nullprep <- nullModelTestPrep(nullmod)

	### create a matrix of genetic variants to test.
	geno <- matrix(rbinom(200*n, size = 2, prob = 0.2), nrow = n, ncol = 200)

	test.wald <- testGenoSingleVar(nullprep, G = geno, E = NULL, test = c("Wald"), GxE.return.cov = FALSE)
	
	# compare to weighted least squares using weighted lm:
	
	res.lm <- test.wald
	for (i in 1:ncol(geno)){
		lm.temp <- lm(y ~ -1 + X + geno[,i], weights = c(rep(1/nullmod$varComp[1], n/2), 1/rep(nullmod$varComp[1], n/2)))
		res.lm[i,] <- summary(lm.temp)$coef[4,]
	}
	
	expect_true(all(abs(res.lm$Est - test.wald$Est ) < 1e-8))
	expect_true(all(abs(res.lm$Wald.Stat - test.wald$Wald.Stat ) < 1e-8))
})
