context("check variant set association tests")

.testNullPrep <- function(n, MM=FALSE) {
	X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
	y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))
	
	group.idx <- list(G1 = c(1:(n/2)), G2 = c((n/2 + 1):n))

        if (MM) {
            cor.mat <- matrix(rnorm(n*n, sd = 0.05),n,n)
            cor.mat <- crossprod(cor.mat)
            covMatList <- list(A = cor.mat)
        } else {
            covMatList <- NULL
        }
	
	nullmod <- fitNullModel(y, X, covMatList = covMatList, group.idx = group.idx, verbose=FALSE)
	
	nullModelTestPrep(nullmod)
}


test_that("SKAT with rho=1 matches burden", {
	n <- 100

	### create a matrix of genetic variants to test.
	geno <- matrix(rbinom(200*n, size = 2, prob = 0.2), nrow = n, ncol = 200)
        weights <- rep(1, ncol(geno))

        ## mixed model
	nullprep <- .testNullPrep(n, MM=TRUE)
        skat <- .testVariantSetSKAT(nullprep, geno, weights, rho=1, pval.method="davies")
        burden <- .testVariantSetBurden(nullprep, geno, weights, burden.test="Score")
        expect_true(abs(skat$pval_1 - burden$Score.pval) < 0.01)
        
        ## basic
	nullprep <- .testNullPrep(n, MM=FALSE)
        skat <- .testVariantSetSKAT(nullprep, geno, weights, rho=1, pval.method="davies")
        burden <- .testVariantSetBurden(nullprep, geno, weights, burden.test="Score")
        expect_true(abs(skat$pval_1 - burden$Score.pval) < 0.01)
})
