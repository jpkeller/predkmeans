context("mlogit")

test_that("demo", {
	expect_equal(nchar("ab"), 2)
})


Y <- cbind(rep(c(0, 1), each=20), rep(c(1, 0), each=20))
set.seed(10)
X <- cbind(1, rnorm(40))


# Check that Loglike functions agree
llrcpp <- logLikeMultinomial_Rcpp(beta=0, Y=Y, X=X)
ll <- logLikeMultinomial(beta=0, Y=Y, X=X)
test_that("LogLike Fns agree", {
	expect_equal(llrcpp, ll)
})

# Test output class of mlogit
mm <- mlogit(Y=Y, X=X, verbose=0)
test_that("mlogit output class", {
	expect_is(mm, "mlogit")
})
