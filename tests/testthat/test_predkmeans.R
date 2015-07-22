context("predkmeans")


set.seed(10)
Y <- cbind(1:20, rnorm(20))
X <- cbind(1, 1:20)


# Check predkmeans warnings
test_that("Predkmeans Warnings",{
	expect_error(predkmeans(X=Y, R=X), "\"K\" is missing")
	expect_error(predkmeans(X=Y, R=X, K=25), "more clusters than observations")
})




#predkmeans(X=Y, R=X, K=25)