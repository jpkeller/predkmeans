#####################################
## Functions for CV of predkmeans  ##
#####################################
##
## This file contains:
##
##	fitCV()
##	predictCV()


##
##' @name fitCV
## 
##' @title Cross-validation of Predictive K-means Clustering
##
##' @description Performs cross-validation of predictive-kmeans on a dataset.
##
## Inputs:
##' @param X Outcome data
##' @param R Covariates
##' @param K Number of clusters
##' @param cv.groups A list or matrix providing the cross-validation
##'		groups for splitting the data.  Alternatively, a single
##'		number giving the number of groups into which the
##'		data are randomly split. A value of '0' implies leave-one-out.
##'		Defaults to 10.
##'	@param doKmeans Should kmeans also be done, in addition to predictive
##'		kmeans? Note: If sigma2=0 and sigma2fixed=TRUE, then doKmeans is
##'		is set to TRUE, and no predictive k-means is done.
##'	@param scale Should the outcomes be re-scaled within each training
##'		group?
##' @param mu = starting values for centers. If NULL (default), 
##'  	then value is chosen according to \code{muStart}.
##' @param muStart Character string indicating how initial value
##'  	of mu should be selected. Only used if mu=NULL.  Possible
##' 	values are "random" (default) or "kmeans".
##' @param sigma2 starting value of sigma2 (ADD SPECIAL CASE OF 0 here)
##' @param sigma2fixed Logical indicating whether sigma2
##'		should be held fixed.  If FALSE, then
##'		sigma2 is estimated using Maximum Likelihood.
##' @param maxitEM Maximum number of EM iterations for
##' 	finding the Mixture of Experts solution
##' @param tol convergence criterion
##' @param maxitMlogit Maximum number of iterations in the=
##'   	mlogit optimization (nested within EM algorithm)
##' @param muRestart  Gives max number of attempts at picking
##'    	starting values. Only used when muStart='random'.
##'		If selected starting values for mu are constant
##'     within each cluster, then the starting values
##'     are re-selected up to muRestart times.
##' @param  convEM  controls the measure of convergence for the 
##'     EM algorithm.  Should be one of "mu", "gamma", or
##'	    "both".  Defaults to "both."  The EM algorithm 
##'     stops when the Frobenius norm of the change in mu,
##'		the change in gamma, or the change in mu and the change in gamma
##'		is less than 'tol'.
##' @param  verbose numeric vector indicating how much output to produce
##' @param nStarts number of times to perform EM algorithm 
##' @param returnAll A list containing all \code{nStarts} solutions is
##'		included in the output.
##' @param ... Additional arguments passed to \code{\link{mlogit}}
##
##' @details The algorithm for sovling the mixture of Experts model is 
##'		based upon the approach presented by Jordan and Jacobs (1994).  If \code{sigma2} is 0 and \code{sigm2fixed} is TRUE, then kmeans clustering is done instead.
##
##' @family 'predkmeans methods'
##' @seealso \code{\link{predictML.predkmeans}}
##' @references Jordan, M. and R. Jacobs (1994). Hierarchical mixtures of
##'		experts and the EM algorithm. \emph{Neural computation 6}(2),
##'		 181-214.
##
##' @author Joshua Keller
##' @export
##'	@return  An object of class \code{predkmeans}, containing the following elements:
##' \item{res.best}{A list containing the results from the best-fitting solution to the Mixture of Experts problem. (TO BE COMPLETED)}
##' \item{center}{Matrix of cluster centers}
##' \item{cluster}{Vector of cluster labels assigned to observations}
##' \item{K}{Number of clusters}
##'	\item{sigam2}{Final value of sigma^2.}
##' \item{wSS}{Mean within-cluster sum-of-squares}
##' \item{sigma2fixed}{Logical indicator of whether sigma2 was held fixed}
##list(mu=mu, gamma=gamma, sigma2=sigma2, conv=converged, objective= obj, iter=iter, mfit=mfit[c("beta", "fitted01", "fitted", "res.best", "status")])
#		Author: J. Keller
#		Original Date: October 2015
#
fitCV <- function(X, R, K, cv.groups=10, doKmeans=FALSE, scale=TRUE, covarnames=NULL, PCA=FALSE, pcacovarnames=colnames(R), npca=5, TPRS=TRUE,TPRSdf=5, TPRSxname="x", TPRSyname="y",sigma2=0,  sigma2fixed=FALSE, 
 mu=NULL, muStart=c("kmeans","random"), maxitEM=100, tol=1e-5, convEM=c("both", "mu", "gamma"), nStarts=1, maxitMlogit=500,verbose=0, muRestart=1000, returnAll=FALSE, ...) {
	
	fncall <- match.call()
	
	
if(is.null(K)){
	stop("K not provided.  Please provide.")
}

# Checks for cv.groups.
# List, matrix, numeric.	
# Assume list for now.

ids <- rownames(X)

if(sigma2==0 & sigma2fixed){
	doPredKmeans <- FALSE
	doKmeans <- TRUE
} else {
	doPredKmeans <- TRUE
}

setup <- vector("list", length(cv.groups))
if (doKmeans) kres <- vector("list", length(cv.groups))
if (doPredKmeans) pres <- vector("list", length(cv.groups))

for (i in 1:length(cv.groups)){
	test.set <- cv.groups[[i]]
	training.set <- setdiff(ids, test.set)

	training.data <- X[training.set, , drop=FALSE]
	test.data <- X[test.set, , drop=FALSE]
    if (scale){
        training.data <- scale(training.data)
        test.data <- scale(test.data, center=attributes(training.data)[["scaled:center"]], scale=attributes(training.data)[["scaled:scale"]])
        scale.center <- attributes(training.data)[["scaled:center"]]
        scale.scale <- attributes(training.data)[["scaled:scale"]]
    } else {
        scale.center <- NULL
        scale.scale <- NULL
    }      

	# Create Training and Test Covariates
	training.covars <- R[training.set, covarnames, drop=FALSE]
	test.covars <- R[test.set, covarnames, drop=FALSE]
	if (PCA){
		pcafit <- get_PCA_matrix(data=R[training.set, pcacovarnames, drop=FALSE], ncomps=npca, matrixonly=FALSE)
		training.covars <- cbind(training.covars, pcafit$X)
		test.covars <- cbind(test.covars , predict(pcafit$pca, newdata= R[test.set,])[, 1:npca, drop=FALSE])
	}
	if (TPRS){
	    tprs.train <- get_TPRS_modelmatrix(data=cbind(training.covars, R[training.set,]), TPRSdf=TPRSdf, covarnames=colnames(training.covars), xname=TPRSxname, yname=TPRSyname, matrixonly=FALSE)
	    training.covars <- tprs.train$X
	    test.covars <- predict(tprs.train$gamfit, newdata=cbind(R[test.set,], test.covars), type="lpmatrix")
	} else{
		training.covars <- cbind(1, training.covars)	
		names(training.covars)[1] <- "(Intercept)"
		training.covars <- as.matrix(training.covars)
		test.covars <- cbind(1, test.covars)	
		names(test.covars)[1] <- "(Intercept)"
		test.covars <- as.matrix(test.covars)
	}

	if (doKmeans){
		kres.i <- predkmeans(X=training.data, R=training.covars, K=K, nStarts=nStarts, sigma2=0, sigma2fixed=TRUE, ...)
		kres[[i]] <- kres.i
		kres[[i]]$test.set <- test.set
		kres[[i]]$training.set <- training.set
		kres[[i]]$scale.center <- scale.center
		kres[[i]]$scale.scale <- scale.scale		
	} else {
		kres <- NULL
	}

	if (doPredKmeans) {
		pres.i <- predkmeans(X=training.data, R=training.covars, K=K, sigma2=sigma2, sigma2fixed=sigma2fixed, nStarts=nStarts,...)
		pres[[i]] <- pres.i
		pres[[i]]$test.set <- test.set
		pres[[i]]$training.set <- training.set	
		pres[[i]]$scale.center <- scale.center
		pres[[i]]$scale.scale <- scale.scale
	} else {
		pres <- NULL
	}
	
	setup[[i]] <- list(training.set=training.set, test.set=test.set, scale.center=scale.center, scale.scale=scale.scale)
}
		
return(list(call=fncall, kres=kres, pres=pres, setup=setup))
}	


predictCV <- function(fitCV, X, R, method=c("ML", "MixExp", "SVM"), doMetrics =FALSE, covarnames=NULL, PCA=FALSE, pcacovarnames=colnames(R), npca=5, TPRS=TRUE,TPRSdf=5, TPRSxname="x", TPRSyname="y", ...) {

predpres <- vector("list", length(fitCV$setup))
metrics <- predpres
setup <- fitCV$setup 
for (i in 1:length(fitCV$setup)){
	test.set <- setup[[i]]$test.set
	training.set <- setup[[i]]$training.set

	training.data <- X[training.set, , drop=FALSE]
	test.data <- X[test.set, , drop=FALSE]
    if (!is.null(setup[[i]]$scale.center)){
        test.data <- scale(test.data, center=setup[[i]]$scale.center, scale=setup[[i]]$scale.scale)
    }
    
   	# Create Training and Test Covariates
	training.covars <- R[training.set, covarnames, drop=FALSE]
	test.covars <- R[test.set, covarnames, drop=FALSE]
	if (PCA){
		pcafit <- get_PCA_matrix(data=R[training.set, pcacovarnames, drop=FALSE], ncomps=npca, matrixonly=FALSE)
		training.covars <- cbind(training.covars, pcafit$X)
		test.covars <- cbind(test.covars , predict(pcafit$pca, newdata= R[test.set,])[, 1:npca, drop=FALSE])
	}
	if (TPRS){
	    tprs.train <- get_TPRS_modelmatrix(data=cbind(training.covars, R[training.set,]), TPRSdf=TPRSdf, covarnames=colnames(training.covars), xname=TPRSxname, yname=TPRSyname, matrixonly=FALSE)
	    training.covars <- tprs.train$X
	    test.covars <- predict(tprs.train$gamfit, newdata=cbind(R[test.set,], test.covars), type="lpmatrix")
	} else{
		training.covars <- cbind(1, training.covars)	
		names(training.covars)[1] <- "(Intercept)"
		training.covars <- as.matrix(training.covars)
		test.covars <- cbind(1, test.covars)	
		names(test.covars)[1] <- "(Intercept)"
		test.covars <- as.matrix(test.covars)
	}  

	predfn <- switch(method, ML=predictML.predkmeans, MixExp=predictMixExp.predkmeans, SVM=predictSVM.predkmeans)
	
	predpres[[i]] <- predfn(object=fitCV$pres[[i]], Rstar=training.covars, R=test.covars,...)	
	if (doMetrics)	metrics[[i]] <- predictionMetrics(centers=fitCV$pres[[i]]$centers, cluster.pred=predpres[[i]]$test.pred,X=test.data, labels=FALSE)
}	

if(doMetrics){
	return(list(out=predpres, metrics=metrics))
} else {
	return(predpres)
}
}	




get_PCA_matrix <- function(data, ncomps, covarnames=colnames(data), center=TRUE, scale=TRUE, matrixonly=TRUE){
	pca <- prcomp(data[,covarnames, drop=FALSE], center=center, scale=scale)
	X <- pca$x[,1:ncomps, drop=FALSE]
	if (matrixonly) {
		return(X)
	} else {
		return(list(X=X, pca=pca))
	}
}


# Note k=df+1
get_TPRS_modelmatrix <- function(data, TPRSdf=5, covarnames=NULL, xname="lambert_x", yname="lambert_y", matrixonly=TRUE, TPRSfx=TRUE){
	require(mgcv)	

	TPRSk <- TPRSdf + 1
	# Create the formula
	f <- ""
	if (length(covarnames)>0){
		f <- paste("+ ", paste0(covarnames, collapse=" + "))
	}
	f <- formula(paste0(xname,"~s(", xname,", ", yname,", fx=TPRSfx, k=TPRSk)", f))
	# Fit a GAM to get the model matrix
	gamfit <- gam(f, data=as.data.frame(data))
	X  <- model.matrix(gamfit)

	if(matrixonly){
		return(X)
	} else {
		return(list(X=X, gamfit=gamfit))
	}
}
