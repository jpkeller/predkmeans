#####################################
## Functions for CV of predkmeans  ##
#####################################
##
## This file contains:
##
##	predkmeansCVest()
##	predkmeansCVpred()
##  print.predkmeansCVest()
##  print.predkmeansCVpred()

##
##' @name predkmeansCVest
## 
##' @title Cross-validation of Predictive K-means Clustering
##
##' @description Performs cross-validation of predictive-kmeans on a dataset.
##
## Inputs:
##' @param X Outcome data
##' @param R Covariates. Coerced to data frame.
##' @param K Number of clusters
##' @param cv.groups A list or matrix providing the cross-validation
##'		groups for splitting the data.  Alternatively, a single
##'		number giving the number of groups into which the
##'		data are randomly split. A value of '0' implies leave-one-out.
##'		Defaults to 10.
##'	@param scale Should the outcomes be re-scaled within each training
##'		group?
##'	@param covarnames Names of covariates to be included directly.
##' @param PCA Logical indicator for whether PCA components should be computed
##'		from R.
##'	@param PCAcontrol Arguments passed to \code{\link{createPCAmodelmatrix}}. This includes \code{ncomps}.
##' @param TPRS Logical indicator for whether thin-plate regression
##'		splines should be created and added to covariates.
##'	@param TPRScontrol Arguments passed to \code{\link{createTPRSmodelmatrix}}. This includes \code{df}.
##' @param sigma2 starting value of sigma2. Note: If sigma2=0 and 
##'		sigma2fixed=TRUE, then regular k-means is done in place of  
##'		predictive k-means.
##' @param sigma2fixed Logical indicating whether sigma2
##'		should be held fixed.  If FALSE, then
##'		sigma2 is estimated using Maximum Likelihood.
##' @param returnAll A list containing all \code{nStarts} solutions is
##'		included in the output.
##' @param ... Additional arguments passed to \code{\link{predkmeans}}
##
##' @details To be added....
##
##' @family 'predkmeans methods'
##' @seealso \code{\link{predkmeans}}, \code{\link{predkmeansCVpred}}, \code{\link{createPCAmodelmatrix}}, \code{\link{createTPRSmodelmatrix}}
##
##' @author Joshua Keller
##' @export
#
#		Author: J. Keller
#		Original Date: October 2015
#
predkmeansCVest <- function(X, R, K, cv.groups=10, sigma2=0,  sigma2fixed=FALSE, scale=TRUE, covarnames=NULL, PCA=FALSE, PCAcontrol=list(covarnames=colnames(R), ncomps=5), TPRS=TRUE,TPRScontrol=list(df=5, xname="x", yname="y"), returnAll=FALSE, ...){ 
	
	R <- as.data.frame(R)
	fncall <- match.call()
	if(is.null(K)){
		stop("K not provided.  Please provide.")
	}
	
	# Fill in defaults for TPRScontrol
	TPRScontrol.default <- eval(formals(predkmeansCVest)$TPRScontrol)
	for (i in names(TPRScontrol.default)){
		if (i %in% names(TPRScontrol)) next;
		TPRScontrol[[i]] <- 	TPRScontrol.default[[i]]
	}
	# Fill in defaults for PCAcontrol
	PCAcontrol.default <- eval(formals(predkmeansCVest)$PCAcontrol)
	for (i in names(TPRScontrol.default)){
		if (i %in% names(PCAcontrol)) next;
		PCAcontrol[[i]] <- 	PCAcontrol.default[[i]]
	}
		
# Checks for cv.groups.
# List, matrix, numeric.	
# Assume list for now.	
if (!is.list(cv.groups)){
	stop("cv.groups must be list. Other formats not yet implemented")
}

ids <- rownames(X)

setup <- vector("list", length(cv.groups))
pkm <- vector("list", length(cv.groups))

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
		pcafit <- do.call(createPCAmodelmatrix, args=c(list(data=R[training.set, ,drop=FALSE], matrixonly=FALSE), PCAcontrol))
		training.covars <- cbind(training.covars, pcafit$X)
		test.covars <- cbind(test.covars , predict(pcafit$pca, newdata= R[test.set,])[, 1:ncol(pcafit$X), drop=FALSE])
	}
	if (TPRS){
	    tprs.train <- do.call(createTPRSmodelmatrix, args=c(list(data=cbind(training.covars, R[training.set,]),  covarnames=colnames(training.covars),  matrixonly=FALSE), TPRScontrol))
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
		
	pres.i <- predkmeans(X=training.data, R=training.covars, K=K, sigma2=sigma2, sigma2fixed=sigma2fixed,...)
	pkm[[i]] <- pres.i
	pkm[[i]]$test.set <- test.set
	pkm[[i]]$training.set <- training.set	
	pkm[[i]]$scale.center <- scale.center
	pkm[[i]]$scale.scale <- scale.scale

	setup[[i]] <- list(training.set=training.set, test.set=test.set, scale.center=scale.center, scale.scale=scale.scale)
}
		
out <- list(call=fncall, pkm=pkm, setup=setup, PCA=PCA, PCAcontrol=PCAcontrol, TPRS=TPRS, TPRScontrol=TPRScontrol, covarnames=covarnames, X=X, R=R, K=K)
class(out) <- "predkmeansCVest"
return(out)
}	


##
##' @name predkmeansCVpred
## 
##' @title Prediction from Cross-validation of Predictive K-means Clustering
##
##' @description Performs cross-validation of predictive-kmeans on a dataset.
##
##' @param object A \code{predkmeansCVest} object. See \code{\link{predkmeansCVest}}.
##'	@param X Matrix of observations
##' @param R matrix of covariates.
##' @param method Character string indicating which prediciton method should be used. Optins are \code{ML}, \code{MixExp}, and \code{SVM}. See \code{\link{predictML}} for more information.
##' @param ... Additional arguments passed to the prediction method.

##' @seealso \code{\link{predkmeansCVest}}, \code{\link{predictML}}
##
##' @author Joshua Keller
##' @export
##
predkmeansCVpred <- function(object, X=object$X, R=object$R, method=c("ML", "MixExp", "SVM"),  ...) {

	if(!inherits(object, "predkmeansCVest")) stop("problem!")
method <- match.arg(method)

pkm <- vector("list", length(object$setup))
metrics <- vector("list", length(object$setup))

for (i in 1:length(object$setup)){
	test.set <- object$setup[[i]]$test.set
	training.set <- object$setup[[i]]$training.set

	training.data <- X[training.set, , drop=FALSE]
	test.data <- X[test.set, , drop=FALSE]
    if (!is.null(object$setup[[i]]$scale.center)){
        test.data <- scale(test.data, center=object$setup[[i]]$scale.center, scale=object$setup[[i]]$scale.scale)
    }
    
   	# Create Training and Test Covariates
	training.covars <- R[training.set, object$covarnames, drop=FALSE]
	test.covars <- R[test.set, object$covarnames, drop=FALSE]
	if (object$PCA){
		pcafit <- do.call(createPCAmodelmatrix, args=c(list(data=R[training.set, ,drop=FALSE], matrixonly=FALSE), object$PCAcontrol))
		training.covars <- cbind(training.covars, pcafit$X)
		test.covars <- cbind(test.covars , predict(pcafit$pca, newdata= R[test.set,])[, 1:ncol(pcafit$X), drop=FALSE])
	}
	if (object$TPRS){
	  	tprs.train <- do.call(createTPRSmodelmatrix, args=c(list(data=cbind(training.covars, R[training.set,]),  covarnames=colnames(training.covars),  matrixonly=FALSE), object$TPRScontrol))
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
	
	

	predobj <- object$pkm[[i]]
	pkm[[i]] <- predfn(object= predobj, Rstar=training.covars, R=test.covars,...)

	metrics[[i]] <- unlist(predictionMetrics(centers= predobj $centers, cluster.pred=pkm[[i]]$test.pred, X=test.data, labels=FALSE))
}

	out <- list(res=pkm, method=method)
	out$metrics <- simplify2array(metrics)

class(out) <- "predkmeansCVpred"
return(out)
}	


##' @title Print details for class \code{predkmeansCVest}
##' @description \code{\link[base:print]{print}} method for class \code{predkmeansCVest}.
##' @param x object of class \code{predkmeansCVest}
##' @param ... Ignored additional arguments.
##' @export
##' @family 'predkmeansCVest methods'
print.predkmeansCVest <- function(x, ...){
	if(class(x)!="predkmeansCVest"){
		stop("x must be of class predkmeansCVest.")
	}
	cat("Cross-validation fit for predictive k-means object with\n" )
	cat("    ", x$K, "Clusters\n")
	cat("    ", length(x$setup), "CV Groups\n")
	cat("Model has:\n")
	if(x$PCA) cat("    ", x$PCAcontrol$npca, "PCA components\n")
	if(x$TPRS) cat("    ", x$TPRScontrol$df, "df TPRS\n")
	invisible(x)
}##print.predkmeansCVest()


##' @title Print details for class \code{predkmeansCVpred}
##' @description \code{\link[base:print]{print}} method for class \code{predkmeansCVpred}.
##' @param x object of class \code{predkmeansCVpred}
##' @param ... Ignored additional arguments.
##' @export
##' @family 'predkmeansCVpred methods'
print.predkmeansCVpred <- function(x, ...){
	if(class(x)!="predkmeansCVpred"){
		stop("x must be of class predkmeansCVpred.")
	}
	cat("Cross-validation predictions for predictive k-means object.\n" )
	cat("Predictions computed for", ifelse(x$type=="predkmeans", "predictive", "regular"), "k-means centers,\n")
	cat("using ", x$method, ".\n", sep="")
	invisible(x)
}##print.predkmeansCVpred()


##' @title Compute summary details for class \code{predkmeansCVpred}
##' @description \code{\link[base:print]{summary}} method for class \code{predkmeansCVpred}.
##' @param object object of class \code{predkmeansCVpred}
##' @param ... Ignored additional arguments.
##' @export
##' @family predkmeans methods
summary.predkmeansCVpred <- function(object, ...){
	class(object) <- "summary.predkmeansCVpred"
	object
}##summary.predkmeansCVpred()

##' @export
print.summary.predkmeansCVpred <- function(x, ...){
	cat("Cross-validation predictions for predictive k-means object.\n" )
	cat("Predictions computed for", ifelse(x$type=="predkmeans", "predictive", "regular"), "k-means centers,\n")
	cat("using ", x$method, ".\n", sep="")
	cat("Prediction metrics are:\n")
	print(rowMeans(x$metrics))
}


