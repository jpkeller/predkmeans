##################################################
## Functions for predicting cluster membership  ##
##################################################
##
## This file contains:
##
##	predictML.predkmeans()
##	predictSVM.predkmeans()
##
##
##


##' @name predictML.predkmeans
##' @aliases predictML.predkmeans predictML 
##' @title Prediction of Cluster Membership
##
##' @description Predicts cluster membership using either multinomial logistic
#'    regression or SVMs.
# 
#' @details Function for predicting cluster membership in clusters identified by k-means or predictive k-means using multinomial logistic regression or support vector machines (SVMs).  For multinomial logitic regression, parameter estimation is handled by \code{mlogit}.  For SVMs, the svm is fit using \code{best.svm} in \code{e1071} package.
#
# INPUT:
#' @param   object A predkmeans object, from which the cluster centers will be extracted.
##' @param 	centers	Matrix of cluster centers, assumed to be K-by-p
##' @param 	K Number of clusters
##' @param	R 	matrix of covariates for observations to be predicted at.  
##' @param	Rstar matrix of covariates at training locations
#' @param	X 	Optional matrix of observations at prediction locations. If provided, measures of prediction performance will be reported.
#' @param Xstar matrix of observation at training locations.  Either this or \code{tr.assign} is required.	
#' @param  tr.assign   vector of cluster assignments at training locations
#' @param muStart starting value for cluster centers in mlogit optimization (IDEA: change to pull from predkmeans object?).  If not provided, starting values are selected randomly. 
##' @param maxitMlogit Maximum number of iterations for \code{mlogit} in prediction
##' @param verbose integer indicating amount of output to be displayed
##' @param nMlogitStarts number of mlogit starts to use in estimation of parameters
##' @param mlogit.control list of control parameters to be passes to \code{mlogit}
#' @seealso \code{\link{mlogit}}, \code{\link{predkmeans}}
#' @author Joshua Keller
#
# Output
#' @return A list containing some or all of the following elements:
#'	\item{tr.assign}{Cluster assignments at training locations}
#' 	\item{mlfit}{A subset of the mlogit object returned by the function of that name}
#' \item{beta}{beta}
#' \item{test.pred}{test.pred}
#' \item{test.assign}{test.assign}
#' \item{pred.acc}{pred.acc}
##' \item{MSEP}{MSEP}
##' \item{wSS.subj}{wSS.subj}
##' @export
predictML.predkmeans <- function(object=NULL, centers=object$centers, K=nrow(centers), R,  Rstar, X=NULL, Xstar=NULL, tr.assign=object$cluster, muStart ="random", maxitMlogit =500, verbose=1, nMlogitStarts=1,  mlogit.control=list(suppressFittedWarning=TRUE)){

if (!is.null(object)){
	if(!inherits(object, c("predkmeans"))){
		stop("Input 'object' must be of class predkmeans.")
	}
}
if (K > nrow(centers)){
	stop("K must be less than or equal to number of centers provided")
}
if (K < nrow(centers)){
	warning("K is less than number of centers provided. Using first K centers.")
	centers <- centers[1:K, ]
}
if (is.null(tr.assign)) {
	if (is.null(Xstar)){
		stop("Requires either assigned clusters 'tr.assign' or data 'Xstar' for training observations")
	} else {
		tr.assign <- assignCluster(Xstar, centers= centers)
	}
}

mm.tr.assign <- model.matrix(~0 + factor(tr.assign, levels=1:K)) 
colnames(mm.tr.assign) <- paste0("C", 1:K)
if (muStart=="random"){
	beta.init <- matrix(rnorm((K-1)*ncol(Rstar)* nMlogitStarts), nrow=(K-1)*ncol(Rstar))
} else {
	beta.init <- muStart
}
mlfit <- do.call(mlogit, args=c(list(Y= mm.tr.assign, X=Rstar, beta=beta.init, iterlim=maxitMlogit, verbose=verbose), mlogit.control))
beta <- mlfit$beta
test.pred <- get.cluster.pred(R= R, beta=beta)

if (!is.null(X) && nrow(X)>0){
	test.assign <- assignCluster(X, centers= centers)
	pred.acc <- mean(test.pred == test.assign)
	MSEP <- sum((centers[test.pred,] -  centers[test.assign,])^2)/nrow(X)
	wSS.subj <- sum((X - centers[test.assign,])^2)/nrow(X)
} else {
	test.assign <- NULL
	pred.acc <- NULL
	pred.SS <- NULL
	MSEP <- NULL
	wSS.subj <- NULL
}

out <- list(tr.assign= tr.assign, mlfit=mlfit[c("beta", "fitted", "res.best", "status")], beta=beta, test.pred = test.pred , test.assign=test.assign, pred.acc=pred.acc, MSEP=MSEP, wSS.subj=wSS.subj)

return(out)
}


##' @rdname predictML.predkmeans
##' @aliases predictSVM.predkmeans predictSVM 
##
##' @param svm.control list of options for \code{best.svm} 
##
##' @export
# Function for predicting cluster membership in clusters
# identified by k-means or predictive k-means
# using SVMs
#
# INPUT:
#	centers -- Matrix of cluster centers, assumed to be K-by-p
#	R -- matrix of covariates for observations to be predicted at.  
#	X -- matrix of observations at prediction locations (optional -- include if 
#			prediction metrics are desired)
#	
#	Rstar -- matrix of covariates at training locations
#   tr.assign -- Assignments of training data
#   
predictSVM.predkmeans <- function(object=NULL,
centers=object$centers, R,  Rstar, K=nrow(centers), X=NULL, Xstar=NULL, tr.assign =object$cluster,  svm.control=list(gamma=c(1/(2:1), 2), cost=seq(20, 100, by=20))){

if(!requireNamespace("e1071", quietly=TRUE)){
	stop("e1071 is required for SVM prediction.  Please install it.", call.=FALSE)
}


if (!is.null(object)){
	if(!inherits(object, c("predkmeans"))){
		stop("Input 'object' must be of class predkmeans.")
	}
}
if (K > nrow(centers)){
	stop("K must be less than or equal to number of centers provided")
}
if (K < nrow(centers)){
	warning("K is less than number of centers provided. Using first K centers.")
	centers <- centers[1:K, ]
}
if (is.null(tr.assign)) {
	if (is.null(Xstar)){
		stop("Requires either assigned clusters 'tr.assign' or data 'Xstar' for training observations")
	} else {
		tr.assign <- assignCluster(Xstar, centers= centers)
	}
}

if (!all(colnames(R) %in% colnames(Rstar))) {
	stop("Covariate mismatch")
}

	if ("(Intercept)" %in% colnames(Rstar)) {
		Rstar <- Rstar[,-which(colnames(Rstar)=="(Intercept)")]
		R <- R[,-which(colnames(R)=="(Intercept)")]
	}

	svm.model <- do.call(e1071::best.svm, args=c(list(x= Rstar, y=factor(tr.assign, levels=1:K)), svm.control))
	# Save space, this is a large call, since it explicitly includes the entire variable structure/objects
	svm.model$call <- NULL 
	test.pred <- predict(svm.model, newdata= R)

if (!is.null(X) && nrow(X)>0){
	test.assign <- assignCluster(X, centers= centers)
	pred.acc <- mean(test.pred == test.assign)
	MSEP <- sum((centers[test.pred,] -  centers[test.assign,])^2)/nrow(X)
	wSS.subj <- sum((X - centers[test.assign,])^2)/nrow(X)
} else {
	test.assign <- NULL
	pred.acc <- NULL
	pred.SS <- NULL
	MSEP <- NULL
	wSS.subj <- NULL
}

out <- list(tr.assign= tr.assign, svm.model=svm.model, test.pred = test.pred , test.assign=test.assign, pred.acc=pred.acc, MSEP=MSEP, wSS.subj=wSS.subj)

return(out)
}
