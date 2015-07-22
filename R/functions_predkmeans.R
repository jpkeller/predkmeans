###############################
## Functions for predkmeans  ##
###############################
##
## This file contains:
##
##	predkmeans()
##	print.predkmeans()
##	summary.predkmeans()
##	print.summary.predkmeans()
##	getExpMahal()
##  getH()
##	getMu()
##	getSigma2()


##
##' @name predkmeans
## 
##' @title Predictive K-means Clustering
##
##' @description Uses a Mixture-of-experts algorithm to find 
##' cluster centers that are influenced by prediction covariates
##
## Inputs:
##' @param X Outcome data
##' @param R Covariates
##' @param K Number of clusters
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
##' @param maxitMlogit Maximum number of iterations in the
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
##' @references Jordan, M. and R. Jacobs (1994). Hierarchical mixtures of
##'		experts and the EM algorithm. \emph{Neural computation 6}(2),
##'		 181-214.
##
##' @author Joshua Keller
##
##'	@return To be added......
#		Author: J. Keller
#		Original Date: February 2014
#		Last modified: June 2015
# 
#	Changelog:
#		16 June 2015 -- Changed to update new function names as part
#						of creating package.
#		15 July 2015 -- Changing to account for the fact that
#						mlogit now beta has a column of zeros
#						......still need to do.
#		17 July 2015 -- mixExp
#		20 July 2015 -- changed function name to predkmeans
predkmeans <- function(X, R, K=3, mu=NULL, muStart=c("kmeans","random"), sigma2=0,  sigma2fixed=FALSE,maxitEM=100, tol=1e-5, convEM=c("both", "mu", "gamma"), nStarts=1, maxitMlogit=500,verbose=0, muRestart=1000, returnAll=FALSE, ...) {
	
	# Check that arguments are correct
	muStart <- match.arg(muStart)
	convEM <- match.arg(convEM)
	
	# Check dimensions
	n <- nrow(X)
	if (n!=nrow(R)){
		stop("Number of outcome observations does not match number of covariate observations.")
	}
		
	# If sigma2 is 0 and fixed, then return kmeans result
	if (sigma2==0 && sigma2fixed){
		
			res.best <- kmeans(x=X, centers=K, nstart=nStarts)
	centers <- res.best$centers
	tr.assign <- res.best$cluster
	wSS <- res.best$tot.withinss
	out <- list(res.best= res.best, centers=centers, cluster=tr.assign, K=K, sigma2=0, wSS=wSS, sigma2fixed= sigma2fixed)
		class(out) <- "predkmeans"
	return(out)
	}
	
	
	d <- ncol(X) # Dimension of outcome
	if (is.null(colnames(X))) colnames(X) <- paste0("X", 1:d)
	p <- ncol(R)
	mu_orig <- mu # Store original value of mu
	results <- vector("list", nStarts)
	for (i in 1:nStarts){
		if (verbose>0) cat("Fit", i, "of", nStarts, "\n")
		
		# Initialize mu
		mu <- mu_orig
		if (is.null(mu)){
			if (muStart=="random"){
				all_same <- TRUE
				all_same_iter <- 1
				while(all_same){
					all_same_iter <- all_same_iter + 1
					if (all_same_iter > muRestart){
						stop(paste0("Picked ", muRestart, " random sets of starting centers and all were equal.  Check whether data is identical.  If necessary, change seed or increase  'all_same_iter_limit'."))
					}
					mu <- X[sample(1:nrow(X), size=K), ]
		 			all_same <- all(apply(mu, 2, function(x) diff(range(x))<.Machine$double.eps ^ 0.5))
				}
			} else if (muStart=="kmeans"){
				mu <- kmeans(X, centers=K)$centers
			}
			rownames(mu) <- 1:K	
		} # if (is.mull(mu))

		# Compute initial values for h (Step E-1)
		h0 <- assignCluster(X, centers=mu)
		h0 <- model.matrix(~0 + factor(h0, levels=1:K))
		colnames(h0) <- paste0("C", 1:K)
		
		# Compute mu, gamma, and sigma2 (Step M-1)
		mu <- getMu(X, h0)
		if(!sigma2fixed){
			sigma2 <- getSigma2(X, mu, h0)
		}
		mfit0 <- mlogit(h0, R, add.intercept=F, verbose=verbose>3, suppressFittedWarning=TRUE, iterlim=maxitMlogit)
		gamma <- mfit0$beta # cbind(C1=0, mfit0$beta)
		
		converged <- 0
		iter <- 0
		convCheck <- .Machine$integer.max	
		if (verbose>2) cat("Mixture of Experts code running with\n", K, " experts each finding a mean in ", d, " variables.\nThe gating networks each have ", p, " parameters\n", sep="")							
		# EM algorithm
		while (!converged){
			iter <- iter + 1
			if (iter> maxitEM) {
				if (verbose>1) cat("Stopping: Reached maximum iterations\n")
				converged <- 9
				break
			}
			
			# E-step
			h <- getH(X, R, gamma, mu, sigma2=sigma2)
			# M-step, Part 1a: mu
			mu_new <- getMu(X, h)
			# M-step, Part 1b: sigma2
			if(!sigma2fixed){
				sigma2 <- getSigma2(X, mu_new, h)
			}
			# M-step, Part 2: gamma
			mfit <- mlogit(Y=h, X=R, beta=as.vector(gamma[,-1]), add.intercept=F, iterlim= maxitMlogit, verbose=verbose>3, checkY=FALSE, ...)	
			# Quesiton: Should I stop freeze the parameters for a certain cluster, if 
			# I have achieved all 0/1's?
			if (!mfit$res.best$conv) {
				if (verbose>1) message("EM-algorithm: mlogit optimization error in updating gamma.")
				converged <- 8
				break
			} else {
				gamma_new <- mfit$beta # cbind(0, mfit$beta)
			}
			# Check for convergence		
			convCheck <- switch(convEM, mu=norm(mu_new - mu, type="F") <tol,
				gamma= norm(gamma_new - gamma, type="F") <tol,
				both=norm(mu_new - mu, type="F") <tol && norm(gamma_new - gamma, type="F") < tol)
	
			 # Force at least three iterations?
			if(convCheck && iter>2) {  
				converged <- 1
				convCheck <- norm(mu_new - mu, type="F")
			}
			gamma <- gamma_new
			mu <- mu_new
		
		}	## while (!converged)								

		# Compute the objective function value
		#Uproxy <- exp(R %*% gamma)/rowSums(exp(R %*% gamma))
		Uproxy <- R %*% gamma
		Uproxy <- apply(Uproxy, 2, function(x) 1/rowSums(exp(Uproxy-x)))
		
		# Objective is the log likelihood
		obj <- 0
		for (kk in 1:nrow(X)){
			q <- sweep(mu, 2, X[kk,])
			obj <- obj + log(sum(sigma2^(-d/2)*exp(-1/(2*sigma2)*rowSums(q*q)) * Uproxy[kk,]))
		}				
				
		results[[i]] <- list(mu=mu, gamma=gamma, sigma2=sigma2, conv=converged, objective= obj, iter=iter, mfit=mfit[c("beta", "fitted01", "fitted", "res.best", "status")])			
			
		if (verbose>1) {
			cat("EM Algorithm completed\n")
			cat("Convergence code:", converged, "\n")
			cat("Number of EM iterations:", iter, "\n")
			cat("Norm of last mu update:", convCheck, "\n")
			cat("Value of objective:", obj, "\n")
			cat("Sigma^2:", sigma2, "\n\n")
		}
	}

	best.fit <- which.max(lapply(results, function(x) x$objective))	
	# if(!results[[best.fit]]$mfit$status$conv){
		# warning("Best Fit has mlogit optimization error.")
	# }
	res.best <- results[[best.fit]]
	centers <- res.best $mu
	tr.assign <- assignCluster(X, centers= centers)
	wSS <- sum((X -  centers[tr.assign,])^2)
	out <- list(res.best= res.best, centers=centers, cluster=tr.assign, K=K, sigma2= res.best$sigma2, wSS=wSS, sigma2fixed= sigma2fixed, h=h, gamma=gamma)

	if (returnAll){
		out <- c(out, res.all=list(results), best.fit=best.fit)
	}
	class(out) <- "predkmeans"
	
	out
}## predkmeans()




###############################
## S3Methods for predkmeans  ##
###############################

##' @title Print details for class \code{predkmeans}
##' @description \code{\link[base:print]{print}} method for class \code{predkmeans}.
##' @param x object of class \code{predkmeans}
##' @param ... Ignored additional arguments.
##' @export
##' @family predkmeans methods
print.predkmeans <- function(x, ...){
	if(class(x)!="predkmeans"){
		stop("x must be of class predkmeans.")
	}
	cat("Predictive k-means object with\n" )
	cat("    ", x$K, "Clusters\n")
	cat("    ", ncol(x$centers), "Variables\n")
#	cat("Convergence status: ", x$res.best$conv)
	invisible(x)
}##print.predkmeans()


##' @title Compute summary details for class \code{predkmeans}
##' @description \code{\link[base:print]{summary}} method for class \code{predkmeans}.
##' @param object object of class \code{predkmeans}
##' @param ... Ignored additional arguments.
##' @export
##' @family predkmeans methods
summary.predkmeans <- function(object, ...){
	class(object) <- "summary.predkmeans"
	object
}##summary.predkmeans()

##' @export
print.summary.predkmeans <- function(x, ...){
	cat("Predictive k-means object with\n" )
	cat("    ", x$K, "Clusters\n")
	cat("    ", ncol(x$centers), "Variables\n")
	cat("Convergence status: ", ifelse(class(x$res.best)=="kmeans", 1, x$res.best$conv), "\n")
	cat("Sigma^2 = ", round(x$sigma2, 2), " (Fixed = ", x$sigma2fixed, ")\n", sep="")
	cat("Within-cluster Sum-of-Squares (wSS) = ", round(x$wSS, 2), "\n")
	cat("Cluster centers are:\n")
	print(x$centers)
}





###############################################
# For each row Xi of X, this computes 
#      exp(-1/(2*sigma2) * (Xi - mu)^T(Xi - mu))
# This is computed for the matrix X at once, rather than row by row
getExpMahal <- function(X, mu, sigma2) {
	q <- sweep(X, 2, mu)
	exp(-1/(2*sigma2)*rowSums(q*q))
}

###############################################
# Computes h, the expected value of z (the latent variable
# indicating cluster membership) conditional on the outcomes
# X and the current parameter values of gamma and mu
#
# Input:
#      X -- Matrix of outcome (e.g., AQS) data.  This should have n rows of
#             observations and d columns of observed variables (e.g., pollutants)
#      R -- Matrix of covariate information (e.g., GIS variables).  This should
#             have n rows of observations and p number of covariate columns.
#             This matrix is taken as-is, so if an intercept is desired it should already
#             be included as a column of R.
#      gamma -- Matrix of parameter estimates for the multinomial logistic
#                        regression of the geographic covariates.  Should contain
#                        p rows of parameters for each of K columns (corresponding
#                        to the cluster centers
#      mu -- Matrix of cluster centers.  Should have d columns for each of 
#                K rows
#             
getH <- function(X, R, gamma, mu, sigma2=sigma2){
	Px <- apply(mu, 1, function(m) getExpMahal(X, m, sigma2=sigma2))
#	Rg <- exp(R %*% gamma)/rowSums(exp(R %*% gamma))
	Rg <- R %*% gamma	
	Rg <- apply(Rg, 2, function(x) 1/rowSums(exp(Rg-x)))
	num <- Px*Rg
	h <- num/rowSums(num)
	h
}

##############################################
# Computes the weighted centers of clusters for the EM
# implementation of the Mixture of Experts problem.
getMu <- function(X, h){
	mu <- apply(h, 2, function(hk) colSums(hk*X)/sum(hk))
	t(mu)
}

getSigma2 <- function(X, mu, h){
	K <- ncol(h)
	p <- ncol(mu)
	mahal <- apply(mu, 1, function(m) rowSums(sweep(X, 2, m)^2))
	num <- sum(h*mahal)
	denom <- p*sum(h)
	num/denom
}

