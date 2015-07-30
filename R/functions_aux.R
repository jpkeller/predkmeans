
get.cluster.pred <- function(R, beta){
	eta.pred <- R %*% beta
	mu.pred <- exp(eta.pred)/rowSums(exp(eta.pred))
	cluster.pred <- apply(mu.pred, 1, which.max)
	return(cluster.pred)
}


##' @name assignCluster
##' @title Make Cluster Assignments
##' @description Assigns observation to the nearest cluster center, using squared Euclidean distance.
##
##'	@param X matrix of observations
##' @param centers matrix of cluster centers
##' @return A vector of cluster labels
##' @author Joshua Keller
##' @export
assignCluster <- function(X, centers){
	k <- nrow(centers)
    cl_dist <- matrix(NA, nrow=dim(X)[1], ncol=k)
    ones <- rep(1, dim(X)[2])
    for (cl in 1:k){		
        C <- matrix(data=unlist(rep(centers[cl,], each=dim(X)[1])), nrow=dim(X)[1], byrow=F, dimnames=dimnames(X))
        cl_dist[,cl] <- (X-C)^2 %*% ones
    }
    labels <- apply(cl_dist, 1, which.min)
    names(labels) <- rownames(X)
    return(labels)
}





# @title Get PCA comps for regression
# @export
get_PCA_matrix <- function(data, ncomps, covarnames=colnames(data), center=TRUE, scale=TRUE, matrixonly=TRUE){
	pca <- prcomp(data[,covarnames, drop=FALSE], center=center, scale=scale)
	X <- pca$x[,1:ncomps, drop=FALSE]
	if (matrixonly) {
		return(X)
	} else {
		return(list(X=X, pca=pca))
	}
}


# @title Get matrix of TPRS for regression
# @export
get_TPRS_modelmatrix <- function(data, TPRSk =5, covarnames=NULL, xname="lambert_x", yname="lambert_y", matrixonly=TRUE, TPRSfx=TRUE){
if(!requireNamespace("mgcv", quietly=TRUE)){
	stop("mgcv is required to create TPRS objects.  Please install it.", call.=FALSE)
}	
	
	# Create the formula
	f <- ""
	if (length(covarnames)>0){
		f <- paste("+ ", paste0(covarnames, collapse=" + "))
	}
	f <- formula(paste0(xname,"~s(", xname,", ", yname,", fx=TPRSfx, k=TPRSk)", f))
	# Fit a GAM to get the model matrix
	gamfit <- mgcv::gam(f, data=as.data.frame(data))
	X  <- model.matrix(gamfit)

	if(matrixonly){
		return(X)
	} else {
		return(list(X=X, gamfit=gamfit))
	}
}


# Function for creating CV groups
create_cv_groups <- function(mons, folds=length(mons)){
    nmons <- length(mons)
    folds <- min(folds, nmons)
    cv.groups <- trunc(folds * 1:nmons / (nmons+1)) + 1
    cv.groups <- cv.groups[order(runif(nmons, 0, 1))]
    cv.list <- list()
    for (i in 1:folds)
    {
        cv.list[[i]] <- mons[cv.groups==i]
    }
    return(cv.list)
}
