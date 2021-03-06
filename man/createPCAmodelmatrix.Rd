% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_aux.R
\name{createPCAmodelmatrix}
\alias{createPCAmodelmatrix}
\title{Create Principal Component Analysis (PCA) scores matrix}
\usage{
createPCAmodelmatrix(
  data,
  ncomps,
  covarnames = colnames(data),
  center = TRUE,
  scale = TRUE,
  matrixonly = TRUE
)
}
\arguments{
\item{data}{Matrix or data frame of data}

\item{ncomps}{Number of PCA components to return.}

\item{covarnames}{Names of variables or column numbers in \code{data}
on which the PCA is to be run.}

\item{center}{Logical indicator of whether \code{data} should be centered. Passed to \code{\link{prcomp}}.}

\item{scale}{Logical indicator of whether \code{data} should be scaled. Passed to \code{\link{prcomp}}.}

\item{matrixonly}{Logical indicator of whether only the model matrix should
be returned, or the full output from \code{\link{prcomp}}.}
}
\value{
If \code{matrixonly=TRUE}, a matrix of PCA scores. Otherwise a list containing two elements: \code{X},  a matrix of scores, and \code{pca}, the output from \code{prcomp}.
}
\description{
Wrapper function for creating PCA scores to be used
	in a regression analysis.
}
\details{
This is a wrapper around \code{\link{prcomp}}, which does 
		the necessary computation.
}
\examples{
n <- 100
d <- 15
X <- matrix(rnorm(n*d), ncol=d, nrow=n)
X <- as.data.frame(X)
mx <- createPCAmodelmatrix(data=X, ncomps=2)
}
\seealso{
\code{\link{createTPRSmodelmatrix}}, \code{\link{predkmeansCVest}}
}
\author{
Joshua Keller
}
