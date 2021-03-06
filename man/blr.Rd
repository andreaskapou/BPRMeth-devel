% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blr.R
\name{blr}
\alias{blr}
\title{Fitting linear models using Basis Functions}
\usage{
blr(x, y, basis, lambda = 0, return.all = TRUE)
}
\arguments{
\item{x}{The observations.}

\item{y}{The response.}

\item{basis}{Basis function object e.g. \code{\link{create_rbf_object}}.}

\item{lambda}{Optional parameter for performing ridge regression.}

\item{return.all}{Optional logical, indicating if all the metrics should be
computed (mainly for efficiency).}
}
\value{
An object of class "basis_lm" is a list containing the following
 components:
\itemize{
 \item \code{coefficients} a named vector of coefficients.
 \item \code{residuals} the residuals, that is response minus fitted values.
 \item \code{fitted.values} the fitted mean values.
 \item \code{df.residuals} the residual degrees of freedom.
 \item \code{sigma} the standard deviation of the residuals.
 \item \code{vcov} The covariance matrix.
 \item \code{basis} The basis object used.
 \item \code{lambda} The regularization parameter.
 \item \code{call} the matched call.
}
}
\description{
\code{blr} is used to fit linear models using basis functions such as
Radial Basis Functions (RBFs), Fourier and Polynomial basis functions.
}
\seealso{
\code{\link{predict.blr}}, \code{\link{create_fourier_object}},
 \code{\link{create_rbf_object}}, \code{\link{summary.blr}}
}
\author{
C.A.Kapourani \email{C.A.Kapourani@ed.ac.uk}
}
