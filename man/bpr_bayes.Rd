% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bpr_bayes.R
\name{bpr_bayes}
\alias{bpr_bayes}
\alias{bpr_bayes}
\alias{bpr_bayes.list}
\alias{bpr_bayes.matrix}
\title{Gibbs sampling approach for the BPR model}
\usage{
bpr_bayes(x, ...)

\method{bpr_bayes}{list}(x, w_mle = NULL, basis = NULL,
  fit_feature = NULL, cpg_dens_feat = FALSE, w_0_mean = NULL,
  w_0_cov = NULL, gibbs_nsim = 20, gibbs_burn_in = 10,
  keep_gibbs_draws = FALSE, is_parallel = TRUE, no_cores = NULL, ...)

\method{bpr_bayes}{matrix}(x, w_mle = NULL, basis = NULL,
  fit_feature = NULL, cpg_dens_feat = FALSE, w_0_mean = NULL,
  w_0_cov = NULL, gibbs_nsim = 20, gibbs_burn_in = 10,
  keep_gibbs_draws = FALSE, ...)
}
\arguments{
\item{x}{The input object, either a \code{\link[base]{matrix}} or a
\code{\link[base]{list}}.}

\item{...}{Additional parameters.}

\item{w_mle}{A vector of parameters (i.e. coefficients of the basis
functions) containing the MLE estimates.}

\item{basis}{A 'basis' object. E.g. see \code{\link{create_rbf_object}}.}

\item{fit_feature}{Return additional feature on how well the profile fits the
methylation data. Either NULL for ignoring this feature or one of the
following: 1) "RMSE" for returning the fit of the profile using the RMSE as
measure of error or 2) "NLL" for returning the fit of the profile using the
Negative Log Likelihood as measure of error.}

\item{cpg_dens_feat}{Logical, whether to return an additional feature for the
CpG density across the promoter region.}

\item{w_0_mean}{The prior mean hyperparameter for w}

\item{w_0_cov}{The prior covariance hyperparameter for w}

\item{gibbs_nsim}{Optional argument giving the number of simulations of the
Gibbs sampler.}

\item{gibbs_burn_in}{Optional argument giving the burn in period of the Gibbs
sampler.}

\item{keep_gibbs_draws}{Logical indicating if we should keep the whole MCMC
chain for further analysis.}

\item{is_parallel}{Logical, indicating if code should be run in parallel.}

\item{no_cores}{Number of cores to be used, default is max_no_cores - 2.}
}
\value{
Depending on the input object \code{x}: \itemize{\item{If \code{x} is
  a \code{\link[base]{list}}:}  An object containing the following elements:
  \itemize{ \item{ \code{W_opt}: An Nx(M+1) matrix with the optimized
  parameter values. Each row of the matrix corresponds to each element of the
  list x. The columns are of the same length as the parameter vector w (i.e.
  number of basis functions). } \item{ \code{Mus}: An N x M matrix with the
  RBF centers if basis object is \code{\link{create_rbf_object}}, otherwise
  NULL.} \item{ \code{basis}: The basis object. } \item{ \code{w}: The
  initial values of the parameters w. } } \item{If \code{x} is a
  \code{\link[base]{matrix}}:} An object containing the following elements:
  \itemize{ \item{ \code{w_opt}: Optimized values for the coefficient vector
  w. The length of the result is the same as the length of the vector w. }
  \item{ \code{basis}: The basis object. } }}
}
\description{
The function \code{bpr_bayes} computes the posterior of the BPR
  model using auxiliary variable approach. Since we cannot compute the
  posterior analytically, a Gibbs sampling scheme is used.
}
\examples{
data <- meth_data
out_opt <- bpr_bayes(x = data, is_parallel = FALSE)

ex_data <- meth_data
basis <- create_rbf_object(M=3)
out_opt <- bpr_bayes(x = ex_data, is_parallel = FALSE, basis = basis)

basis <- create_rbf_object(M=2)
w <- c(0.1, 0.1, 0.1)
w_0_mean <- rep(0, length(w))
w_0_cov <- diag(10, length(w))
data <- meth_data[[1]]
out_opt <- bpr_bayes(x = data, w_mle = w, w_0_mean = w_0_mean,
                     w_0_cov = w_0_cov, basis = basis)

basis <- create_rbf_object(M=0)
w <- c(0.1)
w_0_mean <- rep(0, length(w))
w_0_cov <- diag(10, length(w))
data <- meth_data[[1]]
out_opt <- bpr_bayes(x = data, w_mle = w, w_0_mean = w_0_mean,
                     w_0_cov = w_0_cov, basis = basis)

}
\seealso{
\code{\link{create_basis}}, \code{\link{eval_functions}}
}
\author{
C.A.Kapourani \email{C.A.Kapourani@ed.ac.uk}
}
