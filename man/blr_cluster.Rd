% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blr_cluster.R
\name{blr_cluster}
\alias{blr_cluster}
\title{Cluster methylation profiles with Gaussian noise}
\usage{
blr_cluster(x, K = 3, pi_k = NULL, w = NULL, basis = NULL, s2 = NULL,
  em_max_iter = 100, epsilon_conv = 1e-04, lambda = 1/10,
  is_verbose = FALSE)
}
\arguments{
\item{x}{The Gaussian distributed observations, which has to be a list of
elements of length N, where each element is an L x 2 matrix of
observations, where 1st column contains the locations and the 2nd column
contains the methylation levels.}

\item{K}{Integer denoting the number of clusters K.}

\item{pi_k}{Vector of length K, denoting the mixing proportions.}

\item{w}{A MxK matrix, where each column consists of the basis function
coefficients for each corresponding cluster.}

\item{basis}{A 'basis' object. E.g. see \code{\link{create_rbf_object}}.}

\item{s2}{Vector of initial linear regression variances for each cluster.}

\item{em_max_iter}{Integer denoting the maximum number of EM iterations.}

\item{epsilon_conv}{Numeric denoting the convergence parameter for EM.}

\item{lambda}{The complexity penalty coefficient for ridge regression.}

\item{is_verbose}{Logical, print results during EM iterations.}
}
\value{
A 'blr_cluster' object which, in addition to the input
  parameters, consists of the following variables: \itemize{
  \item{\code{pi_k}: Fitted mixing proportions.} \item{\code{w}: A MxK matrix
  with the fitted coefficients of the basis functions for each cluster k.}
  \item{\code{NLL}: The Negative Log Likelihood after the EM algorithm has
  finished.} \item{\code{post_prob}: Posterior probabilities of each promoter
  region belonging to each cluster.} \item{\code{labels}: Hard clustering
  assignments of each observation/promoter region.} \item{\code{BIC}:
  Bayesian Information Criterion metric.} \item{\code{AIC}: Akaike
  Information Criterion metric.} \item{\code{ICL}: Integrated Complete
  Likelihood criterion metric.} }
}
\description{
\code{blr_cluster} is a wrapper function that clusters methylation
profiles using the EM algorithm. Initially, it performs parameter checking,
and initializes main parameters, such as mixing proportions, basis function
coefficients, then the EM algorithm is applied and finally model selection
metrics are calculated, such as BIC and AIC.
}
\examples{
my_clust <- blr_cluster(x = lm_data, em_max_iter = 100, is_verbose = TRUE)

}
\author{
C.A.Kapourani \email{C.A.Kapourani@ed.ac.uk}
}
