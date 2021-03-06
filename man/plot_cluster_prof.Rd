% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_cluster_prof}
\alias{plot_cluster_prof}
\title{Plot of clustered methylation profiles}
\usage{
plot_cluster_prof(bpr_cluster_obj,
  main_lab = "Clustered methylation profiles")
}
\arguments{
\item{bpr_cluster_obj}{The output of the \code{bpr_cluster_wrap} function.}

\item{main_lab}{The title of the plot}
}
\value{
The figure to be plotted in the device.
}
\description{
\code{plot_cluster_prof} creates a plot of cluster methylation profiles,
where each colour denotes a different cluster.
}
\examples{
# Cluster methylation profiles using 4 RBFs
obs <- meth_data
basis <- create_rbf_object(M = 4)
res   <- bpr_cluster_wrap(x = obs, K = 3, em_max_iter = 5, opt_itnmax = 4,
                          init_opt_itnmax = 5, is_parallel = FALSE)

# Create the plot
plot_cluster_prof(bpr_cluster_obj = res)

}
\seealso{
\code{\link{plot_scatter_gex}},
  \code{\link{plot_fitted_profiles}}, \code{\link{boxplot_cluster_gex}}
}
\author{
C.A.Kapourani \email{C.A.Kapourani@ed.ac.uk}
}
