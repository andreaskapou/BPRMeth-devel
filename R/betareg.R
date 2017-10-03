#' @name betareg_model
#' @rdname betareg_model
#' @aliases betareg
#'
#' @title Compute the Beta regression model
#'
#' @description These functions evaluate the Beta regression model likelihood
#'   and gradient.There are also functions to compute the sum of Beta regression
#'   likelihoods and weighted sum of BPR likelihoods. They are written in C++
#'   for efficiency (not yet!!).
#' @section Mathematical formula: The Beta distributed Probit Regression log
#'   likelihood function is computed by the following formula: \deqn{log p(y |
#'   x, w) = \sum_{l=1}^{L} log Beta(y_{l} | t_{l}, \Phi(w^{T}h(x_{l})))} where
#'   h(x_{l}) are the basis functions, and Beta is reparametrized to contain
#'   mean and dispersion parameters.
#'
#' @param w A vector of parameters (i.e. coefficients of the basis functions)
#' @param H The \code{L x M} matrix design matrix, where L is the number of
#'   observations and M the number of basis functions.
#' @param data An \code{L x 2} matrix containing in the 1st column are the
#'   observations, in the 2nd column are the proportions. Each row corresponds
#'   to each row of the design matrix.
#' @param x A list of elements of length N, where each element is an L x 2
#'   matrix of observations, where 1st column contains the locations. The 2nd
#'   column contains the proportions.
#' @param des_mat A list of length N, where each element contains the \code{L x
#'   M} design matrices, where L is the number of observations and M the number
#'   of basis functions.
#' @param post_prob A vector of length N containing the posterior probabilities
#'   for each element of list x, respectively.
#' @param lambda The complexity penalty coefficient for penalized regression.
#' @param is_NLL Logical, indicating if the Negative Log Likelihood should be
#'   returned.
#'
#' @return Either the Beta regression log likelihood or the gradient.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{eval_functions}}, \code{\link{betareg_optimize}}
NULL


#' @rdname betareg_model
#'
#' @export
betareg_likelihood <- function(w, H, data, lambda = 1/2, is_NLL = FALSE){
    # Extract mean and dispersion parameter
    assertthat::assert_that(is.matrix(data))
    y <- data[, 2]
    g <- data[, 3]

    y[which(y > (1 - 1e-15))] <- 1 - 1e-15
    y[which(y < 1e-15)] <- 1e-15

    # Predictions of the target variables
    # Compute the cdf of N(0,1) distribution (i.e. probit function)
    Phi <- pnorm(H %*% w)

    # In extreme cases where probit is 0 or 1, subtract a tiny number
    # so we can evaluate the log(0) when computing the Binomial
    Phi[which(Phi > (1 - 1e-15))] <- 1 - 1e-15
    Phi[which(Phi < 1e-15)] <- 1e-15

    # Compute the log likelihood
    res <- c(sum(lgamma(g) - lgamma(Phi*g) - lgamma((1-Phi)*g) + (Phi*g - 1)*log(y) +
                                ((1-Phi)*g - 1)*log(1-y)) - lambda * t(w) %*% w)

    # If we required the Negative Log Likelihood
    if (is_NLL){
        res <- (-1) * res
    }
    return(res)
}


#' @rdname betareg_model
#'
#' @export
betareg_gradient <- function(w, H, data, lambda = 1/2, is_NLL = FALSE){
    # Extract mean and dispersion parameter
    assert_that(is.matrix(data))
    y <- data[, 2]
    g <- data[, 3]

    y[which(y > (1 - 1e-15))] <- 1 - 1e-15
    y[which(y < 1e-15)] <- 1e-15

    # Predictions of the target variables
    y_hat <- as.vector(H %*% w)
    # Compute the cdf of N(0,1) distribution (i.e. probit function)
    Phi <- pnorm(y_hat)

    # In extreme cases where probit is 0 or 1, subtract a tiny number
    # so we can evaluate the log(0) when computing the Binomial
    Phi[which(Phi > (1 - 1e-15))] <- 1 - 1e-15
    Phi[which(Phi < 1e-15)] <- 1e-15

    # Compute the density of a N(0,1) distribution
    N <- dnorm(y_hat)
    N[which(N < 1e-15)] <- 1e-15

    # Compute the gradient vector w.r.t the coefficients w
    gr <- c(t(N * g * (log(y) - log(1-y) - digamma(Phi * g) + digamma((1-Phi)*g))) %*% H - 2 * lambda * w)
#
#     gr <- vector(mode = "numeric", length = length(w))
#     for (i in 1:NROW(data)){
#         gr <- gr + (N[i] * g[i] * (log(y[i]) - log(1-y[i]) - digamma(Phi[i] * g[i]) + digamma((1-Phi[i])*g[i]))) * H[i, ]
#     }
#     gr <- gr -  2 * lambda * w
    # If we required the Negative Log Likelihood
    if (is_NLL){
        gr <- (-1) * gr
    }
    return(gr)
}

# (Internal) Sum of weighted BPR log likelihoods
#
# \code{sum_weighted_bpr_lik} computes the sum of the BPR log likelihoods for
# each elements of x, and then weights them by the corresponding posterior
# probabilities. This function is mainly used for the M-Step of the EM algorithm
# \code{\link{bpr_EM}}.
#
# @param w A vector of parameters (i.e. coefficients of the basis functions)
# @param x A list of elements of length N, where each element is an L x 3 matrix
#   of observations, where 1st column contains the locations. The 2nd and 3rd
#   columns contain the total trials and number of successes at the
#   corresponding locations, repsectively.
# @param des_mat A list of length N, where each element contains the \code{L x
#   M} design matrices, where L is the number of observations and M the number
#   of basis functions.
# @param post_prob A vector of length N containing the posterior probabilities
#   fore each element of list x, respectively.
# @param lambda The complexity penalty coefficient for penalized regression.
# @param is_NLL Logical, indicating if the Negative Log Likelihood should be
#   returned.
#
# @return The weighted sum of BPR log likelihoods
#
# @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#
# @seealso \code{\link{bpr_likelihood}}, \code{\link{bpr_EM}}
#

#' @rdname betareg_model
#'
#' @export
sum_weighted_betareg_lik <- function(w, x, des_mat, post_prob, lambda = 1/2,
                                     is_NLL = TRUE){
    N <- length(x)

    # TODO: Create tests
    # For each element in x, evaluate the BPR log likelihood
    res <- vapply(X   = 1:N,
                  FUN = function(y) betareg_likelihood(w = w,
                                                       H = des_mat[[y]],
                                                       data = x[[y]],
                                                       lambda = lambda,
                                                       is_NLL = is_NLL),
                  FUN.VALUE = numeric(1),
                  USE.NAMES = FALSE)

    # Return the dot product of the result and the posterior probabilities
    return(post_prob %*% res)
}


# (Internal) Sum of weighted gradients of the BPR log likelihood
#
# \code{sum_weighted_bpr_grad} computes the sum of the gradients of BPR log
# likelihood for each elements of x, and then weights them by the corresponding
# posterior probabilities. This function is mainly used for the M-Step of the EM
# algorithm \code{\link{bpr_EM}}.
#
# @inheritParams sum_weighted_bpr_lik
#
# @return A vector with weighted sum of the gradients of BPR log likelihood.
#
# @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#
# @seealso \code{\link{bpr_gradient}}, \code{\link{bpr_EM}}
#

#' @rdname betareg_model
#'
#' @export
sum_weighted_betareg_grad <- function(w, x, des_mat, post_prob, lambda = 1/2,
                                      is_NLL = TRUE){
    N <- length(x)

    # TODO: Create tests
    # For each element in x, evaluate the gradient of the BPR log likelihood
    res <- vapply(X   = 1:N,
                  FUN = function(y) betareg_gradient(w = w,
                                                     H = des_mat[[y]],
                                                     data = x[[y]],
                                                     lambda = lambda,
                                                     is_NLL = is_NLL),
                  FUN.VALUE = numeric(length(w)),
                  USE.NAMES = FALSE)

    # Return the dot product of the result and the posterior probabilities
    return(post_prob %*% t(res))
}
