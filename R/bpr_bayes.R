#' @name bpr_bayes
#' @rdname bpr_bayes
#' @aliases bpr_bayes
#'
#' @title Gibbs sampling approach for the BPR model
#'
#' @description The function \code{bpr_bayes} computes the posterior of the BPR
#'   model using auxiliary variable approach. Since we cannot compute the
#'   posterior analytically, a Gibbs sampling scheme is used.
#'
#' @param x The input object, either a \code{\link[base]{matrix}} or a
#'   \code{\link[base]{list}}.
#' @param ... Additional parameters.
#' @param w_mle A vector of parameters (i.e. coefficients of the basis
#'   functions) containing the MLE estimates.
#' @param basis A 'basis' object. E.g. see \code{\link{create_rbf_object}}.
#' @param fit_feature Return additional feature on how well the profile fits the
#'   methylation data. Either NULL for ignoring this feature or one of the
#'   following: 1) "RMSE" for returning the fit of the profile using the RMSE as
#'   measure of error or 2) "NLL" for returning the fit of the profile using the
#'   Negative Log Likelihood as measure of error.
#' @param cpg_dens_feat Logical, whether to return an additional feature for the
#'   CpG density across the promoter region.
#' @param w_0_mean The prior mean hyperparameter for w
#' @param w_0_cov The prior covariance hyperparameter for w
#' @param gibbs_nsim Optional argument giving the number of simulations of the
#'   Gibbs sampler.
#' @param gibbs_burn_in Optional argument giving the burn in period of the Gibbs
#'   sampler.
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores Number of cores to be used, default is max_no_cores - 2.
#'
#' @return Depending on the input object \code{x}: \itemize{\item{If \code{x} is
#'   a \code{\link[base]{list}}:}  An object containing the following elements:
#'   \itemize{ \item{ \code{W_opt}: An Nx(M+1) matrix with the optimized
#'   parameter values. Each row of the matrix corresponds to each element of the
#'   list x. The columns are of the same length as the parameter vector w (i.e.
#'   number of basis functions). } \item{ \code{Mus}: An N x M matrix with the
#'   RBF centers if basis object is \code{\link{create_rbf_object}}, otherwise
#'   NULL.} \item{ \code{basis}: The basis object. } \item{ \code{w}: The
#'   initial values of the parameters w. } } \item{If \code{x} is a
#'   \code{\link[base]{matrix}}:} An object containing the following elements:
#'   \itemize{ \item{ \code{w_opt}: Optimized values for the coefficient vector
#'   w. The length of the result is the same as the length of the vector w. }
#'   \item{ \code{basis}: The basis object. } }}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{create_basis}}, \code{\link{eval_functions}}
NULL

#' @rdname bpr_bayes
#'
#' @examples
#' data <- meth_data
#' out_opt <- bpr_bayes(x = data, is_parallel = FALSE)
#'
#' @export
bpr_bayes <- function(x, ...){
    UseMethod("bpr_bayes")
}


# Default function for the generic function 'bpr_bayes'
bpr_bayes.default <- function(x, ...){
    stop("Object x should be either matrix or list!")
}


#' @rdname bpr_bayes
#'
#' @examples
#' ex_data <- meth_data
#' basis <- create_rbf_object(M=3)
#' out_opt <- bpr_bayes(x = ex_data, is_parallel = FALSE, basis = basis)
#'
#' @export
bpr_bayes.list <- function(x, w_mle = NULL, basis = NULL, fit_feature = NULL,
                           cpg_dens_feat = FALSE, w_0_mean = NULL,
                           w_0_cov = NULL, gibbs_nsim = 20, gibbs_burn_in = 10,
                           is_parallel = TRUE, no_cores = NULL, ...){
    # Check that x is a list object
    assertthat::assert_that(is.list(x))

    # Extract number of observations
    N <- length(x)
    assertthat::assert_that(N > 0)

    # Perform checks for initial parameter values
    out <- .do_checks_bpr_bayes(w = w_mle, basis = basis)
    w_mle <- out$w
    basis <- out$basis

    # Number of coefficients
    M <- basis$M + 1

    if (is.vector(w_mle)){
        w_mle <- matrix(w_mle, ncol = M, nrow = N, byrow = TRUE)
    }

    if (is.null(w_0_mean)){
        w_0_mean <- rep(0, M)
    }else{
        if (length(w_0_mean) != M){
            message("Initializing prior mean due to length differences.")
            w_0_mean <- rep(0, M)
        }
    }
    if (is.null(w_0_cov)){
        w_0_cov <- diag(2, M)
    }else{
        if (NROW(w_0_cov) != M){
            message("Initializing prior mean due to length differences.")
            w_0_cov <- diag(2, M)
        }
    }

    # Initialize so the CMD check on R passes without NOTES
    i <- 0

    # If parallel mode is ON
    if (is_parallel){
        # If number of cores is not given
        if (is.null(no_cores)){
            no_cores <- parallel::detectCores() - 2
        }else{
            if (no_cores >= parallel::detectCores()){
                no_cores <- parallel::detectCores() - 1
            }
        }
        if (is.na(no_cores)){
            no_cores <- 2
        }
        # Create cluster object
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)

        # Parallel optimization for each element of x, i.e. for each region i.
        res <- foreach::"%dopar%"(obj = foreach::foreach(i = 1:N),
                                  ex  = {
                  out_opt <- bpr_bayes.matrix(x           = x[[i]],
                                              w_mle       = w_mle[i, ],
                                              basis       = basis,
                                              fit_feature = fit_feature,
                                              cpg_dens_feat = cpg_dens_feat,
                                              w_0_mean    = w_0_mean,
                                              w_0_cov     = w_0_cov,
                                              gibbs_nsim  = gibbs_nsim,
                                              gibbs_burn_in = gibbs_burn_in)
                                  })
        # Stop parallel execution
        parallel::stopCluster(cl)
    }else{
        # Sequential optimization for each element of x, i.e. for each region i.
        res <- foreach::"%do%"(obj = foreach::foreach(i = 1:N),
                               ex  = {
               out_opt <- bpr_bayes.matrix(x           = x[[i]],
                                           w_mle       = w_mle[i, ],
                                           basis       = basis,
                                           fit_feature = fit_feature,
                                           cpg_dens_feat = cpg_dens_feat,
                                           w_0_mean    = w_0_mean,
                                           w_0_cov     = w_0_cov,
                                           gibbs_nsim  = gibbs_nsim,
                                           gibbs_burn_in = gibbs_burn_in)
                               })
    }

    message("Creating summary statistics from Gibbs chain...")
    # Matrix for storing optimized coefficients
    W_opt <- sapply(res, function(x) x$w_opt)
    W_var <- sapply(res, function(x) x$w_var)
    if (is.matrix(W_opt)){
        W_opt <- t(W_opt)
        W_var <- t(W_var)
    }else{
        W_opt <- as.matrix(W_opt)
        W_var <- as.matrix(W_var)
    }
    colnames(W_opt) <- paste("w", seq(1, NCOL(W_opt)), sep = "")
    colnames(W_var) <- paste("w", seq(1, NCOL(W_var)), sep = "")

    # Matrix for storing the centers of RBFs if object class is 'rbf'
    Mus <- NULL
    if (methods::is(basis, "rbf")){
        if (is.null(basis$mus)){
            Mus <- sapply(lapply(res, function(x) x$basis), function(y) y$mus)
            if (is.matrix(Mus)){
                Mus <- t(Mus)
            }else{
                Mus <- as.matrix(Mus)
            }
            colnames(Mus) <- paste("mu", seq(1, NCOL(Mus)), sep = "")
        }
    }

    return(list(W_opt = W_opt,
                W_var = W_var,
                Mus = Mus,
                basis = basis))
}


#' @rdname bpr_bayes
#'
#' @examples
#' basis <- create_rbf_object(M=2)
#' w <- c(0.1, 0.1, 0.1)
#' w_0_mean <- rep(0, length(w))
#' w_0_cov <- diag(10, length(w))
#' data <- meth_data[[1]]
#' out_opt <- bpr_bayes(x = data, w_mle = w, w_0_mean = w_0_mean,
#'                      w_0_cov = w_0_cov, basis = basis)
#'
#' basis <- create_rbf_object(M=0)
#' w <- c(0.1)
#' w_0_mean <- rep(0, length(w))
#' w_0_cov <- diag(10, length(w))
#' data <- meth_data[[1]]
#' out_opt <- bpr_bayes(x = data, w_mle = w, w_0_mean = w_0_mean,
#'                      w_0_cov = w_0_cov, basis = basis)
#'
#' @importFrom stats optim
#' @importFrom truncnorm rtruncnorm
#' @importFrom mvtnorm rmvnorm
#'
#' @export
bpr_bayes.matrix <- function(x, w_mle = NULL, basis = NULL, fit_feature = NULL,
                             cpg_dens_feat = FALSE, w_0_mean = NULL,
                             w_0_cov = NULL, gibbs_nsim = 20,
                             gibbs_burn_in = 10, ...){

    # Vector for storing CpG locations relative to TSS
    obs <- as.vector(x[ ,1])

    # Create design matrix H
    des_mat <- .design_matrix(x = basis, obs = obs)
    H       <- des_mat$H
    basis   <- des_mat$basis

    # Number of coefficients
    M <- basis$M + 1

    # Conjugate prior on the coefficients \w ~ N(w_0_mean, w_0_cov)
    if (is.null(w_0_mean)){
        w_0_mean <- rep(0, M)
    }else{
        if (length(w_0_mean) != M){
            w_0_mean <- rep(0, M)
        }
    }
    if (is.null(w_0_cov)){
        w_0_cov <- diag(2, M)
    }else{
        if (NROW(w_0_cov) != M){
            w_0_cov <- diag(2, M)
        }
    }
    # Initialize regression coefficients
    if (is.null(w_mle)){
        w_mle <- rep(0.5, M)
    }

    # If we have Binomial distributed data
    if (NCOL(x) == 3){
        # Total number of reads for each CpG
        N_i <- as.vector(x[, 2])
        # Corresponding number of methylated reads for each CpG
        m_i <- as.vector(x[, 3])

        # Sum of total trials for each observation i
        N <- sum(N_i)

        # Create extended vector y of length (N x 1)
        y <- vector(mode = "integer")
        for (i in 1:NROW(x)){
            y <- c(y, rep(1, m_i[i]), rep(0, N_i[i] - m_i[i]))
        }

        # Create extended design matrix Xx of dimension (N x M)
        H <- as.matrix(H[rep(1:NROW(H), N_i), ])
    }else{ # If we have Bernoulli distributed data
        N <- length(obs)
        y <- x[, 2]
    }

    N1  <- sum(y)  # Number of successes
    N0  <- N - N1  # Number of failures

    # ---------------------------------
    # Gibbs sampling algorithm
    # ---------------------------------
    w_chain <- .gibbs_bpr(H  = H,
                          y  = y,
                          N  = N,
                          N0 = N0,
                          N1 = N1,
                          w_mle      = w_mle,
                          w_0_mean   = w_0_mean,
                          w_0_cov    = w_0_cov,
                          gibbs_nsim = gibbs_nsim)

    if (M == 1){
        w_opt <- mean(w_chain[-(1:gibbs_burn_in)])
        w_var <- var(w_chain[-(1:gibbs_burn_in)])
    }else{
        w_opt <- colMeans(w_chain[-(1:gibbs_burn_in), ])
        w_var <- apply(w_chain[-(1:gibbs_burn_in), ], 2, var)
    }

    # If we need to add the goodness of fit to the data as feature
    if (!is.null(fit_feature)){
        if (identical(fit_feature, "NLL")){
            fit <- .bpr_likelihood(w = w_opt,
                                   H = H,
                                   data = x[ ,2:3],
                                   is_NLL = TRUE)
        }else if (identical(fit_feature, "RMSE")){
            # Predictions of the target variables
            f_pred <- as.vector(pnorm(H %*% w_opt))
            f_true <- x[ ,3] / x[ ,2]
            fit <- sqrt(mean( (f_pred - f_true) ^ 2) )
        }
        w_opt <- c(w_opt, fit)
    }

    # Add as feature the CpG density in the promoter region
    if (cpg_dens_feat){
        w_opt <- c(w_opt, length(obs))
    }

    return(list(w_opt = w_opt,
                w_var = w_var,
                basis = basis))
}


# Internal function to make all the appropriate type checks.
.do_checks_bpr_bayes <- function(w = NULL, basis = NULL){
    if (is.null(basis)){
        basis <- create_rbf_object(M = 3)
    }
    if (is.null(w)){
        w <- rep(0.5, basis$M + 1)
    }
    if (is.matrix(w)){
        if (length(w[1,]) != (basis$M + 1) ){
            stop("Coefficients vector should be M+1!")
        }
    }else{
        if (length(w) != (basis$M + 1) ){
            stop("Coefficients vector should be M+1, M: number of basis functions!")
        }
    }
    return(list(w = w, basis = basis))
}
