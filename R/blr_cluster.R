#' Cluster methylation profiles with Gaussian noise
#'
#' \code{blr_cluster} is a wrapper function that clusters methylation
#' profiles using the EM algorithm. Initially, it performs parameter checking,
#' and initializes main parameters, such as mixing proportions, basis function
#' coefficients, then the EM algorithm is applied and finally model selection
#' metrics are calculated, such as BIC and AIC.
#'
#' @param x The Gaussian distributed observations, which has to be a list of
#'   elements of length N, where each element is an L x 2 matrix of
#'   observations, where 1st column contains the locations and the 2nd column
#'   contains the methylation levels.
#' @param K Integer denoting the number of clusters K.
#' @param pi_k Vector of length K, denoting the mixing proportions.
#' @param w A MxK matrix, where each column consists of the basis function
#'   coefficients for each corresponding cluster.
#' @param basis A 'basis' object. E.g. see \code{\link{create_rbf_object}}.
#' @param s2 Vector of initial linear regression variances for each cluster.
#' @param em_max_iter Integer denoting the maximum number of EM iterations.
#' @param epsilon_conv Numeric denoting the convergence parameter for EM.
#' @param lambda The complexity penalty coefficient for ridge regression.
#' @param is_verbose Logical, print results during EM iterations.
#'
#' @return A 'blr_cluster' object which, in addition to the input
#'   parameters, consists of the following variables: \itemize{
#'   \item{\code{pi_k}: Fitted mixing proportions.} \item{\code{w}: A MxK matrix
#'   with the fitted coefficients of the basis functions for each cluster k.}
#'   \item{\code{NLL}: The Negative Log Likelihood after the EM algorithm has
#'   finished.} \item{\code{post_prob}: Posterior probabilities of each promoter
#'   region belonging to each cluster.} \item{\code{labels}: Hard clustering
#'   assignments of each observation/promoter region.} \item{\code{BIC}:
#'   Bayesian Information Criterion metric.} \item{\code{AIC}: Akaike
#'   Information Criterion metric.} \item{\code{ICL}: Integrated Complete
#'   Likelihood criterion metric.} }
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @examples
#' ex_data <- meth_data
#' my_clust <- blr_cluster(x = lm_data, em_max_iter = 100, is_verbose = TRUE)
#'
#' @export
blr_cluster <- function(x, K = 3, pi_k = NULL, w = NULL, basis = NULL,
                             s2 = NULL, em_max_iter = 100, epsilon_conv = 1e-04,
                             lambda = 1/10, is_verbose = FALSE){
    # Check that x is a list object
    assertthat::assert_that(is.list(x))

    # Extract number of observations
    N <- length(x)
    assertthat::assert_that(N > 0)

    # Perform checks for initial parameter values
    out <- .do_blr_EM_checks(x = x,
                            K = K,
                            pi_k = pi_k,
                            w = w,
                            s2 = s2,
                            basis = basis,
                            lambda = lambda)
    w     <- out$w
    basis <- out$basis
    pi_k  <- out$pi_k
    s2    <- out$s2

    # Store weighted PDFs
    weighted_pdf <- matrix(0, nrow = N, ncol = K)
    # Initialize and store NLL for each EM iteration
    NLL <- c(1e+40)

    # Create design matrix for each observation
    des_mat <- lapply(X = x, FUN = function(y) design_matrix(x = basis, obs = y[, 1])$H)

    # Precompute total observations for faster implementation
    L_n <- vector("numeric", length = N)
    for (n in 1:N){
        L_n[n] <- NROW(x[[n]])
    }

    # Apply EM algorithm to cluster similar methylation profiles
    message("Clustering methylation profiles via EM ...\n")

    # Run EM algorithm until convergence
    for (t in 1:em_max_iter){
        #
        # E-Step -----------------------------------------------
        for (k in 1:K){
            # For each element in x, evaluate the BPR log likelihood
            weighted_pdf[, k] <- vapply(X = 1:N,
                                        FUN = function(y)
                                            sum(dnorm(x = x[[y]][,2], mean = des_mat[[y]] %*% w[, k],
                                                      sd = sqrt(s2[k]), log = TRUE)) - lambda * t(w[, k]) %*% w[, k],
                                        FUN.VALUE = numeric(1),
                                        USE.NAMES = FALSE
                                        )
            weighted_pdf[, k] <- log(pi_k[k]) + weighted_pdf[, k]
        }

        # Calculate probs using the logSumExp trick for numerical stability
        Z <- apply(weighted_pdf, 1, .log_sum_exp)
        # Get actual posterior probabilities, i.e. responsibilities
        post_prob <- exp(weighted_pdf - Z)
        # Evaluate and store the NLL
        NLL  <- c(NLL, (-1) * sum(Z))

        #
        # M-Step -----------------------------------------------
        #
        # Compute sum of posterior probabilities for each cluster
        N_k <- colSums(post_prob)
        # Update mixing proportions for each cluster
        pi_k <- N_k / N

        for (k in 1:K){
            # Update basis function coefficient vector w for each cluster
            tmp_HH <- matrix(0, ncol = basis$M +1, nrow = basis$M + 1)
            tmp_Hy <- vector("numeric", length = basis$M +1)
            for (n in 1:N){
                tmp_HH <- tmp_HH + crossprod(des_mat[[n]]) * post_prob[n, k]
                tmp_Hy <- tmp_Hy + crossprod(des_mat[[n]], x[[n]][,2]) * post_prob[n, k]
            }
            w[, k] <- solve(tmp_HH + lambda) %*% tmp_Hy

            # Update variance of regression model for each cluster
            tmp <- sum(vapply(X = 1:N,
                          FUN = function(y)
                              crossprod(x[[y]][,2] - des_mat[[y]] %*% w[, k]) * post_prob[y, k],
                          FUN.VALUE = numeric(1),  USE.NAMES = FALSE ))
            s2[k] <- tmp / post_prob[, k] %*% L_n
        }

        if (is_verbose){
            cat("It:\t", t, "\tNLL:\t", NLL[t + 1],
                "\tNLL_diff:\t", NLL[t] - NLL[t + 1], "\n")
        }
        if (NLL[t + 1] > NLL[t]){
            message("Negative Log Likelihood increases - Stopping EM!\n")
            break
        }
        # Check for convergence
        if (NLL[t] - NLL[t + 1] < epsilon_conv){
            break
        }
    }
    # Check if EM converged in the given maximum iterations
    if (t == em_max_iter){
        warning("EM did not converge with the given maximum iterations!\n")
    }

    bpr_cluster <- list(K = K,
                        N = N,
                        w = w,
                        s2 = s2,
                        pi_k = pi_k,
                        lambda = lambda,
                        em_max_iter = em_max_iter,
                        NLL = NLL,
                        basis = basis,
                        post_prob = post_prob)
    message("Finished clustering!\n\n")

    # Add names to the estimated parameters for clarity
    names(bpr_cluster$pi_k) <- paste0("clust", 1:K)
    colnames(bpr_cluster$w) <- paste0("clust", 1:K)

    # Get hard cluster assignments for each observation
    bpr_cluster$labels <- apply(X      = bpr_cluster$post_prob,
                                MARGIN = 1,
                                FUN    = function(x)
                                    which(x == max(x, na.rm = TRUE)))

    # Perform model selection
    total_params <- (K - 1) + K * NROW(w)

    # Bayesian Information Criterion
    bpr_cluster$BIC <- 2 * utils::tail(bpr_cluster$NLL, n = 1) +
        total_params * log(N)

    # Akaike Iformation Criterion
    bpr_cluster$AIC <- 2 * utils::tail(bpr_cluster$NLL, n = 1) +
        2 * total_params

    # Integrated Complete Likelihood criterion
    entropy <- (-1) * sum(bpr_cluster$post_prob * log(bpr_cluster$post_prob),
                          na.rm = TRUE)
    bpr_cluster$ICL <- bpr_cluster$BIC + entropy

    class(bpr_cluster) <- "blr_cluster"

    return(bpr_cluster)
}


# Internal function to make all the appropriate type checks.
.do_blr_EM_checks <- function(x, K = 2, pi_k = NULL,  w = NULL, s2 = NULL, basis = NULL,
                          lambda = 1/10){
    if (is.null(basis)){
        basis <- create_rbf_object(M = 3)
    }
    if (is.null(w)){
        w <- rep(0.5, basis$M + 1)

        # Keep only the optimized coefficients
        W_opt <- matrix(0, ncol = length(w), nrow = length(x))
        for (i in 1:length(x)){
            W_opt[i, ] <- blr(x = x[[i]][, 1],
                                   y = x[[i]][, 2],
                                   basis = basis,
                                   lambda = lambda,
                                   return.all = FALSE)$coefficients
        }
        # Use Kmeans with random starts
        cl <- stats::kmeans(W_opt, K, nstart = 25)
        # Get the mixture components
        C_n <- cl$cluster
        # Mean for each cluster
        w <- t(cl$centers)

        # Mixing proportions
        if (is.null(pi_k)){
            N <- length(x)
            pi_k <- as.vector(table(C_n) / N)
        }
    }
    if (is.null(pi_k)){
        pi_k <- rep(1 / K, K)
    }
    if (is.null(s2)){
        s2 <- rep(0.2, K)
    }
    if (NROW(w) != (basis$M + 1) ){
        stop("Coefficient vector should be M+1, M: number of basis functions!")
    }
    return(list(w = w, basis = basis, pi_k = pi_k, s2 = s2))
}
