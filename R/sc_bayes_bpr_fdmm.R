#' Gibbs sampling algorithm for sc-Bayesian BPR finite mixture model
#'
#' \code{sc_bayes_bpr_fdmm} implements the Gibbs sampling algorithm for
#' performing clustering of single cells based on their DNA methylation
#' profiles, where the observation model is the Bernoulli distributed Probit
#' Regression likelihood.
#'
#' @param x A list of length I, where I are the total number of cells. Each
#'   element of the list contains another list of length N, where N is the total
#'   number of genomic regions. Each element of the inner list is an L x 2
#'   matrix of observations, where 1st column contains the locations and the 2nd
#'   column contains the methylation level of the corresponding CpGs.
#' @param K Integer denoting the number of clusters K.
#' @param pi_k Vector of length K, denoting the mixing proportions.
#' @param w A N x M x K array, where each column contains the basis function
#'   coefficients for the corresponding cluster.
#' @param basis A 'basis' object. E.g. see \code{\link{create_rbf_object}}
#' @param w_0_mean The prior mean hyperparameter for w
#' @param w_0_cov The prior covariance hyperparameter for w
#' @param dir_a The Dirichlet concentration parameter, prior over pi_k
#' @param lambda The complexity penalty coefficient for penalized regression.
#' @param gibbs_nsim Argument giving the number of simulations of the Gibbs
#'   sampler.
#' @param gibbs_burn_in Argument giving the burn in period of the Gibbs sampler.
#' @param inner_gibbs Logical, indicating if we should perform Gibbs sampling to
#'   sample from the augmented BPR model.
#' @param gibbs_inner_nsim Number of inner Gibbs simulations.
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores Number of cores to be used, default is max_no_cores - 1.
#' @param is_verbose Logical, print results during EM iterations
#'
#' @importFrom stats rmultinom rnorm
#' @importFrom MCMCpack rdirichlet
#' @importFrom truncnorm rtruncnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
sc_bayes_bpr_fdmm <- function(x, K = 2, pi_k = rep(1/K, K), w = NULL, basis = NULL, w_0_mean = NULL, w_0_cov = NULL,
                              dir_a = rep(1/K, K), lambda = 1/2, gibbs_nsim = 5000, gibbs_burn_in = 1000,
                              inner_gibbs = FALSE, gibbs_inner_nsim = 50, is_parallel = TRUE, no_cores = NULL, is_verbose = TRUE){

    # Check that x is a list object
    assertthat::assert_that(is.list(x))
    assertthat::assert_that(is.list(x[[1]]))

    if (is_parallel){ # If parallel mode is ON
        # If number of cores is not given
        if (is.null(no_cores)){ no_cores <- parallel::detectCores() - 2
        }else{
            if (no_cores >= parallel::detectCores()){ no_cores <- parallel::detectCores() - 1 }
        }
        if (is.na(no_cores)){ no_cores <- 2 }
        if (no_cores > K){ no_cores <- K }
        # Create cluster object
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
    }

    I <- length(x)      # Extract number of cells
    N <- length(x[[1]]) # Extract number of promoter regions
    if (is.null(basis)){ basis <- create_rbf_object(M = 3) }
    M <- basis$M + 1    # Number of coefficient parameters

    # Initialize priors over the parameters
    if (is.null(w_0_mean)){ w_0_mean <- rep(0, M) }
    if (is.null(w_0_cov)){ w_0_cov <- diag(3, M) }

    prec_0 <- solve(w_0_cov)          # Invert covariance matrix to get the precision matrix
    w_0_prec_0 <- prec_0 %*% w_0_mean # Compute product of prior mean and prior precision matrix

    # Matrices / Lists for storing results
    w_pdf        <- matrix(0, nrow = I, ncol = K)  # Store weighted PDFs
    post_prob    <- matrix(0, nrow = I, ncol = K)  # Hold responsibilities
    C            <- matrix(0, nrow = I, ncol = K)  # Mixture components
    C_prev       <- C                              # Keep previous components
    C_matrix     <- matrix(0, nrow = I, ncol = K)  # Total mixture components
    NLL          <- vector(mode = "numeric", length = gibbs_nsim)
    NLL[1]       <- 1e+100

    H = y = z = V <- list()
    for (k in 1:K){
        H[[k]] <- vector("list", N) # List of concatenated design matrices
        y[[k]] <- vector("list", N) # List of observed methylation data
        z[[k]] <- vector("list", N) # List of auxiliary latent variables
        V[[k]] <- vector("list", N) # List of posterior variances
    }
    len_y <- matrix(0, nrow = K, ncol = N) # Total CpG observations per region
    sum_y <- matrix(0, nrow = K, ncol = N) # Total methylated CpGs per region

    # Store mixing proportions draws
    pi_draws <- matrix(NA_real_, nrow = gibbs_nsim, ncol = K)
    pi_draws[1, ] <- pi_k

    # Store BPR coefficient draws
    w_draws <- array(data = 0, dim = c(gibbs_nsim - gibbs_burn_in, N, M , K))
    # TODO: Initialize w in a sensible way
    if (is.null(w)){ w <- array(data = 0, dim = c(N, M, K)) }

    ind <- list()     # Keep a list of promoter regions with CpG coverage
    des_mat <- list() # Create design matrix for each cell for each promoter region
    for (i in 1:I){
        des_mat[[i]] <- vector(mode = "list", length = N)
        ind[[i]]     <- which(!is.na(x[[i]]))
        if (is_parallel){  # TODO: Make this function faster?
            des_mat[[i]][ind[[i]]] <- parallel::mclapply(X = x[[i]][ind[[i]]], FUN = function(y)
                design_matrix(x = basis, obs = y[, 1])$H, mc.cores = no_cores)
        }else{
            des_mat[[i]][ind[[i]]] <- lapply(X = x[[i]][ind[[i]]], FUN = function(y) design_matrix(x = basis, obs = y[, 1])$H)
        }
    }
    if (is_parallel){ # Stop parallel execution
        parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
    }

    message("Starting Gibbs sampling...")
    # Show progress bar
    pb <- txtProgressBar(min = 1, max = gibbs_nsim, style = 3)
    ##---------------------------------
    # Start Gibbs sampling
    ##---------------------------------
    for (t in 2:gibbs_nsim){
        ## ---------------------------------------------------------------
        # Compute weighted pdfs for each cluster
        for (k in 1:K){
            # Apply to each cell and only to regions with CpG coverage
            w_pdf[, k] <- log(pi_k[k]) + vapply(X = 1:I, FUN = function(i) sum(vapply(X = ind[[i]], FUN = function(y)
                bpr_likelihood(w = w[y, , k], H = des_mat[[i]][[y]], data = x[[i]][[y]], lambda = lambda, is_NLL = FALSE),
                FUN.VALUE = numeric(1), USE.NAMES = FALSE)), FUN.VALUE = numeric(1), USE.NAMES = FALSE)
        }
        # Use the logSumExp trick for numerical stability
        Z <- apply(w_pdf, 1, .log_sum_exp)
        # Get actual posterior probabilities, i.e. responsibilities
        post_prob <- exp(w_pdf - Z)
        NLL[t] <- -sum(Z)   # Evaluate NLL

        ## -------------------------------------------------------------------
        # Draw mixture components for ith simulation
        # Sample one point from a Multinomial i.e. ~ Discrete
        for (i in 1:I){ C[i, ] <- rmultinom(n = 1, size = 1, post_prob[i, ]) }

        ## -------------------------------------------------------------------
        # Check for empty clusters
        cell_changed <- c()
        for (k in 1:K){
            Cn_k <- colSums(C)
            if (Cn_k[k] == 0){
                celli <- sample(which(C[, which.max(Cn_k)] == 1), 1)
                while(celli %in% cell_changed){ celli <- sample(which(C[, which.max(Cn_k)] == 1), 1)}
                cell_changed <- c(cell_changed, celli)
                C[celli, ] <- 0; C[celli, k] <- 1
            }
        }
        # TODO: Should we keep all data
        if (t > gibbs_burn_in){ C_matrix <- C_matrix + C }

        ## -------------------------------------------------------------------
        # Update mixing proportions using updated cluster component counts
        pi_k <- as.vector(rdirichlet(n = 1, alpha = dir_a + colSums(C)))
        pi_draws[t, ] <- pi_k

        # Matrix to keep promoters with no CpG coverage
        empty_region <- matrix(0, nrow = N, ncol = K)
        for (k in 1:K){
            # Which cells are assigned to cluster k
            C_idx <- which(C[, k] == 1)
            # TODO: Handle cases when we have empty clusters...
            if (length(C_idx) == 0){ message("Warning: Empty cluster..."); next }
            # Check if current clusters ids are not equal to previous ones
            if (!identical(C[, k], C_prev[, k])){
                message(t, ": Not identical in cluster ", k)
                # Iterate over each promoter region
                for (n in 1:N){
                    # Initialize empty vector for observed methylation data
                    yy <- vector(mode = "integer")
                    # Concatenate the nth promoter from all cells in cluster k
                    tmp <- lapply(des_mat, "[[", n)[C_idx]
                    tmp <- do.call(rbind, tmp[!is.na(tmp)])
                    # TODO: Check when we have empty promoters....
                    if (is.null(tmp)) {
                        H[[k]][[n]] <- NA
                        empty_region[n, k] <- 1
                    }else{
                        H[[k]][[n]] <- tmp
                        # Obtain the corresponding methylation levels
                        for (cell in C_idx){
                            obs <- x[[cell]][[n]]
                            if (length(obs) > 1){ yy <- c(yy, obs[, 2]) }
                        }
                        # Precompute for faster computations
                        len_y[k, n] <- length(yy)
                        sum_y[k, n] <- sum(yy)
                        y[[k]][[n]] <- yy
                        z[[k]][[n]] <- rep(NA_real_, len_y[k, n])
                        # Compute posterior variance of w_nk
                        V[[k]][[n]] <- solve(prec_0 + crossprod(H[[k]][[n]], H[[k]][[n]]))
                    }
                }
            }
            for (n in 1:N){
                # In case we have no CpG data in this promoter
                if (is.vector(H[[k]][[n]])){ next }
                # Perform Gibbs sampling on the augmented BPR model
                if (inner_gibbs & t > 50){
                    w_inner <- matrix(0, nrow = gibbs_inner_nsim, ncol = M)
                    w_inner[1, ] <- w[n, , k]
                    for (tt in 2:gibbs_inner_nsim){
                        # Update Mean of z
                        mu_z <- H[[k]][[n]] %*% w_inner[tt - 1, ]
                        # Draw latent variable z from z | w, y, X
                        if (sum_y[k, n] == 0){
                            z[[k]][[n]] <- rtruncnorm(len_y[k, n], mean = mu_z, sd = 1, a = -Inf, b = 0)
                        }else if (sum_y[k, n] == len_y[k, n]){
                            z[[k]][[n]] <- rtruncnorm(len_y[k, n], mean = mu_z, sd = 1, a = 0, b = Inf)
                        }else{
                            z[[k]][[n]][y[[k]][[n]] == 1] <- rtruncnorm(sum_y[k, n], mean = mu_z[y[[k]][[n]] == 1], sd = 1, a = 0, b = Inf)
                            z[[k]][[n]][y[[k]][[n]] == 0] <- rtruncnorm(len_y[k, n] - sum_y[k, n], mean = mu_z[y[[k]][[n]] == 0], sd = 1, a = -Inf, b = 0)
                        }
                        # Compute posterior mean of w
                        Mu <- V[[k]][[n]] %*% (w_0_prec_0 + crossprod(H[[k]][[n]], z[[k]][[n]]))
                        # Draw variable \w from its full conditional: \w | z, X
                        if (M == 1){ w_inner[tt, ] <- c(rnorm(n = 1, mean = Mu, sd = V[[k]][[n]])) }
                        else{ w_inner[tt, ] <- c(rmvnorm(n = 1, mean = Mu, sigma = V[[k]][[n]])) }
                    }
                    if (M == 1){ w[n, , k] <- mean(w_inner[-(1:(gibbs_inner_nsim/2)), ]) }
                    else{ w[n, , k] <- colMeans(w_inner[-(1:(gibbs_inner_nsim/2)), ]) }
                }else{
                    # Update Mean of z
                    mu_z <- H[[k]][[n]] %*% w[n, , k]
                    # Draw latent variable z from its full conditional: z | w, y, X
                    if (sum_y[k, n] == 0){
                        z[[k]][[n]] <- rtruncnorm(len_y[k, n], mean = mu_z, sd = 1, a = -Inf, b = 0)
                    }else if (sum_y[k, n] == len_y[k, n]){
                        z[[k]][[n]] <- rtruncnorm(len_y[k, n], mean = mu_z, sd = 1, a = 0, b = Inf)
                    }else{
                        z[[k]][[n]][y[[k]][[n]] == 1] <- rtruncnorm(sum_y[k, n], mean = mu_z[y[[k]][[n]] == 1], sd = 1, a = 0, b = Inf)
                        z[[k]][[n]][y[[k]][[n]] == 0] <- rtruncnorm(len_y[k, n] - sum_y[k, n], mean = mu_z[y[[k]][[n]] == 0], sd = 1, a = -Inf, b = 0)
                    }
                    # Compute posterior mean of w
                    Mu <- V[[k]][[n]] %*% (w_0_prec_0 + crossprod(H[[k]][[n]], z[[k]][[n]]))
                    # Draw variable \w from its full conditional: \w | z, X
                    if (M == 1){ w[n, , k] <- c(rnorm(n = 1, mean = Mu, sd = V[[k]][[n]])) }
                    else{ w[n, , k] <- c(rmvnorm(n = 1, mean = Mu, sigma = V[[k]][[n]])) }
                }
            }
        }

        # For each empty promoter region, take the methyation profile of the
        # promoter regions that belong to another cluster
        for (n in 1:N){
            clust_empty_ind <- which(empty_region[n, ] == 1)
            # No empty promoter regions
            if (length(clust_empty_ind) == 0){ next }
            # Case that should never happen with the right preprocessing step
            else if (length(clust_empty_ind) == K){
                for (k in 1:K){ w[n, , k] <- c(rmvnorm(1, w_0_mean, w_0_cov)) }
            }else{
                # TODO: Perform a better imputation approach...
                cover_ind <- which(empty_region[n, ] == 0)
                # Randomly choose a cluster to obtain the methylation profiles
                k_imp <- sample(length(cover_ind), 1)
                for (k in seq_along(clust_empty_ind)){
                    w[n, , clust_empty_ind[k]] <- w[n, , k_imp]
                }
            }
        }
        C_prev <- C            # Make current cluster indices same as previous
        if (t > gibbs_burn_in){w_draws[t - gibbs_burn_in, , ,] <- w}   # Store the coefficient draw
        setTxtProgressBar(pb,t)
    }
    close(pb)
    message("Finished Gibbs sampling...")

    ##-----------------------------------------------
    message("Computing summary statistics...")
    # Compute summary statistics from Gibbs simulation
    if (K == 1){ pi_post <- mean(pi_draws[gibbs_burn_in:gibbs_nsim, ]) }
    else{ pi_post <- colMeans(pi_draws[gibbs_burn_in:gibbs_nsim, ]) }
    C_post <- C_matrix / (gibbs_nsim - gibbs_burn_in)
    w_post <- array(0, dim = c(N, M, K))
    for(k in 1:K){ w_post[, , k] <- colSums(w_draws[, , , k]) / (gibbs_nsim - gibbs_burn_in)  }

    # Object to keep input data
    dat <- list(K = K, N = N, I = I, M = M, basis = basis, dir_a = dir_a, lambda = lambda,
                w_0_mean = w_0_mean, w_0_cov = w_0_cov, gibbs_nsim = gibbs_nsim, gibbs_burn_in = gibbs_burn_in)
    # Object to hold all the Gibbs draws
    draws <- list(pi = pi_draws, w = w_draws, C = C_matrix, NLL = NLL)
    # Object to hold the summaries for the parameters
    summary <- list(pi = pi_post, w = w_post, C = C_post)
    # Create sc_bayes_bpr_fdmm object
    obj <- structure(list(summary = summary, draws = draws, dat = dat),
                     class = "sc_bayes_bpr_fdmm")
    return(obj)
}
