#' Cluster single cells based on methylation profiles
#'
#' \code{sc_bpr_cluster_wrap} is a wrapper function that clusters single-cells
#' based on their DNA methylation profiles using the EM algorithm, where the
#' observation model is the Binomial/Bernoulli distributed Probit Regression
#' likelihood. Initially, it performs parameter checking, runs a 'mini' EM to
#' fnd the optimal starting parameter values, and then the EM algorithm is
#' applied and finally model selection metrics are calculated, such as BIC and
#' AIC.
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
#' @param em_max_iter Integer denoting the maximum number of EM iterations.
#' @param epsilon_conv Numeric denoting the convergence parameter for EM.
#' @param use_kmeans Logical, use k-means for initializing centres or randmoly
#'   picking a point a cluster centre.
#' @param em_init_nstart Number of EM random starts for finding optimal
#'   likelihood.
#' @param em_init_max_iter Maximum number of EM iterations for the 'small' init
#'   EM.
#' @param init_opt_itnmax Optimization iterations for obtaining the initial EM
#'   parameter values.
#' @param is_verbose Logical, print results during EM iterations
#' @inheritParams bpr_optimize
#'
#' @return A 'sc_bpr_cluster' object which, in addition to the input parameters,
#'   consists of the following variables: \itemize{ \item{\code{pi_k}: Fitted
#'   mixing proportions.} \item{\code{w}: A N x M x K array matrix with the
#'   fitted coefficients of the basis functions for each cluster k and region
#'   n.} \item{\code{NLL}: The Negative Log Likelihood after the EM algorithm
#'   has finished.} \item{\code{post_prob}: Posterior probabilities of each cell
#'   belonging to each cluster.} \item{\code{labels}: Hard clustering
#'   assignments of each cell.} \item{\code{BIC}: Bayesian Information Criterion
#'   metric.} \item{\code{AIC}: Akaike Information Criterion metric.}
#'   \item{\code{ICL}: Integrated Complete Likelihood criterion metric.} }
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
sc_bpr_cluster_wrap <- function(x, K = 2, pi_k = NULL, w = NULL, basis = NULL,
                                lambda = 1/8, em_max_iter = 100, epsilon_conv = 1e-05,
                                use_kmeans = TRUE, em_init_nstart = 10,
                                em_init_max_iter = 10, opt_method = "CG", opt_itnmax = 50,
                                init_opt_itnmax = 100, is_parallel = TRUE,
                                no_cores = NULL, is_verbose = FALSE){

    # Check that x is a list object
    assertthat::assert_that(is.list(x))
    assertthat::assert_that(is.list(x[[1]]))

    I <- length(x)      # Extract number of cells
    N <- length(x[[1]]) # Extract number of promoter regions
    assertthat::assert_that(N > 0)

    # If parallel mode is ON
    if (is_parallel){
        # If number of cores is not given
        if (is.null(no_cores)){ no_cores <- parallel::detectCores() - 2 }
        else{
            if (no_cores >= parallel::detectCores()){
                no_cores <- parallel::detectCores() - 1
            }
        }
        if (is.na(no_cores)){ no_cores <- 2 }
        # Create cluster object
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
    }

    if (is.null(basis)){ basis <- create_rbf_object(M = 3) }

    ind <- list() # Keep a list of promoter regions with CpG coverage
    H <- list()   # Create design matrix for each cell for each genomic region
    for (i in 1:I){
        H[[i]]   <- vector(mode = "list", length = N)
        H[[i]]   <- lapply(H[[i]], function(x) x = NA)
        ind[[i]] <- which(!is.na(x[[i]]))
        if (is_parallel){
            H[[i]][ind[[i]]] <- parallel::mclapply(X = x[[i]][ind[[i]]], FUN = function(y)
                design_matrix(x = basis, obs = y[, 1])$H, mc.cores = no_cores)
        }else{
            H[[i]][ind[[i]]] <- lapply(X = x[[i]][ind[[i]]], FUN = function(y)
                design_matrix(x = basis, obs = y[, 1])$H)
        }
    }
    if (is_parallel){ # Stop parallel execution
        parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
    }

    # Apply EM algorithm to cluster similar cells based on methylation profiles
    if (is_verbose) { message("Running mini EM ...\n") }
    # Perform checks for initial parameter values
    out <- .do_scEM_checks(x = x, H = H, reg_ind = ind, K = K, pi_k = pi_k, w = w,
                           basis = basis, lambda = lambda, em_init_nstart = em_init_nstart,
                           em_init_max_iter = em_init_max_iter, epsilon_conv = epsilon_conv,
                           opt_method = opt_method, opt_itnmax = opt_itnmax,
                           init_opt_itnmax = init_opt_itnmax,
                           is_parallel = is_parallel, no_cores = no_cores,
                           is_verbose = is_verbose)

    w    <- out$w        # Initial weights for each cluster and cell region
    pi_k <- out$pi_k     # Initial mixing proportions
    M    <- basis$M + 1  # Number of coefficient parameters

    # Apply EM algorithm to cluster similar cells based on methylation profiles
    if (is_verbose) { message("Clustering cells based on methylation profiles via EM ...\n") }
    scbpr_cluster <- .scbpr_EM(x = x, H = H, reg_ind = ind, K = K, pi_k = pi_k,
                               w = w, basis = basis, lambda = lambda,
                               em_max_iter = em_max_iter, epsilon_conv = epsilon_conv,
                               opt_method = opt_method, opt_itnmax = opt_itnmax,
                               is_parallel = is_parallel, no_cores = no_cores,
                               is_verbose = is_verbose)
    if (is_verbose) { message("Finished clustering!\n\n") }

    # Add names to the estimated parameters for clarity
    names(scbpr_cluster$pi_k) <- paste0("clust", 1:K)
    ###colnames(bpr_cluster$w) <- paste0("clust", 1:K)

    # Get hard cluster assignments for each observation
    scbpr_cluster$labels <- apply(X = scbpr_cluster$post_prob, MARGIN = 1,
                                  FUN = function(x) which(x == max(x, na.rm = TRUE)))
    # Perform model selection
    total_params <- (K - 1) + K * M * N
    # Bayesian Information Criterion
    scbpr_cluster$BIC <- 2 * utils::tail(scbpr_cluster$NLL, n = 1) + total_params * log(I)
    # Akaike Iformation Criterion
    scbpr_cluster$AIC <- 2 * utils::tail(scbpr_cluster$NLL, n = 1) + 2 * total_params
    # Integrated Complete Likelihood criterion
    entropy <- (-1) * sum(scbpr_cluster$post_prob * log(scbpr_cluster$post_prob), na.rm = TRUE)
    scbpr_cluster$ICL <- scbpr_cluster$BIC + entropy
    # Store initial max optimization iterations
    scbpr_cluster$init_opt_itnmax <- init_opt_itnmax
    class(scbpr_cluster) <- "scbpr_cluster"

    return(scbpr_cluster)
}


#'
#' EM algorithm
#'
.scbpr_EM <- function(x, H, reg_ind, K = 2, pi_k, w, basis, lambda = 1/6,
                      em_max_iter = 100, epsilon_conv = 1e-05, opt_method = "CG",
                      opt_itnmax = 50, is_parallel = TRUE, no_cores = NULL,
                      is_verbose = FALSE){

    I <- length(x)      # Number of cells
    N <- length(x[[1]]) # Number of regions
    M <- basis$M + 1    # Number of basis functions
    NLL <- 1e+100       # Initialize and store NLL for each EM iteration
    n <- 0

    # Matrices / Lists for storing results
    w_pdf     <- matrix(0, nrow = I, ncol = K)  # Store weighted PDFs
    post_prob <- matrix(0, nrow = I, ncol = K)  # Hold responsibilities
    w_tmp     <- array(data = 0, dim = c(N, M, K))

    # Run EM algorithm until convergence
    for (t in 1:em_max_iter){
        # TODO: Handle empty clusters!!!
        # TODO: Handle empty clusters!!!
        ## ---------------------------------------------------------------
        # Compute weighted pdfs for each cluster
        for (k in 1:K){
            # Apply to each cell and only to regions with CpG coverage
            w_pdf[, k] <- log(pi_k[k]) + vapply(X = 1:I, FUN = function(i) sum(vapply(X = reg_ind[[i]], FUN = function(y)
                bpr_likelihood(w = w[y, , k], H = H[[i]][[y]], data = x[[i]][[y]], lambda = lambda, is_NLL = FALSE),
                FUN.VALUE = numeric(1), USE.NAMES = FALSE)), FUN.VALUE = numeric(1), USE.NAMES = FALSE)
        }
        # Use the logSumExp trick for numerical stability
        Z <- apply(w_pdf, 1, .log_sum_exp)
        # Get actual posterior probabilities, i.e. responsibilities
        post_prob <- exp(w_pdf - Z)
        NLL  <- c(NLL, (-1) * sum(Z))   # Evaluate NLL

        #
        # M-Step -----------------------------------------------
        #
        # Compute sum of posterior probabilities for each cluster
        I_k <- colSums(post_prob)
        # Update mixing proportions for each cluster
        pi_k <- I_k / I

        # Update basis function coefficient matrix w
        # If parallel mode is ON
        if (is_parallel){
            # Create cluster object
            cl <- parallel::makeCluster(no_cores)
            doParallel::registerDoParallel(cl)

            # Parallel optimization for each region n
            res_out <- foreach::"%dopar%"(obj = foreach::foreach(n = 1:N),
                                       ex  = {
                               out <- optim_regions(x = lapply(x, "[[", n),
                                                    H = lapply(H, "[[", n),
                                                    w = w[n, , ],
                                                    K = K,
                                                    opt_method = opt_method,
                                                    opt_itnmax = opt_itnmax,
                                                    post_prob = post_prob,
                                                    lambda = lambda)
                                       })
            for (k in 1:K){
                tmp <- sapply(res_out, function(x) x[, k])
                if (is.matrix(tmp)){ w_tmp[, , k] <- t(tmp) }
                else{ w_tmp[, 1, k] <- tmp }
            }
            w <- w_tmp

            parallel::stopCluster(cl)
            doParallel::stopImplicitCluster()
        }else{
            # Sequential optimization for each region n
            res_out <- foreach::"%do%"(obj = foreach::foreach(n = 1:N),
                                          ex  = {
                              out <- optim_regions(x = lapply(x, "[[", n),
                                                   H = lapply(H, "[[", n),
                                                   w = w[n, , ],
                                                   K = K,
                                                   opt_method = opt_method,
                                                   opt_itnmax = opt_itnmax,
                                                   post_prob = post_prob,
                                                   lambda = lambda)
                                          })
            for (k in 1:K){
                tmp <- sapply(res_out, function(x) x[, k])
                if (is.matrix(tmp)){ w_tmp[, , k] <- t(tmp) }
                else{ w_tmp[, 1, k] <- tmp }
            }
            w <- w_tmp
        }

        if (is_verbose){
            cat("\r", "It:\t", t, "\tNLL:\t", NLL[t + 1], "\tNLL_diff:\t", NLL[t] - NLL[t + 1])
        }
        if (NLL[t + 1] > NLL[t]){
            message("Negative Log Likelihood increases - Stopping EM!\n")
            break
        }
        # Check for convergence
        if (NLL[t] - NLL[t + 1] < epsilon_conv){ break }
    }

    # Check if EM converged in the given maximum iterations
    if (t == em_max_iter){
        warning("EM did not converge with the given maximum iterations!\n")
    }

    obj <- structure(list(K = K, N = N, w = w, pi_k = pi_k, lambda = lambda,
                          em_max_iter = em_max_iter, opt_method = opt_method,
                          opt_itnmax = opt_itnmax, NLL = NLL,
                          basis = basis, post_prob = post_prob),
                     class = "scbpr_EM")
    return(obj)
}

#
# Optimize a promoter regions across cells, which are weighted by the
# responsibilities of belonging to each cluster.
#
optim_regions <- function(x, H, w, K, opt_method = opt_method, opt_itnmax, post_prob, lambda){
    covered_ind <- which(!is.na(H))
    if (is.vector(w)){ w <- matrix(w, ncol = K) }
    for (k in 1:K){  # For each cluster k
        # TODO: How to handle empty regions???
        w[, k] <- optim(par       = w[, k],
                        fn        = sum_weighted_bpr_lik,
                        gr        = sum_weighted_bpr_grad,
                        method    = opt_method,
                        control   = list(maxit = opt_itnmax),
                        x         = x[covered_ind],
                        des_mat   = H[covered_ind],
                        post_prob = post_prob[covered_ind, k],
                        lambda    = lambda,
                        is_NLL    = TRUE)$par
    }
    return(w)
}

# Internal function to make all the appropriate type checks.
.do_scEM_checks <- function(x, H, reg_ind, K, pi_k = NULL, w = NULL, basis,
                            lambda = 1/6, use_kmeans = TRUE, em_init_nstart = 10,
                            em_init_max_iter = 10, epsilon_conv = 1e-05,
                            opt_method = "CG", opt_itnmax = 50, init_opt_itnmax = 100,
                            is_parallel = TRUE, no_cores = NULL, is_verbose = TRUE){
    I <- length(x)
    N <- length(x[[1]])
    M <- basis$M + 1

    if (is.null(w)){
        ww <- array(data = rnorm(N*M*I, 0, 0.01), dim = c(N, M, I))
        w_init <- rep(0.5, M)
        for (i in 1:I){
            # Compute regression coefficients using MLE
            ww[reg_ind[[i]], ,i] <- bpr_optim_fast(x = x[[i]][reg_ind[[i]]],
                                                   H = H[[i]][reg_ind[[i]]],
                                                   w = w_init,
                                                   lambda = lambda,
                                                   opt_method = opt_method,
                                                   opt_itnmax = init_opt_itnmax,
                                                   is_parallel = is_parallel,
                                                   no_cores = no_cores)$W_opt
        }

        # Transform to long format to perform k-means
        W_opt <- matrix(0, nrow = I, ncol = N * M)
        for (i in 1:I){ W_opt[i, ] <- as.vector(ww[,,i]) }
        w <- array(data = 0, dim = c(N, M, K))
        NLL_prev <- 1e+120
        optimal_w = optimal_pi_k <- NULL

        # Run 'mini' EM algorithm to find optimal starting points
        for (t in 1:em_init_nstart){
            if (use_kmeans){
                # Use Kmeans with random starts
                cl <- stats::kmeans(W_opt, K, nstart = 1)
                # Get the mixture components
                C_n <- cl$cluster
                # TODO: Check that k-means does not return empty clusters..
                # Sample randomly one point from each cluster as initial centre
                for (k in 1:K){ w[, ,k] <- ww[, , sample(which(C_n == k), 1)] }
                # Mixing proportions
                if (is.null(pi_k)){ pi_k <- as.vector(table(C_n) / I ) }
            }else{
                w <- array(data = ww[, ,sample(I, K)], dim = c(N, M, K))
                if (is.null(pi_k)){ pi_k <- rep(1/K, K) }
            }
            # Run mini EM
            em <- .scbpr_EM(x = x, H = H, reg_ind = reg_ind, K = K, pi_k = pi_k,
                            w = w, basis = basis, lambda = lambda,
                            em_max_iter = em_init_max_iter, epsilon_conv = epsilon_conv,
                            opt_method = opt_method,
                            opt_itnmax = opt_itnmax, is_parallel = is_parallel,
                            no_cores = no_cores, is_verbose = is_verbose)

            # Check if NLL is lower and keep the optimal params
            NLL_cur <- utils::tail(em$NLL, n = 1)
            if (NLL_cur < NLL_prev){
                # TODO:: Store optimal pi_k from EM
                optimal_w <- w
                optimal_pi_k <- pi_k
                NLL_prev <- NLL_cur
            }
        }
    }
    if (is.null(pi_k)){ optimal_pi_k <- rep(1 / K, K) }
    if (length(w[1,,1]) != (basis$M + 1) ){
        stop("Coefficient vector should be M+1, M: number of basis functions!")
    }
    return(list(w = optimal_w, basis = basis, pi_k = optimal_pi_k))
}
