#' Perform likelihood ratio test between methylation samples
#'
#' \code{bpr_diff_test_wrap} is a function that wraps all the necessary
#' subroutines for performing statistical testing between methylation samples
#' using the likelihood ratio test.
#'
#' @param x The binomial distributed observations. A list containing two lists
#'   for control and treatment samples. Each list has elements of length N,
#'   where each element is an L x 3 matrix of observations, where 1st column
#'   contains the locations. The 2nd and 3rd columns contain the total reads and
#'   number of successes at the corresponding locations, repsectively. See
#'   \code{\link{process_haib_caltech_wrap}} on a possible way to get this data
#'   structure.
#' @param w Optional vector of initial parameter / coefficient values.
#' @param basis Optional basis function object, default is an 'rbf' object, see
#'   \code{\link{create_rbf_object}}.
#' @inheritParams bpr_optimize
#'
#' @return A 'bpr_test' object which, in addition to the input
#'   parameters, consists of the following variables: \itemize{ \item{
#'   \code{W_opt}: An Nx(2M+2) matrix with the optimized parameter values. Each
#'   row of the matrix corresponds to the concatenated coefficients of the
#'   methylation profiles from both samples. The columns are of the same length
#'   as the concatenated parameter vector [w_contr, w_treat] (i.e. number of
#'   basis functions). } \item{ \code{Mus}: A list containing two matrices of
#'   size N x M with the RBF centers for each sample, if basis object is
#'   \code{\link{create_rbf_object}}, otherwise NULL.} \item{train}: The
#'   training data. \item{test}: The test data. \item \code{gex_model}: The
#'   fitted regression model. \item \code{train_pred} The predicted values for
#'   the training data. \item \code{test_pred} The predicted values for the test
#'   data. \item \code{train_errors}: The training error metrics. \item
#'   \code{test_errors}: The test error metrics.}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{bpr_optimize}}, \code{\link{create_basis}},
#'   \code{\link{eval_functions}}, \code{\link{train_model_gex}},
#'   \code{\link{predict_model_gex}}
#'
#' @export
bpr_diff_test_wrap <- function(x, w = NULL, basis = NULL,
                               opt_method = "CG", opt_itnmax = 100,
                               is_parallel = TRUE, no_cores = NULL){

    # Check that x is a list object
    assertthat::assert_that(is.list(x))

    # Learn methylation profiles for control samples
    message("Learning control methylation profiles ...\n")
    out_contr_opt <- bpr_optim(x           = x$control,
                               w           = w,
                               basis       = basis,
                               fit_feature = "NLL",
                               cpg_dens_feat = FALSE,
                               opt_method  = opt_method,
                               opt_itnmax  = opt_itnmax,
                               is_parallel = is_parallel,
                               no_cores    = no_cores)

    # Learn methylation profiles for treatment samples
    message("Learning treatment methylation profiles ...\n")
    out_treat_opt <- bpr_optim(x           = x$treatment,
                               w           = w,
                               basis       = basis,
                               fit_feature = "NLL",
                               cpg_dens_feat = FALSE,
                               opt_method  = opt_method,
                               opt_itnmax  = opt_itnmax,
                               is_parallel = is_parallel,
                               no_cores    = no_cores)

    # Number of basis functions + bias term
    params <- basis$M + 1

    ##-------------------------------------
    #      Alternative hypothesis
    ##-------------------------------------
    # Obtain the NLL of the control (alternative)
    nll_contr_alt <- out_contr_opt$W_opt[, params + 1]
    # Obtain the NLL of the treatment (alternative)
    nll_treat_alt <- out_treat_opt$W_opt[, params + 1]
    # NLL for alternative hypothesis
    nll_alt <- nll_contr_alt + nll_treat_alt


    ##-------------------------------------
    #      NULL hypothesis
    ##-------------------------------------
    # Obtain the NLL of the control (null)
    nll_contr_null <- out_contr_opt$W_opt[, params + 1]
    # Obtain the NLL of the treatment (null)
    nll_treat_null <- vector(mode = "numeric", length = length(x$treatment))
    for (i in 1:length(x$treatment)){
        # Create design matrix H
        H <- .design_matrix.rbf(x   = basis,
                                obs = x$treatment[[i]][,1])$H

        # Evaluate the likelihood under control parameters
        nll_treat_null[i] <- bpr_likelihood(w = out_contr_opt$W_opt[i, 1:params],
                                            H = H,
                                            data = x$treatment[[i]],
                                            is_NLL = TRUE)
    }
    # NLL for null hypothesis
    nll_null <- nll_contr_null + nll_treat_null

    ##------------------------------------
    # Compute log likelihood ratio test
    ##------------------------------------
    lr_test <- 2 * (nll_null - nll_alt)

    message("Done!\n\n")

    # Create 'bpr_predict' object
    obj <- structure(list(basis        = out_treat_opt$basis,
                          opt_method   = opt_method,
                          opt_itnmax   = opt_itnmax,
                          W_opt_contr  = out_contr_opt$W_opt,
                          W_opt_treat  = out_treat_opt$W_opt,
                          Mus          = list(control = out_contr_opt$Mus,
                                              treatment = out_treat_opt$Mus),
                          nll_contr_alt  = nll_contr_alt,
                          nll_treat_alt  = nll_treat_alt,
                          nll_contr_null = nll_contr_null,
                          nll_treat_null = nll_treat_null,
                          nll_alt = nll_alt,
                          nll_null = nll_null,
                          lr_test = lr_test),
                     class = "bpr_test")
    return(obj)
}


