#' Predict differential gene expression from differential methylation profiles
#'
#' \code{bpr_diff_test_data_wrap} is a function that wraps all the necessary
#' subroutines for performing prediction of differential gene expression levels.
#' Initially, it optimizes the parameters of the basis functions so as to learn
#' the methylation profiles for the control and the treatment samples Then, the
#' two learned methylation profiles are concatenated to keep all coefficients
#' for both profiles. Then the learned parameters / coefficients of the basis
#' functions are given as input features for performing regression in order to
#' predict the corresponding differential (log2 fold-change) gene expression
#' levels.
#'
#' @param formula An object of class \code{\link[stats]{formula}}, e.g. see
#'   \code{\link[stats]{lm}} function. If NULL, the simple linear regression
#'   model is used.
#' @param x The binomial distributed observations. A list containing two lists
#'   for control and treatment samples. Each list has elements of length N,
#'   where each element is an L x 3 matrix of observations, where 1st column
#'   contains the locations. The 2nd and 3rd columns contain the total reads and
#'   number of successes at the corresponding locations, repsectively. See
#'   \code{\link{process_haib_caltech_wrap}} on a possible way to get this data
#'   structure.
#' @param y Corresponding gene expression data. A list containing two vectors
#'   for control and treatment samples.
#' @param model_name A string denoting the regression model. Currently,
#'   available models are: \code{"svm"}, \code{"randomForest"}, \code{"rlm"},
#'   \code{"mars"} and \code{"lm"}.
#' @param w Optional vector of initial parameter / coefficient values.
#' @param basis Optional basis function object, default is an 'rbf' object, see
#'   \code{\link{create_rbf_object}}.
#' @param train_ind Optional vector containing the indices for the train set.
#' @param train_perc Optional parameter for defining the percentage of the
#'   dataset to be used for training set, the remaining will be the test set.
#' @param is_summary Logical, print the summary statistics.
#' @inheritParams bpr_optimize
#'
#' @return A 'bpr_diff_predict' object which, in addition to the input
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
bpr_diff_test_data_wrap <- function(formula = NULL, x, y, model_name = "svm",
                               w = NULL, basis = NULL, train_ind = NULL,
                               train_perc = 0.7, fit_feature = "RMSE",
                               cpg_dens_feat = TRUE, opt_method = "CG",
                               opt_itnmax = 100, is_parallel = TRUE,
                               no_cores = NULL, is_summary = TRUE){

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
    yy <- list()
    for (i in 1:length(x$treatment)){
        yy[[i]] <- rbind(x$control[[i]], x$treatment[[i]])
        Order <- base::order(yy[[i]][,1])
        yy[[i]] <- yy[[i]][Order, ]
    }

    # Learn methylation profiles for control samples
    message("Learning control methylation profiles ...\n")
    out_opt_conc <- bpr_optim(x           = yy,
                              w           = w,
                              basis       = basis,
                              fit_feature = "NLL",
                              cpg_dens_feat = FALSE,
                              opt_method  = opt_method,
                              opt_itnmax  = opt_itnmax,
                              is_parallel = is_parallel,
                              no_cores    = no_cores)

    # NLL for null hypothesis
    nll_null_conc <- out_opt_conc$W_opt[, params + 1]

    ##------------------------------------
    # Compute log likelihood ratio test
    ##------------------------------------
    lr_test_conc <- 2 * (nll_null_conc - nll_alt)

    message("Done!\n\n")

    # Create 'bpr_predict' object
    obj <- structure(list(formula      = formula,
                          model_name   = model_name,
                          basis        = out_treat_opt$basis,
                          #train_ind    = dataset$train_ind,
                          train_perc   = train_perc,
                          fit_feature  = fit_feature,
                          cpg_dens_feat = cpg_dens_feat,
                          opt_method   = opt_method,
                          opt_itnmax   = opt_itnmax,
                          #W_opt_conc   = W_diff,
                          W_opt_contr  = out_contr_opt$W_opt,
                          W_opt_treat  = out_treat_opt$W_opt,
                          Mus          = list(control = out_contr_opt$Mus,
                                              treatment = out_treat_opt$Mus)
                          #train        = dataset$train,
                          #test         = dataset$test,
                          #gex_model    = train_model$gex_model,
                          #train_pred   = train_model$train_pred,
                          #test_pred    = predictions$test_pred,
                          #train_errors = train_model$train_errors,
                          #test_errors  = predictions$test_errors,
                          #train_conf_mat = train_model$confusion_matrix,
                          #test_conf_mat  = predictions$confusion_matrix
                          ),
                     class = "bpr_diff_classify")
    return(obj)














    my_bpr_likelihood <- function(w, H, data, is_NLL = FALSE){
        total <- data[, 1]
        succ  <- data[, 2]

        # Predictions of the target variables
        # Compute the cdf of N(0,1) distribution (i.e. probit function)
        Phi <- pnorm(H %*% w)

        # In extreme cases where probit is 0 or 1, subtract a tiny number
        # so we can evaluate the log(0) when computing the Binomial
        Phi[which(Phi > (1 - 1e-289))] <- 1 - 1e-289
        Phi[which(Phi < 1e-289)] <- 1e-289

        # Compute the log likelihood
        res <- sum(dbinom(x = succ, size = total, prob = Phi, log = TRUE))

        # If we required the Negative Log Likelihood
        if (is_NLL){
            res <- (-1) * res
        }
        return(res)
    }









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
    nll_contr_alt <- vector(mode = "numeric", length = length(x$treatment))
    nll_treat_alt <- vector(mode = "numeric", length = length(x$treatment))
    for (i in 1:length(x$treatment)){
        # Create design matrix H
        H <- .design_matrix.rbf(x   = basis,
                            obs = x$control[[i]][,1])$H
        # Evaluate the likelihood under control parameters
        nll_contr_alt[i] <- my_bpr_likelihood(w = out_contr_opt$W_opt[i, 1:params],
                                              H = H,
                                              data = x$control[[i]][,2:3],
                                              is_NLL = TRUE)

        # Create design matrix H
        H <- .design_matrix.rbf(x   = basis,
                            obs = x$treatment[[i]][,1])$H
        # Evaluate the likelihood under control parameters
        nll_treat_alt[i] <- my_bpr_likelihood(w = out_treat_opt$W_opt[i, 1:params],
                                              H = H,
                                              data = x$treatment[[i]][,2:3],
                                              is_NLL = TRUE)

    }
    # NLL for alternative hypothesis
    nll_alt <- nll_contr_alt + nll_treat_alt


    ##-------------------------------------
    #      NULL hypothesis
    ##-------------------------------------
    nll_contr_null <- nll_contr_alt
    nll_treat_null <- vector(mode = "numeric", length = length(x$treatment))
    for (i in 1:length(x$treatment)){
        # Create design matrix H
        H <- .design_matrix.rbf(x   = basis,
                            obs = x$treatment[[i]][,1])$H

        # Evaluate the likelihood under control parameters
        nll_treat_null[i] <- my_bpr_likelihood(w = out_contr_opt$W_opt[i, 1:params],
                                               H = H,
                                               data = x$treatment[[i]][,2:3],
                                               is_NLL = TRUE)
    }
    # NLL for null hypothesis
    nll_null <- nll_contr_null + nll_treat_null

    ##------------------------------------
    # Compute log likelihood ratio test
    ##------------------------------------
    lr_test <- 2 * (nll_null - nll_alt)













    ##-------------------------------------
    #      NULL hypothesis
    ##-------------------------------------
    yy <- list()
    for (i in 1:length(x$treatment)){
        yy[[i]] <- rbind(x$control[[i]], x$treatment[[i]])
        Order <- base::order(yy[[i]][,1])
        yy[[i]] <- yy[[i]][Order, ]
    }

    # Learn methylation profiles for control samples
    message("Learning control methylation profiles ...\n")
    out_opt_conc <- bpr_optim(x           = yy,
                              w           = w,
                              basis       = basis,
                              fit_feature = "NLL",
                              cpg_dens_feat = FALSE,
                              opt_method  = opt_method,
                              opt_itnmax  = opt_itnmax,
                              is_parallel = is_parallel,
                              no_cores    = no_cores)

    # NLL for null hypothesis
    nll_null_conc <- vector(mode = "numeric", length = length(x$treatment))
    for (i in 1:length(yy)){
        # Create design matrix H
        H <- .design_matrix.rbf(x   = basis,
                            obs = yy[[i]][,1])$H

        # Evaluate the likelihood under control parameters
        nll_null_conc[i] <- my_bpr_likelihood(w = out_opt_conc$W_opt[i, 1:params],
                                              H = H,
                                              data = yy[[i]][,2:3],
                                              is_NLL = TRUE)
    }

    ##------------------------------------
    # Compute log likelihood ratio test
    ##------------------------------------
    lr_test_conc <- 2 * (nll_null_conc - nll_alt)






    ##-------------------------------------
    #      NULL hypothesis
    ##-------------------------------------
    zz <- list()
    for (i in 1:length(x$treatment)){
        zz[[i]] <- cbind(x$control[[i]][,1],
                         x$control[[i]][,2] + x$treatment[[i]][,2],
                         x$control[[i]][,3] + x$treatment[[i]][,3])
    }

    # Learn methylation profiles for control samples
    message("Learning control methylation profiles ...\n")
    out_opt_sum <- bpr_optim(x           = zz,
                             w           = w,
                             basis       = basis,
                             fit_feature = "NLL",
                             cpg_dens_feat = FALSE,
                             opt_method  = opt_method,
                             opt_itnmax  = opt_itnmax,
                             is_parallel = is_parallel,
                             no_cores    = no_cores)

    # NLL for null hypothesis
    nll_sum_contr_null <- vector(mode = "numeric", length = length(x$treatment))
    nll_sum_treat_null <- vector(mode = "numeric", length = length(x$treatment))
    for (i in 1:length(zz)){
        # Create design matrix H
        H <- .design_matrix.rbf(x   = basis,
                            obs = zz[[i]][,1])$H

        # Evaluate the likelihood under control parameters
        nll_sum_contr_null[i] <- my_bpr_likelihood(w = out_opt_sum$W_opt[i, 1:params],
                                                   H = H,
                                                   data = x$control[[i]][,2:3],
                                                   is_NLL = TRUE)

        nll_sum_treat_null[i] <- my_bpr_likelihood(w = out_opt_sum$W_opt[i, 1:params],
                                                   H = H,
                                                   data = x$treatment[[i]][,2:3],
                                                   is_NLL = TRUE)
    }
    # NLL for null hypothesis
    nll_sum_null <- nll_sum_contr_null + nll_sum_treat_null

    ##------------------------------------
    # Compute log likelihood ratio test
    ##------------------------------------
    lr_test_sum <- 2 * (nll_sum_null - nll_alt)



}
