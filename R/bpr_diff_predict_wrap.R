#' Predict differential gene expression from differential methylation profiles
#'
#' \code{bpr_diff_predict_wrap} is a function that wraps all the necessary
#' subroutines for performing prediction of differential gene expression levels.
#' Initially, it optimizes the parameters of the basis functions so as to learn
#' the methylation profiles for the control and the treatment samples Then, the
#' two learned methylation profiles are concatenated to keep all coefficients
#' for both profiles. Then the learned parameters / coefficients of the basis
#' functions are given as input features for performing regression in order to
#' predict the corresponding differential (fold-change) gene expression levels.
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
#' @examples
#' obs <- meth_data
#' y   <- gex_data
#' basis <- create_rbf_object(M = 5)
#' out   <- bpr_predict_wrap(x = obs, y = y, basis = basis,
#'                           is_parallel = FALSE, opt_itnmax = 10)
#'
#' @export
bpr_diff_predict_wrap <- function(formula = NULL, x, y, model_name = "svm",
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
                             fit_feature = fit_feature,
                             cpg_dens_feat = cpg_dens_feat,
                             opt_method  = opt_method,
                             opt_itnmax  = opt_itnmax,
                             is_parallel = is_parallel,
                             no_cores    = no_cores)

  # Learn methylation profiles for treatment samples
  message("Learning treatment methylation profiles ...\n")
  out_treat_opt <- bpr_optim(x           = x$treatment,
                             w           = w,
                             basis       = basis,
                             fit_feature = fit_feature,
                             cpg_dens_feat = cpg_dens_feat,
                             opt_method  = opt_method,
                             opt_itnmax  = opt_itnmax,
                             is_parallel = is_parallel,
                             no_cores    = no_cores)

  # Compute fold change gene expression levels
  y_diff <- gtools::foldchange(y$control, y$treatment)

  # Concatenate coefficients from both samples
  W_diff <- cbind(out_contr_opt$W_opt, out_treat_opt$W_opt)
  colnames(W_diff) <- NULL

  # Create training and test sets
  message("Partitioning to test and train data ...\n")
  dataset <- .partition_data(x          = W_diff,
                             y          = y_diff,
                             train_ind  = train_ind,
                             train_perc = train_perc)

  # Train regression model from methylation profiles
  message("Training linear regression model ...\n")
  train_model <- train_model_gex(formula    = formula,
                                 model_name = model_name,
                                 train      = dataset$train,
                                 is_summary = is_summary)

  # Predict gene expression from methylation profiles
  message("Making predictions ...\n")
  predictions <- predict_model_gex(model      = train_model$gex_model,
                                   test       = dataset$test,
                                   is_summary = is_summary)
  message("Done!\n\n")

  # Create 'bpr_predict' object
  obj <- structure(list(formula      = formula,
                        model_name   = model_name,
                        basis        = out_treat_opt$basis,
                        train_ind    = dataset$train_ind,
                        train_perc   = train_perc,
                        fit_feature  = fit_feature,
                        cpg_dens_feat = cpg_dens_feat,
                        opt_method   = opt_method,
                        opt_itnmax   = opt_itnmax,
                        W_opt_conc   = W_diff,
                        W_opt_contr  = out_contr_opt$W_opt,
                        W_opt_treat  = out_treat_opt$W_opt,
                        Mus          = list(control = out_contr_opt$Mus,
                                            treatment = out_treat_opt$Mus),
                        train        = dataset$train,
                        test         = dataset$test,
                        gex_model    = train_model$gex_model,
                        train_pred   = train_model$train_pred,
                        test_pred    = predictions$test_pred,
                        train_errors = train_model$train_errors,
                        test_errors  = predictions$test_errors),
                   class = "bpr_diff_predict")
  return(obj)
}
