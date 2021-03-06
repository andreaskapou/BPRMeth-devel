#' Fitting linear models using Basis Functions
#'
#' \code{blr} is used to fit linear models using basis functions such as
#' Radial Basis Functions (RBFs), Fourier and Polynomial basis functions.
#'
#' @param x The observations.
#' @param y The response.
#' @param basis Basis function object e.g. \code{\link{create_rbf_object}}.
#' @param lambda Optional parameter for performing ridge regression.
#' @param return.all Optional logical, indicating if all the metrics should be
#'  computed (mainly for efficiency).
#'
#' @return An object of class "basis_lm" is a list containing the following
#'  components:
#' \itemize{
#'  \item \code{coefficients} a named vector of coefficients.
#'  \item \code{residuals} the residuals, that is response minus fitted values.
#'  \item \code{fitted.values} the fitted mean values.
#'  \item \code{df.residuals} the residual degrees of freedom.
#'  \item \code{sigma} the standard deviation of the residuals.
#'  \item \code{vcov} The covariance matrix.
#'  \item \code{basis} The basis object used.
#'  \item \code{lambda} The regularization parameter.
#'  \item \code{call} the matched call.
#' }
#'
#' @seealso \code{\link{predict.blr}}, \code{\link{create_fourier_object}},
#'  \code{\link{create_rbf_object}}, \code{\link{summary.blr}}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
blr <- function(x, y, basis, lambda = 0, return.all = TRUE){
  est <- list(basis = basis, lambda = lambda)
  # Create the design matrix
  des_mat <- design_matrix(x = basis, obs = x)
  H <- des_mat$H
  if (lambda == 0){
    # Compute QR decomposition of H
    qx <- qr(H)
    #  Compute (H'H)^(-1)H'y
    est$coefficients <- solve.qr(qx, y)
  }else{
    I <- diag(1, NCOL(H))  # Identity matrix
    # TODO: Should we do this or not for the bias term??
    # I[1,1]  <- 1e-10  # Do not change the intercept coefficient
    # Compute QR decomposition
    qx <- qr(lambda * I + t(H) %*% H)
    # Compute (lambda * I + H'H)^(-1)H'y
    est$coefficients <- solve.qr(qx, t(H) %*% y)
  }

  # Standard deviation of residuals
  est$sigma <- sqrt(sum(est$residuals ^ 2) / est$df.residuals)
  if (return.all){
    # Fitted values
    est$fitted.values <- as.vector(H %*% est$coefficients)
    # Residuals
    est$residuals <- y - est$fitted.values
    # Degrees of freedom
    est$df.residuals <- NROW(H) - NCOL(H)
    # Compute sigma^2 * (x'x)^(-1)
    est$vcov <- est$sigma ^ 2 * chol2inv(qx$qr)
    colnames(est$vcov) <- rownames(est$vcov) <- colnames(H)
  }
  est$call <- match.call()
  class(est) <- "blr"

  return(est)
}


#' Make predictions using the Basis Linear Model
#'
#' S3 method for class blr, that makes predictions for \code{newdata} using
#' the Basis Linear Model. This funciton is similar to \code{\link[stats]{lm}},
#' however currently it does not provide the \code{\link[stats]{formula}}
#' functionality.
#'
#' @param object Object of class \code{\link{blr}}.
#' @param newdata An optional data frame in which to look for variables with
#'  which to predict. If omitted, the fitted values are used.
#' @param ... Optional additional parameters.
#'
#' @return \code{predict.blr} produces a vector of predictions.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
predict.blr <- function(object, newdata = NULL, ...){
  if(is.null(newdata)){
    y <- stats::fitted(object)
  }
  else{
    # Create the design matrix
    H <- design_matrix(x = object$basis, obs = newdata)$H
    y <- as.vector(H %*% stats::coef(object))
  }
  return(y)
}


#' Print the output of the Basis Linear Model
#'
#' S3 method for class basis_lm, that prints the output of the Basis Linear Model
#' in a similar way to \code{\link[stats]{lm}}.
#'
#' @param x Object of class \code{\link{blr}}.
#' @param ... Optional additional parameters.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
print.blr <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(as.vector(x$coefficients))
}


#' Summary output of the Basis Linear Model
#'
#' S3 method for class blr, that computes the summary of the Basis Linear
#' Model in a similar way to \code{\link[stats]{lm}}.
#'
#' @param object Object of class \code{\link{blr}}.
#' @param ... Optional additional parameters.
#'
#' @return A summary.blr object
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @method summary blr
#' @export
summary.blr <- function(object, ...){
  se <- sqrt(diag(object$vcov))
  tval <- stats::coef(object) / se

  TAB <- cbind(Estimate = stats::coef(object),
               StdErr = se,
               t.value = tval,
               p.value = 2 * stats::pt(-abs(tval), df = object$df))
  colnames(TAB) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
  res <- structure(list(call = object$call,
                        coefficients = TAB),
                   class = "summary.blr")
  return(res)
}


#' Print summary output of the Basis Linear Model
#'
#' S3 method for class summary.blr, that prints the summary of the Basis
#' Linear Model in a similar way to \code{\link[stats]{lm}}.
#'
#' @param x Object of class \code{\link{summary.blr}}.
#' @param ... Optional additional parameters.
#'
#' @export
print.summary.blr <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
}
