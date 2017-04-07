#' @name design_matrix
#' @rdname design_matrix
#' @aliases designmatrix, des_matrix, des_mat
#'
#' @title Generic function for creating design matrices
#'
#' @description These functions call the appropriate methods depending on the
#'   class of the object \code{x} to create RBF, polynomial or Fourier design
#'   matrices.
#'
#' @param x A basis function object.
#' @param obs A vector of observations.
#' @param ... Additional parameters.
#'
#' @return A design matrix object
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{create_polynomial_object}},
#'   \code{\link{create_rbf_object}},  \code{\link{create_fourier_object}}
#'
#' @export
NULL

#' @rdname design_matrix
#'
#' @examples
#' obj <- create_polynomial_object(M=2)
#' obs <- c(0,.2,.5)
#' polyn <- design_matrix(obj, obs)
#'
#' #----------------
#'
#' obj <- create_rbf_object(M=2)
#' obs <- c(0,.2,.5)
#' rbf <- design_matrix(obj, obs)
#'
#' @export
design_matrix <- function(x, ...){
    UseMethod("design_matrix")
}


#' @rdname design_matrix
#'
design_matrix.default <- function(x, ...){
    stop(paste("Object type '", class(x), "' is not implemented.", sep = ""))
}


#' @rdname design_matrix
#'
#' @export
design_matrix.polynomial <- function(x, obs, ...){
    assertthat::assert_that(methods::is(x, "polynomial"))
    assertthat::assert_that(is.vector(obs))

    N <- length(obs)  # Length of the dataset
    H <- matrix(1, nrow = N, ncol = x$M + 1)
    if (x$M > 0){
        for (j in 1:x$M){
            H[, j + 1] <- .polynomial_basis(obs, j)  # Compute X^(j)
        }
    }
    return(list(H = H, basis = x))
}


#' @rdname design_matrix
#'
#' @export
design_matrix.rbf <- function(x, obs, ...){
    assertthat::assert_that(methods::is(x, "rbf"))
    assertthat::assert_that(is.vector(obs))

    N   <- length(obs)  # Length of the dataset
    if (x$M == 0){
        H <- matrix(1, nrow = N, ncol = 1)
        x$mus <- 0
    }else{
        if (is.null(x$mus)){
            if (x$eq_spaced_mus){
                x$mus <- vector(mode = "numeric", x$M)
                if (! x$whole_region){
                    # TODO: Should this be deleted?
                    for (i in 1:x$M){
                        x$mus[i] <- i * (max(obs) - min(obs)) /
                            (x$M + 1) + min(obs)
                    }
                }
            }else{
                repeat{
                    # TODO: Should this be deleted?
                    km <- stats::kmeans(obs, x$M, iter.max = 30, nstart = 10)
                    if (min(km$size) > 0){
                        break  # Only accept non-empty clusters
                    }
                }
                x$mus <- km$centers  # RBF centers
            }
        }
        # Convert the 'obs' vector to an N x 1 dimensional matrix
        obs <- as.matrix(obs)
        H <- matrix(1, nrow = N, ncol = x$M + 1)
        for (j in 1:x$M){
            H[, j + 1] <- apply(obs, 1, .rbf_basis,
                                mus = x$mus[j], gamma = x$gamma)
        }
    }
    return(list(H = H, basis = x))
}


#' @rdname design_matrix
#'
#' @export
design_matrix.fourier <- function(x, obs, ...){
    assertthat::assert_that(methods::is(x, "fourier"))
    assertthat::assert_that(is.vector(obs))
    # Similar implementation to the FDA package.
    # Compute base frequency
    omega <- 2 * pi / x$period
    # Initialize design matrix
    H <- matrix(1 / sqrt(2), nrow = length(obs), ncol = x$M + 1)
    if (x$M > 1){
        # Get the basis
        j <- seq(2, x$M, 2)
        k <- j / 2
        # Compute outer product
        evals <- outer(omega * obs, k)
        # Make the sine evaluations
        H[, j] <- sin(evals)
        # Make the cosine evaluations
        H[, j + 1] <- cos(evals)
    }
    # Divide by this constant to get actual magnitude
    # TODO: Understand better this step
    H <- H / sqrt(x$period / 2)
    return(list(H = H, basis = x))
}
