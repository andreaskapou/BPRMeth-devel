create_meth_data <- function(N = 300, pi.c = c(0.45, 0.35, 0.2), max_L = 25,
                            xmin = -100, xmax=100, fmin = -1, fmax = 1){
  set.seed(3)
  # Create a list to store data for each methylation region
  X       <- list()
  # A vector for storing corresponding gene expression data
  Y       <- vector(mode = "numeric", length = N)

  # For each of the N objects
  for (i in 1:N){
    # L is the number of CpGs found in the ith region
    L <- rbinom(n = 1, size = max_L, prob = .8)
    X[[i]] <- matrix(0, nrow = L, ncol = 3)
    # Randomly sample locations for the CpGs
    obs <- sort(sample(xmin:xmax, L))
    # Scale them, so the data lie in the (fmin, fmax) range
    X[[i]][ ,1] <- BPRMeth:::.minmax_scaling(data = obs,
                                  xmin = xmin,
                                  xmax = xmax,
                                  fmin = fmin,
                                  fmax = fmax)

    if (i < N * pi.c[1]){   # First methylation profile
      lb <- round(L / 4)

      X[[i]][1:lb,2] <- rbinom(lb, 20, .9)
      repeat{
        X[[i]][1:lb,3] <- rbinom(lb, 14, .9)
        if(all(X[[i]][1:lb,2] > X[[i]][1:lb,3]))
          break
      }

      X[[i]][(lb + 1):L,2] <- rbinom(L - lb, 20, .9)
      repeat{
        X[[i]][(lb + 1):L,3] <- rbinom(L - lb, 2, .9)
        if (all(X[[i]][(lb + 1):L,2] > X[[i]][(lb + 1):L,3]))
          break
      }
      Y[i] <- rpois(1, lambda=200)
    }else if (i < (N * pi.c[2] + N * pi.c[1])){ # Second methylation profile
      lb <- round(L / 1.5)

      X[[i]][1:lb,2] <- rbinom(lb, 20, .9)
      repeat{
        X[[i]][1:lb,3] <- rbinom(lb, 2, .8)
        if(all(X[[i]][1:lb,2] > X[[i]][1:lb,3]))
          break
      }

      X[[i]][(lb + 1):L,2] <- rbinom(L - lb, 20, .9)
      repeat{
        X[[i]][(lb + 1):L,3] <- rbinom(L-lb, 14, .9)
        if (all(X[[i]][(lb + 1):L,2] > X[[i]][(lb + 1):L,3]))
          break
      }
      Y[i] <- rpois(1, lambda=100)
    }else{                  # Third methylation profile
      lb <- round(L / 2.5)
      mb <- round(L / 3.5)

      X[[i]][1:lb,2] <- rbinom(lb, 20, .9)
      repeat{
        X[[i]][1:lb,3] <- rbinom(lb, 2, .9)
        if(all(X[[i]][1:lb,2] > X[[i]][1:lb,3]))
          break
      }

      X[[i]][(lb + 1):(lb + mb),2] <- rbinom(mb, 20, .9)
      repeat{
        X[[i]][(lb + 1):(lb + mb),3] <- rbinom(mb, 14, .9)
        if (all(X[[i]][(lb + 1):(lb + mb),2] > X[[i]][(lb + 1):(lb + mb),3]))
          break
      }

      X[[i]][(lb + 1 + mb):L,2] <- rbinom(L - mb - lb, 20, .9)
      repeat{
        X[[i]][(lb + 1 + mb):L,3] <- rbinom(L - mb - lb, 2, .9)
        if (all(X[[i]][(lb + 1 + mb):L,2] > X[[i]][(lb + 1 + mb):L],3))
          break
      }
      Y[i] <- rpois(1, lambda=50)
    }
  }
  return(list(X = X, Y = Y))
}

set.seed(1)
bpr <- create_meth_data(N=600)
meth_data <- bpr$X
gex_data <- bpr$Y
devtools::use_data(meth_data, gex_data, overwrite = TRUE)




create_lm_data <- function(N = 300, pi.c = c(0.45, 0.35, 0.2), max_L = 25,
                           true_w = list(w1 = c(-1.1, -2, 2.7, -2),
                                         w2 = c(2, 1, -2, 1),
                                         w3 = c(0.4, -0.7, 1.7, 2.8)),
                             xmin = -100, xmax=100, fmin = -1, fmax = 1){
    set.seed(123)
    # Create a list to store data for each methylation region
    X       <- list()

    # For each of the N objects
    for (i in 1:N){
        # L is the number of CpGs found in the ith region
        L <- rbinom(n = 1, size = max_L, prob = .8)
        X[[i]] <- matrix(0, nrow = L, ncol = 2)
        # Randomly sample locations for the CpGs
        obs <- sort(sample(xmin:xmax, L))
        # Scale them, so the data lie in the (fmin, fmax) range
        X[[i]][, 1] <- BPRMeth:::.minmax_scaling(data = obs,
                                                 xmin = xmin,
                                                 xmax = xmax,
                                                 fmin = fmin,
                                                 fmax = fmax)

        # Choose function to generate the data
        if (i < N * pi.c[1])
            cl_label <- 1
        else if (i < (N * pi.c[2] + N * pi.c[1]))
            cl_label <- 2
        else
            cl_label <- 3
        # Randomly sample a cluster to generate data from
        # cl_label <- sample(x = 3, size = 1, prob = pi.c)

        # Set corresponding weights
        w <- true_w[[cl_label]]
        # Create basis function object
        basis <- create_rbf_object(M = 3)
        # Create design matrix
        H <- BPRMeth::design_matrix(basis, X[[i]][ ,1])$H
        # Compute mean of data
        mu <- H %*% w
        # Randomly generate data from this function
        X[[i]][, 2] <- rnorm(L, mu, 0.1)
    }
    return(X)
}

set.seed(1)
lm_out <- create_lm_data(N=600)
lm_data <- lm_out
devtools::use_data(lm_data, overwrite = TRUE)





create_meth_data_2 <- function(N = 300, pi.c = c(0.45, 0.35, 0.2), max_L = 25,
                             xmin = -100, xmax=100, fmin = -1, fmax = 1){
    set.seed(3)
    # Create a list to store data for each methylation region
    #X       <- list()
    # A vector for storing corresponding gene expression data
    Y       <- vector(mode = "numeric", length = N)

    X <- list(control = list(), treatment = list())

    L <- list()
    for (i in 1:N){
        # L is the number of CpGs found in the ith region
        L[[i]] <- rbinom(n = 1, size = max_L, prob = .8)
    }
    # For each of the N objects
    for (i in 1:N){
        # L is the number of CpGs found in the ith region
        #L <- rbinom(n = 1, size = max_L, prob = .8)
        X$control[[i]] <- matrix(0, nrow = L[[i]], ncol = 3)
        # Randomly sample locations for the CpGs
        obs <- sort(sample(xmin:xmax, L[[i]]))
        # Scale them, so the data lie in the (fmin, fmax) range
        X$control[[i]][ ,1] <- BPRMeth:::.minmax_scaling(data = obs,
                                                         xmin = xmin,
                                                         xmax = xmax,
                                                         fmin = fmin,
                                                         fmax = fmax)
        X$treatment[[i]] <- X$control[[i]]
    }
    for (i in 1:N){
        if (i < N * pi.c[1]){   # First methylation profile
            lb <- round(L[[i]] / 4)

            X$control[[i]][1:lb,2] <- rbinom(lb, 20, .9)
            repeat{
                X$control[[i]][1:lb,3] <- rbinom(lb, 14, .9)
                if(all(X$control[[i]][1:lb,2] > X$control[[i]][1:lb,3]))
                    break
            }

            X$control[[i]][(lb + 1):L[[i]],2] <- rbinom(L[[i]] - lb, 20, .9)
            repeat{
                X$control[[i]][(lb + 1):L[[i]],3] <- rbinom(L[[i]] - lb, 2, .9)
                if (all(X$control[[i]][(lb + 1):L[[i]],2] > X$control[[i]][(lb + 1):L[[i]],3]))
                    break
            }
            Y[i] <- rpois(1, lambda=200)
        }else if (i < (N * pi.c[2] + N * pi.c[1])){ # Second methylation profile
            lb <- round(L[[i]] / 1.5)

            X$control[[i]][1:lb,2] <- rbinom(lb, 20, .9)
            repeat{
                X$control[[i]][1:lb,3] <- rbinom(lb, 2, .8)
                if(all(X$control[[i]][1:lb,2] > X$control[[i]][1:lb,3]))
                    break
            }

            X$control[[i]][(lb + 1):L[[i]],2] <- rbinom(L[[i]] - lb, 20, .9)
            repeat{
                X$control[[i]][(lb + 1):L[[i]],3] <- rbinom(L[[i]]-lb, 14, .9)
                if (all(X$control[[i]][(lb + 1):L[[i]],2] > X$control[[i]][(lb + 1):L[[i]],3]))
                    break
            }
            Y[i] <- rpois(1, lambda=100)
        }else{                  # Third methylation profile
            lb <- round(L[[i]] / 2.5)
            mb <- round(L[[i]] / 3.5)

            X$control[[i]][1:lb,2] <- rbinom(lb, 20, .9)
            repeat{
                X$control[[i]][1:lb,3] <- rbinom(lb, 2, .9)
                if(all(X$control[[i]][1:lb,2] > X$control[[i]][1:lb,3]))
                    break
            }

            X$control[[i]][(lb + 1):(lb + mb),2] <- rbinom(mb, 20, .9)
            repeat{
                X$control[[i]][(lb + 1):(lb + mb),3] <- rbinom(mb, 14, .9)
                if (all(X$control[[i]][(lb + 1):(lb + mb),2] > X$control[[i]][(lb + 1):(lb + mb),3]))
                    break
            }

            X$control[[i]][(lb + 1 + mb):L[[i]],2] <- rbinom(L[[i]] - mb - lb, 20, .9)
            repeat{
                X$control[[i]][(lb + 1 + mb):L[[i]],3] <- rbinom(L[[i]] - mb - lb, 2, .9)
                if (all(X$control[[i]][(lb + 1 + mb):L[[i]],2] > X$control[[i]][(lb + 1 + mb):L[[i]]],3))
                    break
            }
            Y[i] <- rpois(1, lambda=50)
        }
    }


    for (i in 1:N){
        if (i < N * pi.c[1]){   # First methylation profile
            lb <- round(L[[i]] / 4)

            X$treatment[[i]][1:lb,2] <- rbinom(lb, 20, .9)
            repeat{
                X$treatment[[i]][1:lb,3] <- rbinom(lb, 14, .9)
                if(all(X$treatment[[i]][1:lb,2] > X$treatment[[i]][1:lb,3]))
                    break
            }

            X$treatment[[i]][(lb + 1):L[[i]],2] <- rbinom(L[[i]] - lb, 20, .9)
            repeat{
                X$treatment[[i]][(lb + 1):L[[i]],3] <- rbinom(L[[i]] - lb, 2, .9)
                if (all(X$treatment[[i]][(lb + 1):L[[i]],2] > X$treatment[[i]][(lb + 1):L[[i]],3]))
                    break
            }
            Y[i] <- rpois(1, lambda=200)
        }else if (i < (N * pi.c[2] + N * pi.c[1])){ # Second methylation profile
            lb <- round(L[[i]] / 1.5)

            X$treatment[[i]][1:lb,2] <- rbinom(lb, 20, .9)
            repeat{
                X$treatment[[i]][1:lb,3] <- rbinom(lb, 2, .8)
                if(all(X$treatment[[i]][1:lb,2] > X$treatment[[i]][1:lb,3]))
                    break
            }

            X$treatment[[i]][(lb + 1):L[[i]],2] <- rbinom(L[[i]] - lb, 20, .9)
            repeat{
                X$treatment[[i]][(lb + 1):L[[i]],3] <- rbinom(L[[i]]-lb, 14, .9)
                if (all(X$treatment[[i]][(lb + 1):L[[i]],2] > X$treatment[[i]][(lb + 1):L[[i]],3]))
                    break
            }
            Y[i] <- rpois(1, lambda=100)
        }else{                  # Third methylation profile
            lb <- round(L[[i]] / 2.5)
            mb <- round(L[[i]] / 3.5)

            X$treatment[[i]][1:lb,2] <- rbinom(lb, 20, .9)
            repeat{
                X$treatment[[i]][1:lb,3] <- rbinom(lb, 2, .9)
                if(all(X$treatment[[i]][1:lb,2] > X$treatment[[i]][1:lb,3]))
                    break
            }

            X$treatment[[i]][(lb + 1):(lb + mb),2] <- rbinom(mb, 20, .9)
            repeat{
                X$treatment[[i]][(lb + 1):(lb + mb),3] <- rbinom(mb, 14, .9)
                if (all(X$treatment[[i]][(lb + 1):(lb + mb),2] > X$treatment[[i]][(lb + 1):(lb + mb),3]))
                    break
            }

            X$treatment[[i]][(lb + 1 + mb):L[[i]],2] <- rbinom(L[[i]] - mb - lb, 20, .9)
            repeat{
                X$treatment[[i]][(lb + 1 + mb):L[[i]],3] <- rbinom(L[[i]] - mb - lb, 2, .9)
                if (all(X$treatment[[i]][(lb + 1 + mb):L[[i]],2] > X$treatment[[i]][(lb + 1 + mb):L[[i]]],3))
                    break
            }
            Y[i] <- rpois(1, lambda=50)
        }
    }

    return(list(X = X, Y = Y))
}




create_meth_data_3 <- function(N = 300, pi.c = c(0.45, 0.35, 0.2), max_L = 25,
                               xmin = -100, xmax=100, fmin = -1, fmax = 1){
    set.seed(3)
    # Create a list to store data for each methylation region
    #X       <- list()
    # A vector for storing corresponding gene expression data
    Y       <- vector(mode = "numeric", length = N)

    X <- list(control = list(), treatment = list())

    L <- list()
    for (i in 1:N){
        # L is the number of CpGs found in the ith region
        L[[i]] <- rbinom(n = 1, size = max_L, prob = .8)
    }
    # For each of the N objects
    for (i in 1:N){
        # L is the number of CpGs found in the ith region
        #L <- rbinom(n = 1, size = max_L, prob = .8)
        X$control[[i]] <- matrix(0, nrow = L[[i]], ncol = 3)
        # Randomly sample locations for the CpGs
        obs <- sort(sample(xmin:xmax, L[[i]]))
        # Scale them, so the data lie in the (fmin, fmax) range
        X$control[[i]][ ,1] <- BPRMeth:::.minmax_scaling(data = obs,
                                                         xmin = xmin,
                                                         xmax = xmax,
                                                         fmin = fmin,
                                                         fmax = fmax)
        X$treatment[[i]] <- X$control[[i]]
    }
    for (i in 1:N){
        if (i < N * pi.c[1]){   # First methylation profile
            lb <- round(L[[i]] / 4)

            X$control[[i]][1:lb,2] <- rbinom(lb, 20, .9)
            repeat{
                X$control[[i]][1:lb,3] <- rbinom(lb, 14, .9)
                if(all(X$control[[i]][1:lb,2] > X$control[[i]][1:lb,3]))
                    break
            }

            X$control[[i]][(lb + 1):L[[i]],2] <- rbinom(L[[i]] - lb, 20, .9)
            repeat{
                X$control[[i]][(lb + 1):L[[i]],3] <- rbinom(L[[i]] - lb, 2, .9)
                if (all(X$control[[i]][(lb + 1):L[[i]],2] > X$control[[i]][(lb + 1):L[[i]],3]))
                    break
            }
            Y[i] <- rpois(1, lambda=200)
        }else if (i < (N * pi.c[2] + N * pi.c[1])){ # Second methylation profile
            lb <- round(L[[i]] / 1.5)

            X$control[[i]][1:lb,2] <- rbinom(lb, 20, .9)
            repeat{
                X$control[[i]][1:lb,3] <- rbinom(lb, 2, .8)
                if(all(X$control[[i]][1:lb,2] > X$control[[i]][1:lb,3]))
                    break
            }

            X$control[[i]][(lb + 1):L[[i]],2] <- rbinom(L[[i]] - lb, 20, .9)
            repeat{
                X$control[[i]][(lb + 1):L[[i]],3] <- rbinom(L[[i]]-lb, 14, .9)
                if (all(X$control[[i]][(lb + 1):L[[i]],2] > X$control[[i]][(lb + 1):L[[i]],3]))
                    break
            }
            Y[i] <- rpois(1, lambda=100)
        }else{                  # Third methylation profile
            lb <- round(L[[i]] / 2.5)
            mb <- round(L[[i]] / 3.5)

            X$control[[i]][1:lb,2] <- rbinom(lb, 20, .9)
            repeat{
                X$control[[i]][1:lb,3] <- rbinom(lb, 2, .9)
                if(all(X$control[[i]][1:lb,2] > X$control[[i]][1:lb,3]))
                    break
            }

            X$control[[i]][(lb + 1):(lb + mb),2] <- rbinom(mb, 20, .9)
            repeat{
                X$control[[i]][(lb + 1):(lb + mb),3] <- rbinom(mb, 14, .9)
                if (all(X$control[[i]][(lb + 1):(lb + mb),2] > X$control[[i]][(lb + 1):(lb + mb),3]))
                    break
            }

            X$control[[i]][(lb + 1 + mb):L[[i]],2] <- rbinom(L[[i]] - mb - lb, 20, .9)
            repeat{
                X$control[[i]][(lb + 1 + mb):L[[i]],3] <- rbinom(L[[i]] - mb - lb, 2, .9)
                if (all(X$control[[i]][(lb + 1 + mb):L[[i]],2] > X$control[[i]][(lb + 1 + mb):L[[i]]],3))
                    break
            }
            Y[i] <- rpois(1, lambda=50)
        }
    }


    for (i in 1:N){
        if (i < N * pi.c[1]){   # First methylation profile
            lb <- round(L[[i]] / 4)

            X$treatment[[i]][1:lb,2] <- rbinom(lb, 20, .9)
            repeat{
                X$treatment[[i]][1:lb,3] <- rbinom(lb, 14, .9)
                if(all(X$treatment[[i]][1:lb,2] > X$treatment[[i]][1:lb,3]))
                    break
            }

            X$treatment[[i]][(lb + 1):L[[i]],2] <- rbinom(L[[i]] - lb, 20, .9)
            repeat{
                X$treatment[[i]][(lb + 1):L[[i]],3] <- rbinom(L[[i]] - lb, 8, .9)
                if (all(X$treatment[[i]][(lb + 1):L[[i]],2] > X$treatment[[i]][(lb + 1):L[[i]],3]))
                    break
            }
            Y[i] <- rpois(1, lambda=200)
        }else if (i < (N * pi.c[2] + N * pi.c[1])){ # Second methylation profile
            lb <- round(L[[i]] / 1.5)

            X$treatment[[i]][1:lb,2] <- rbinom(lb, 20, .9)
            repeat{
                X$treatment[[i]][1:lb,3] <- rbinom(lb, 2, .8)
                if(all(X$treatment[[i]][1:lb,2] > X$treatment[[i]][1:lb,3]))
                    break
            }

            X$treatment[[i]][(lb + 1):L[[i]],2] <- rbinom(L[[i]] - lb, 20, .9)
            repeat{
                X$treatment[[i]][(lb + 1):L[[i]],3] <- rbinom(L[[i]]-lb, 14, .9)
                if (all(X$treatment[[i]][(lb + 1):L[[i]],2] > X$treatment[[i]][(lb + 1):L[[i]],3]))
                    break
            }
            Y[i] <- rpois(1, lambda=100)
        }else{                  # Third methylation profile
            lb <- round(L[[i]] / 2.5)
            mb <- round(L[[i]] / 3.5)

            X$treatment[[i]][1:lb,2] <- rbinom(lb, 20, .9)
            repeat{
                X$treatment[[i]][1:lb,3] <- rbinom(lb, 2, .9)
                if(all(X$treatment[[i]][1:lb,2] > X$treatment[[i]][1:lb,3]))
                    break
            }

            X$treatment[[i]][(lb + 1):(lb + mb),2] <- rbinom(mb, 20, .9)
            repeat{
                X$treatment[[i]][(lb + 1):(lb + mb),3] <- rbinom(mb, 14, .9)
                if (all(X$treatment[[i]][(lb + 1):(lb + mb),2] > X$treatment[[i]][(lb + 1):(lb + mb),3]))
                    break
            }

            X$treatment[[i]][(lb + 1 + mb):L[[i]],2] <- rbinom(L[[i]] - mb - lb, 20, .9)
            repeat{
                X$treatment[[i]][(lb + 1 + mb):L[[i]],3] <- rbinom(L[[i]] - mb - lb, 2, .9)
                if (all(X$treatment[[i]][(lb + 1 + mb):L[[i]],2] > X$treatment[[i]][(lb + 1 + mb):L[[i]]],3))
                    break
            }
            Y[i] <- rpois(1, lambda=50)
        }
    }

    return(list(X = X, Y = Y))
}
