#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @rdname bpr_model
//'
//' @export
// [[Rcpp::export]]
double bpr_likelihood(const arma::vec& w,
                      const arma::mat& H,
                      const arma::mat& data,
                      const double lambda,
                      const bool is_NLL){

    int nrows = data.n_rows;
    int ncols = data.n_cols;
    double lik = 0.0;

    // Predictions of the target variables
    Rcpp::NumericVector g = Rcpp::wrap(H * w);
    // Compute the cdf of N(0,1) distribution (i.e. probit function)
    Rcpp::NumericVector Phi = Rcpp::pnorm(g);

    for (int i = 0; i < nrows; i++){
        // In extreme cases where probit is 0 or 1, subtract a tiny number
        // so we can evaluate the log(0) when computing the likelihood
        if (Phi[i] > (1 - 1e-15)){
            Phi[i] = 1 - 1e-15;
        }else if (Phi[i] < 1e-15){
            Phi[i] = 1e-15;
        }
        // Compute the log likelihood
        if (ncols == 3){ // If Binomial distributed data
            lik += R::dbinom(data(i, 2), data(i, 1), Phi[i], true);
        }else{           // If Bernoulli distributed data
            lik += R::dbinom(data(i, 1), 1, Phi[i], true);
        }
    }

    // Compute ridge regression likelihood
    lik = lik - lambda * arma::as_scalar(w.t() * w);
    // If we require the Negative Log Likelihood
    if (is_NLL == true){
        lik = (-1) * lik;
    }
    return lik;
}


//' @rdname bpr_model
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector bpr_gradient(const arma::vec& w,
                                 const arma::mat& H,
                                 const arma::mat& data,
                                 const double lambda,
                                 const bool is_NLL){

    int nrows = data.n_rows;
    int ncols = data.n_cols;
    int M = w.size();

    // Predictions of the target variables
    Rcpp::NumericVector g = Rcpp::wrap(H * w);
    // Compute the cdf of N(0,1) distribution (i.e. probit function)
    Rcpp::NumericVector Phi = Rcpp::pnorm(g);
    // Compute the density of a N(0,1) distribution
    Rcpp::NumericVector N = Rcpp::dnorm(g);
    Rcpp::NumericVector gr(M);

    for (int i = 0; i < nrows; i++){
        // In extreme cases where probit is 0 or 1, subtract a tiny number
        // so we can evaluate the log(0) when computing the likelihood
        if (Phi[i] > (1 - 1e-15)){
            Phi[i] = 1 - 1e-15;
        }else if (Phi[i] < 1e-15){
            Phi[i] = 1e-15;
        }
        if (N[i] < 1e-15){
            N[i] = 1e-15;
        }
        // Compute the gradient vector w.r.t the coefficients w
        if (ncols == 3){  // If Binomial distributed data
            for (int m = 0; m < M; m++){
                gr[m] += N[i] * H(i, m) * (data(i, 2) / Phi[i] - (data(i, 1) -
                    data(i, 2)) / (1 - Phi[i]) );
            }
        }else{            // If Bernoulli distributed data
            for (int m = 0; m < M; m++){
                gr[m] += N[i] * H(i, m) * (data(i, 1) / Phi[i] -
                    (1 - data(i, 1)) / (1 - Phi[i]) );
            }
        }
    }
    for (int m = 0; m < M; m++){
        // Compute ridge regression likelihood
        gr[m] -= 2 * lambda * w[m];
        // If we require the Negative Log Likelihood
        if (is_NLL == true){
            gr[m] *= -1;
        }
    }
    return gr;
}


//' @rdname bpr_model
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector bpr_lik_region(const arma::vec& w,
                                   const Rcpp::List& x,
                                   const Rcpp::List& des_mat,
                                   const double lambda,
                                   const bool is_NLL){
    // Number of regions
    int N = x.size();
    Rcpp::NumericVector res(N);
    for (int i = 0; i < N; i++){
        // Extract observations in each region
        arma::mat data = x[i];
        // Extract deisgn matrix of each region
        arma::mat H = des_mat[i];
//         Rcpp::List des_obj = des_mat[i];
//         arma::mat H = as<arma::mat>(des_obj["H"]);
        // Compute BPR likelihood
        res[i] = bpr_likelihood(w, H, data, lambda, is_NLL);
    }
    return res;
}


//' @rdname bpr_model
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix bpr_lik_resp(const arma::mat& w,
                                 const Rcpp::List& x,
                                 const Rcpp::List& des_mat,
                                 const arma::vec pi_k,
                                 const double lambda,
                                 const bool is_NLL){
    // Number of regions
    int N = x.size();
    int K = w.n_cols;
    int k;
    Rcpp::NumericMatrix res(N, K);
    for (int i = 0; i < N; i++){
        // Extract observations in each region
        arma::mat data = x[i];
        // Extract deisgn matrix of each region
        arma::mat H = des_mat[i];
        // Compute BPR likelihood
        for (k = 0; k < K; k++){
            arma::vec ww = w.col(k);
            res(i, k) = pi_k[k] + bpr_likelihood(ww, H, data, lambda, is_NLL);
        }
    }
    return res;
}


//' @rdname bpr_model
//'
//' @export
// [[Rcpp::export]]
double sum_weighted_bpr_lik(const arma::vec& w,
                            const Rcpp::List& x,
                            const Rcpp::List& des_mat,
                            const arma::vec& post_prob,
                            const double lambda,
                            const bool is_NLL){
    // Number of regions
    int N = x.size();
    Rcpp::NumericVector res(N);
    for (int i = 0; i < N; i++){
        // Extract observations in each region
        arma::mat data = x[i];
        // Extract deisgn matrix of each region
        arma::mat H = des_mat[i];
        // Compute BPR likelihood
        res[i] = bpr_likelihood(w, H, data, lambda, is_NLL);
    }
    // Inner product with the weight vector of posterior probabilities
    arma::vec lik = as<arma::vec>(res);
    return arma::as_scalar(post_prob.t() * lik);
}


//' @rdname bpr_model
//'
//' @export
// [[Rcpp::export]]
arma::rowvec sum_weighted_bpr_grad(const arma::vec& w,
                                   const Rcpp::List& x,
                                   const Rcpp::List& des_mat,
                                   const arma::vec& post_prob,
                                   const double lambda,
                                   const bool is_NLL){
    // Number of regions
    int N = x.size();
    // Number of basis functions
    int M = w.size();
    Rcpp::NumericMatrix res(N, M);
    for (int i = 0; i < N; i++){
        // Extract observations in each region
        arma::mat data = x[i];
        // Extract deisgn matrix of each region
        arma::mat H = des_mat[i];
        // Compute the gradient of BPR model
        res(i, _) = bpr_gradient(w, H, data, lambda, is_NLL);
    }
    // Inner product with the weight vector of posterior probabilities
    arma::mat lik = as<arma::mat>(res);
    arma::rowvec w_lik = post_prob.t() * lik;
    return w_lik;
}
