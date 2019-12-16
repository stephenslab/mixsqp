// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>

// This depends statement is needed to tell R where to find the
// additional header files.
//
// [[Rcpp::depends(RcppArmadillo)]]
//

// FUNCTION DECLARATIONS
// ---------------------
void get_q(const arma::mat& L, arma::mat& Q, int m, int l,
           int power_iteration_parameter);

// FUNCTION DEFINITIONS
// --------------------
// [[Rcpp::export]]
Rcpp::List randomized_svd(const arma::mat& L, int k,
                          int oversampling_parameter,
                          int power_iteration_parameter,
                          double relative_tolerance,
                          std::string sampling_distribution,
                          std::string return_type,
                          bool rsvd_verbose) {
  
  // -------------------------------------------------------
  // STEP 1: miscellaneous things
  // -------------------------------------------------------
  // set seed
  // please set seed at the R level via set.seed()
  
  // set size
  int n         = L.n_rows;
  int m         = L.n_cols;
  int l         = k + oversampling_parameter;
  
  // check if rsvd is necessary
  if (n < m) {
    // if n < m, rsvd is not necessary for our purpose
    return Rcpp::List::create(Rcpp::Named("status") = "fail",
                              Rcpp::Named("reason") = "n < m");
  }
  if (m < k) {
    // decrease k
    return Rcpp::List::create(Rcpp::Named("status") = "fail",
                              Rcpp::Named("reason") = "m < k");
  }
  if (m < l) {
    // decrease oversampling_paramaeter
    return Rcpp::List::create(Rcpp::Named("status") = "fail",
                              Rcpp::Named("reason") = "m < k + p");
  }
  // -------------------------------------------------------
  
  
  // -------------------------------------------------------
  // STEP 2: matrix factorization routines and report first results
  // -------------------------------------------------------
  // run randomized matrix factorization routine
  arma::mat Q, U, V;
  get_q(L, Q, m, l, power_iteration_parameter);
  arma::vec s;
  arma::svd_econ(U,s,V,Q.t() * L, "both", "std");
  if (k > s.n_elem) {
    return Rcpp::List::create(Rcpp::Named("status") = "fail",
                              Rcpp::Named("reason") = "L has rank < k");
  }
  // -------------------------------------------------------
  // report k and rtol
  double smax   = s[0];
  double smin   = s[k-1];
  if (rsvd_verbose) {
    Rprintf("with oversampling  ---------   k: %3d,    rtol: %0.3e\n", l, s[l-1] / smax);
    Rprintf("with requested k   ---------   k: %3d,    rtol: %0.3e\n", k, smin / smax);
  }
  // -------------------------------------------------------
  
  // -------------------------------------------------------
  // STEP 3: truncate if necessary and return results
  // -------------------------------------------------------
  // First, when the obtained relative tolerance is smaller than the requested relative tolerance,
  // we will do truncation and return results
  // -------------------------------------------------------
  int k_new     = 0;
  for (; k_new < k; ++k_new) {
    if (s[k_new] / smax < relative_tolerance) break;
  }
  smin          = s[k_new - 1];
  
  // report k and rtol
  if (rsvd_verbose) {
    Rprintf("with truncated k   ---------   k: %3d,    rtol: %0.3e\n", k_new, smin / smax);
    if (k_new == k) {
      Rprintf("please increase k if relative tolerance obtained is not satisfactory\n");
    }
  }
  
  // return results
  if (return_type == std::string("svd")) {
    return Rcpp::List::create(Rcpp::Named("s") = s.subvec(0,k_new - 1),
                              Rcpp::Named("U") = Q * U.cols(0,k_new - 1),
                              Rcpp::Named("V") = V.cols(0,k_new - 1),
                              Rcpp::Named("status") = "truncated");
  } else if (return_type == std::string("qb")) {
    U.each_row()       %= s.t(); 
    return Rcpp::List::create(Rcpp::Named("Q") = Q,
                              Rcpp::Named("Bt") = V.cols(0,k_new - 1) * (U.cols(0,k_new - 1)).t(),
                              Rcpp::Named("status") = "truncated");
  } else {
    return Rcpp::List::create(Rcpp::Named("status") = "fail",
                              Rcpp::Named("reason") = "invalid return type");
  }
  
  
  /*if (smin / smax < relative_tolerance) {
    int k_new     = 0;
    for (; k_new < k; k_new++) {
      if (s[k_new] / smax < relative_tolerance) break;
    }
    smin          = s[k_new];
    
    // report k and rtol
    if (rsvd_verbose) {
      Rprintf("k: %3d,    rtol: %0.3e\n", k_new + 1, smin / smax);
    }
    
    // return results
    if (return_type == std::string("svd")) {
      return Rcpp::List::create(Rcpp::Named("s") = s.subvec(0,k_new),
                                Rcpp::Named("U") = Q * U.cols(0,k_new),
                                Rcpp::Named("V") = V.cols(0,k_new),
                                Rcpp::Named("status") = "truncated");
    } else if (return_type == std::string("qb")) {
      U.each_row()       %= s.t(); 
      return Rcpp::List::create(Rcpp::Named("Q") = Q,
                                Rcpp::Named("B") = U.cols(0,k_new) * (V.cols(0,k_new)).t(),
                                Rcpp::Named("status") = "truncated");
    } else {
      return Rcpp::List::create(Rcpp::Named("status") = "fail",
                                Rcpp::Named("reason") = "invalid return type");
    }
  // -------------------------------------------------------
  // Second, when the obtained relative tolerance is not smaller than the requested relative tolerance,
  // we will do truncation and return results.
  // Also, we will warn that the result is not satisfactory in terms of relative tolerance
  // -------------------------------------------------------
  } else {
    // report k and rtol
    if (rsvd_verbose) {
      Rprintf("not truncated due to small relative tolerance requested\n");
      Rprintf("please increase k if relative tolerance obtained is not satisfactory\n");
    }
    // return results
    if (return_type == std::string("svd")) {
      return Rcpp::List::create(Rcpp::Named("s") = s,
                                Rcpp::Named("U") = Q * U,
                                Rcpp::Named("V") = V,
                                Rcpp::Named("status") = "not truncated");
    } else if (return_type == std::string("qb")) {
      U.each_row()       %= s.t(); 
      return Rcpp::List::create(Rcpp::Named("Q") = Q,
                                Rcpp::Named("B") = U * V.t(),
                                Rcpp::Named("status") = "not truncated");
    } else {
      return Rcpp::List::create(Rcpp::Named("status") = "fail",
                                Rcpp::Named("reason") = "invalid return type");
    }
  } */
}

void get_q(const arma::mat& L, arma::mat& Q, int m, int l,
           int power_iteration_parameter) {
  
  // construct random projection matrix
  arma::mat M   = arma::randn<arma::mat>(m,l);
  
  // construct projected matrix
  arma::mat Y   = L * M;
  
  // predefine
  arma::mat U(l,l);
  arma::mat Z(m,l);
  
  // subspace iteration to find better projection
  if (power_iteration_parameter > 0) {
    arma::qr_econ(Q,U,Y);
    Z = L.t() * Q;
    arma::qr_econ(M,U,Z);
    Y = L * M;
  }
  arma::qr_econ(Q,U,Y);
}