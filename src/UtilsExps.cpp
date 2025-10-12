#include <Rcpp.h>
#include "pfm_utils.hpp"

/*
 * R wrapper for pfmutils::DiscMat2PFM
 *
 * Input:
 *   x : NumericMatrix (n_samples rows, n_vars columns). NA (NaN) values are allowed;
 *       rows with any NA are ignored internally.
 *
 * Output:
 *   NumericVector with 'dim' attribute set to the per-variable category counts (dims).
 *   This R vector behaves as an n-dimensional array of probabilities (sum ~= 1).
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppDiscMat2PFM(const Rcpp::NumericMatrix& x, double eps = 1e-14) {
  // Convert R matrix (rows=samples, cols=vars) to std::vector<std::vector<double>>
  const int nrow = x.nrow();
  const int ncol = x.ncol();
  if (nrow <= 0 || ncol <= 0) {
    Rcpp::stop("RcppDiscMat2PFM: input matrix must have positive dimensions.");
  }

  std::vector<std::vector<double>> cpp_x;
  cpp_x.resize(static_cast<std::size_t>(ncol));
  for (int j = 0; j < ncol; ++j) {
    cpp_x[j].reserve(static_cast<std::size_t>(nrow));
    for (int i = 0; i < nrow; ++i) {
      cpp_x[j].push_back(static_cast<double>(x(i, j)));
    }
  }

  // Call core C++ function
  std::pair<std::vector<double>, std::vector<int>> res;
  try {
    res = pfmutils::DiscMat2PFM(cpp_x, eps);
  } catch (const std::exception& e) {
    Rcpp::stop(std::string("RcppDiscMat2PFM: Exception: ") + e.what());
  }

  const std::vector<double>& hist = res.first;
  const std::vector<int>& dims = res.second;

  // Convert hist to R NumericVector and set dim attribute
  Rcpp::NumericVector out(hist.begin(), hist.end());
  Rcpp::IntegerVector rdim(dims.begin(), dims.end());
  out.attr("dim") = rdim;

  return out;
}
