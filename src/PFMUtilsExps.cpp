#include <Rcpp.h>
#include "pfm_utils.hpp"

/*
 * R wrapper for pfmutils::DiscMat2PFM
 *
 * Description
 * ------------
 * Converts a discretized numeric matrix into an n-dimensional joint probability
 * frequency matrix (PFM). Each column of the input matrix is treated as a discrete
 * variable whose categories may be non-contiguous numeric labels.
 *
 * The function automatically:
 *   1. Removes rows containing any NA/NaN values.
 *   2. Maps each variable’s unique category values to contiguous integer indices [0..K-1].
 *   3. Accumulates joint frequencies across all valid samples into an n-dimensional
 *      histogram (array), normalized to sum to 1.
 *
 * Input
 * -----
 *   x   : NumericMatrix (n_samples rows, n_vars columns). NA (NaN) values are allowed;
 *         rows with any NA are ignored internally.
 *   eps : (optional) double, a small constant added to empty cells before normalization
 *         to avoid zero probabilities. Default = 1e-14.
 *
 * Output
 * ------
 *   NumericVector with a 'dim' attribute corresponding to the number of unique
 *   categories per variable. The resulting object behaves as an n-dimensional
 *   probability array in R, summing to approximately 1.
 *
 * Details
 * --------
 * - Category mappings are determined independently for each variable based on
 *   the sorted order of unique observed values.
 * - This wrapper calls the C++ function pfmutils::DiscMat2PFM implemented in
 *   `pfmutils_discmat2pfm.hpp`.
 * - The returned array is consistent with R’s column-major ordering; i.e.,
 *   indexing via `pfm[i1, i2, ...]` in R corresponds to the correct cell.
 *
 * Example
 * -------
 * @examples
 * # Example usage in R:
 * mat <- matrix(c(
 *   1, 2,
 *   1, 3,
 *   2, 2,
 *   3, 1
 * ), ncol = 2, byrow = TRUE)
 *
 * # Compute the probability frequency matrix
 * pfm <- RcppDiscMat2PFM(mat)
 *
 * # Check the dimension (number of categories per variable)
 * dim(pfm)
 *
 * # Verify it sums to 1
 * sum(pfm)
 *
 * # Inspect the probability array
 * print(pfm)
 *
 * # Expected output (approximate, depending on unique values):
 * # , , Var2 = 1
 * #     [,1] [,2] [,3]
 * # [1,] 0.25 0.00 0.00
 * # [2,] 0.00 0.00 0.25
 * # [3,] 0.00 0.00 0.00
 * #
 * # , , Var2 = 2
 * # [1,] 0.25 0.00 0.00
 * # [2,] 0.00 0.25 0.00
 * # [3,] 0.00 0.00 0.00
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
