#include <Rcpp.h>
#include <limits>

/**
 * Rcpp wrapper for generating lagged versions of multiple time series variables.
 *
 * This function applies a time-lag transformation to each column (variable)
 * of the input matrix. For each variable, it shifts the values downward by
 * `lagNum` steps, filling the initial entries with NaN. The output matrix
 * has the same dimensions as the input.
 *
 * Parameters
 * ----------
 * @param mat NumericMatrix representing multiple time series variables.
 *            Each column corresponds to one variable, and each row is a time step.
 * @param lagNum Integer specifying the lag number (default = 1).
 *
 * Returns
 * -------
 * NumericMatrix of the same size as the input, where each column contains the
 * lagged values of the corresponding variable. The first `lagNum` rows of each
 * column are filled with NaN due to lack of previous values.
 *
 * Example
 * -------
 * @examples
 * # Example usage in R:
 * # mat <- matrix(1:9, ncol = 3)
 * # result <- RcppGenTSLagMulti(mat, 1)
 * # print(result)
 * #
 * # Expected output:
 * #      [,1] [,2] [,3]
 * # [1,]  NaN  NaN  NaN
 * # [2,]    1    4    7
 * # [3,]    2    5    8
 *
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppGenTSLagMulti(const Rcpp::NumericMatrix& mat,
                                      int lagNum = 1) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();

  // Initialize result matrix with NaN
  Rcpp::NumericMatrix lagged(nrow, ncol);
  double nan = std::numeric_limits<double>::quiet_NaN();

  // Fill result with NaN
  std::fill(lagged.begin(), lagged.end(), nan);

  // Apply lag to each column independently
  for (int j = 0; j < ncol; ++j) {
    for (int i = lagNum; i < nrow; ++i) {
      lagged(i, j) = mat(i - lagNum, j);
    }
  }

  return lagged;
}
