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
 * # result <- RcppGenTSLagMultiSingle(mat, 1)
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
Rcpp::NumericMatrix RcppGenTSLagMultiSingle(const Rcpp::NumericMatrix& mat,
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

/**
 * @title Rcpp wrapper for generating lagged versions of multiple time series variables
 *
 * @description
 * This function applies a time-lag transformation to each column (variable)
 * of the input matrix. Unlike the single-lag version, it allows specifying
 * a different lag order for each variable using the `lagNums` vector.
 *
 * For each variable, the function shifts its values downward by the corresponding
 * lag amount, filling the initial entries with NaN. The output matrix has the same
 * dimensions as the input.
 *
 * @param mat NumericMatrix representing multiple time series variables.
 *            Each column corresponds to one variable, and each row is a time step.
 * @param lagNums Integer vector specifying the lag number for each variable (column of `mat`).
 *                Must be the same length as the number of columns in `mat`.
 *
 * @return NumericMatrix of the same size as the input, where each column contains
 *         the lagged values of the corresponding variable. The first `lagNums[j]`
 *         rows of each column are filled with NaN due to lack of previous values.
 *
 * @examples
 * \dontrun{
 * # Example usage in R:
 * mat <- matrix(1:9, ncol = 3)
 * lagNums <- c(1, 2, 3)
 * result <- RcppGenTSLagMulti(mat, lagNums)
 * print(result)
 *
 * # Expected output (approximately):
 * #      [,1] [,2] [,3]
 * # [1,]  NaN  NaN  NaN
 * # [2,]    1  NaN  NaN
 * # [3,]    2    1  NaN
 * # [4,]    3    2    1
 * }
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppGenTSLagMulti(const Rcpp::NumericMatrix& mat,
                                      const Rcpp::IntegerVector& lagNums) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();

  // --- Validate inputs ---
  if (nrow == 0 || ncol == 0) {
    Rcpp::stop("Input matrix 'mat' must not be empty.");
  }
  if (lagNums.size() != ncol) {
    Rcpp::stop("Length of 'lagNums' must match the number of columns in 'mat'.");
  }
  for (int j = 0; j < ncol; ++j) {
    if (lagNums[j] < 0) {
      Rcpp::stop("All elements of 'lagNums' must be non-negative.");
    }
  }

  // --- Initialize result matrix with NaN ---
  Rcpp::NumericMatrix lagged(nrow, ncol);
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::fill(lagged.begin(), lagged.end(), nan);

  // --- Apply variable-specific lag ---
  for (int j = 0; j < ncol; ++j) {
    int lag = lagNums[j];
    if (lag == 0) {
      // No lag: copy original column
      for (int i = 0; i < nrow; ++i) {
        lagged(i, j) = mat(i, j);
      }
    } else if (lag < nrow) {
      // Shift values downward by lag positions
      for (int i = lag; i < nrow; ++i) {
        lagged(i, j) = mat(i - lag, j);
      }
    }
    // If lag >= nrow → entire column stays NaN
  }

  return lagged;
}
