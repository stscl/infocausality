#include <Rcpp.h>
#include <spatial_lagging.hpp>


/**
 * Rcpp wrapper for computing the mean of lagged values for a lattice structure at a specified lag number.
 *
 * @param vec Numeric vector representing the spatial data for each unit.
 * @param nb List of integer vectors representing the neighborhood structure.
 *           Each element contains indices of immediate neighbors for a spatial unit.
 * @param lagNum Integer specifying the number of lag steps to compute.
 *
 * @return Numeric vector where each element represents the mean of lagged values
 *         for the corresponding spatial unit at the specified lag number.
 *         Returns NaN for spatial units with no valid neighbors.
 *
 * @examples
 * # Example usage in R:
 * # vec <- c(1.0, 2.0, 3.0, 4.0, 5.0)
 * # nb <- list(c(2, 3), c(1, 4), c(1, 5), c(2), c(3))
 * # result <- RcppGenLatticeLag(vec, nb, 1)
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppGenLatticeLag(Rcpp::NumericVector vec,
                                      Rcpp::List nb,
                                      int lagNum) {
  // Convert R numeric vector to std::vector<double>
  std::vector<double> cpp_vec = Rcpp::as<std::vector<double>>(vec);

  // Get the number of elements in the nb object
  int n = nb.size();

  // Create a std::vector<std::vector<int>> to store the result
  std::vector<std::vector<int>> cpp_nb(n);

  // Iterate over each element in the nb object
  for (int i = 0; i < n; ++i) {
    // Get the current element (should be an integer vector)
    Rcpp::IntegerVector current_nb = nb[i];

    // Create a vector<int> to store the current subset of elements
    std::vector<int> current_subset;

    // Iterate over each element in the current subset
    for (int j = 0; j < current_nb.size(); ++j) {
      // Subtract one from each element to convert from R's 1-based indexing to C++'s 0-based indexing
      current_subset.push_back(current_nb[j] - 1);
    }

    // Add the current subset to the result
    cpp_nb[i] = current_subset;
  }

  // Call the C++ function
  std::vector<double> result = SpatialLagging::GenLatticeLag(cpp_vec, cpp_nb, lagNum);

  // Convert result back to R numeric vector
  return Rcpp::wrap(result);
}

/**
 * Rcpp wrapper for computing the mean of lagged values for a grid structure at a specified lag number.
 *
 * @param mat Numeric matrix representing the grid data.
 * @param lagNum Integer specifying the number of steps to lag when considering
 *               neighbors in the Moore neighborhood.
 *
 * @return Numeric vector where each element represents the mean of lagged values
 *         for the corresponding grid cell at the specified lag number.
 *         The results are arranged in row-major order (same as the input grid).
 *         Returns NaN for grid cells with no valid neighbors.
 *
 * @examples
 * # Example usage in R:
 * # mat <- matrix(1:9, nrow = 3, ncol = 3)
 * # result <- RcppGenGridLag(mat, 1)
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppGenGridLag(Rcpp::NumericMatrix mat,
                                   int lagNum) {
  // Convert R matrix to std::vector<std::vector<double>>
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  std::vector<std::vector<double>> cpp_mat(nrow, std::vector<double>(ncol));

  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      cpp_mat[i][j] = mat(i, j);
    }
  }

  // Call the C++ function
  std::vector<double> result = SpatialLagging::GenGridLag(cpp_mat, lagNum);

  // Convert result back to R numeric vector
  return Rcpp::wrap(result);
}
