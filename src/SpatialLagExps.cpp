#include <Rcpp.h>
#include <spatial_lagging.hpp>

/**
 * @title Rcpp wrapper for single-variable lattice lag computation
 *
 * @description
 * This function provides an Rcpp interface to the C++ core function `GenLatticeLagUni()`,
 * which computes lagged mean values for a lattice-based spatial structure defined by
 * explicit neighbor relationships. The computation aggregates values from neighbors
 * up to a specified lag order and calculates the mean while ignoring missing (NaN) values.
 *
 * The function converts an input R numeric vector and neighborhood list into C++ standard
 * containers, applies zero-based indexing adjustment, and invokes the underlying C++ routine.
 *
 * @param vec Numeric vector representing spatial observations for each lattice unit.
 * @param nb List of integer vectors specifying the neighborhood structure.
 *           Each element contains indices of immediate neighbors (1-based in R, converted to 0-based in C++).
 * @param lagNum Integer specifying the lag order (non-negative) to compute.
 *
 * @return Numeric vector where each element represents the mean of lagged values
 *         for the corresponding spatial unit at the specified lag order.
 *         Returns NaN for spatial units with no valid neighbors.
 *
 * @examples
 * \dontrun{
 * # Example usage in R:
 * vec <- c(1.0, 2.0, 3.0, 4.0, 5.0)
 * nb <- list(c(2, 3), c(1, 4), c(1, 5), c(2), c(3))
 * lagNum <- 1
 * result <- RcppGenLatticeLagUni(vec, nb, lagNum)
 * print(result)
 * }
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppGenLatticeLagUni(const Rcpp::NumericVector& vec,
                                         const Rcpp::List& nb,
                                         int lagNum = 1) {
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
  std::vector<double> result = SpatialLagging::GenLatticeLagUni(cpp_vec, cpp_nb, lagNum);

  // Convert result back to R numeric vector
  return Rcpp::wrap(result);
}

/**
 * @title Rcpp wrapper for multi-variable lattice lag computation
 *
 * @description
 * This function provides an Rcpp interface to the C++ function `GenLatticeLagMulti()`,
 * which computes lagged mean values for multiple spatial variables defined on the same lattice structure.
 * Each column in the input matrix is treated as a separate variable.
 *
 * For each variable, the function computes lagged means based on the specified neighborhood structure
 * and lag order (number of steps away in the adjacency graph).
 *
 * @param vecs Numeric matrix where each column represents one spatial variable.
 *             Rows correspond to spatial units (must match the neighborhood size).
 * @param nb List of integer vectors representing the neighborhood structure.
 *           Each element contains indices (1-based in R) of immediate neighbors for each spatial unit.
 * @param lagNum Integer specifying the lag order (non-negative).
 *
 * @return Numeric matrix where each column corresponds to the lagged mean values
 *         for the respective input variable. The dimensions match the input `vecs`.
 *         If all neighbor values are invalid (NaN) for a unit, the result is NaN.
 *
 * @examples
 * \dontrun{
 * # Example usage in R:
 * vecs <- matrix(c(1,2,3,4,5,
 *                  2,3,4,5,6), ncol = 2)
 * nb <- list(c(2,3), c(1,4), c(1,5), c(2), c(3))
 * result <- RcppGenLatticeLagMulti(vecs, nb, 1)
 * print(result)
 * }
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppGenLatticeLagMulti(const Rcpp::NumericMatrix& vecs,
                                           const Rcpp::List& nb,
                                           int lagNum = 1) {
  // --- Validate inputs ---
  if (vecs.nrow() == 0 || vecs.ncol() == 0) {
    Rcpp::stop("Input matrix 'vecs' must not be empty.");
  }
  if (lagNum < 0) {
    Rcpp::stop("Parameter 'lagNum' must be non-negative.");
  }

  // --- Convert nb from R list to std::vector<std::vector<int>> ---
  int n_units = nb.size();
  std::vector<std::vector<int>> cpp_nb(n_units);

  for (int i = 0; i < n_units; ++i) {
    Rcpp::IntegerVector current_nb = nb[i];
    std::vector<int> neighbors;
    neighbors.reserve(current_nb.size());

    for (int j = 0; j < current_nb.size(); ++j) {
      // Convert from 1-based (R) to 0-based (C++)
      neighbors.push_back(current_nb[j] - 1);
    }

    cpp_nb[i] = std::move(neighbors);
  }

  // --- Convert each column of vecs into std::vector<double> ---
  int n_rows = vecs.nrow();
  int n_cols = vecs.ncol();

  std::vector<std::vector<double>> cpp_vecs;
  cpp_vecs.reserve(n_cols);

  for (int j = 0; j < n_cols; ++j) {
    std::vector<double> col(n_rows);
    for (int i = 0; i < n_rows; ++i) {
      col[i] = vecs(i, j);
    }
    cpp_vecs.emplace_back(std::move(col));
  }

  // --- Call the core C++ computation function ---
  std::vector<std::vector<double>> result_cpp =
    SpatialLagging::GenLatticeLagMulti(cpp_vecs, cpp_nb, lagNum);

  // --- Validate output size ---
  if (result_cpp.empty()) {
    Rcpp::stop("Computation returned empty result.");
  }

  // --- Convert back to Rcpp::NumericMatrix ---
  Rcpp::NumericMatrix result(n_rows, n_cols);

  for (int j = 0; j < n_cols; ++j) {
    if (j >= static_cast<int>(result_cpp.size())) {
      Rcpp::stop("Result dimension mismatch between input and output.");
    }
    const auto& col = result_cpp[j];
    if (static_cast<int>(col.size()) != n_rows) {
      Rcpp::stop("Output column length mismatch for variable %d.", j + 1);
    }
    for (int i = 0; i < n_rows; ++i) {
      result(i, j) = col[i];
    }
  }

  return result;
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
 * # result <- RcppGenGridLagUni(mat, 1)
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppGenGridLagUni(const Rcpp::NumericMatrix& mat,
                                      int lagNum = 1) {
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
  std::vector<double> result = SpatialLagging::GenGridLagUni(cpp_mat, lagNum);

  // Convert result back to R numeric vector
  return Rcpp::wrap(result);
}

/**
 * @title Rcpp wrapper for multi-variable grid lag computation (Moore neighborhood)
 *
 * @description
 * This function provides an Rcpp interface to the C++ function `GenGridLagMulti()`,
 * which computes lagged mean values for multiple 2D grid variables using the Moore neighborhood
 * (also known as queen's case). Each column in the input matrix represents one spatial variable,
 * stored in row-major (flattened) order.
 *
 * The function reconstructs each variable’s 2D grid from the flattened column vectors,
 * computes lagged mean values for all variables simultaneously, and returns the results
 * as a numeric matrix with the same structure as the input.
 *
 * @param mat Numeric matrix where each column represents a spatial variable,
 *            stored in row-major flattened order (each row corresponds to one grid cell).
 * @param nrow Integer specifying the number of rows in each 2D grid.
 *             The number of columns in the grid is inferred from `mat.nrow() / nrow`.
 * @param lagNum Integer specifying the lag distance in the Moore neighborhood (non-negative).
 *
 * @return Numeric matrix of the same dimensions as `mat`, where each column represents
 *         the lagged mean values for the corresponding input variable.
 *         NaN values are returned for cells without valid neighbors.
 *
 * @examples
 * \dontrun{
 * # Example usage in R:
 * nrow <- 3
 * mat <- matrix(c(1:9, 9:1), ncol = 2)
 * lagNum <- 1
 * result <- RcppGenGridLagMulti(mat, nrow, lagNum)
 * print(result)
 * }
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppGenGridLagMulti(const Rcpp::NumericMatrix& mat,
                                        int nrow,
                                        int lagNum = 1) {
  // --- Validate inputs ---
  if (mat.nrow() == 0 || mat.ncol() == 0) {
    Rcpp::stop("Input matrix 'mat' must not be empty.");
  }
  if (lagNum < 0) {
    Rcpp::stop("Parameter 'lagNum' must be non-negative.");
  }
  if (nrow <= 0) {
    Rcpp::stop("Parameter 'nrow' must be positive.");
  }

  const int n_cells = mat.nrow();
  const int n_vars = mat.ncol();

  // Derive the number of columns in each grid
  if (n_cells % nrow != 0) {
    Rcpp::stop("nrow does not evenly divide the number of rows in 'mat'.");
  }
  const int ncol = n_cells / nrow;

  // --- Convert input matrix to 3D std::vector: arr[var][row][col] ---
  std::vector<std::vector<std::vector<double>>> cpp_arr;
  cpp_arr.reserve(n_vars);

  for (int v = 0; v < n_vars; ++v) {
    std::vector<std::vector<double>> grid(nrow, std::vector<double>(ncol));
    const Rcpp::NumericVector col_data = mat.column(v);

    for (int i = 0; i < nrow; ++i) {
      for (int j = 0; j < ncol; ++j) {
        // row-major reconstruction
        grid[i][j] = col_data[i * ncol + j];
      }
    }

    cpp_arr.emplace_back(std::move(grid));
  }

  // --- Call core computation function ---
  std::vector<std::vector<double>> result_cpp =
    SpatialLagging::GenGridLagMulti(cpp_arr, lagNum);

  // --- Validate output size ---
  if (result_cpp.size() != static_cast<size_t>(n_vars)) {
    Rcpp::stop("Mismatch between input and output variable counts.");
  }

  // --- Convert back to Rcpp::NumericMatrix ---
  Rcpp::NumericMatrix result(n_cells, n_vars);

  for (int v = 0; v < n_vars; ++v) {
    const auto& flat_vec = result_cpp[v];
    if (static_cast<int>(flat_vec.size()) != n_cells) {
      Rcpp::stop("Output vector size mismatch for variable %d.", v + 1);
    }

    for (int i = 0; i < n_cells; ++i) {
      result(i, v) = flat_vec[i];
    }
  }

  return result;
}
