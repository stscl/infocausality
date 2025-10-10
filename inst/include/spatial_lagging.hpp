#include <vector>
#include <numeric>
#include <cmath>
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

/**
 * @namespace SpatialLagging
 * @brief This namespace contains functions for computing lagged values and neighbors
 *        for lattice and grid structures. It encapsulates related functionalities
 *        to avoid naming conflicts.
 */
namespace SpatialLagging {

/**
 * Computes the lagged neighbors for a lattice structure up to a specified lag number.
 * This function recursively expands the neighbors at each lag step, starting with direct neighbors
 * (lag 0), and including neighbors from previous lags, until reaching the specified lag number.
 *
 * For lagNum = 0, each spatial unit's neighbor is itself.
 * For lagNum >= 1, the function accumulates neighbors from all previous lags and deduplicates the results.
 * Empty results are filled with `std::numeric_limits<int>::min()` to indicate no neighbors.
 *
 * Parameters:
 *   spNeighbor - A 2D vector where each element contains indices of immediate neighbors for each spatial unit.
 *   lagNum     - The number of lag steps to compute (must be non-negative).
 *
 * Returns:
 *   A 2D vector where each element represents the list of lagged neighbors for a spatial unit.
 */
std::vector<std::vector<int>> CppLaggedNeighbor4Lattice(const std::vector<std::vector<int>>& spNeighbor,
                                                        int lagNum) {
  // Handle negative lagNum: return empty vector
  if (lagNum < 0) {
    return {};
  }

  // If lagNum is 0, return a vector of indices
  if (lagNum == 0) {
    std::vector<std::vector<int>> result;
    for (size_t i = 0; i < spNeighbor.size(); ++i) {
      result.push_back({static_cast<int>(i)});
    }
    return result;
  }

  // // Handle lagNum=1: return the immediate neighbors directly
  // if (lagNum == 1) {
  //   return spNeighbor;
  // }

  // Recursively compute results for lagNum-1
  std::vector<std::vector<int>> prevResult = CppLaggedNeighbor4Lattice(spNeighbor, lagNum - 1);
  std::vector<std::vector<int>> currentResult;

  int n = spNeighbor.size();
  // Process each spatial unit to compute current lagNum's neighbors
  for (int i = 0; i < n; ++i) {
    // Check if prevResult[i] size is equal to n
    if (prevResult[i].size() == spNeighbor.size()) {
      currentResult.push_back(prevResult[i]);
      continue; // Skip further processing for this index
    }

    std::unordered_set<int> mergedSet;

    // Add previous lag results (lag from 0 to lagNum-1)
    for (int elem : prevResult[i]) {
      if (elem != std::numeric_limits<int>::min()) {
        mergedSet.insert(elem);
      }
    }

    // Collect new elements from neighbors of previous lag's results
    std::unordered_set<int> newElements;
    for (int j : prevResult[i]) {
      // Skip invalid indices and placeholder min value
      if (j == std::numeric_limits<int>::min() || j < 0 || j >= n) {
        continue;
      }
      // Aggregate neighbors of j from spNeighbor
      for (int k : spNeighbor[j]) {
        newElements.insert(k);
      }
    }

    // Merge new elements into the set
    for (int elem : newElements) {
      mergedSet.insert(elem);
    }

    // Convert set to sorted vector and deduplicate
    std::vector<int> vec(mergedSet.begin(), mergedSet.end());
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());

    // Handle empty result by filling with min value
    if (vec.empty()) {
      vec.push_back(std::numeric_limits<int>::min());
    }

    currentResult.push_back(vec);
  }

  return currentResult;
}

/**
 * Computes the lagged values for a given vector based on the neighborhood structure and lag number.
 * This function first determines the lagged neighbors for each spatial unit using
 * the `CppLaggedNeighbor4Lattice` function. If `lagNum > 0`, it removes duplicate indices that
 * appeared in previous lag levels to ensure each lag level captures only new neighbors.
 *
 * For each spatial unit, the function extracts values from `vec` corresponding to the computed
 * lagged neighbors. If no valid neighbors exist, the function fills the result with `NaN`.
 *
 * Parameters:
 *   vec    - A vector of double values representing the spatial data for each unit.
 *   nb     - A 2D vector where each row contains indices of immediate neighbors in the lattice.
 *   lagNum - The number of lag steps to compute (must be non-negative).
 *
 * Returns:
 *   A 2D vector where each element contains the lagged values corresponding to the computed
 *   lagged neighbors for each spatial unit.
 */
std::vector<std::vector<double>> CppLaggedVal4Lattice(const std::vector<double>& vec,
                                                      const std::vector<std::vector<int>>& nb,
                                                      int lagNum) {
  int n = vec.size();

  // Compute the lagged neighbors using the provided function
  std::vector<std::vector<int>> laggedNeighbors = CppLaggedNeighbor4Lattice(nb, lagNum);
  // Remove duplicates with previous lagNum (if lagNum > 0)
  if (lagNum > 0) {
    std::vector<std::vector<int>> prevLaggedResults = CppLaggedNeighbor4Lattice(nb, lagNum - 1);
    for (int i = 0; i < n; ++i) {
      // Convert previous lagged results to a set for fast lookup
      std::unordered_set<int> prevSet(prevLaggedResults[i].begin(), prevLaggedResults[i].end());

      // Remove duplicates from current lagged results
      std::vector<int> newIndices;
      for (int index : laggedNeighbors[i]) {
        if (prevSet.find(index) == prevSet.end()) {
          newIndices.push_back(index);
        }
      }

      // If the new indices are empty, set it to a special value (e.g., std::numeric_limits<int>::min())
      if (newIndices.empty()) {
        newIndices.push_back(std::numeric_limits<int>::min());
      }

      // Update the lagged results
      laggedNeighbors[i] = newIndices;
    }
  }

  // Initialize the result vector with the same number of rows as the lagged neighbors
  std::vector<std::vector<double>> result(laggedNeighbors.size());

  // Iterate over each point in the lattice
  for (size_t i = 0; i < laggedNeighbors.size(); ++i) {
    // Initialize the lagged values for the current point
    std::vector<double> laggedValues;

    if (laggedNeighbors[i].size() == 1 && laggedNeighbors[i][0] == std::numeric_limits<int>::min()){
      // If the index is out of bounds, push a default value (e.g., nan)
      laggedValues.push_back(std::numeric_limits<double>::quiet_NaN());
    } else {
      // Iterate over each neighbor index and extract the corresponding value from `vec`
      for (int neighborIndex : laggedNeighbors[i]) {
        // Check if the neighbor index is valid
        if (neighborIndex >= 0 && neighborIndex < n) {
          laggedValues.push_back(vec[neighborIndex]);
        } else {
          // If the index is out of bounds, push a default value (e.g., nan)
          laggedValues.push_back(std::numeric_limits<double>::quiet_NaN());
        }
      }
    }

    // Add the lagged values to the result
    result[i] = laggedValues;
  }

  return result;
}

/**
 * Computes lagged mean values for lattice structures using neighbor relationships up to a specified lag number.
 * For each spatial unit, this function aggregates values from its lagged neighbors (up to lagNum steps away),
 * calculates the mean, and ignores NaN values. If all neighbor values are invalid (NaN), the result is NaN.
 *
 * Parameters:
 *   vec     - Input vector of spatial observations.
 *   nb      - Immediate neighbor indices for lattice elements (2D vector).
 *   lagNum  - Non-negative integer specifying the lag order.
 *
 * Returns:
 *   A vector where each element represents the lagged mean value for the corresponding spatial unit.
 */
std::vector<double> GenLatticeLagUni(
    const std::vector<double>& vec,
    const std::vector<std::vector<int>>& nb,
    int lagNum
) {
  if (vec.empty() || lagNum < 0) {
    return {};
  }

  // Get lagged values via existing utility
  std::vector<std::vector<double>> laggedValues = CppLaggedVal4Lattice(vec, nb, lagNum);
  std::vector<double> laggedMeans;

  for (const auto& neighborVals : laggedValues) {
    double sum = 0.0;
    int count = 0;

    for (double val : neighborVals) {
      if (!std::isnan(val)) {
        sum += val;
        ++count;
      }
    }

    if (count > 0) {
      laggedMeans.push_back(sum / count);
    } else {
      laggedMeans.push_back(std::numeric_limits<double>::quiet_NaN());
    }
  }

  return laggedMeans;
}

/**
 * Computes lagged mean values for multiple variables in lattice structures.
 * This function processes a matrix of spatial variables and computes lagged means
 * for each variable using neighbor relationships up to a specified lag number.
 *
 * The function leverages existing infrastructure for single-variable lag computation
 * and efficiently applies it to multiple variables with minimal overhead.
 *
 * Parameters:
 *   mat     - Input matrix where each row represents a spatial variable (vector of doubles).
 *   nb      - Immediate neighbor indices for lattice elements (2D vector).
 *   lagNum  - Non-negative integer specifying the lag order.
 *
 * Returns:
 *   A matrix where each row contains lagged mean values for the corresponding input variable.
 */
std::vector<std::vector<double>> GenLatticeLagMultiSingle(
    const std::vector<std::vector<double>>& mat,
    const std::vector<std::vector<int>>& nb,
    int lagNum) {

  // Early return for invalid inputs
  if (mat.empty() || lagNum < 0) {
    return {};
  }

  const size_t numVars = mat.size();
  const size_t numUnits = mat[0].size();

  // // Validate input matrix consistency
  // for (size_t i = 1; i < numVars; ++i) {
  //   if (mat[i].size() != numUnits) {
  //     return {}; // Inconsistent spatial unit counts
  //   }
  // }

  // Precompute lagged neighbors once for efficiency
  std::vector<std::vector<int>> laggedNeighbors;
  if (lagNum > 0) {
    laggedNeighbors = CppLaggedNeighbor4Lattice(nb, lagNum);

    // Remove duplicates from previous lag levels if needed
    if (lagNum > 0) {
      std::vector<std::vector<int>> prevLaggedResults = CppLaggedNeighbor4Lattice(nb, lagNum - 1);
      for (size_t i = 0; i < laggedNeighbors.size(); ++i) {
        std::unordered_set<int> prevSet(prevLaggedResults[i].begin(), prevLaggedResults[i].end());
        std::vector<int> newIndices;

        for (int index : laggedNeighbors[i]) {
          if (prevSet.find(index) == prevSet.end()) {
            newIndices.push_back(index);
          }
        }

        if (newIndices.empty()) {
          newIndices.push_back(std::numeric_limits<int>::min());
        }

        laggedNeighbors[i] = newIndices;
      }
    }
  }

  // Initialize result matrix
  std::vector<std::vector<double>> result(numVars, std::vector<double>(numUnits));

  // Process each variable in parallel (if needed) or sequentially
  for (size_t varIdx = 0; varIdx < numVars; ++varIdx) {
    const std::vector<double>& currentVar = mat[varIdx];

    // Handle lagNum = 0 case efficiently
    if (lagNum == 0) {
      for (size_t unitIdx = 0; unitIdx < numUnits; ++unitIdx) {
        result[varIdx][unitIdx] = currentVar[unitIdx];
      }
      continue;
    }

    // Compute lagged means for current variable
    for (size_t unitIdx = 0; unitIdx < numUnits; ++unitIdx) {
      const std::vector<int>& neighbors = laggedNeighbors[unitIdx];
      double sum = 0.0;
      int validCount = 0;

      // Special case: no valid neighbors
      if (neighbors.size() == 1 && neighbors[0] == std::numeric_limits<int>::min()) {
        result[varIdx][unitIdx] = std::numeric_limits<double>::quiet_NaN();
        continue;
      }

      // Aggregate values from valid neighbors
      for (int neighborIdx : neighbors) {
        if (neighborIdx >= 0 && neighborIdx < static_cast<int>(numUnits)) {
          double value = currentVar[neighborIdx];
          if (!std::isnan(value)) {
            sum += value;
            ++validCount;
          }
        }
      }

      // Compute mean or set to NaN if no valid neighbors
      if (validCount > 0) {
        result[varIdx][unitIdx] = sum / validCount;
      } else {
        result[varIdx][unitIdx] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  return result;
}

/**
 * Computes lagged mean values for multiple spatial variables on a lattice structure,
 * where each variable can have its own distinct lag order.
 *
 * This function generalizes the GenLatticeLagMultiSingle by accepting a vector
 * of lag orders (lagNums), enabling per-variable lag computation. It ensures
 * computational efficiency by only computing neighbor sets for the unique lag values
 * present in lagNums, and their immediate lower levels (lagNum - 1), to derive
 * "pure" (non-overlapping) neighbor sets for each lag order.
 *
 * Core idea:
 *   1. Identify all unique lag levels from lagNums.
 *   2. For each unique lag, compute its neighbors using CppLaggedNeighbor4Lattice,
 *      as well as those for (lag - 1) if needed.
 *   3. Subtract (lag - 1) neighbors from lag neighbors to form pure lag neighbors.
 *   4. For each variable, retrieve the corresponding pure neighbor set for its lagNum
 *      and compute the mean of its neighbors' values.
 *
 * Parameters:
 *   mat      - Input matrix where each row represents a spatial variable
 *              and each column represents a spatial unit.
 *   nb       - A 2D vector of immediate neighbors for each spatial unit.
 *   lagNums  - A vector of lag orders corresponding to each variable in mat.
 *
 * Returns:
 *   A matrix (same shape as mat) where each row contains the lagged mean values
 *   for the corresponding variable. Missing or invalid results are filled with NaN.
 */
std::vector<std::vector<double>> GenLatticeLagMulti(
    const std::vector<std::vector<double>>& mat,
    const std::vector<std::vector<int>>& nb,
    const std::vector<int>& lagNums) {

  // ==== Input validation ====
  if (mat.empty() || lagNums.empty() || mat.size() != lagNums.size()) {
    return {};
  }

  const size_t numVars = mat.size();
  const size_t numUnits = mat[0].size();

  // ==== Identify unique lag values (only those present in lagNums) ====
  std::unordered_set<int> uniqueLagsSet;
  for (int lag : lagNums) {
    if (lag >= 0) uniqueLagsSet.insert(lag);
  }
  std::vector<int> uniqueLags(uniqueLagsSet.begin(), uniqueLagsSet.end());
  std::sort(uniqueLags.begin(), uniqueLags.end());

  // ==== Precompute required neighbor sets (lag and lag-1 for each unique lag) ====
  std::unordered_map<int, std::vector<std::vector<int>>> laggedNeighborMap;
  for (int lag : uniqueLags) {
    // Compute lagNum neighbors
    laggedNeighborMap[lag] = CppLaggedNeighbor4Lattice(nb, lag);

    // Compute lagNum - 1 neighbors if needed and not already computed
    if (lag > 0 && laggedNeighborMap.find(lag - 1) == laggedNeighborMap.end()) {
      laggedNeighborMap[lag - 1] = CppLaggedNeighbor4Lattice(nb, lag - 1);
    }
  }

  // ==== Derive pure neighbor sets only for lag values appearing in lagNums ====
  std::unordered_map<int, std::vector<std::vector<int>>> pureLagNeighborMap;

  for (int lag : uniqueLags) {
    if (lag == 0) {
      // Lag 0 → the unit itself
      pureLagNeighborMap[lag] = laggedNeighborMap[lag];
      continue;
    }

    const auto& curr = laggedNeighborMap[lag];
    const auto& prev = laggedNeighborMap[lag - 1];
    std::vector<std::vector<int>> pure(curr.size());

    for (size_t i = 0; i < curr.size(); ++i) {
      std::unordered_set<int> prevSet(prev[i].begin(), prev[i].end());
      std::vector<int> diff;

      for (int idx : curr[i]) {
        if (prevSet.find(idx) == prevSet.end()) {
          diff.push_back(idx);
        }
      }

      if (diff.empty()) {
        diff.push_back(std::numeric_limits<int>::min());
      }
      pure[i] = diff;
    }

    pureLagNeighborMap[lag] = pure;
  }

  // ==== Compute lagged mean values for each variable ====
  std::vector<std::vector<double>> result(numVars, std::vector<double>(numUnits));

  for (size_t varIdx = 0; varIdx < numVars; ++varIdx) {
    const std::vector<double>& values = mat[varIdx];
    int lag = lagNums[varIdx];

    // Handle invalid lag
    if (lag < 0) {
      std::fill(result[varIdx].begin(), result[varIdx].end(),
                std::numeric_limits<double>::quiet_NaN());
      continue;
    }

    const auto& lagNeighbors = pureLagNeighborMap[lag];

    for (size_t unitIdx = 0; unitIdx < numUnits; ++unitIdx) {
      const std::vector<int>& neighbors = lagNeighbors[unitIdx];

      // If no valid neighbors, set NaN
      if (neighbors.size() == 1 && neighbors[0] == std::numeric_limits<int>::min()) {
        result[varIdx][unitIdx] = std::numeric_limits<double>::quiet_NaN();
        continue;
      }

      double sum = 0.0;
      int validCount = 0;

      for (int nbIdx : neighbors) {
        if (nbIdx >= 0 && nbIdx < static_cast<int>(numUnits)) {
          double val = values[nbIdx];
          if (!std::isnan(val)) {
            sum += val;
            ++validCount;
          }
        }
      }

      result[varIdx][unitIdx] = (validCount > 0)
        ? (sum / validCount)
        : std::numeric_limits<double>::quiet_NaN();
    }
  }

  return result;
}

/**
 * Computes the lagged values for each element in a grid matrix based on a specified lag number and Moore neighborhood.
 * For each element in the matrix, the function calculates the values of its neighbors at a specified lag distance
 * in each of the 8 directions of the Moore neighborhood. If a neighbor is out of bounds, it is assigned a NaN value.
 *
 * Parameters:
 *   mat    - A 2D vector representing the grid data.
 *   lagNum - The number of steps to lag when considering the neighbors in the Moore neighborhood.
 *
 * Returns:
 *   A 2D vector containing the lagged values for each element in the grid, arranged by the specified lag number.
 *   If a neighbor is out of bounds, it is filled with NaN.
 *
 * Note:
 *   The return value for each element is the lagged value of the neighbors, not the index of the neighbor.
 */
std::vector<std::vector<double>> CppLaggedVal4Grid(
    const std::vector<std::vector<double>>& mat,
    int lagNum
) {
  // Validate input
  if (mat.empty() || mat[0].empty() || lagNum < 0) {
    return {};
  }

  const int rows = mat.size();
  const int cols = mat[0].size();
  const int numCells = rows * cols;
  const int numNeighbors = 8 * lagNum;

  // If lagNum is 0, return the current values of gird row by row
  if (lagNum == 0) {
    std::vector<std::vector<double>> result;
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        result.push_back({mat[i][j]});
      }
    }
    return result;
  }

  // Generate all valid offsets for the given lagNum (Queen's case)
  std::vector<std::pair<int, int>> offsets;
  for (int dx = -lagNum; dx <= lagNum; ++dx) {
    for (int dy = -lagNum; dy <= lagNum; ++dy) {
      if (std::max(std::abs(dx), std::abs(dy)) == lagNum) {
        offsets.emplace_back(dx, dy);
      }
    }
  }

  // Initialize result with NaN
  std::vector<std::vector<double>> result(
      numCells,
      std::vector<double>(numNeighbors, std::numeric_limits<double>::quiet_NaN())
  );

  // Populate neighbor values
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      const int cellIndex = i * cols + j;
      for (size_t k = 0; k < offsets.size(); ++k) {
        const auto& [dx, dy] = offsets[k];
        const int ni = i + dx;
        const int nj = j + dy;
        if (ni >= 0 && ni < rows && nj >= 0 && nj < cols) {
          result[cellIndex][k] = mat[ni][nj];
        }
        // Else remains NaN
      }
    }
  }

  return result;
}

/**
 * Computes lagged mean values for grid structures using Moore neighborhood (queen's case) at a specified lag distance.
 * For each grid cell, this function aggregates neighbor values at exactly lagNum steps away in all 8 directions,
 * calculates the mean, and ignores NaN values. If all neighbor values are invalid (NaN), the result is NaN.
 *
 * Parameters:
 *   mat     - Input grid (2D matrix) of spatial observations.
 *   lagNum  - Non-negative integer specifying the lag distance.
 *
 * Returns:
 *   A vector (flattened grid) where each element represents the lagged mean value for the corresponding grid cell.
 */
std::vector<double> GenGridLagUni(
    const std::vector<std::vector<double>>& mat,
    int lagNum
) {
  if (mat.empty() || mat[0].empty() || lagNum < 0) {
    return {};
  }

  // Get lagged values via existing utility
  std::vector<std::vector<double>> laggedValues = CppLaggedVal4Grid(mat, lagNum);
  std::vector<double> laggedMeans;

  for (const auto& neighborVals : laggedValues) {
    double sum = 0.0;
    int count = 0;

    for (double val : neighborVals) {
      if (!std::isnan(val)) {
        sum += val;
        ++count;
      }
    }

    if (count > 0) {
      laggedMeans.push_back(sum / count);
    } else {
      laggedMeans.push_back(std::numeric_limits<double>::quiet_NaN());
    }
  }

  return laggedMeans;
}

/**
 * Computes lagged mean values for multiple grid variables using Moore neighborhood (queen's case).
 * This function processes a 3D array of grid data where each element represents a separate grid variable,
 * and computes lagged means for all variables at the specified lag order simultaneously.
 *
 *
 * Parameters:
 *   arr     - 3D vector input where arr[varIndex] is a 2D grid matrix (same as GenGridLagUni input).
 *   lagNum  - Non-negative integer specifying the lag distance in Moore neighborhood.
 *
 * Returns:
 *   A 2D vector where result[varIndex] contains flattened lagged mean values for the corresponding
 *   input variable, arranged in row-major order (consistent with GenGridLagUni output format).
 *
 * Note:
 *   - All input grids must have identical dimensions
 *   - NaN values in input data are ignored during mean calculation
 *   - Returns empty vector if inputs are invalid or dimensions mismatch
 */
std::vector<std::vector<double>> GenGridLagMultiSingle(
    const std::vector<std::vector<std::vector<double>>>& arr,
    int lagNum) {

  // Validate input parameters
  if (arr.empty() || lagNum < 0) {
    return {};
  }

  // Reserve space for results: one row per subset
  std::vector<std::vector<double>> results;
  results.reserve(arr.size());

  // Compute lagged means for each subset independently
  for (const auto& mat : arr) {
    // Each mat is a 2D grid; GenGridLagUni returns a 1D vector of length = rows*cols
    std::vector<double> subsetLagMeans = GenGridLagUni(mat, lagNum);
    results.emplace_back(std::move(subsetLagMeans));
  }

  return results;
}

/**
 * @title Multi-grid lagged mean computation with variable lag orders
 *
 * @description
 * Computes lagged mean values for multiple grid variables (2D matrices) using
 * the Moore neighborhood (queen’s case), allowing each variable to have a distinct lag order.
 *
 * Each grid is processed independently using its specified lag order.
 * For each grid cell, neighbors exactly `lagNum` steps away (in all 8 directions)
 * are averaged to compute the lagged mean. NaN values are ignored in mean calculation.
 *
 * @param arr
 *   A 3D vector representing multiple grid variables.
 *   - arr[v] is a 2D vector (rows × cols) of spatial observations for variable v.
 *   - All grids must have the same dimensions.
 *
 * @param lagNums
 *   A vector of integers specifying the lag order for each variable.
 *   - Must be the same length as `arr`.
 *   - Each lag value must be non-negative.
 *
 * @return
 *   A 2D vector where `result[v]` contains flattened lagged mean values
 *   for the v-th grid variable (in row-major order).
 *   - If inputs are invalid or inconsistent, returns an empty vector.
 *   - If all neighbor values for a cell are NaN, the output for that cell is NaN.
 *
 * @note
 *   - Each grid uses its own lag order (flexible lag specification).
 *   - This function generalizes `GenGridLagMultiSingle()` to support per-variable lag orders.
 *   - All grids must be rectangular and share identical dimensions.
 */
std::vector<std::vector<double>> GenGridLagMulti(
    const std::vector<std::vector<std::vector<double>>>& arr,
    const std::vector<int>& lagNums
) {
  // --- Validate inputs ---
  if (arr.empty()) {
    return {};
  }
  if (arr.size() != lagNums.size()) {
    return {};
  }

  // // --- Validate grid dimensions ---
  // const int rows = arr[0].size();
  // const int cols = arr[0].empty() ? 0 : arr[0][0].size();
  // if (rows == 0 || cols == 0) {
  //   return {};
  // }
  //
  // for (size_t v = 1; v < arr.size(); ++v) {
  //   if (arr[v].size() != static_cast<size_t>(rows) ||
  //       (rows > 0 && arr[v][0].size() != static_cast<size_t>(cols))) {
  //     // Inconsistent grid dimensions across variables
  //     return {};
  //   }
  // }

  // --- Compute lagged means per variable ---
  std::vector<std::vector<double>> results;
  results.reserve(arr.size());

  for (size_t v = 0; v < arr.size(); ++v) {
    int lag = lagNums[v];
    if (lag < 0) {
      // Invalid lag, append empty vector
      results.emplace_back();
      continue;
    }

    // Compute lagged means for this variable
    std::vector<double> laggedMeans = GenGridLagUni(arr[v], lag);
    results.emplace_back(std::move(laggedMeans));
  }

  return results;
}

} // namespace SpatialLagging
