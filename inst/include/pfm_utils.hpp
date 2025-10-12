#ifndef PFM_UTILS_HPP
#define PFM_UTILS_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <utility>
#include <unordered_map>
#include <stdexcept>
#include <cstddef>

/**
 * @brief Header-only utilities to convert discretized matrix (per-variable columns)
 *        into a joint probability frequency matrix (PFM).
 *
 * Core function:
 *   std::pair<std::vector<double>, std::vector<int>> DiscMat2PFM(
 *       const std::vector<std::vector<double>>& x,
 *       double eps = 1e-14);
 *
 * Namespace: pfmutils
 *
 * Behaviour summary:
 *  - `x` is outer: variables (nvars), inner: samples (nsamples).
 *  - Rows with any NaN (std::isnan) are skipped.
 *  - For each variable, unique category values are discovered and mapped to
 *    contiguous indices [0..K-1] deterministically (sorted order).
 *  - The joint counts are accumulated into a flattened vector in row-major order
 *    where the first variable is the fastest-changing dimension (compatible with R array dim ordering).
 *  - A small floor `eps` is applied to zero-count cells to avoid exact zeros (optional).
 *  - The histogram is normalized to sum to 1 and returned as probabilities.
 */

namespace pfmutils {

inline std::pair<std::vector<double>, std::vector<int>>
  DiscMat2PFM(const std::vector<std::vector<double>>& x, double eps = 1e-14) {
    // Validate shape
    const std::size_t nvars = x.size();
    if (nvars == 0) {
      throw std::invalid_argument("DiscMat2PFM: input must contain at least one variable.");
    }
    const std::size_t nsamples = x[0].size();
    if (nsamples == 0) {
      throw std::invalid_argument("DiscMat2PFM: each variable must contain at least one sample.");
    }
    for (std::size_t j = 1; j < nvars; ++j) {
      if (x[j].size() != nsamples) {
        throw std::invalid_argument("DiscMat2PFM: all variables must have the same number of samples.");
      }
    }

    // Determine valid sample indices (skip rows where any variable is NaN)
    std::vector<std::size_t> valid_idx;
    valid_idx.reserve(nsamples);
    for (std::size_t i = 0; i < nsamples; ++i) {
      bool ok = true;
      for (std::size_t j = 0; j < nvars; ++j) {
        if (std::isnan(x[j][i])) { ok = false; break; }
      }
      if (ok) valid_idx.push_back(i);
    }
    if (valid_idx.empty()) {
      throw std::invalid_argument("DiscMat2PFM: no valid rows after removing NaNs.");
    }

    // For each variable (column), compute sorted unique values and a map value->index
    std::vector<std::vector<double>> unique_vals(nvars);
    std::vector<std::unordered_map<double, int>> maps(nvars);
    std::vector<int> dims(nvars);

    for (std::size_t j = 0; j < nvars; ++j) {
      // collect values for valid rows
      unique_vals[j].reserve(valid_idx.size());
      for (std::size_t idx : valid_idx) unique_vals[j].push_back(x[j][idx]);

      // sort and unique
      std::sort(unique_vals[j].begin(), unique_vals[j].end());
      auto it = std::unique(unique_vals[j].begin(), unique_vals[j].end());
      unique_vals[j].erase(it, unique_vals[j].end());

      // build map value -> contiguous index
      maps[j].reserve(unique_vals[j].size() * 2 + 1);
      for (int k = 0; k < static_cast<int>(unique_vals[j].size()); ++k) {
        maps[j].emplace(unique_vals[j][k], k);
      }

      dims[j] = static_cast<int>(unique_vals[j].size());
      if (dims[j] <= 0) {
        throw std::runtime_error("DiscMat2PFM: variable has zero unique categories.");
      }
    }

    // Optionally you could check all dims equal; here we do not require equality,
    // but you can enable a check if needed:
    // bool all_equal = std::all_of(dims.begin()+1, dims.end(), [&](int v){ return v == dims[0]; });

    // Compute total number of cells and row-major strides (first var fastest)
    std::size_t total_cells = 1;
    for (std::size_t j = 0; j < nvars; ++j) {
      total_cells *= static_cast<std::size_t>(dims[j]);
    }
    if (total_cells == 0) {
      throw std::runtime_error("DiscMat2PFM: total joint state space size is zero.");
    }

    std::vector<std::size_t> strides(nvars);
    strides[0] = 1;
    for (std::size_t j = 1; j < nvars; ++j) {
      strides[j] = strides[j - 1] * static_cast<std::size_t>(dims[j - 1]);
    }

    // Allocate histogram and count occurrences
    std::vector<double> hist(total_cells, 0.0);
    for (std::size_t s = 0; s < valid_idx.size(); ++s) {
      std::size_t i = valid_idx[s];
      std::size_t lin = 0;
      for (std::size_t j = 0; j < nvars; ++j) {
        double v = x[j][i];
        auto it_map = maps[j].find(v);
        if (it_map == maps[j].end()) {
          // Should not happen because maps were built from valid rows
          throw std::runtime_error("DiscMat2PFM: value not found in mapping (unexpected).");
        }
        lin += static_cast<std::size_t>(it_map->second) * strides[j];
      }
      hist[lin] += 1.0;
    }

    // Apply small eps floor and normalize
    double sum = 0.0;
    for (std::size_t k = 0; k < total_cells; ++k) {
      if (hist[k] <= eps) hist[k] = eps;
      sum += hist[k];
    }
    if (!(sum > 0.0)) {
      throw std::runtime_error("DiscMat2PFM: invalid histogram sum.");
    }
    for (std::size_t k = 0; k < total_cells; ++k) hist[k] /= sum;

    return std::make_pair(std::move(hist), std::move(dims));
  }

} // namespace pfmutils

#endif // PFM_UTILS_HPP
