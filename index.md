# infocausality

[![infocausality website:
https://stscl.github.io/infocausality/](reference/figures/infocausality.png)](https://stscl.github.io/infocausality/)

*Information-Theoretic Measure of Causality*

`infocausality` is an R package for information-theoretic causal
analysis.  
It quantifies temporal and spatial causality through information flow,
and decomposes it into unique, redundant, and synergistic components.
The package provides native support for `data.frame`, `sf`, and
`SpatRaster` objects, offering a unified interface for both time-series
and spatial cross-sectional causal analysis.

> *Refer to the package documentation
> <https://stscl.github.io/infocausality/> for more detailed
> information.*

> ⚠️ **Note**: The SURD (Synergistic-Unique-Redundant Decomposition)
> core computations in `infocausality` are executed via Python bindings.
> A pure C++ implementation with improved performance and easier
> deployment has been developed in the
> [`infoxtr`](https://github.com/stscl/infoxtr) package. For new
> projects, we recommend using `infoxtr` as a drop-in replacement.

## Installation

- Install from [CRAN](https://CRAN.R-project.org/package=infocausality)
  with:

``` r
install.packages("infocausality", dep = TRUE)
```

- Install binary version from
  [R-universe](https://stscl.r-universe.dev/infocausality) with:

``` r
install.packages("infocausality",
                 repos = c("https://stscl.r-universe.dev",
                           "https://cloud.r-project.org"),
                 dep = TRUE)
```

- Install from source code on
  [GitHub](https://github.com/stscl/infocausality) with:

``` r
if (!requireNamespace("devtools")) {
    install.packages("devtools")
}
devtools::install_github("stscl/infocausality",
                         build_vignettes = TRUE,
                         dep = TRUE)
```
