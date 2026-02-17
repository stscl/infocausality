# synergistic-unique-redundant decomposition of causality

synergistic-unique-redundant decomposition of causality

## Usage

``` r
# S4 method for class 'data.frame'
surd(
  data,
  target,
  agents,
  lag = 1,
  bin = 5,
  max.combs = NULL,
  cores = 1,
  backend = "threading"
)

# S4 method for class 'sf'
surd(
  data,
  target,
  agents,
  lag = 1,
  bin = 5,
  max.combs = NULL,
  cores = 1,
  backend = "threading",
  nb = NULL
)

# S4 method for class 'SpatRaster'
surd(
  data,
  target,
  agents,
  lag = 1,
  bin = 5,
  max.combs = NULL,
  cores = 1,
  backend = "threading"
)
```

## Arguments

- data:

  observation data.

- target:

  name of the target variable.

- agents:

  names of agent variables.

- lag:

  (optional) lag order.

- bin:

  (optional) number of discretization bins.

- max.combs:

  (optional) maximum combination order. If `NULL`, the standard SURD
  decomposition is applied.

- cores:

  (optional) number of cores for parallel computation.

- backend:

  (optional) `Joblib` backend: `loky`, `threading`, or
  `multiprocessing`.

- nb:

  (optional) neighbours list.

## Value

A list.

- unique:

  Unique information contributions per variable.

- synergistic:

  Synergistic information components by agent combinations.

- redundant:

  Redundant information shared by agent subsets.

- mutual_info:

  Mutual information measures for each combination.

- info_leak:

  Information leak ratio.

## References

Martinez-Sanchez, A., Arranz, G. & Lozano-Duran, A. Decomposing
causality into its synergistic, unique, and redundant components. Nat
Commun 15, 9296 (2024).

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
# \donttest{
tryCatch(
  surd(columbus, "hoval", c("inc", "crime")),
  error = \(e) message("Skipping Python-dependent example: ", e$message)
)
#> Downloading uv...
#> Done!
# }
```
