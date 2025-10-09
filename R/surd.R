.surd_ts = \(data, target, agents, lag = 1, bin = 5, max.combs = NULL, cores = 1, backend = "threading"){
  obs = cbind(
    data[,target,drop = TRUE],
    RcppGenTSLagMulti(as.matrix(data[,agents,drop = FALSE]),lag)
    )
  utils_run_surd(obs, bin, max.combs, cores)
}

.surd_lattice = \(data, target, agents, lag = 1, bin = 5, max.combs = NULL, cores = 1, backend = "threading", nb = NULL){
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  data = sf::st_drop_geometry(data)
  obs = cbind(
    data[,target,drop = TRUE],
    RcppGenLatticeLagMulti(as.matrix(data[,agents,drop = FALSE]),nb,lag)
    )
  utils_run_surd(obs, bin, max.combs, cores)
}

.surd_grid = \(data, target, agents, lag = 1, bin = 5, max.combs = NULL, cores = 1, backend = "threading"){
  obs = cbind(
    terra::values(data[[target]],mat = TRUE,na.rm = FALSE),
    RcppGenGridLagMulti(terra::values(data[[agents]],mat = TRUE,na.rm = FALSE),
                        terra::nrow(data), lag)
    )
  utils_run_surd(obs, bin, max.combs, cores)
}

#' synergistic-unique-redundant decomposition of causality
#'
#' @param data observation data.
#' @param target name of the target variable.
#' @param agents names of agent variables.
#' @param lag (optional) lag order.
#' @param bin (optional) number of discretization bins.
#' @param max.combs (optional) maximum combination order. If `NULL`, the standard SURD decomposition is applied.
#' @param cores (optional) number of cores for parallel computation.
#' @param backend (optional) `Joblib` backend: `loky`, `threading`, or `multiprocessing`.
#' @param nb (optional) neighbours list.
#'
#' @return A list.
#' \describe{
#'   \item{synergistic}{Synergistic information components by agent combinations.}
#'   \item{unique}{Unique information contributions per variable.}
#'   \item{redundant}{Redundant information shared by agent subsets.}
#'   \item{mutual_info}{Mutual information measures for each combination.}
#'   \item{info_leak}{Information leak ratio.}
#' }
#'
#' @export
#' @name surd
#' @aliases surd,data.frame-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' surd(columbus,"hoval",c("inc","crime"))
#'
methods::setMethod("surd", "data.frame", .surd_ts)

#' @rdname surd
methods::setMethod("surd", "sf", .surd_lattice)

#' @rdname surd
methods::setMethod("surd", "SpatRaster", .surd_grid)
