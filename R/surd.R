.surd_lattice = \(data, target, agents, lag = 1, bin = 5, max.combs = NULL, cores = 1, nb = NULL){
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  data = sf::st_drop_geometry(data)
  obs = cbind(data[,target,drop = TRUE],
              RcppGenLatticeLagMulti(as.matrix(data[,agents,drop = FALSE])),nb,lag)

}

.surd_grid = \(data, target, agents, lag = 1, bin = 5, max.combs = NULL, cores = 1){
  obs = cbind(terra::values(data[[target]], mat = TRUE, na.rm = FALSE),
              RcppGenGridLagMulti(terra::values(data[[agents]], mat = TRUE, na.rm = FALSE)),
                                  terra::nrow(data), lag)

}

.surd_ts = \(data, target, agents, lag = 1, bin = 5, max.combs = NULL, cores = 1){
  obs = cbind(data[,target,drop = TRUE],
              RcppGenTSLagMulti(as.matrix(data[,agents,drop = FALSE])),lag)

}
