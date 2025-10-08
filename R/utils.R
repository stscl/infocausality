utils_source_python = \(oauth_func_path) {
  module_name = gsub("\\.py$", "", basename(oauth_func_path))
  module_path = dirname(oauth_func_path)
  reticulate::import_from_path(module_name, path = module_path, convert = TRUE)
}

utils_run_surd = \(obs, bin = 5L, max.combs = NULL, cores = 1, backend = "threading"){
  ic = utils_source_python(system.file("python", "InfoCausality.py",
                                       package = "infocausality"))
  ict = ic$InfoCausality(obs, nbins = as.integer(bin))
  invisible(ict$surd(max.combs,cores,backend))
}
