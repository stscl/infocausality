.onLoad = function(...) {
  loadNamespace("sf")
  loadNamespace("terra")
  reticulate::py_require(c("numpy","joblib"))
}
