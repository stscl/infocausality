.onLoad = function(...) {
  loadNamespace("sf")
  loadNamespace("terra")
  reticulate::py_require(c("numpy","joblib"))
}


.onAttach = function(...) {
  packageStartupMessage(" Note: infocausality is now superseded by the 'infoxtr' package.")
}
