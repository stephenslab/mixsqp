# This is called when the package is first loaded.
.onLoad <- function (libname, pkgname) {
  options(mixsqp.debug.mode = FALSE)
  options(mixsqp.debug.file = "mixsqp.RData")
  invisible()
}
