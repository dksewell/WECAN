.onLoad <- function(libname, pkgname) {
  if (nzchar(system.file("libs/x64", "WECAN.dll", package = "WECAN"))) {
    library.dynam("WECAN", pkgname, libname)
  }
}
