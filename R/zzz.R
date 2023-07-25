#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {
    if (!interactive()) {
        return()
    }

    packageStartupMessage(
        sprintf("Banksy version %s", packageVersion(pkg = "Banksy"))
    )
}
