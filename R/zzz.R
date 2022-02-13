
#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {
    
    if (!interactive()) return()
    
    packageStartupMessage(
        sprintf('BANKSY version %s', packageVersion(pkg = 'Banksy'))
    )
    
}
