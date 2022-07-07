##' @import DAPARdata
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("\nThis is the 'DAPAR' version ",
    utils::packageVersion("DAPAR"), ".\n\n",
    "  To get started, visit\n",
    "    http://www.prostar-proteomics.org/\n\n"))
}
