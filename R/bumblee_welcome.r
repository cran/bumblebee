.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to bumblebee \n see https://magosil86.github.io/bumblebee/ for examples and documentation")
}


.onLoad <- function(libname, pkgname) {
  op <- options()
  op.bumblebee <- list(
    bumblebee.path = "",
    bumblebee.install.args = "",
    bumblebee.name = "Lerato E. Magosi",
    bumblebee.desc.author = '"Lerato E. Magosi <magosil86@gmail.com> [aut, cre]"',
    bumblebee.desc.license = "MIT",
    bumblebee.desc.suggests = NULL,
    bumblebee.desc = list()
  )
  toset <- !(names(op.bumblebee) %in% names(op))
  if(any(toset)) options(op.bumblebee[toset])

  invisible()
}
