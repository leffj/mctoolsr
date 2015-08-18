# devtools::use_package('tools')
# devtools::use_package('biom')
# devtools::use_package('reshape2')
# devtools::use_package('dplyr')
# devtools::use_package('ggplot2')
# devtools::use_package('vegan')
# devtools::use_package('nlme')
# devtools::use_package('VennDiagram', 'Suggests')

.onAttach = function(libname, pkgname) {
  packageStartupMessage("You're loading mctoolsr. Direct inquiries to:
                        'https://github.com/leffj/mctoolsr'")
}

.onLoad = function(libname, pkgname) {
  op = options()
  op.devtools = list(
    mctoolsr.path = "",
    mctoolsr.install.args = "",
    mctoolsr.name = "mctoolsr",
    mctoolsr.desc.author = '"Jon Leff <jonathan.leff@colorado.edu> [aut, cre]"',
    mctoolsr.desc.license = "GPL-3",
    mctoolsr.desc.suggests = NULL,
    mctoolsr.desc = list()
  )
  toset = !(names(op.devtools) %in% names(op))
  if(any(toset)) options(op.devtools[toset])
  
  invisible()
}