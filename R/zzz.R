.pkgenv <- new.env(parent = emptyenv())

.onAttach <- function(...){
  #date <- date()
  #x <- regexpr("[0-9]{4}", date)
  #this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)

  #this.version = utils::packageVersion("fastqq")

  packageStartupMessage("** ----------------------------------------------------------------- **")
  packageStartupMessage("** fastqq")
  packageStartupMessage("**  - Faster QQ plots for large number of points")
  packageStartupMessage("**")
  #packageStartupMessage("** Version    : ",this.version,"       (",this.year,")",sep="")
  packageStartupMessage("** Maintainer : Gudmundur Einarsson  (gudmundur.einarsson.phd@gmail.com)")
  packageStartupMessage("**")
  packageStartupMessage("** Please file issues on github.com/gumeo/fastqq.")
  packageStartupMessage("** If you use this package, please cite it.")
  packageStartupMessage("** ----------------------------------------------------------------- **")
}

.onUnload <- function (libpath) {
  library.dynam.unload("fastqq", libpath)
}
