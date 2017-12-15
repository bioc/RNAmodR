#' @include RNAmod.R
NULL



#' @name mod
#' 
#' @title mod
#'
#' @return
#' @export
#'
#' @examples
setClass("mod",
         contains = "VIRTUAL")

setMethod(
  f = "show", 
  signature = signature(object = "mod"),
  definition = function(object) {
    
  }
)

#' @rdname mod
#'
#' @description 
#' \code{mod_m7G}
#'
#' @return
#' @export
#'
#' @examples
setClass("mod_m7G",
         contains = "mod")

#' @rdname mod
#'
#' @description 
#' \code{mod_D}
#'
#' @return
#' @export
#'
#' @examples
setClass("mod_D",
         contains = "mod")

#' @rdname mod
#'
#' @description 
#' \code{mod_m3C}
#'
#' @return
#' @export
#'
#' @examples
setClass("mod_m3C",
         contains = "mod")