#' @include RNAmod.R
NULL

#' @rdname parseMod
#' 
#' @description 
#' \code{mod_m7G}
#' 
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "parseMod",
  signature = signature(object = "mod_m7G"),
  definition = function(object,
                        data) {
    
    browser()
    
    
    print("m7G")
  }
)

#' @rdname parseMod
#' 
#' @description 
#' \code{mod_D}
#'
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "parseMod",
  signature = signature(object = "mod_D"),
  definition = function(object,
                        data) {
    
    
    print("D")
    
    
  }
)