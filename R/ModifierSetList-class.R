#' @include RNAmodR.R
#' @include ModifierSet-class.R
NULL

#' @name ModifierSetList
#' 
#' @title ModifierSetList
#' @description 
#' title
NULL

#' @rdname ModifierSetList
#' @export
setClass("ModifierSetList",
         contains = c("VIRTUAL",
                      "SimpleList"))