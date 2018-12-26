#' @include RNAmodR.R
#' @include Modifier-class.R
NULL

#' @name ModifierSet
#' 
#' @title ModifierSet
#' @description 
#' title
NULL

#' @rdname ModifierSet
#' @export
setClass("ModifierSet",
         contains = c("VIRTUAL",
                      "SimpleList"))