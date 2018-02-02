#' @include RNAmodR.R
NULL

# these are used just internally and can be switched to different package

.is_object_class <- function(x, .xname, className){
  if( class(x) != className ){
    return(false("%s is not of an object of class '%s'.", .xname, className))
  }
  if( !validObject(x) ){
    return(false("%s is not a valid object of class '%s'.", .xname, className))
  }
  TRUE
}

# GRanges

is_GRanges <- function(x, .xname = assertive.base::get_name_in_parent(x)){
  .is_object_class(x, .xname, "GRanges")
}

assert_is_GRanges <- function(x, 
                              severity = getOption("assertive.severity","stop")){
  assert_engine(is_GRanges, 
                x, 
                .xname = assertive.base::get_name_in_parent(x),
                severity = severity)
}

assert_all_are_GRanges <- function(x, 
                                   severity = getOption("assertive.severity","stop")){
  .xname <- assertive.base::get_name_in_parent(x)
  msg <- gettextf("%s are not all all GRanges objects.", .xname)
  assert_engine(is_GRanges, 
                x, 
                .xname = .xname, 
                msg = msg,
                severity = severity)
}

# FaFile

is_FaFile <- function(x, .xname = assertive.base::get_name_in_parent(x)){
  .is_object_class(x, .xname, "FaFile")
}

assert_is_FaFile <- function(x, 
                             severity = getOption("assertive.severity","stop")){
  assert_engine(is_FaFile, 
                x, 
                .xname = assertive.base::get_name_in_parent(x),
                severity = severity)
}

assert_all_are_FaFile <- function(x, 
                                  severity = getOption("assertive.severity","stop")){
  .xname <- assertive.base::get_name_in_parent(x)
  msg <- gettextf("%s are not all all FaFile objects.", .xname)
  assert_engine(is_FaFile, 
                x, 
                .xname = .xname, 
                msg = msg,
                severity = severity)
}



#' @name is_TxDb
#' 
#' @title assertive: TxDb
#' 
#' @description 
#' \code{is_TxDb} Check 
#' if input is a TxDb object
#' 
#' @param x object to test
#' @param .xname name of the object in parent env
#' @param severity assertions only: severity of failing assertion
#'
#' @return logical
#' @export
is_TxDb <- function(x, .xname = assertive.base::get_name_in_parent(x)){
  .is_object_class(x, .xname, "TxDb")
}

# SummarizedExperiment

#' @name is_SummarizedExperiment
#' 
#' @title assertive: SummarizedExperiment
#' 
#' @description 
#' \code{is_SummarizedExperiment}/\code{assert_is_SummarizedExperiment}: Check 
#' if input is a SummarizedExperiment object
#' \code{assert_all_are_SummarizedExperiment}: Check if input is a list of 
#' SummarizedExperiment objects
#' 
#' @param x object to test
#' @param .xname name of the object in parent env
#' @param severity assertions only: severity of failing assertion
#'
#' @return logical
#' @export
is_SummarizedExperiment <- function(x, .xname = assertive.base::get_name_in_parent(x)){
  .is_object_class(x, .xname, "SummarizedExperiment")
}

#' @rdname is_SummarizedExperiment
#' 
#' @export
#' @return invisible logical
assert_is_SummarizedExperiment <- function(x, 
                             severity = getOption("assertive.severity","stop")){
  assert_engine(is_SummarizedExperiment, 
                x, 
                .xname = assertive.base::get_name_in_parent(x),
                severity = severity)
}

#' @rdname is_SummarizedExperiment
#' 
#' @export
#' @return invisible logical
assert_all_are_SummarizedExperiment <- function(x, 
                                  severity = getOption("assertive.severity","stop")){
  .xname <- assertive.base::get_name_in_parent(x)
  msg <- gettextf("%s are not all all SummarizedExperiment objects.", .xname)
  assert_engine(is_SummarizedExperiment, 
                x, 
                .xname = .xname, 
                msg = msg,
                severity = severity)
}