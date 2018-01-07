#' @include RNAmod.R
NULL

# these are used just internally and can be switched to different package

is_object_class <- function(x, .xname, className){
  if( class(x) != className ){
    return(false("%s is not of an object of class '%s'.", .xname, className))
  }
  if( !validObject(x) ){
    return(false("%s is not a valid object of class '%s'.", .xname, className))
  }
  TRUE
}

# GRanges

is_GRanges <- function(x, .xname = get_name_in_parents(x)){
  is_object_class(x, .xname, "GRanges")
}

assert_is_GRanges <- function(x, 
                              severity = getOption("assertive.severity","stop")){
  assert_engine(is_GRanges, 
                x, 
                .xname = get_name_in_parents(x),
                severity = severity)
}

assert_all_are_GRanges <- function(x, 
                                   severity = getOption("assertive.severity","stop")){
  .xname <- get_name_in_parents(x)
  msg <- gettextf("%s are not all all GRanges objects.", .xname)
  assert_engine(is_GRanges, 
                x, 
                .xname = .xname, 
                msg = msg,
                severity = severity)
}

# FaFile

is_FaFile <- function(x, .xname = get_name_in_parents(x)){
  is_object_class(x, .xname, "FaFile")
}

assert_is_FaFile <- function(x, 
                              severity = getOption("assertive.severity","stop")){
  assert_engine(is_GRanges, 
                x, 
                .xname = get_name_in_parents(x),
                severity = severity)
}

assert_all_are_FaFile <- function(x, 
                                   severity = getOption("assertive.severity","stop")){
  .xname <- get_name_in_parents(x)
  msg <- gettextf("%s are not all all GRanges objects.", .xname)
  assert_engine(is_GRanges, 
                x, 
                .xname = .xname, 
                msg = msg,
                severity = severity)
}