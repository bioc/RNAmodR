#' @include RNAmodR.R
NULL

# maybe use the SYSArgs from systemPipeR in the future as parental class

#' @name RNAmodR-args-class
#' 
#' @title RNAmodR arguments class
#'
#' @description 
#' \code{RNAmodRargs}
#' 
NULL

#' @rdname RNAmodR-args-class
#'
#' @slot files input files for arguments
#' @slot conditions identifier for sample condition. Treated or Control
#' @slot args a data.frame containing the arguments
#'
#' @return a RNAmodRargs object
#' 
#' @export
setClass("RNAmodRargs",
         slots = c(
           files = "character",
           conditions = "character",
           args = "data.frame"),
         prototype = list(
           args = data.frame()
         )
)

#' @rdname RNAmodR-args-class
#'
#' @importFrom utils read.delim
#'
#' @export
RNAmodRargs <- function(param,
                        files,
                        conditions){
  # input check
  assertive::assert_all_are_existing_files(param)
  assertive::assert_all_are_existing_files(files)
  if(length(files) != length(conditions)){
    stop("Number of files and conditions do not match",
         call. = FALSE)
  }
  .check_sample_conditions(conditions)
  # create class
  return(new("RNAmodRargs",
             param,
             files,
             conditions))
}
setMethod(
  f = "initialize", 
  signature = signature(.Object = "RNAmodRargs"),
  definition = function(.Object, 
                        param,
                        files,
                        conditions) {
    # assigning data to slots
    .Object@args <- read.delim(param, comment.char = "#")
    .Object@files <- files
    .Object@conditions <- conditions
    # validation step
    
    #
    return(.Object)
  }
)