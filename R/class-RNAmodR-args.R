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
           args = data.frame(),
           files = NULL,
           conditions = NULL
         )
)

#' @rdname RNAmodR-args-class
#' 
#' @param param parameter file(s) to be loaded
#' 
#' @importFrom utils read.delim
#'
#' @export
RNAmodRargs <- function(param,
                        files,
                        conditions){
  # input check
  assertive::assert_all_are_existing_files(param)
  if(!missing(files) && !missing(conditions)){
    if(length(files) != length(conditions)){
      stop("Number of files and conditions do not match",
           call. = FALSE)
    }
    assertive::assert_all_are_existing_files(files)
    .check_sample_conditions(conditions)
  }
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
    .Object@args <- do.call(rbind,lapply(param, 
                                         read.delim, 
                                         comment.char = "#",
                                         stringsAsFactors = FALSE))
    if(!missing(files) && 
       !missing(conditions) && 
       length(files) > 0 && 
       length(conditions) > 0){
      .Object@files <- files
      .Object@conditions <- conditions
    }
    # validation step
    if(length(RNAMODR_PARAM_COL) != length(colnames(.Object@args)) |
      !all(colnames(.Object@args) %in% RNAMODR_PARAM_COL)){
      stop("Parameter file contains invalid columns. Allowed colums: '",
           paste(RNAMODR_PARAM_COL, collapse = "','"),
           "'.")
    }
    modifications <- unique(.Object@args$Identifier)
    for(i in seq_along(modifications)){
      mod <- modifications[i]
      assertive::assert_all_are_numeric_strings(getParam(.Object,
                                                         mod,
                                                         "map_quality"))
    }
    #
    return(.Object)
  }
)