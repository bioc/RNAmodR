#' @include RNAmod.R
NULL

#' @describeIn getExperimentData
#' 
#' @param .Object a RNAmod object
#' @param number optional number to return results for one experiment only
#'
#' @return a list containg the information on the experiment or a data.frame with all the information
#' 
#' @export
#'
#' @examples
#' data("RNAmod_test", package = "RNAmod")
#' 
#' # All experiments
#' getSampleData(experiment)
#' # One experiment as list
#' getSampleData(experiment, 1)
#' # Two experiments as data.frame
#' getSampleData(experiment, c(1,2))
setMethod(
  f = "getExperimentData", 
  signature = signature(.Object = "RNAmod", 
                        number = "ANY"),
  definition = function(.Object, 
                        number){
    return(.get_experiment_data(.Object,number))
  }
)

.get_experiment_data <- function(.Object,number){
  if( missing(number) || assertive::is_identical_to_false(number) ) {
    return(.Object@.dataSamples)
  }
  if( all(vapply(number,assertive::is_a_number,logical(1)))) {
    df <- .Object@.dataSamples
    if (nrow(df[df$ExperimentNo %in% number, ]) == 1) {
      ret <- as.list(df[df$ExperimentNo == number,])
      ret["ExperimentNo"] <- as.numeric(ret["ExperimentNo"])
      ret[!(names(ret) == "ExperimentNo")] <- 
        as.character(ret[!(names(ret) == "ExperimentNo")])
      return(ret)
    }
    if (nrow(df[df$ExperimentNo %in% number, ]) > 1) {
      ret <- df[df$ExperimentNo %in% number,]
      ret$ExperimentNo <- as.numeric(ret$ExperimentNo)
      return(ret)
    }
  }
  stop("Something went wrong. Invalid experiment number!", call. = FALSE)
}