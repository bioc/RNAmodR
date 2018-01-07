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


#' @describeIn getSummarizedExperiment
#' 
#' @param .Object an RNAmod object 
#' @param number a number defining the experiment. For 
#' getSummarizedExperiments() more than number can be defined as numeric vector.
#' @param modification name of modification, one or more character 
#'
#' @return \code{getSummarizedExperiment()}: SummarizedExperiment
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' getSummarizedExperiment(experiment, 1)
#' getSummarizedExperiments(experiment, c(1, 2))
#' }
setMethod(
  f = "getSummarizedExperiment", 
  signature = signature(.Object = "RNAmod", 
                        number = "numeric",
                        modification = "character"),
  definition = function(.Object, 
                        number,
                        modification) {
    assertive::assert_all_are_non_empty_character(modification)
    if( !assertive::is_scalar(number)) 
    {
      stop(paste0("Multiple numbers provided. only one excepted"))
    }
    experiment <- getExperimentData(.Object, number)
    
    if( assertive::is_a_bool(experiment)) if( assertive::is_false(experiment) ) 
    {
      stop(paste0("Incorrect experiment identifier given: ",number))
    }
    return(.loadSummarizedExperiment(.Object, 
                                     experiment, 
                                     modification, 
                                     failOnNonExist = TRUE))
  }
)
#' @describeIn getSummarizedExperiment
#'
#' @return \code{getSummarizedExperiments()}: list of SummarizedExperiments
#' 
#' @export
setMethod(
  f = "getSummarizedExperiments", 
  signature = signature(.Object = "RNAmod",
                        number = "numeric",
                        modification = "character"),
  definition = function(.Object, 
                        number,
                        modification) {
    if( length(number) == 0 ) stop("no experiment number given")
    assertive::assert_all_are_whole_numbers(number)
    assertive::assert_all_are_non_empty_character(modification)
    
    ses <- lapply(number, function(x,.Object)
      {
        getSummarizedExperiment(.Object, 
                                x)
      },
      .Object)
    names(ses) <- number
    return(ses)
  }
)


#' @title .loadSummarizedExperiment
#' 
#' @description
#' Loads saved SummarizedExperiment from file
#' 
#' @param .Object an RNAmod object 
#' @param experiment a list containing data for one experiment
#' @param modification name of modification, one or more character 
#' @param failOnNonExist a logical value, whether a error should be thrown, if 
#' the file does not exist or just return FALSE for further evaluation
#'
#' @return SummarizedExperiment
setMethod(
  f = ".loadSummarizedExperiment", 
  signature = signature(.Object = "RNAmod", 
                        experiment = "list",
                        modification = "character"),
  definition = function(.Object, 
                        experiment,
                        modification,
                        failOnNonExist) {
    assertive::assert_all_are_non_empty_character(modification)
    assertive::assert_is_a_bool(failOnNonExist)
    
    # Setting se default
    se <- NULL
    
    folder <- fileName <- paste0(getOutputFolder(.Object),
                                 "SE/",
                                 "RNAmod_",
                                 unique(experiment["SampleName"]),
                                 "_")
    fileNames <- paste0(folder,
                        modification,
                        ".RData")
    
    ses <- lapply(fileNames,function(fileName){
      if(assertive::is_existing_file(fileName)) {
        load(fileName)
        return(se)
      } else {
        if( failOnNonExist ) {
          stop("SummarizedExperiment file can not be loaded, since it does not ",
               "exist at the expected location:\n",
               fileName,
               call. = FALSE)
        }
      }
      return(FALSE)
    })
    
    if( any(vapply(ses)) ){
      
    }
    
    se <- .merge_se_for_modifications(ses,modification)
    return(se)
  }
)

.merge_se_for_modifications <- function(ses,modification){
  return(ses[[1]])
}


#' @title .saveSummarizedExperiments
#' 
#' @description
#' Saves SummarizedExperiment to file
#'
#' @param .Object an RNAmod object 
#' @param se a RangedSummarizedExperiment object 
#' @param experiment a list containing data for one experiment
#' @param modification name of modification, one or more character 
#'
#' @return SummarizedExperiment.
#' 
#' @import SummarizedExperiment
setMethod(
  f = ".saveSummarizedExperiments", 
  signature = signature(.Object = "RNAmod", 
                        se = "RangedSummarizedExperiment", 
                        experiment = "list",
                        modification = "character"),
  definition = function(.Object, 
                        se, 
                        experiment,
                        modification) {
    
    folder <- fileName <- paste0(getOutputFolder(.Object),
                                 "SE/")
    if(!assertive::is_dir(folder)){
      dir.create(folder, recursive = TRUE)
    }
    fileNames <- paste0(folder,
                        "RNAmod_",
                        unique(experiment["SampleName"]),
                        "_",
                        modification,
                        ".RData")
    
    .save_single_ses(.split_se_for_each_modification(se,modification), fileNames)
    return(se)
  }
)

.split_se_for_each_modification <- function(se,modification){
  return(list(rep(se,length(modification))))
}

.save_single_ses <- function(ses,fileNames){
  for(i in seq_along(ses)){
    se <- ses[[i]]
    save(se, file = fileNames[[i]])
  }
}


#' @describeIn setSummarizedExperiment
#'
#' @param .Object a RNAmod object 
#' @param se a RangedSummarizedExperiment object
#' @param number a number defining the experiment
#' @param modification name of modification, one or more character 
#' 
#' @return the RpfSummarizedExperiment
#' @export
#'
#' @examples
#' \donttest{
#' setSummarizedExperiment(experiment, se, 1, c("m7G"))
#' }
setMethod(
  f = "setSummarizedExperiment", 
  signature = signature(.Object = "RNAmod", 
                        se = "RangedSummarizedExperiment", 
                        number = "numeric",
                        modification = "character"),
  definition = function(.Object, 
                        se, 
                        number,
                        modification) {
    assertive::assert_all_are_non_empty_character(modification)
    
    experiment <- getExperimentData(.Object, number)
    
    if( assertive::is_a_bool(experiment)) {
      if( assertive::is_false(experiment) ) {
        stop("Incorrect experiment identifier given: ",
             number,
             call. = FALSE)
      }
    }
    
    se <- .saveSummarizedExperiments(.Object, se, experiment, modification)
    return(invisible(se))
  }
)

#' 
#' 
#' #' @rdname getGff
#' #' 
#' #' @title getGff
#' #'
#' #' @param .Object RNAmod. 
#' #' @param number numeric. 
#' #' @param modification character. 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' setMethod(
#'   f = "getGff", 
#'   signature = signature(.Object = "RNAmod", 
#'                         number = "numeric",
#'                         modification = "character"),
#'   definition = function(.Object, 
#'                         number,
#'                         modification) {
#'     
#'   }
#' )




#' @title .saveGff
#' 
#' @description
#' Saves GRanges object as gff3 file file
#'
#' @param .Object an RNAmod object 
#' @param gff a GRanges object 
#' @param experiment a list containing data for one experiment
#' @param modification name of modification, one or more character 
#'
#' @return GRanges object
setMethod(
  f = ".saveGff", 
  signature = signature(.Object = "RNAmod", 
                        gff = "GRanges", 
                        experiment = "list",
                        modification = "character"),
  definition = function(.Object, 
                        gff, 
                        experiment,
                        modification) {
    
    folder <- fileName <- paste0(getOutputFolder(.Object),
                                 "gff/")
    if(!assertive::is_dir(folder)){
      dir.create(folder, recursive = TRUE)
    }
    fileNames <- paste0(folder,
                        "RNAmod_",
                        unique(experiment["SampleName"]),
                        "_",
                        modification,
                        ".gff3")
    .save_single_gffs(.split_gff_for_each_modification(gff,modification),
                      fileNames)
    return(gff)
  }
)

.split_gff_for_each_modification <- function(gff,modification){
  return(list(rep(gff,length(modification))))
}


.save_single_gffs <- function(gffs,fileNames){
  for(i in seq_along(gffs)){
    gff <- gffs[[i]]
    rtracklayer::export.gff3(gff, con = fileNames[[i]])
  }
}

#' @rdname setGff
#'
#' @title setGff
#'
#' @param .Object RNAmod.
#' @param gff GRanges.
#' @param number numeric.
#' @param modification character.
#'
#' @return GRanges object
#' @export
#'
#' @examples
#' \donttest{
#' setGff(experiment, gff, 1, c("m7G"))
#' }
setMethod(
  f = "setGff",
  signature = signature(.Object = "RNAmod",
                        gff = "GRanges",
                        number = "numeric",
                        modification = "character"),
  definition = function(.Object,
                        gff,
                        number,
                        modification) {
    assertive::assert_all_are_non_empty_character(modification)
    
    experiment <- getExperimentData(.Object, number)
    
    if( assertive::is_a_bool(experiment)) {
      if( assertive::is_false(experiment) ) {
        stop("Incorrect experiment identifier given: ",
             number,
             call. = FALSE)
      }
    }
    
    gff <- .saveGff(.Object, gff, experiment, modification)
    return(invisible(gff))
  }
)

# load classes for modification analysis
.load_mod_classes <- function(modifications){
  
  modClasses <- vector(mode = "list", length = length(modifications))
  for(i in seq_along(modifications)){
    className <- paste0("mod_",modifications[[i]])
    
    # try to create modification detection classes
    tryCatch(
      class <- new(className),
      error = function(e) stop("Class for detecting ",
                               modifications[[1]],
                               " does not exist (",className,").",
                               call. = FALSE)
    )
    if( !existsMethod("convertReadsToPositions",signature(class(class),
                                                          "numeric",
                                                          "GRanges",
                                                          "DataFrame") ) )
      stop("Function convertReadsToPositions() not defined for ",class(class))
    if( !existsMethod("parseMod",signature(class(class),
                                           "GRanges",
                                           "FaFile",
                                           "list") ) )
      stop("Function parseMod() not defined for ",class(class))
    if( !existsMethod("mergePositionsOfReplicates",signature(class(class),
                                                             "GRanges",
                                                             "FaFile",
                                                             "list") ) )
      stop("Function mergePositionsOfReplicates() not defined for ",class(class))
    modClasses[[i]] <- class
  }
  names(modClasses) <- modifications
  
  return(modClasses)
}