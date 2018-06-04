#' @include class-RNAmodR.R
#' @include class-RNAmodR-args.R
NULL


# getters and setters ----------------------------------------------------------

#' @rdname RNAmodR-Accessors
#' @export
setMethod(
  f = "getInputFolder",
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object){
    return(paste0(.Object@.wd, .Object@.inputFolder))
  }
)
#' @rdname RNAmodR-Accessors
#' @export
setMethod(
  f = "getOutputFolder",
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object){
    return(paste0(.Object@.wd, .Object@.outputFolder))
  }
)
#' @rdname RNAmodR-Accessors
#' @export
setMethod(
  f = "getExperimentName",
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object){
    return(.Object@experimentName)
  }
)
#' @rdname RNAmodR-Accessors
#' @export
setMethod(
  f = "getGFFFile",
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object){
    return(.Object@.dataGFF)
  }
)
#' @rdname RNAmodR-Accessors
#' @export
setMethod(
  f = "getFastaFile",
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object){
    return(.Object@.dataFasta)
  }
)

#' @rdname getExperimentData
#' 
#' @title getExperimentData
#' 
#' @description 
#' returns setup data for the experiments based on the input
#' layout csv-file.
#' 
#' @param x a RNAmodR object
#' @param number optional number to return results for one experiment only
#'
#' @return a list containg the information on the experiment or a data.frame 
#' with all the information
#' 
#' @export
#'
#' @examples
#' \donttest{
#' data("RNAmodR_test", package = "RNAmodR")
#' 
#' # All experiments
#' getSampleData(experiment)
#' # One experiment as list
#' getSampleData(experiment, 1)
#' # Two experiments as data.frame
#' getSampleData(experiment, c(1,2))
#' }
setMethod(
  f = "getExperimentData", 
  signature = signature(x = "RNAmodR", 
                        number = "ANY"),
  definition = function(x, 
                        number){
    return(.get_experiment_data(x,
                                number))
  }
)

.get_experiment_data <- function(x,
                                 number){
  if( missing(number) || assertive::is_identical_to_false(number) ) {
    return(x@.dataSamples)
  }
  if( all(vapply(number,assertive::is_a_number,logical(1)))) {
    df <- x@.dataSamples
    if (nrow(df[df$ExperimentNo %in% number, ]) == 1) {
      ret <- as.list(df[df$ExperimentNo == number,])
      ret[!(names(ret) == "ExperimentNo")] <- 
        as.character(ret[!(names(ret) == "ExperimentNo")])
      return(ret)
    }
    if (nrow(df[df$ExperimentNo %in% number, ]) > 1) {
      ret <- df[df$ExperimentNo %in% number,]
      return(ret)
    }
  }
  stop("Something went wrong. Invalid experiment number!", call. = FALSE)
}


#' @rdname saveScores
#' @aliases loadScores saveScores
#'
#' @title Saves and loads scores
#'
#' @param x RNAmodR
#' @param grl GRangesList with scores per gene
#' @param geneNames optional: gene names for which scores should be loaded
#'
#' @return GRangesList
#' 
#' @importFrom rtracklayer export.gff3
#' 
#' @export
#'
#' @examples
#' \donttest{
#' saveScores(experiment, grl)
#' loadScores(experiment, "RNA28S5")
#' }
setMethod(
  f = "saveScores",
  signature = signature(x = "RNAmodR",
                        number = "numeric",
                        scores = "list"),
  definition = function(x,
                        number,
                        scores) {
    experiment <- getExperimentData(.Object, number)
    if( assertive::is_a_bool(experiment)) {
      if( assertive::is_false(experiment) ) {
        stop("Incorrect experiment identifier given: ",
             number,
             call. = FALSE)
      }
    }
    folder <- paste0(getOutputFolder(x),
                     "Scores/",
                     experiment$ExperimentName,
                     "/")
    if(!assertive::is_dir(folder)){
      dir.create(folder, recursive = TRUE)
    }
    fileName <- paste0(folder,
                       "/",
                       experiment$ExperimentName,
                       "_scores.gff3")
    
    
    rtracklayer::export.gff3(unlist(grl), con = fileName)
    for(i in seq_along(scores)){
      fileName <- paste0(folder,
                         "/",
                         experiment$ExperimentName,
                         "_scores_",
                         names(scores[i]),
                         ".gff3")
      rtracklayer::export.gff3(scores[[i]], con = fileName)
      # export as csv
      
    }
    return(invisible(scores))
  }
)

#' @rdname saveScores
#' @importFrom rtracklayer import.gff3
#' @export
setMethod(
  f = "loadScores",
  signature = signature(x = "RNAmodR",
                        number = "numeric"),
  definition = function(x,
                        number,
                        geneNames) {
    experiment <- getExperimentData(.Object, number)
    if( assertive::is_a_bool(experiment)) {
      if( assertive::is_false(experiment) ) {
        stop("Incorrect experiment identifier given: ",
             number,
             call. = FALSE)
      }
    }
    folder <- paste0(getOutputFolder(x),
                     "Scores/",
                     experiment$ExperimentName,
                     "/")
    fileName <- paste0(folder,
                       "/",
                       experiment$ExperimentName,
                       "_scores.gff3")
    gr <- rtracklayer::import.gff3(con = fileName)
    grl <- split(gr,mcols(gr$ID))
    if(is.null(geneNames)){
      return(grl)
    }
    grl <- grl[names(grl) %in% geneNames]
    if(length(grl) == 0){
      stop("No scores found for gene names given.",
           call. = FALSE)
    }
    return(grl)
  }
)


#' @rdname saveModifications
#' @aliases loadModifications saveModifications
#'
#' @title Saves and loads identified modifications
#'
#' @param x RNAmodR
#' @param grl GRangesList with modifications found per gene
#' @param geneNames optional: gene names for which modifications should be loaded
#'
#' @return GRangesList
#' 
#' @importFrom rtracklayer export.gff3
#' 
#' @export
#'
#' @examples
#' \donttest{
#' saveModifications(experiment, grl)
#' loadModifications(experiment, "RNA28S5")
#' }
setMethod(
  f = "saveModifications",
  signature = signature(x = "RNAmodR",
                        number = "numeric",
                        modifications = "list"),
  definition = function(x,
                        number,
                        modifications) {
    experiment <- getExperimentData(.Object, number)
    if( assertive::is_a_bool(experiment)) {
      if( assertive::is_false(experiment) ) {
        stop("Incorrect experiment identifier given: ",
             number,
             call. = FALSE)
      }
    }
    folder <- paste0(getOutputFolder(x),
                     "Modifications/",
                     experiment$ExperimentName,
                     "/")
    if(!assertive::is_dir(folder)){
      dir.create(folder, recursive = TRUE)
    }
    fileName <- paste0(folder,
                       "/",
                       experiment$ExperimentName,
                       "_modifications.gff3")
    rtracklayer::export.gff3(unlist(grl), con = fileName)
    for(i in seq_along(modifications)){
      fileName <- paste0(folder,
                         "/",
                         experiment$ExperimentName,
                         "_modifications_",
                         names(modifications[i]),
                         ".gff3")
      rtracklayer::export.gff3(modifications[[i]], con = fileName)
      # export as csv
      
    }
    return(invisible(grl))
  }
)

#' @rdname saveScores
#' @importFrom rtracklayer import.gff3
#' @export
setMethod(
  f = "loadModifications",
  signature = signature(x = "RNAmodR",
                        number = "numeric"),
  definition = function(x,
                        number,
                        geneNames) {
    experiment <- getExperimentData(.Object, number)
    if( assertive::is_a_bool(experiment)) {
      if( assertive::is_false(experiment) ) {
        stop("Incorrect experiment identifier given: ",
             number,
             call. = FALSE)
      }
    }
    folder <- paste0(getOutputFolder(x),
                     "Modifications/",
                     experiment$ExperimentName,
                     "/")
    fileName <- paste0(folder,
                       "/",
                       experiment$ExperimentName,
                       "_modifications.gff3")
    gr <- rtracklayer::import.gff3(con = fileName)
    grl <- split(gr,mcols(gr$ID))
    if(is.null(geneNames)){
      return(grl)
    }
    grl <- grl[names(grl) %in% geneNames]
    if(length(grl) == 0){
      stop("No modifications found for gene names given.",
           call. = FALSE)
    }
    return(grl)
  }
)












#' 
#' 
#' #' @rdname getGffResult
#' #'
#' #' @title getGffResult
#' #'
#' #' @description 
#' #' With this function the results of the parseForModification are
#' #' returned, if they are saved on disk. Any combination of modifications can be
#' #' requested, but only found results for the given modification are returned.
#' #' The results contain the local coordinates per transscript.
#' #'
#' #' @param .Object RNAmod.
#' #' @param number numeric.
#' #' @param modification character.
#' #' @param genomicCoordinates logical value, whether to return results with local
#' #' transcript coordinates or genomic coordinates. (default = FALSE)
#' #'
#' #' @return the found modifications as a GRanges object in the gff3 format.
#' #' @export
#' #' 
#' #' @importFrom GenomeInfoDb seqnames
#' #' @importFrom rtracklayer import.gff3
#' #' @importFrom IRanges ranges
#' #' @importFrom S4Vectors mcols
#' #' @importFrom BiocGenerics strand
#' #'
#' #' @examples
#' #' \donttest{
#' #' getGffResult(mod,1,"m7G")
#' #' }
#' setMethod(
#'   f = "getGffResult",
#'   signature = signature(.Object = "RNAmodR",
#'                         number = "numeric"),
#'   definition = function(.Object,
#'                         number,
#'                         genomicCoordinates) {
#'     # Input check
#'     assertive::assert_is_a_bool(genomicCoordinates)
#'     # get experiment data
#'     experiment <- getExperimentData(.Object, number)
#'     .check_for_experiment(experiment, number)
#'     # get file names
#'     fileNames <- .get_gff_filenames(.Object,
#'                                     unique(experiment$SampleName),
#'                                     unique(unlist(experiment$Modifications)))
#'     fileNames <- fileNames[vapply(fileNames, 
#'                                   assertive::is_existing_file,
#'                                   logical(1))]
#'     if(length(fileNames) == 0){
#'       stop("No result gff files found.",
#'            call. = FALSE)
#'     }
#'     gffs <- lapply(fileNames, rtracklayer::import.gff3)
#'     return(.merge_GRanges(gffs))
#'   }
#' )
#' 
#' .get_gff_filenames <- function(.Object,
#'                                sampleName,
#'                                modification){
#'   folder <- paste0(getOutputFolder(.Object),
#'                    "gff/")
#'   paste0(folder,
#'          "RNAmodR_",
#'          sampleName,
#'          "_",
#'          modification,
#'          ".gff3")
#' }
#' 
#' .merge_GRanges <- function(gffs){
#'   # If only one GRanges, return directly
#'   if(length(gffs) == 1) return(gffs[[1]])
#'   # get chromosome names
#'   chrom_names <- lapply(gffs, function(g){as.character(GenomeInfoDb::seqnames(g))})
#'   chrom_names <- unlist(chrom_names)
#'   # get ranges
#'   ranges <- do.call(c, lapply(gffs, IRanges::ranges))
#'   # get strand information
#'   strand <- unlist(lapply(gffs, function(g){as.character(BiocGenerics::strand(g))}))
#'   # get mcols and check identity
#'   mcols <- lapply(gffs, S4Vectors::mcols)
#'   if(!.ident(lapply(mcols, colnames)))
#'     stop("gff files are not compatible. Incompatible metadata columns.",
#'          call. = FALSE)
#'   mcols <- do.call(rbind, lapply(gffs, S4Vectors::mcols))
#'   # check overall compatibility
#'   lengths <- append(lapply(list(chrom_names,
#'                          ranges,
#'                          strand), length),
#'                     list(nrow(mcols)))
#'   if(!.ident(lengths))
#'     stop("gff files are not compatible. Incompatible data.",
#'          call. = FALSE)
#'   # create GRanges
#'   gff <- GenomicRanges::GRanges(S4Vectors::Rle(chrom_names),
#'                                 ranges = ranges,
#'                                 strand = strand,
#'                                 mcols)
#'   return(gff)
#' }
#' 
#' #' @title .saveGffResult
#' #' 
#' #' @description
#' #' Saves GRanges object as gff3 file file
#' #'
#' #' @param .Object an RNAmod object 
#' #' @param gff a GRanges object 
#' #' @param experiment a list containing data for one experiment
#' #' @param modification name of modification, one or more character 
#' #'
#' #' @return GRanges object
#' #' 
#' #' @importFrom rtracklayer export.gff3
#' setMethod(
#'   f = ".saveGffResult", 
#'   signature = signature(.Object = "RNAmodR", 
#'                         gff = "GRanges", 
#'                         experiment = "DataFrame",
#'                         modification = "character"),
#'   definition = function(.Object, 
#'                         gff, 
#'                         experiment,
#'                         modification) {
#'     fileNames <- .get_gff_filenames(.Object,
#'                                     unique(experiment$SampleName),
#'                                     modification)
#'     folder <- unique(dirname(fileNames))
#'     if(!assertive::is_dir(folder)){
#'       dir.create(folder, recursive = TRUE)
#'     }
#'     .save_single_gffs(.split_gff_for_each_modification(gff,modification),
#'                       fileNames)
#'     return(gff)
#'   }
#' )
#' 
#' .split_gff_for_each_modification <- function(gff,modification){
#'   # If only one GRanges, return directly
#'   if(length(modification) == 1) return(list(gff))
#'   #
#'   gffs <- lapply(modification, function(type){
#'     gff[S4Vectors::mcols(gff)$RNAmodR_type == type]
#'   })
#'   return(gffs)
#' }
#' 
#' 
#' .save_single_gffs <- function(gffs,fileNames){
#'   for(i in seq_along(gffs)){
#'     gff <- gffs[[i]]
#'     if(!is.null(gff) && length(gff) > 0){
#'       rtracklayer::export.gff3(gff, con = fileNames[[i]])
#'     }
#'   }
#' }
#' 
#' #' @rdname setGffResult
#' #'
#' #' @title setGffResult
#' #'
#' #' @param .Object RNAmod.
#' #' @param gff GRanges.
#' #' @param number numeric.
#' #' @param modification character.
#' #'
#' #' @return GRanges object
#' #' @export
#' #'
#' #' @examples
#' #' \donttest{
#' #' setGffResult(experiment, gff, 1, c("m7G"))
#' #' }
#' setMethod(
#'   f = "setGffResult",
#'   signature = signature(.Object = "RNAmodR",
#'                         gff = "GRanges",
#'                         number = "numeric"),
#'   definition = function(.Object,
#'                         gff,
#'                         number) {
#'     experiment <- getExperimentData(.Object, number)
#'     
#'     if( assertive::is_a_bool(experiment)) {
#'       if( assertive::is_false(experiment) ) {
#'         stop("Incorrect experiment identifier given: ",
#'              number,
#'              call. = FALSE)
#'       }
#'     }
#'     
#'     gff <- .saveGffResult(.Object, 
#'                           gff, 
#'                           experiment, 
#'                           unique(unlist(experiment$Modifications)))
#'     return(invisible(gff))
#'   }
#' )
#' 
#' 
#' 
