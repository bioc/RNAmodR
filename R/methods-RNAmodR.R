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
#' @param .Object a RNAmod object
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
  signature = signature(.Object = "RNAmodR", 
                        number = "ANY"),
  definition = function(.Object, 
                        number){
    return(.get_experiment_data(.Object,
                                number))
  }
)

.get_experiment_data <- function(.Object,
                                 number){
  if( missing(number) || assertive::is_identical_to_false(number) ) {
    return(.Object@.dataSamples)
  }
  if( all(vapply(number,assertive::is_a_number,logical(1)))) {
    df <- .Object@.dataSamples
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


#' @rdname getGffResult
#'
#' @title getGffResult
#'
#' @description 
#' With this function the results of the parseForModification are
#' returned, if they are saved on disk. Any combination of modifications can be
#' requested, but only found results for the given modification are returned.
#' The results contain the local coordinates per transscript.
#'
#' @param .Object RNAmod.
#' @param number numeric.
#' @param modification character.
#' @param genomicCoordinates logical value, whether to return results with local
#' transcript coordinates or genomic coordinates. (default = FALSE)
#'
#' @return the found modifications as a GRanges object in the gff3 format.
#' @export
#' 
#' @importFrom GenomeInfoDb seqnames
#' @importFrom rtracklayer import.gff3
#' @importFrom IRanges ranges
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics strand
#'
#' @examples
#' \donttest{
#' getGffResult(mod,1,"m7G")
#' }
setMethod(
  f = "getGffResult",
  signature = signature(.Object = "RNAmodR",
                        number = "numeric"),
  definition = function(.Object,
                        number,
                        genomicCoordinates) {
    # Input check
    assertive::assert_is_a_bool(genomicCoordinates)
    # get experiment data
    experiment <- getExperimentData(.Object, number)
    .check_for_experiment(experiment, number)
    # get file names
    fileNames <- .get_gff_filenames(.Object,
                                    unique(experiment$SampleName),
                                    unique(unlist(experiment$Modifications)))
    fileNames <- fileNames[vapply(fileNames, 
                                  assertive::is_existing_file,
                                  logical(1))]
    if(length(fileNames) == 0){
      stop("No result gff files found.",
           call. = FALSE)
    }
    gffs <- lapply(fileNames, rtracklayer::import.gff3)
    return(.merge_GRanges(gffs))
  }
)

.get_gff_filenames <- function(.Object,
                               sampleName,
                               modification){
  folder <- paste0(getOutputFolder(.Object),
                   "gff/")
  paste0(folder,
         "RNAmodR_",
         sampleName,
         "_",
         modification,
         ".gff3")
}

.merge_GRanges <- function(gffs){
  # If only one GRanges, return directly
  if(length(gffs) == 1) return(gffs[[1]])
  # get chromosome names
  chrom_names <- lapply(gffs, function(g){as.character(GenomeInfoDb::seqnames(g))})
  chrom_names <- unlist(chrom_names)
  # get ranges
  ranges <- do.call(c, lapply(gffs, IRanges::ranges))
  # get strand information
  strand <- unlist(lapply(gffs, function(g){as.character(BiocGenerics::strand(g))}))
  # get mcols and check identity
  mcols <- lapply(gffs, S4Vectors::mcols)
  if(!.ident(lapply(mcols, colnames)))
    stop("gff files are not compatible. Incompatible metadata columns.",
         call. = FALSE)
  mcols <- do.call(rbind, lapply(gffs, S4Vectors::mcols))
  # check overall compatibility
  lengths <- append(lapply(list(chrom_names,
                         ranges,
                         strand), length),
                    list(nrow(mcols)))
  if(!.ident(lengths))
    stop("gff files are not compatible. Incompatible data.",
         call. = FALSE)
  # create GRanges
  gff <- GenomicRanges::GRanges(S4Vectors::Rle(chrom_names),
                                ranges = ranges,
                                strand = strand,
                                mcols)
  return(gff)
}

#' @title .saveGffResult
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
#' 
#' @importFrom rtracklayer export.gff3
setMethod(
  f = ".saveGffResult", 
  signature = signature(.Object = "RNAmodR", 
                        gff = "GRanges", 
                        experiment = "DataFrame",
                        modification = "character"),
  definition = function(.Object, 
                        gff, 
                        experiment,
                        modification) {
    fileNames <- .get_gff_filenames(.Object,
                                    unique(experiment$SampleName),
                                    modification)
    folder <- unique(dirname(fileNames))
    if(!assertive::is_dir(folder)){
      dir.create(folder, recursive = TRUE)
    }
    .save_single_gffs(.split_gff_for_each_modification(gff,modification),
                      fileNames)
    return(gff)
  }
)

.split_gff_for_each_modification <- function(gff,modification){
  # If only one GRanges, return directly
  if(length(modification) == 1) return(list(gff))
  #
  gffs <- lapply(modification, function(type){
    gff[S4Vectors::mcols(gff)$RNAmodR_type == type]
  })
  return(gffs)
}


.save_single_gffs <- function(gffs,fileNames){
  for(i in seq_along(gffs)){
    gff <- gffs[[i]]
    if(!is.null(gff) && length(gff) > 0){
      rtracklayer::export.gff3(gff, con = fileNames[[i]])
    }
  }
}

#' @rdname setGffResult
#'
#' @title setGffResult
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
#' setGffResult(experiment, gff, 1, c("m7G"))
#' }
setMethod(
  f = "setGffResult",
  signature = signature(.Object = "RNAmodR",
                        gff = "GRanges",
                        number = "numeric"),
  definition = function(.Object,
                        gff,
                        number) {
    experiment <- getExperimentData(.Object, number)
    
    if( assertive::is_a_bool(experiment)) {
      if( assertive::is_false(experiment) ) {
        stop("Incorrect experiment identifier given: ",
             number,
             call. = FALSE)
      }
    }
    
    gff <- .saveGffResult(.Object, 
                          gff, 
                          experiment, 
                          unique(unlist(experiment$Modifications)))
    return(invisible(gff))
  }
)
