#' @include AllGenerics.R
#' @include AllGenerics-export.R
NULL

#' @name RNAmodR-class
#' 
#' @title RNAmodR Class
#' 
#' @description 
#' \code{RNAmodR}
NULL

#' @title RNAmodR-class
#' 
#' @rdname RNAmodR-class 
#'
#' @description 
#' title
#'
#' @return a RNAmodR object
#' @export
#' 
#' @import Rsamtools
#' @importClassesFrom GenomicRanges GRanges
#'
#' @examples
#' \donttest{
#' unzip(system.file("extdata", package = "RNAmodR", file = "RNAmodR.zip" ))
#' mod <- RNAmodR("test",
#'                "test_layout.csv",
#'                "test_gff.gff3",
#'                "test_masked.fasta")
#' }
setClass("RNAmodR",
         slots = c(
           experimentName = "character",
           inputFile = "character",
           inputFasta = "character",
           .dataSamples = "DataFrame",
           .dataFasta = "FaFile",
           .txdb = "TxDb",
           .wd = "character",
           .inputFolder = "character",
           .outputFolder = "character"
         ),
         prototype=list(
           .dataSamples = S4Vectors::DataFrame(),
           .wd = "",
           .inputFolder = "data/",
           .outputFolder = "results/"
         )
)


#' @rdname RNAmodR-class
#'
#' @importFrom utils read.csv
#'
#' @export
RNAmodR <- function(experimentName, 
                    experimentLayout, 
                    experimentTxDb, 
                    experimentFasta){
  
  # Input check
  if(missing(experimentName)) stop("Experimental name not provided",
                                   call. = FALSE)
  if(missing(experimentLayout)) stop("Experimental layout not provided",
                                     call. = FALSE)
  if(missing(experimentTxDb)) stop("experimentTxDb information not provided",
                                   call. = FALSE)
  if(missing(experimentFasta)) 
    stop("experimentFasta information not provided",
         call. = FALSE)
  
  if(!is.character(experimentName)) 
    stop("experimentName can only be assigned a character vector of length 1",
         call. = FALSE)
  if(!is.character(experimentLayout)) 
    stop("experimentLayout can only be assigned a character vector of length 1",
         call. = FALSE)
  if(!is_TxDb(experimentTxDb) & !is.character(experimentTxDb))
    stop("experimentTxDb can only be assigned a TxDb object or a character ",
         "vector of length 1",
         call. = FALSE)
  if(!is.character(experimentFasta)) 
    stop("experimentFasta can only be assigned a character vector of length 1",
         call. = FALSE)
  
  class <- new("RNAmodR", 
               experimentName, 
               experimentLayout, 
               experimentTxDb, 
               experimentFasta)
  return(class)
}

setMethod(
  f = "initialize", 
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object, 
                        experimentName, 
                        experimentLayout, 
                        experimentTxDb, 
                        experimentFasta) {
    requireNamespace("rtracklayer", quietly = TRUE)
    # Save data to slots
    .Object@experimentName <- experimentName
    .Object@.wd <- paste0("./", .Object@experimentName, "/")
    
    .Object@inputFile <- paste0(.Object@.wd, 
                                .Object@.inputFolder, 
                                experimentLayout)
    .Object@inputFasta <- paste0(.Object@.wd, 
                                 .Object@.inputFolder, 
                                 experimentFasta)
    # check for existing and readable files
    assertive::assert_all_are_existing_files(c(.Object@inputFile,
                                               .Object@inputFasta))
    
    # Load fasta and gff file
    .Object@.dataFasta <- Rsamtools::FaFile(.Object@inputFasta)
    Rsamtools::indexFa(.Object@inputFasta)
    # save, load, create txdb object
    .Object@.txdb <- .get_txdb_object(.Object,
                                      experimentTxDb)
    if( !all( names(Rsamtools::scanFa(.Object@.dataFasta)) %in% seqlevels(.Object@.txdb) ) ) {
      stop("Names of sequences in the fasta file do not match the chromosome names in the gff file")
    }
    #load sample data
    samples <- .check_sample_data(utils::read.csv(.Object@inputFile, 
                                                  header = TRUE, 
                                                  stringsAsFactors = FALSE, 
                                                  sep = ";",
                                                  fileEncoding = "UTF-8-BOM"),
                                   getInputFolder(.Object))
    .Object@.dataSamples <- samples
    return(.Object)
  }
)

# save, load, create txdb object
.get_txdb_object <- function(.Object,
                             input){
  if(is_TxDb(input)) return(input)
  fileName <- paste0(.Object@.wd, 
                     .Object@.inputFolder,
                     input)
  if(assertive::is_existing_file(fileName)){
    gff <- rtracklayer::import.gff3(fileName)
    # Fix bullshit annotations
    types <- as.character(gff$type)
    if(any(types == "rRNA_gene")) {
      types[types == "rRNA_gene"] <- "rRNA"
    }
    if(any(types == "tRNA_gene")) {
      types[types == "tRNA_gene"] <- "tRNA"
    }
    gff$type <- types
    rtracklayer::export.gff3(object = gff, 
                             con = fileName)
    return(GenomicFeatures::makeTxDbFromGFF(fileName))
  }
  stop("Could not create TxDb object. No TxDb object given and now file found.",
       call. = FALSE)
}

# Checks experimental data for validity
.check_sample_data <- function(samples,
                               inputFolder){
  samples <- S4Vectors::DataFrame(samples)
  # Names of the first seven cols will be overwritten and are expected to fit
  # to the data type
  cols <- colnames(samples)
  if( length(cols) < length(RNAMODR_DEFAULT_COLNAMES) |
      !all(cols[1:length(RNAMODR_DEFAULT_COLNAMES)] == RNAMODR_DEFAULT_COLNAMES) ) {
    stop("Please check the format of the csv file. It should be ; delimited ",
         "and contain at least the following columns: ",
         paste(RNAMODR_DEFAULT_COLNAMES, collapse = ","),
         call. = FALSE)
  }
  cols[1:length(RNAMODR_DEFAULT_COLNAMES)] <- RNAMODR_DEFAULT_COLNAMES
  colnames(samples) <- cols
  # Check if ExperimentNo and ExperimentName match
  tmpExperimentNo <- as.numeric(samples$ExperimentNo)
  tmpSampleName <- paste0(as.character(samples$SampleName),
                          "_",
                          as.character(samples$Replicate))
  if( !assertive::are_same_length(tmpExperimentNo, tmpSampleName)){
    stop("ExperimentNo and SampleName/Replicate are ambigeous. Check layout ",
         "file for errors in these columns.",
         call. = FALSE)
  }
  # Check inputs
  assertive::assert_all_are_positive(samples$ExperimentNo)
  assertive::assert_all_are_whole_numbers(samples$ExperimentNo)
  assertive::assert_all_are_non_missing_nor_empty_character(samples$SampleName)
  assertive::assert_all_are_non_missing_nor_empty_character(samples$ShortName)
  assertive::assert_all_are_positive(samples$Replicate)
  assertive::assert_all_are_whole_numbers(samples$Replicate)
  assertive::assert_all_are_existing_files(paste0(inputFolder,samples$BamFile))
  assertive::assert_all_are_whole_numbers(samples$MapQuality)
  # check modifications
  samples$Modifications <- lapply(samples$Modifications,
                                  function(str){str_split(str, ",")})
  modifications <- unique(unlist(samples$Modifications))
  modClasses <- .load_mod_classes(modifications)
  # check for valid condition types
  types <- unique(samples$Conditions)
  if(!all(types %in% RNAMODR_SAMPLE_TYPES)){
    stop("Not all condition types are valid. Valid types are:\n",
         paste(RNAMODR_SAMPLE_TYPES, collapse = ","),
         call. = FALSE)
  }
  return(samples)
}

# default class methods --------------------------------------------------------

#' @rdname RNAmodR-class
#' 
#' @param object a RNAmod object.
setMethod(
  f = "show",
  signature = signature(object = "RNAmodR"),
  definition = function(object) {
    print("test")
  }
)


# getters and setters ----------------------------------------------------------

#' @rdname RNAmodR-Accessors
#'
#' @export
setMethod(
  f = "getInputFolder",
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object){
    return(paste0(.Object@.wd, .Object@.inputFolder))
  }
)
#' @rdname RNAmodR-Accessors
#'
#' @export
setMethod(
  f = "getOutputFolder",
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object){
    return(paste0(.Object@.wd, .Object@.outputFolder))
  }
)
#' @rdname RNAmodR-Accessors
#'
#' @export
setMethod(
  f = "getExperimentName",
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object){
    return(.Object@experimentName)
  }
)
#' @rdname RNAmodR-Accessors
#'
#' @export
setMethod(
  f = "getTxDb",
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object){
    return(.Object@.txdb)
  }
)
#' @rdname RNAmodR-Accessors
#'
#' @export
setMethod(
  f = "getFastaFile",
  signature = signature(.Object = "RNAmodR"),
  definition = function(.Object){
    return(.Object@.dataFasta)
  }
)