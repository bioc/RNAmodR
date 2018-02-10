#' @include RNAmodR.R
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


#' @rdname getSummarizedExperiment
#' 
#' @title returns one or more SummarizedExperiment
#' 
#' @description
#' Global access to SummarizedExperiments stored by RpfExperiment. 
#' \code{getSummarizedExperiment()} returns the result of experiment, whereas 
#' \code{getSummarizedExperiments()}, is the vectorized version returning a
#' list of experiment results.
#' 
#' @param .Object an RNAmod object 
#' @param number a number defining the experiment. For 
#' getSummarizedExperiments() more than number can be defined as numeric vector.
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
  signature = signature(.Object = "RNAmodR", 
                        number = "numeric"),
  definition = function(.Object, 
                        number) {
    if( !assertive::is_scalar(number)){
      stop(paste0("Multiple numbers provided. only one excepted.",
                  call. = FALSE))
    }
    experiment <- getExperimentData(.Object, number)
    .check_for_experiment(experiment, number)
    return(.loadSummarizedExperiment(.Object, 
                                     experiment, 
                                     unique(unlist(experiment$Modifications)), 
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
  signature = signature(.Object = "RNAmodR",
                        number = "numeric"),
  definition = function(.Object, 
                        number) {
    if( length(number) == 0 ) stop("no experiment number given")
    assertive::assert_all_are_whole_numbers(number)
    
    ses <- lapply(number, function(x){
        getSummarizedExperiment(.Object, 
                                x)
      })
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
#' 
#' @importFrom S4Vectors SimpleList
setMethod(
  f = ".loadSummarizedExperiment", 
  signature = signature(.Object = "RNAmodR", 
                        experiment = "DataFrame",
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
                                 "RNAmodR_",
                                 unique(experiment$SampleName),
                                 "_")
    fileNames <- paste0(folder,
                        modification,
                        ".RData")
    fileNames <- fileNames[vapply(fileNames, 
                                  assertive::is_existing_file,
                                  logical(1))]
    if(length(fileNames) == 0){
      stop("No result SummarizedExperiment files found.",
           call. = FALSE)
    }
    ses <- lapply(fileNames,function(file){
      load(file)
      return(se)
    })
    return(.merge_se_for_modifications(ses,modification))
  }
)

# check for equality of all elements of a list
.ident <- function(l){
  all(vapply(l, function(x){identical(x,l[[length(l)]])}, logical(1)))
}

.merge_se_for_modifications <- function(ses,modification){
  # If only one SummarizedExperiment, return directly
  if(length(ses) == 1) return(ses[[1]])
  #
  rowColNames <- lapply(ses, function(se){
    colnames(SummarizedExperiment::rowData(se))
    })
  if(!.ident(rowColNames))
    stop("SummarizedExperiments are not compatible. ",
         "Incompatible rowData columns.",
         call. = FALSE)
  if(!.ident(lapply(ses, rownames)))
    stop("SummarizedExperiments are not compatible. ",
         "Incompatible rownames.",
         call. = FALSE)
  colColNames <- lapply(ses, function(se){
    colnames(SummarizedExperiment::colData(se))
    })
  if(!.ident(colColNames))
    stop("SummarizedExperiments are not compatible. ",
         "Incompatible colData columns.",
         call. = FALSE)
  colRowNames <- lapply(ses, function(se){
    rownames(SummarizedExperiment::colData(se))
    })
  if(!.ident(colRowNames))
    stop("SummarizedExperiments are not compatible. ",
         "Incompatible colData rownames.",
         call. = FALSE)
  # get combined assays data
  assayData <- lapply(ses, SummarizedExperiment::assays)
  # get combined mod data
  mods <- lapply(ses, function(se){
    SummarizedExperiment::rowData(se)$mods
  })
  mods <- lapply(seq_along(mods[[1]]), function(i){
    m <- lapply(seq_along(mods), function(j){
      mods[[j]][[i]]
    })
    do.call(rbind, m)
  })
  
  # Construct joind SummarizedExperiment
  se <- ses[[1]]
  SummarizedExperiment::assays(se) <- do.call(c, assayData)
  SummarizedExperiment::rowData(se)$mods <- mods
  return(se)
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
  signature = signature(.Object = "RNAmodR", 
                        se = "RangedSummarizedExperiment", 
                        experiment = "DataFrame",
                        modification = "character"),
  definition = function(.Object, 
                        se, 
                        experiment,
                        modification) {
    folder <- paste0(getOutputFolder(.Object),
                     "SE/")
    if(!assertive::is_dir(folder)){
      dir.create(folder, recursive = TRUE)
    }
    modification <- intersect(modification,names(assays(se)))
    if( length(modification) == 0 ){
      stop("No modification data found, while saving SummarizedExperiment.",
           call. = FALSE)
    }
    fileNames <- paste0(folder,
                        "RNAmodR_",
                        unique(experiment$SampleName),
                        "_",
                        modification,
                        ".RData")
    .save_single_ses(.split_se_for_each_modification(se,modification), 
                     fileNames)
    return(se)
  }
)

.split_se_for_each_modification <- function(se,modification){
  # If only one SummarizedExperiment, return directly
  if(length(modification) == 1) return(list(se))
  #
  ses <- BiocParallel::bplapply(modification, function(type, se){
    tmp <- se
    assays(tmp) <- assays(se)[type]
    rowData(tmp)$mods <- lapply(rowData(tmp)$mods, function(df){
      df[df$RNAmodR_type == type,]
    })
    return(tmp)
  }, se)
  return(ses)
}

.save_single_ses <- function(ses,fileNames){
  for(i in seq_along(ses)){
    se <- ses[[i]]
    save(se, file = fileNames[[i]])
  }
}


#' @rdname setSummarizedExperiment
#' 
#' @title sets a SummarizedExperiment object
#' 
#' @description
#' Saves/overwrites a SummarizedExperiment object for certain experiment and 
#' passes the SummarizedExperiment object on to be saved to disk as .RData file
#' in the \code{results\\SE} folder
#' 
#' @param .Object a RNAmod object 
#' @param se a RangedSummarizedExperiment object
#' @param number a number defining the experiment
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
  signature = signature(.Object = "RNAmodR", 
                        se = "RangedSummarizedExperiment", 
                        number = "numeric"),
  definition = function(.Object, 
                        se, 
                        number) {
    experiment <- getExperimentData(.Object, number)
    .check_for_experiment(experiment, number)
    
    se <- .saveSummarizedExperiments(.Object, 
                                     se, 
                                     experiment,
                                     unique(unlist(experiment$Modifications)))
    return(invisible(se))
  }
)

.check_for_experiment <- function(experiment, 
                                  number){
  if( assertive::is_a_bool(experiment)) {
    if( assertive::is_false(experiment) ) {
      stop("Incorrect experiment identifier given: ",
           number,
           call. = FALSE)
    }
  }
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


# easy access functions to results ---------------------------------------------

#' @rdname seGetModifications
#' @aliases seGetModifications seGetPositions
#' @title Getting Position and Modification data from SummarizedExperiment 
#' objects
#' 
#' @description 
#' \code{seGetModifications} and \code{seGetPositions} returns modification and
#' position data, respectively, from a SummarizedExperiment object. The function
#' expects the data to be stored as provided and defined by the RNAmodR package 
#' as rowData and assay data.
#'
#' @param se a SummarizedExperiment object
#' @param genes gene names used for subsetting the results 
#' @param modifications \code{seGetModifications} modification types used for 
#' subsetting the results 
#'
#' @return
#' DataFrame containing the information on saved modifications
#' @export
#'
#' @examples
#' \donttest{
#' seGetModifications(se)
#' seGetModifications(se,"RDN18-1")
#' seGetModifications(se,"RDN18-1","m7G")
#' }
setMethod(
  f = "seGetModifications",
  signature = signature(se = "SummarizedExperiment"),
  definition = function(se,
                        genes,
                        modifications){
    # check validity of SE
    .check_row_data(se, "mods")
    # get data
    mods <- SummarizedExperiment::rowData(se)$mods
    names(mods) <- rownames(se)
    # subset if defined
    if(!missing(genes)){
      assertive::assert_all_are_non_empty_character(genes)
      mods <- mods[names(mods) %in% genes]
    }
    if(!missing(modifications)){
      assertive::assert_all_are_non_empty_character(modifications)
      modAvail <- unique(unlist(lapply(mods, function(x){
        unique(as.character(x$RNAmodR_type))})))
      if(length(modAvail[modifications %in% modAvail]) == 0){
        stop("None of the modifications '",
             paste(modifications,collapse = "','"),
             "' are available in results. Following modifications are ",
             "present: '",
             paste(modAvail,collapse = "','"),"'",
             call. = FALSE)
      }
      modAvail <- modAvail[modAvail %in% modifications]
      res <- lapply(mods, function(x){
        x[x$RNAmodR_type %in% modAvail,]
      })
      names(res) <- names(mods)
      mods <- res
    }
    mods <- .clean_output(mods, function(x){
      if(nrow(x) == 0) return(TRUE)
      FALSE
    })
    return(mods)
  }
)

.check_row_data <- function(se, name){
  colnames <- colnames(SummarizedExperiment::rowData(se))
  if(!(name %in% colnames)){
    stop("Invalid SummarizedExperiment object. Column ", name, " not found in ",
         "rowData.",
         call. = FALSE)
  }
  return(invisible(TRUE))
}

#' @rdname seGetModifications
#'
#' @return
#' a list of DataFrames containing the position data
#' @export
#'
#' @examples
#' \donttest{
#' seGetPositions(se)
#' seGetPositions(se,"RDN18-1")
#' }
setMethod(
  f = "seGetPositions",
  signature = signature(se = "SummarizedExperiment"),
  definition = function(se,
                        genes){
    # check validity of SE
    .check_row_data(se, "positions")
    # get data
    positions <- SummarizedExperiment::rowData(se)$positions
    # subset if defined
    if(!missing(genes)){
      assertive::assert_all_are_non_empty_character(genes)
      positions <- positions[names(positions) %in% genes]
    }
    positions <- .clean_output2(positions, is.null)
    return(positions)
  }
)
