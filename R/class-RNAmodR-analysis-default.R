#' @include class-RNAmodR-mod-type.R
#' @include class-RNAmodR-analysis-type.R
NULL

RNAMODR_DEFAULT_COVERAGE_MIN <- 50
RNAMODR_DEFAULT_AVR_COVERAGE_MIN <- 10

#' @rdname RNAmodR-analysis-class
#'
#' @description 
#' \code{analysis_default}: the class can be used for analyzing the 5'-end of
#' reads to detect specific post-transcriptional modifications.
#' 
#' @export
#' @importFrom scales scientific
#'
#' @examples
#' \donttest{
#' ad <- new("analysis_default")
#' }
setClass("analysis_default",
         contains = "analysis",
         prototype = list(plotType = "default",
                          dataLabel = "mean(relative arrest rate)",
                          dataFormat = scales::percent))


#' @rdname convertReadsToPositions
#'
#' @description
#' \code{analysis_default}: load the input files, eg. bam files and 
#' pre-processes the information to be used in the modification detection.
#'
#' @return a modifified \code{analysis} class object containing now the read 
#' data required for modification detection.
#' @export
#' 
#' @importFrom IRanges findOverlaps extractList
#' @importFrom S4Vectors split from to
#'
#' @examples
#' \donttest{
#' ad <- new("analysis_default")
#' convertReadsToPositions(ad,
#'                        "test.bam",
#'                        "test.gff",
#'                        test_param)
#' }
setMethod(
  f = "convertReadsToPositions",
  signature = signature(object = "analysis_default",
                        files = "character",
                        conditions = "character",
                        gff = "GRanges",
                        param = "ScanBamParam"),
  definition = function(object,
                        files,
                        conditions,
                        gff,
                        param) {
    object@conditions <- conditions
    # detect modifications in each file
    data <- lapply(files,
                   FUN = .get_transcript_data,
                   gff,
                   param)
    names(data) <- conditions
    data <- data[!is.null(data)]
    if(length(data) == 0){
      stop("No reads detected in any bam file :\n",
           paste(files, collapse = "\n"),
           call. = FALSE)
    }
    # Process only genes found in all treated datasets
    IDs <- lapply(data,names)
    IDs <- IDs[names(IDs) == "Treated"]
    IDs <- Reduce(intersect, IDs)
    res <- lapply(IDs,
                  FUN = .get_data,
                  data = data,
                  conditions = names(data))
    names(res) <- IDs
    res <- res[!vapply(res,is.null,logical(1))]
    object@data <- res
    return(object)
  }
)

# returns the position data for analysis as a list of data per replicate for
# individual transcript
.get_data <- function(ID,
                      data,
                      conditions){
  res <- lapply(data,"[[",ID)
  res <- res[!vapply(res,is.null,logical(1))]
  return(res)
}

# detect modifications in each file
.get_transcript_data <- function(bamFile,
                                 gff,
                                 param){
  bamData <- GenomicAlignments::readGAlignments(bamFile,
                                                param = param)
  # Total counts
  totalCounts <- Rsamtools::idxstatsBam(bamFile,
                                        param = param)
  totalCounts <- sum(totalCounts$mapped)
  # process result and split into chunks based on gff
  IDs <- .get_IDs_from_scanBamParam(param)
  gff_subset <- gff[.get_unique_identifiers(gff) %in% IDs,]
  hits <- IRanges::findOverlaps(bamData,
                                gff_subset)
  bamData <- IRanges::extractList(bamData, 
                                  S4Vectors::split(
                                    S4Vectors::from(hits),
                                    as.factor(S4Vectors::to(hits))))
  names(bamData) <- .get_unique_identifiers(gff_subset)[
    as.numeric(names(bamData))]
  
  # bamData <- bamData[names(bamData) %in% c("RDN18-1") |
  #                      grepl("^t",names(bamData))]
  # bamData <- bamData[names(bamData) %in% c("tE(TTC)B",
  #                                          "tG(GCC)B",
  #                                          "tH(GTG)E1",
  #                                          "tS(CGA)C")]  
  # bamData <- bamData[names(bamData) %in% c("tS(CGA)C")]
  # bamData <- bamData[names(bamData) %in% c("tH(GTG)E1")]
  # bamData <- bamData[names(bamData) %in% c("RDN25-1",
  #                                          "RDN18-1")]
  # bamData <- bamData[names(bamData) %in% c("YMR116C")]
  # bamData <- bamData[names(bamData) %in% c("RDN18-1")]
  # bamData <- bamData[names(bamData) %in% c("YJL047C",
  #                                          "YHR199C-A")]
  # bamData <- bamData[names(bamData) %in% c("YBR041W",
  #                                          "YBR056W",
  #                                          "YAL030W")]
  # bamData <- bamData[names(bamData) %in% c("YBR056W")]
  
  if(length(bamData) == 0){
    warning("No reads detected in bam file '",
            bamFile,
            "'")
    return(NULL)
  }
  # param <- BiocParallel::bpparam()
  # bak_param <- param
  # param$workers <- 2
  # BiocParallel::register(param)
  # transcripts <- mapply(
  transcripts <- BiocParallel::bpmapply(
                                        FUN = .get_position_data_of_transcript,
                                        bamData,
                                        names(bamData),
                                        MoreArgs = list(totalCounts,
                                                        gff),
                                        SIMPLIFY = FALSE)
  # BiocParallel::register(bak_param)
  names(transcripts) <- names(bamData)
  # remove entries for transcript for which position data is sufficient
  transcripts <- transcripts[!vapply(transcripts, is.null, logical(1))]
  if(length(transcripts) == 0) return(NULL)
  return(transcripts)
}

.estimate_number_of_workers <- function(){
  
}


# For each transcript get positional data
# This can be individually done for different modification types
.get_position_data_of_transcript <- function(data,
                                             id,
                                             counts,
                                             gff){
  # debug
  if( getOption("RNAmodR_debug") ){
    message(id)
  }
  # skip if transcript does not have data
  if(length(data) == 0) return(NULL)
  # get ID and GRanges
  gr <- .subset_gff_for_unique_transcript(gff, 
                                          id)
  # get a list of introns and the position which are void
  posToBeRemoved <- .get_intron_positions(gff,
                                          gr$ID)
  # move position based on strand
  data <- data[.is_on_correct_strand(data,.get_unique_strand(gr))]
  # discard reads out of boundaries
  data <- data[BiocGenerics::end(data) <= BiocGenerics::end(gr),]
  data <- data[BiocGenerics::start(data) >= BiocGenerics::start(gr),]
  # get coverage per position
  coverage <- .convert_global_coverage_to_local_coverage(gff,
                                                         gr,
                                                         data,
                                                         posToBeRemoved)
  # If avergae coverage is to low
  if((sum(coverage)/length(coverage)) < RNAMODR_DEFAULT_AVR_COVERAGE_MIN) {
    return(NULL)
  }
  # do position conversion to translate genomic position to local transcript
  # position. take care of introns, etc
  stopsData <- .convert_global_stops_to_local_stops(gff,
                                                    gr,
                                                    data,
                                                    posToBeRemoved)
  if(length(stopsData) != length(coverage)) browser()
  # calculate relative amount of stops per coverage as percent
  posData <- .calc_relative_stop_data(stopsData,
                                      coverage)
  return(list(data = posData,
              coverage = coverage))
}

# if control and treated sample is available
.calc_relative_stop_data <- function(stopsData,
                                     coverage){
  # remove low coverage positions by setting pos data to zero
  toLowCoverage <- as.numeric(names(coverage[coverage < RNAMODR_DEFAULT_COVERAGE_MIN]))
  zeroCoverage <- as.numeric(names(coverage[coverage == 0]))
  # offset to not cover the end of coverage
  zeroCoverage <- zeroCoverage + 1
  zeroCoverage <- zeroCoverage[zeroCoverage <= max(as.numeric(names(coverage)))]
  # avoid deviding by zero
  coverage[toLowCoverage] <- 1
  # calc relative stop data
  res <- stopsData / coverage
  res[toLowCoverage] <- NA
  # this quenches unsepecific results
  res[c(1,zeroCoverage)] <- NA
  res
}

#' @rdname parseMod
#' 
#' @description 
#' \code{analysis_default}: parses the pre-processed data in the data slot
#' to detect modifications
#' 
#' @return a modified \code{analysis} class object containting pre-processed
#' data and found modifications
#' 
#' @export
#' 
#' @importFrom stringr str_locate_all
#'
#' @examples
#' \donttest{
#' ad <- new("analysis_default")
#' convertReadsToPositions(ad,
#'                        "test.bam",
#'                        "test.gff",
#'                        test_param)
#' parseMod(ad,
#'          "test.gff",
#'          "test.fasta",
#'          modclasses)
#' }
setMethod(
  f = "parseMod",
  signature = signature(object = "analysis_default",
                        gff = "GRanges",
                        fafile = "FaFile",
                        modClasses = "list"),
  definition = function(object,
                        gff,
                        fafile,
                        modClasses) {
    # detect modification per transcript
    # res <- mapply(
    res <- BiocParallel::bpmapply(
                                  FUN = .analyze_transcript_prep,
                                  names(object@data),
                                  object@data,
                                  MoreArgs = list(gff = gff,
                                                  fafile = fafile,
                                                  modClasses = modClasses),
                                  SIMPLIFY = FALSE)
    names(res) <- names(object@data)
    res <- res[!vapply(res, is.null, logical(1))]
    # If not results are present return NA instead of NULL
    if(is.null(res)){
      res <- NA
    }
    object@modifications <- res
    return(object)
  }
)

# detect and merge modification positions in one transcript
.analyze_transcript_prep <- function(ID,
                                     data,
                                     gff,
                                     fafile,
                                     modClasses){
  # debug
  if( getOption("RNAmodR_debug") ){
    message(ID)
  }
  # get sequence of transcript and subset gff for single transcript data
  gr <- .subset_gff_for_unique_transcript(gff, ID)
  # get the genomic sequences
  seq <- .get_seq_for_unique_transcript(gr,fafile,ID)
  # get transcript sequence by removing intron sequences
  seq <- .get_transcript_sequence(gff,gr$ID,seq)
  # generate a location vector
  locations <- 1:length(seq)
  names(locations) <- seq
  # analyze the transcript
  # Retrieve modifications positions
  modifications <- lapply(locations,
                          .check_for_modification,
                          modClasses,
                          data,
                          locations)
  # name the locations based on sequence position
  names(modifications) <- paste0(ID,
                                 "_",
                                 names(locations),
                                 "_",
                                 locations)
  modifications <- modifications[!vapply(modifications,is.null,logical(1))]
  if(is.null(modifications) || length(modifications) == 0) return(NULL)
  modifications <- modifications[order(as.numeric(unlist(lapply(modifications, 
                                                                "[[", 
                                                                "location"))))]
  return(modifications)
}

# check for modifications at given position
.check_for_modification <- function(location, 
                                    modClasses,
                                    data,
                                    locations){
  res <- lapply(modClasses, function(class){
    return(checkForModification(class,
                                location = location,
                                locations = locations,
                                data = data))
  })
  res <- res[!vapply(res, is.null, logical(1))]
  if(length(res) != 1) return(NULL)
  # Return data
  return(list(location = location,
              type = res[[1]]$type,
              signal = res[[1]]$signal,
              signal.sd = res[[1]]$signal.sd,
              z = res[[1]]$z,
              nbsamples = res[[1]]$nbsamples))
}

#' @rdname mergePositionsOfReplicates
#'
#' @description
#' \code{analysis_default}
#'
#' @return a modified \code{analysis} class object now containing the aggregated
#' pre-processed data of the input. This data is used for saving.
#' 
#' @importFrom matrixStats rowSds
#' @importFrom stats complete.cases
#' 
#' @export
#'
#' @examples
#' \donttest{
#' load(system.file("data", file = "test_example.RData", package = "RNAmodR"))
#' mergePositionsOfReplicates(ad)
#' }
setMethod(
  f = "mergePositionsOfReplicates",
  signature = signature(object = "analysis_default"),
  definition = function(object) {
    res <- mapply(
    # res <- BiocParallel::bpmapply(
                                  FUN = .merge_positions,
                                  names(object@data),
                                  object@data)
    names(res) <- names(object@data)
    res <- res[!is.null(res)]
    # If not results are present return NA instead of NULL
    if(is.null(res)){
      res <- NA
    }
    object@data <- res
    return(object)
  }
)

# merge positions in one transcript
.merge_positions <- function(ID,data){
  # debug
  if( getOption("RNAmodR_debug") ){
    message(ID)
  }
  # get mean of coverage
  coverage <- lapply(data,"[[","coverage")
  coverage <- split(coverage,names(coverage))
  coverage <- .merge_coverage_data(coverage)
  # get mean of position stop data
  data <- lapply(data,"[[","data")
  positions <- as.numeric(unique(unlist(lapply(data,names))))
  data <- split(data,names(data))
  data <- .merge_position_data(data)
  df <- data.frame(pos = positions,
                   mean = data$mean,
                   sd = data$sd,
                   coverage = coverage$mean,
                   stringsAsFactors = FALSE)
  return(df)
}

# iterates on every position and calculates the difference of the means and sd
.merge_position_data <- function(data){
  treated <- data$Treated
  control <- data$Control
  if(is.null(control)){
    df <- S4Vectors::DataFrame(treated)
    colnames(df) <- names(treated)
    # remove incomplete 
    df[!stats::complete.cases(as.data.frame(df)),] <- 0
    mean <- rowMeans(as.matrix(df))
    sd <- matrixStats::rowSds(as.matrix(df))
  } else {
    df <- S4Vectors::DataFrame(append(treated,control))
    colnames(df) <- c(names(treated),names(control))
    treated <- df[,colnames(df) == "Treated",drop = FALSE]
    treated[!stats::complete.cases(as.data.frame(treated)),] <- 0
    control <- df[,colnames(df) == "Control",drop = FALSE]
    control[!stats::complete.cases(as.data.frame(control)),] <- 0
    mean <- rowMeans(as.matrix(treated)) - 
      rowMeans(as.matrix(control))
    mean[mean < 0] <- 0
    sd <- matrixStats::rowSds(as.matrix(treated)) + 
      matrixStats::rowSds(as.matrix(control))
  }
  return(list(mean = mean,
              sd = sd))
}

# iterates on every position and calculates a mean coverage
.merge_coverage_data <- function(coverage){
  treated <- coverage$Treated
  df <- data.frame(treated)
  colnames(df) <- names(treated)
  df[!complete.cases(df),] <- rep(0, ncol(df))
  mean <- rowMeans(df)
  return(list(mean = mean))
}