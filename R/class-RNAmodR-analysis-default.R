#' @include class-RNAmodR-mod-type.R
#' @include class-RNAmodR-analysis-type.R
NULL

RNAMODR_DEFAULT_FPK_THRESHOLD <- 100
RNAMODR_DEFAULT_COVERAGE_MIN <- 10

#' @name analysis_default
#' 
#' @title the default analysis class for RNAmodR
#'
#' @description 
#' \code{analysis_default}: the class can be used for analyzing the 5'-end of
#' reads to detect specific post-transcriptional modifications.
#' 
#' @export
#'
#' @examples
#' ad <- new("analysis_default")
setClass("analysis_default",
         contains = "analysis")


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
                        gff = "GRanges",
                        param = "ScanBamParam"),
  definition = function(object,
                        files,
                        gff,
                        param) {
    message(Sys.time())
    # detect modifications in each file
    data <- lapply(files,
                   FUN = get_transcript_data,
                   gff,
                   param)
    data <- data[!is.null(data)]
    if(length(data) == 0){
      stop("No reads detected in any bam file :\n",
           paste(files, collapse = "\n"),
           call. = FALSE)
    }
    object@data <- data
    return(object)
  }
)


# detect modifications in each file
get_transcript_data <- function(bamFile,
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
  
  bamData <- bamData[names(bamData) %in% c("RDN18-1") |
                       grepl("^t",names(bamData))]
  # bamData <- bamData[names(bamData) %in% c("RDN18-1",
  #                                          "tE(TTC)B")]
  
  # bamData <- bamData[names(bamData) %in% c("YMR116C")]
  # bamData <- bamData[names(bamData) %in% c("RDN18-1",
  #                                          "YAL030W")]
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
  transcripts <- mapply(
    # transcripts <- BiocParallel::bpmapply(
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
  # do position conversion to translate genomic position to local transcript
  # position. take care of introns, etc
  stopsData <- .convert_global_stops_to_local_stops(gff,
                                                    gr,
                                                    data,
                                                    posToBeRemoved)
  # get coverage per position
  coverage <- .convert_global_coverage_to_local_coverage(gff,
                                                         gr,
                                                         data,
                                                         posToBeRemoved)
  if(length(stopsData) != length(coverage)) browser()
  # calculate relative amount of stops per coverage as percent
  posData <- .calc_relative_stop_data(stopsData,
                                      coverage)
  return(list(data = posData))
}

.calc_relative_stop_data <- function(stopsData,
                                     coverage){
  # remove low coverage positions by setting pos data to zero
  toLowCoverage <- names(coverage[coverage < RNAMODR_DEFAULT_COVERAGE_MIN])
  # add first position since has by definition a value of 1
  toLowCoverage <- c(toLowCoverage,1)
  # avoid deviding by zero
  coverage[toLowCoverage] <- 1
  # calc relative stop data
  res <- stopsData / coverage
  res[toLowCoverage] <- 0
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
    message(Sys.time())
    # browser()
    # Process only genes found in all datasets
    IDs <- lapply(object@data,names)
    IDs <- Reduce(intersect, IDs)
    # detect modification per transcript
    res <- lapply(IDs,
    # res <- BiocParallel::bplapply(IDs,
                                  FUN = .analyze_transcript_prep,
                                  data = object@data,
                                  gff = gff,
                                  fafile = fafile,
                                  modClasses = modClasses)
    names(res) <- IDs
    res <- res[!vapply(res, is.null, logical(1))]
    # If not results are present return NA instead of NULL
    if(is.null(res)){
      res <- NA
    }
    object@modifications <- res
    return(object)
  }
)

.get_single_position_letters <- function(x) {
  x <- as.character(x)
  substring(x, 1:nchar(x), 1:nchar(x))  
}

# detect and merge modification positions in one transcript
.analyze_transcript_prep <- function(ID,
                                     data,
                                     gff,
                                     fafile,
                                     modClasses){
  # debug
  if( getOption("RNAmodR_debug") ){
    .print_transcript_info(paste(ID," prep"), "")
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
  # get data
  data <- .get_data(ID,data)
  # split data
  data <- lapply(data,"[[","data")
  # analyze the transcript
  res <- .analyze_transcript(ID = ID,
                             modClasses = modClasses,
                             data = data,
                             locations = locations)
  res <- res[order(as.numeric(unlist(lapply(res, "[[", "location"))))]
  if(is.null(res)) return(NULL)
  return(res)
}

# returns the position data for analysis as a list of data per replicate for
# individual transcript
.get_data <- function(ID,data){
  res <- lapply(data,"[[",ID)
  return(res)
}

# merge positions in one transcript
.analyze_transcript <- function(ID,
                                modClasses,
                                data,
                                locations){
  if( iterationN > .get_transcript_max_iteration()) return(NULL)
  # debug
  if( getOption("RNAmodR_debug") ){
    .print_transcript_info(paste(ID," detect"), iterationN)
  }
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
  if(length(modifications) == 0) return(NULL)
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
  if(length(res) > 1) return(NULL)
  if(length(res) == 0) return(NULL)
  # Return data
  return(list(location = location,
              type = res[[1]]$type,
              signal = res[[1]]$signal,
              signal.sd = res[[1]]$signal.sd,
              p.value = res[[1]]$p.value,
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
    message(Sys.time())
    # Process only genes found in all datasets
    IDs <- lapply(object@data,names)
    IDs <- Reduce(intersect, IDs)
    res <- lapply(IDs,
    # res <- BiocParallel::bplapply(IDs,
                                  FUN = .merge_positions,
                                  object@data)
    names(res) <- IDs
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
  data <- .get_data(ID,data)
  data <- lapply(data,"[[","data")
  positions <- as.numeric(unique(unlist(lapply(data,names))))
  res <- lapply(positions,
                FUN = .merge_position,
                data)
  df <- data.frame(pos = unlist(lapply(res,"[[","pos")),
                   mean = unlist(lapply(res,"[[","mean")),
                   sd = unlist(lapply(res,"[[","sd")),
                   stringsAsFactors = FALSE)
  return(df)
}

# merge position data for one position
.merge_position <- function(pos,
                            data){
  data <- unlist(lapply(data, "[[",pos))
  data[is.na(data)] <- 0
  return(list(pos = pos,
              n = length(data),
              mean = mean(data),
              sd = stats::sd(data)))
}

