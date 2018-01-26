#' @include class-RNAmod-mod-type.R
#' @include class-RNAmod-analysis-type.R
NULL

#' @rdname mod
#'
#' @description 
#' \code{analysis_default}
#'
#' @return
#' @export
#'
#' @examples
setClass("analysis_default",
         contains = "analysis")


#' @rdname convertReadsToPositions
#'
#' @description
#' \code{analysis_default}: calls the default method
#'
#' @return
#' @export
#' 
#' @importFrom IRanges findOverlaps extractList
#' @importFrom S4Vectors split from to
#'
#' @examples
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
                   FUN = .get_positions,
                   gff,
                   param)
    data <- data[!is.null(data)]
    if(length(data) == 0){
      stop("No reads detected in any bam file for '",
           paste(modifications, collapse = "', '"),
           "' modifications",
           call. = FALSE)
    }
    object@data <- data
    return(object)
  }
)


# detect modifications in each file
.get_positions <- function(bamFile,
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
  
  bamData <- bamData[names(bamData) %in% c("RDN18-1",
                                           "YBR041W",
                                           "YBR056W",
                                           "YAL030W")]
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
  positions <- mapply(
  # positions <- BiocParallel::bpmapply(
                                      FUN = .get_positions_in_transcript,
                                      bamData,
                                      names(bamData),
                                      MoreArgs = list(totalCounts,
                                                      gff),
                                      SIMPLIFY = FALSE)
  names(positions) <- names(bamData)
  return(positions)
}


# For each transcript get positional data
# This can be individually done for different modification types
.get_positions_in_transcript <- function(data,
                                         id,
                                         counts,
                                         gff){
  # skip if transcript does not have data
  if(length(data) == 0) return(NULL)
  # get ID and GRanges
  gr <- .subset_gff_for_unique_transcript(gff, 
                                          id)
  # get the genomic distance
  # also used in .convert_global_to_local_position
  length <- width(gr)
  # do position conversion to translate genomic position to local transcript
  # position. take care of introns, etc
  pos <- .convert_global_to_local_position(gff,gr,data)
  # Normalize counts per positions against million of reads in BamFile
  posData <- table(pos)/(counts/10^6)
  # spread table with zero values to the length of transcript
  posData <- setNames(as.double(unlist(lapply(1:length, 
                                              function(i){
    if(length(posData[names(posData) == i]) == 0) return(0)
    posData[names(posData) == i]
  }))),1:length)
  return(posData)
}


#' @rdname parseMod
#' 
#' @description 
#' \code{analysis_default}
#' 
#' @return
#' @export
#' 
#' @importFrom stringr str_locate_all
#'
#' @examples
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
  substring(x, 1:nchar(x), 1:nchar(x))  
}

# detect and merge modification positions in one transcript
.analyze_transcript_prep <- function(ID,
                                     data,
                                     gff,
                                     fafile,
                                     modClasses){
  # browser()
  # debug
  if( getOption("RNAmod_debug") ){
    .print_transcript_info(paste(ID," prep"), "")
  }
  data <- .get_data(ID,data)
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
  res <- .analyze_transcript(ID = ID,
                             modClasses = modClasses,
                             data = data,
                             locations = locations,
                             iterationN = 1)
  res <- res[order(as.numeric(unlist(lapply(res, "[[", "location"))))]
  if(is.null(res)) return(NULL)
  # recalculate sigma values by masking all positions except the one testing for
  names <- names(res)
  res <- lapply(seq_along(res), function(i){
    x <- data
    testPosition <- res[[i]]
    nonTestPositions <- res[seq_along(res) != i]
    if( length(nonTestPositions) > 0 ){
      locs <- setNames(
        as.numeric(unlist(lapply(nonTestPositions, "[[","location"))),
        unlist(lapply(nonTestPositions, "[[","type")))
      x <- .mask_data(modClasses,
                      x,
                      locs)
    }
    return(checkForModification(modClasses[[testPosition$type]],
                                location = testPosition$location,
                                locations = locations,
                                data = x))
  })
  names(res) <- names
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
                                locations,
                                iterationN){
  if( iterationN > .get_transcript_max_iteration()) return(NULL)
  # debug
  if( getOption("RNAmod_debug") ){
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
  # check if by masking found modifications positions additional positions are
  # picked up - remember the +1 location data
  modLocations <- setNames(
    as.numeric(unlist(lapply(modifications, "[[","location"))),
    unlist(lapply(modifications, "[[","type")))
  # do not take into account positions found as modification
  data <- .mask_data(modClasses,
                     data,
                     modLocations)
  locations <- locations[!(locations %in% modLocations)]
  return(append(modifications,.analyze_transcript(ID = ID,
                                                  modClasses = modClasses,
                                                  data = data,
                                                  locations = locations,
                                                  (iterationN+1))))
}

# mask data for know modification locations
.mask_data <- function(modClasses,
                       data,
                       modLocations){
  # browser()
  # only use class for which data is present
  modTypes <- unique(names(modLocations))
  modClassTypes <- unlist(lapply(modClasses, getModType))
  data <- lapply(data, function(x){
    modClassesSubset <- modClasses[modClassTypes %in% modTypes]
    for(i in seq_along(modClassesSubset)){
      subsetLocations <- 
        modLocations[names(modLocations) == getModType(modClassesSubset[[i]])]
      x <- maskPositionData(modClassesSubset[[i]],
                       x,
                       subsetLocations)
    }
    x
  })
  return(data)
}

# check for modifications at given position
.check_for_modification <- function(location, 
                                    modClasses,
                                    data,
                                    locations){
  
  # if(location == 1575 | location == 599 | location == 1420) { browser() }
  res <- lapply(modClasses, function(class){
    # short cut if amount/properties of data are not sufficient
    if( is.null(preTest(class,
                        location,
                        data,
                        locations))) return(NULL)
    # If potential modification right in front of current location
    if(length(locations[locations == (location-1)]) != 0) {
      locTestPre <- .check_for_modification((location-1),
                                            modClasses,
                                            .mask_data(list(class),
                                                       data, 
                                                       .get_location_vector(location,
                                                                            getModType(class))),
                                            locations[locations != location])
      if(!is.null(locTestPre)){
        # udpate data accordingly
        data <- .mask_data(modClasses,
                           data, 
                           .get_location_vector(locTestPre$location,
                                                locTestPre$type))
      }
    }
    # If potential modification right after of current location
    if(length(locations[locations == (location+1)]) != 0) {
      locTestPost <- .check_for_modification((location+1), 
                                             modClasses,
                                             .mask_data(list(class),
                                                        data, 
                                                        .get_location_vector(location,
                                                                             getModType(class))),
                                             locations[locations != location])
      if(!is.null(locTestPost)){
        # udpate data accordingly
        data <- .mask_data(modClasses,
                           data, 
                           .get_location_vector(locTestPost$location,
                                                locTestPost$type))
      }
    }
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

.get_location_vector <- function(location, 
                                 type){
  setNames(location,
           rep(type,length(location)))
}

#' @rdname mergePositionsOfReplicates
#'
#' @description
#' \code{analysis_default}
#'
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "mergePositionsOfReplicates",
  signature = signature(object = "analysis_default",
                        gff = "GRanges",
                        fafile = "FaFile"),
  definition = function(object,
                        gff,
                        fafile) {
    message(Sys.time())
    # Process only genes found in all datasets
    IDs <- lapply(object@data,names)
    IDs <- Reduce(intersect, IDs)
    # res <- lapply(IDs,
    res <- BiocParallel::bplapply(IDs,
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

# returns the position data for m7G analysis.
# each entry in list a result for a gene of all replicates
.get_default_data <- function(ID,data){
  res <- lapply(data,function(x){
    return(x[[ID]])
  })
  return(res)
}

# merge positions in one transcript
.merge_positions <- function(ID,data){
  data <- .get_data(ID,data)
  positions <- unique(unlist(lapply(data,names)))
  res <- lapply(positions,
                FUN = .merge_position,
                data)
  df <- data.frame(pos = unlist(lapply(res,"[[","pos")),
                   mean = unlist(lapply(res,"[[","mean")),
                   sd = unlist(lapply(res,"[[","sd")))
  return(df)
}

# merge position data for one position
.merge_position <- function(pos,data){
  data <- unlist(lapply(data, "[[",pos))
  data[is.na(data)] <- 0
  return(list(pos = pos,
              mean = mean(data),
              sd = stats::sd(data)))
}

