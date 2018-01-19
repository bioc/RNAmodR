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
    browser()
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
  # Construct Dataframe from scanBam data
  bamData <- Rsamtools::scanBam(bamFile, param=param)
  bamData <- .convert_bam_to_DataFrame(bamData, param=param)
  # Total counts
  totalCounts <- Rsamtools::idxstatsBam(bamFile, param=param)
  totalCounts <- sum(totalCounts$mapped)
  # process result
  bamData <- S4Vectors::split(bamData,
                              bamData$ID)
  # browser()
  # for testing
  # m7G 
  # bamData <- bamData[names(bamData) %in% c("RDN18-1")]
  # bamData <- bamData[names(bamData) %in% c("tC(GCA)B")]
  
  # D 
  # bamData <- bamData[names(bamData) %in% c("tH(GUG)E1")]
  # bamData <- bamData[names(bamData) %in% c("tI(AAU)B")]
  # bamData <- bamData[names(bamData) %in% c("tD(GUC)B")]
  
  # m3C
  # bamData <- bamData[names(bamData) %in% c("tS(CGA)C")]
  # bamData <- bamData[names(bamData) %in% c("RDN25-1")]
  
  # combination
  # bamData <- bamData[names(bamData) %in% c("RDN18-1",
  #                                          "tS(CGA)C",
  #                                          "tC(GCA)B")]
  # if( getOption("RNAmod_debug") ){
  #   bamData <- bamData[names(bamData) %in% getOption("RNAmod_debug_transcripts")]
  # }
  
  if(length(bamData) == 0){
    warning("No reads detected in bam file '",
            bamFile,
            "'")
    return(NULL)
  }
  positions <- lapply(bamData,
                      # res <- BiocParallel::bplapply(bamData,
                      FUN = .get_positions_in_transcript,
                      totalCounts,
                      gff)
  names(positions) <- names(bamData)
  return(positions)
}


# For each transcript get positional data
# This can be individually done for different modification types
.get_positions_in_transcript <- function(data,counts,gff){
  # get ID and GRanges
  gff <- .subset_gff_for_unique_transcript(gff, unique(data$ID))
  # fill up empty positions
  strand <- as.character(strand(gff))
  if(strand == "-"){
    pos <- data$pos + data$qwidth - 1
  } else {
    pos <- data$pos
  }
  # Normalize counts per positions against million of reads in BamFile
  posData <- table(pos)/(counts/10^6)
  # spread table with zero values to the length of transcript
  posData <- setNames(as.double(unlist(lapply(1:width(gff), function(i){
    if(length(posData[names(posData) == i]) == 0) return(0)
    posData[names(posData) == i]
  }))),1:width(gff))
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
                        modifications = "list"),
  definition = function(object,
                        gff,
                        fafile,
                        modifications) {
    browser()
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
                                  mods = modifications)
    names(res) <- IDs
    res <- res[!is.null(res)]
    # If not results are present return NA instead of NULL
    if(is.null(res)){
      res <- NA
    }
    object@modifications <- res
    return(object)
  }
)

# returns the position data for m7G analysis.
# each entry in list a result for a gene of all replicates
.get_data <- function(ID,data){
  res <- lapply(data,"[[",ID)
  return(res)
}

# detect and merge modification positions in one transcript
.analyze_transcript_prep <- function(ID,
                                     data,
                                     gff,
                                     fafile,
                                     mods){
  browser()
  # debug
  if( getOption("RNAmod_debug") ){
    .print_transcript_info(paste(ID," prep"), "")
  }
  data <- .get_data(ID,data)
  # do not take into account position 1
  data <- lapply(data, function(x){
    x[as.numeric(names(x)) == 1] <- 0
    x
  })
  # get sequence of transcript and subset gff for single transcript data
  gff <- .subset_gff_for_unique_transcript(gff, ID)
  seq <- .get_seq_for_unique_transcript(gff,fafile,ID)
  # generate a location vector
  locations <- 1:width(seq)
  # Convert local G position to global positions
  globalLocations <- .convert_local_to_global_locations(gff, locations)
  
  
  browser()
  
  # detect all G positions
  locations <- stringr::str_locate_all(as.character(seq), "G")
  locations <- locations[[1]][,"start"]
  if(length(locations) == 0) return(NULL)
  
  
  
  
  
  res <- .analyze_transcript(ID = ID,
                             mods = mods,
                             data = data,
                             globalLocations = globalLocations,
                             iterationN = 1)
  res <- res[order(as.numeric(unlist(lapply(res, "[[", "location"))))]
  return(res)
}

# merge positions in one transcript
.analyze_transcript <- function(ID,
                                mods,
                                data,
                                globalLocations,
                                iterationN){
  if( iterationN > .get_transcript_max_iteration()) return(NULL)
  # debug
  if( getOption("RNAmod_debug") ){
    .print_transcript_info(paste(ID," detect"), iterationN)
  }
  
  
  
  
  
  
  
  # Retrieve m7G positions
  modifications <- lapply(globalLocations,
                          .check_for_M7G,
                          data,
                          globalLocations)
  if(length(modifications) == 0) return(NULL)
  
  
  
  
  
  
  # name the locations based on sequence position
  names(modifications) <- paste0(ID,"_G_",globalLocations)
  modifications <-  modifications[!vapply(modifications,is.null,logical(1))]
  if(length(modifications) == 0) return(NULL)
  # check if by masking found modifications positions additional positions are
  # picked up - remember the +1 location data
  modLocations <- as.numeric(unlist(lapply(modifications, "[[","location")))
  # do not take into account positions found as modification
  data <- .mask_data(data, modLocations)
  globalLocations <- globalLocations[!(globalLocations %in% modLocations)]
  return(append(modifications,.analyze_M7G_transcript(ID = ID,
                                                      data = data,
                                                      globalLocations = globalLocations,
                                                      (iterationN+1))))
}

# mask data for know modification locations
# extend to N+2,N+3 location?
.mask_data <- function(data, modLocations){
  lapply(data, function(x){
    x[as.numeric(names(x)) %in% (modLocations+1)] <- 
      x[as.numeric(names(x)) %in% (modLocations+1)]*(1-RNAMOD_M7G_ARREST_RATE)
    x
  })
}

# check for m7G at given position
.check_for_M7G <- function(location, 
                           data,
                           locs){
  # short cut if amount of data is not sufficient
  if( is.null(.do_M7G_pretest(location,
                              data))) return(NULL)
  # If potential modification right in front of current location
  if(length(locs[locs == (location-1)]) != 0) {
    locTestPre <- .check_for_M7G((location-1), 
                                 .mask_data(data, location),
                                 locs[locs != location])
    if(!is.null(locTestPre)){
      # udpate data accordingly
      data <- .mask_data(data, (location-1))
    }
  }
  # If potential modification right after of current location
  if(length(locs[locs == (location+1)]) != 0) {
    locTestPost <- .check_for_M7G((location+1), 
                                  .mask_data(data, location),
                                  locs[locs != location])
    if(!is.null(locTestPost)){
      # udpate data accordingly
      data <- .mask_data(data, (location+1))
    }
  }
  # Calculate the arrest rate per position
  arrestData <- lapply(data, .get_arrest_rate)
  # get test result for the current location
  locTest <- .calc_M7G_test_values(location,
                                   data,
                                   arrestData)
  # If insufficient data is present
  if(is.null(locTest)) return(NULL)
  # dynamic threshold based on the noise of the signal (high sd)
  if(!.validate_M7G_pos(RNAMOD_M7G_SIGMA_THRESHOLD, 
                        RNAMOD_M7G_P_THRESHOLD, 
                        locTest$sig.mean, 
                        locTest$p.value) ) {
    # debug
    if( getOption("RNAmod_debug") ){
      .print_location_info(paste(location,"_no"),locs)
    }
    return(NULL)
  }
  # debug
  if( getOption("RNAmod_debug") ){
    .print_location_info(paste(location,"_yes"), locs)
  }
  # Return data
  return(list(location = location,
              signal = locTest$sig.mean,
              signal.sd = locTest$sig.sd,
              p.value = locTest$p.value,
              nbsamples = locTest$n))
}

# check if any data is available to proceed with test
.do_M7G_pretest <- function(location,
                            data){
  # do not take into account position 1
  if(location == 1) return(NULL)
  
  # merge data for positions
  # data on the N+1 location
  testData <- .aggregate_location_data(data, (location+1))
  testData <- testData[testData > 0]
  # base data to compare against
  baseData <- .aggregate_not_location_data(data, (location+1))
  
  # number of replicates
  n <- length(data)
  # if not enough data is present
  if(length(testData) == 0 | 
     length(baseData) < (3*n)) return(NULL)
  return(list(n = n,
              testData = testData,
              baseData = baseData))
}

# test for m7G at current location
.calc_M7G_test_values <- function(location,
                                  data,
                                  arrestData){
  # short cut if amount of data is not sufficient
  pretestData <- .do_M7G_pretest(location,
                                 data)
  if(is.null(pretestData)) return(NULL)
  
  # data from pretest
  testData <- pretestData$testData
  baseData <- pretestData$baseData
  n <- pretestData$n
  
  # data on the arrect direction
  testArrestData <- .aggregate_location_data(arrestData, location)
  # No read arrest detectable
  if( sum(testArrestData) < 0 ) return(NULL)
  # To low arrest detectable
  testArrest <- length(testArrestData[testArrestData >= RNAMOD_M7G_ARREST_RATE])
  if( length(testArrestData) != testArrest ) {
    return(NULL)
  }
  # get test values
  # overall mean and sd
  mean <-  mean(baseData)
  sd <-  stats::sd(baseData)
  # Use the sigma level as value for signal strength
  sig <- (as.numeric(as.character(testData)) - mean) %/% sd
  sig.mean <- mean(sig)
  sig.sd <- stats::sd(sig)
  # Since normality of distribution can not be assumed use the MWU
  # generate p.value for single position
  p.value <- suppressWarnings(wilcox.test(baseData, testData)$p.value)
  return(list(sig = sig,
              sig.mean = sig.mean,
              sig.sd = sig.sd,
              p.value = p.value,
              n = n))
}

# call yes or nor position
.validate_M7G_pos <- function(sig.threshold, 
                              p.threshold, 
                              sig, 
                              p.value){
  ((sig > sig.threshold &&
      p.value <= p.threshold) ||
     (sig > sig.threshold &&
        !.get_use_p()))
}

# create data.frame for sample plotting
.create_M7G_plot_data <- function(testData, baseData, name){
  data.frame(x = c(rep(name,(length(testData) + 
                               length(baseData)))),
             y = c(testData, 
                   baseData),
             group = c(rep("Position",length(testData)), 
                       rep("Baseline", length(baseData))))
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
                        fafile = "FaFile",
                        data = "list"),
  definition = function(object,
                        gff,
                        fafile,
                        data) {
    browser()
    # Process only genes found in all datasets
    IDs <- lapply(data,names)
    IDs <- Reduce(intersect, IDs)
    res <- lapply(IDs,
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
    return(x[[ID]][["default"]])
  })
  return(res)
}


# merge positions in one transcript
.merge_positions <- function(ID,data){
  data <- .get_default_data(ID,data)
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

