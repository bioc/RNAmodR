#' @include class-RNAmod-type.R
NULL

RNAMOD_M7G_ROLLING_MEAN_WINDOW_WIDTH <- 9
RNAMOD_M7G_BASE_WINDOW_WIDTH <- 100
RNAMOD_M7G_P_THRESHOLD <- 0.05
RNAMOD_M7G_SIGMA_THRESHOLD <- 10


#' @rdname mod
#'
#' @description 
#' \code{mod_m7G}
#'
#' @return
#' @export
#'
#' @examples
setClass("mod_m7G",
         contains = "mod",
         prototype = list(modType = "m7G")
)


#' @rdname convertReadsToPositions
#'
#' @description
#' \code{mod_m7G}: calls the default method
#'
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "convertReadsToPositions",
  signature = signature(object = "mod_m7G",
                        counts = "numeric",
                        gff = "GRanges",
                        data = "DataFrame"),
  definition = function(object,
                        counts,
                        gff,
                        data) {
    return(NA)
  }
)


#' @rdname parseMod
#' 
#' @description 
#' \code{mod_m7G}
#' 
#' @return
#' @export
#' 
#' @importFrom stringr str_locate_all
#'
#' @examples
setMethod(
  f = "parseMod",
  signature = signature(object = "mod_m7G",
                        gff = "GRanges",
                        seq = "FaFile",
                        data = "list"),
  definition = function(object,
                        gff,
                        seq,
                        data) {
    # browser()
    # Process only genes found in all datasets
    IDs <- lapply(data,names)
    IDs <- Reduce(intersect, IDs)
    
    # detect modification per transcript
    res <- lapply(IDs,
    # res <- bplapply(IDs,
                  FUN = .analyze_transcript,
                  data,
                  gff,
                  seq)
    names(res) <- IDs
    res <- res[!is.null(res)]
    
    # If not results are present return NA instead of NULL
    if(is.null(res)){
      return(NA)
    }
    return(res)
  }
)


# returns the position data for m7G analysis.
# each entry in list a result for a gene of all replicates
.get_m7G_data <- function(ID,data){
  res <- lapply(data,function(x){
    return(x[[ID]][["default"]])
  })
  return(res)
}


# merge positions in one transcript
.analyze_transcript <- function(ID,data,gff,fafile){
  data <- .get_m7G_data(ID,data)
  
  # get sequence of transcript and subset gff for single transcript data
  gff <- .subset_gff_for_unique_transcript(gff, ID)
  seq <- .get_seq_for_unique_transcript(gff,fafile,ID)
  
  # detect all G positions
  loc <- stringr::str_locate_all(as.character(seq), "G")
  loc <- loc[[1]][,"start"]
  if(length(loc) == 0) return(NULL)
  # remove data on pos 1
  data <- lapply(data, function(x){
    x[as.numeric(names(x)) == 1] <- 0
    x
  })
  
  # Convert local G position to global positions
  locations <- .convert_local_to_global_locations(gff, loc)
  
  # if a transcript is encountered having ID present in the global variable
  # RNAmod_sample_transcript trigger sample plotting to access quality on
  # know locations
  # ID %in% options("RNAmod_sample_transcripts")
  
  name <- NULL
  if(ID %in% options("RNAmod_sample_transcripts")){
    name <- paste0(ID,
                   "_m7G_")
  }
  
  # Calculate a rolling mean to reduce noise artefacts
  # rolling_mean <- function(x,
  #                          n=RNAMOD_M7G_ROLLING_MEAN_WINDOW_WIDTH){
  #   y <- setNames(filter(x,rep(1/n,n), sides=2),names(x))
  #   y[is.na(y)] <- x[is.na(y)]
  #   return(y)
  # }
  # meanData <- lapply(data, rolling_mean)
  # Testing per position will be done against the arrest difference data
  # .get_arrest_diff <- function(i){
  #   x <- data[[i]]
  #   y <- unlist(lapply(seq_along(x), function(j){
  #     x[j]-meanData[[i]][j]
  #   }))
  #   setNames(y,names(x))
  # }
  # data <- lapply(seq_along(data), .get_arrest_diff)
  
  # Calculate the arrest rate per position
  .get_arrest_rate <- function(x){
    y <- unlist(lapply(seq_along(x), function(i){
      a <- x[i]
      b <- x[i+1]
      if(is.na(b) || b == 0) return(-1)
      if( a <= b ) return(1-a/b)
      if(a == 0) return(-1)
      return(-b/a)
    }))
    setNames(y,names(x))
  }
  arrestData <- lapply(data, .get_arrest_rate)
  
  # Retrieve m7G positions
  modifications <- lapply(locations,
                          .check_for_m7G,
                          data,
                          arrestData,
                          locations,
                          name)

  if(length(modifications) == 0) return(NULL)
  # name the locations based on sequence position
  names(modifications) <- paste0("G_",loc)
  modifications <-  modifications[!vapply(modifications,is.null,logical(1))]
  if(length(modifications) == 0) return(NULL)
  return(modifications)
}

.check_for_m7G <- function(location, 
                           data, 
                           arrestData,
                           locs, 
                           name = NULL){
  # if( location == 1575){
  #   browser()
  # }
  
  # number of replicates
  nbSamples <- length(data)
  
  # merge data for positions
  # data on the N+1 location
  testData <- .aggregate_location_data(data, (location+1))
  # data on the arrect direction
  testArrestData <- .aggregate_location_data(arrestData, location)
  # base data to compare against
  # use only G position
  baseData <- .aggregate_not_location_data(lapply(data, function(x){
    x[as.numeric(names(x)) %in% (locs+1)]
  }), (location+1))
  
  locTest <- .calc_m7G_test_values(location,
                                   testData,
                                   baseData,
                                   testArrestData,
                                   nbSamples)
  # If insufficient data is present
  if(is.null(locTest)) return(NULL)
  
  # dynamic threshold based on the noise of the signal (high sd)
  threshold <- RNAMOD_M7G_SIGMA_THRESHOLD + locTest$thresholdAddition
  if(!.validate_m7G_pos(threshold, 
                       RNAMOD_M7G_P_THRESHOLD, 
                       locTest$sig.mean, 
                       locTest$p.value) ) return(NULL)
  
  # browser()
  # if location is among sample location (name is not null)
  # plot the data
  if(!is.null(name)){
    .plot_sample_data(.create_plot_data(testData,
                                        baseData,
                                        paste0(name,location)), 
                      paste0(name,location))
  }
    
  # Return data
  return(list(location = location,
              signal = locTest$sig.mean,
              signal.sd = locTest$sig.sd,
              p.value = locTest$p.value,
              nbsamples = nbSamples))
}

.aggregate_location_data <- function(data, 
                                     location){
  unlist(lapply(data,function(dataPerReplicate){
    dataPerReplicate <- dataPerReplicate[dataPerReplicate > 0]
    return(dataPerReplicate[as.numeric(names(dataPerReplicate)) == location])
  }))
}
.aggregate_not_location_data <- function(data,
                                         location){
  unlist(lapply(data,function(dataPerReplicate){
    dataPerReplicate <- dataPerReplicate[dataPerReplicate > 0]
    return(dataPerReplicate[as.numeric(names(dataPerReplicate)) != location])
  }))
}
.aggregate_area_data <- function(data, 
                                 location, 
                                 width){
  unlist(lapply(data,function(dataPerReplicate){
    return(dataPerReplicate[as.numeric(names(dataPerReplicate)) < (location+width) &
                              as.numeric(names(dataPerReplicate)) > (location-width) &
                              as.numeric(names(dataPerReplicate)) != location])
  }))
}

.calc_m7G_test_values <- function(location,
                                  testData,
                                  baseData,
                                  testArrestData,
                                  n){
  
  if(length(testData) == 0 | 
     length(baseData) < n) return(NULL)
  
  # No read arrest detectable
  if( sum(testArrestData) < 0 ) return(NULL)
  
  # get test values
  # overall mean and sd
  mean <-  mean(baseData)
  sd <-  stats::sd(baseData)
  thresholdAddition <- abs(sd %/% mean)
    
  # Use the sigma level as value for signal strength
  sig <- (as.numeric(as.character(testData)) - mean) %/% sd
  sig.mean <- mean(sig)
  sig.sd <- stats::sd(sig)
  # Since normality of distribution can not be assumed use the MWU
  # generate p.value for single position
  p.value <- suppressWarnings(wilcox.test(baseData, testData)$p.value)
  
  return(list(thresholdAddition = thresholdAddition,
              sig = sig,
              sig.mean = sig.mean,
              sig.sd = sig.sd,
              p.value = p.value))
}

.validate_m7G_pos <- function(sig.threshold, 
                              p.threshold, 
                              sig, 
                              p.value){
  ((sig > sig.threshold &&
      p.value <= p.threshold) ||
     (sig > sig.threshold &&
        !.getUseP()))
}

.create_plot_data <- function(testData, baseData, name){
  data.frame(x = c(rep(name,(length(testData) + 
                               length(baseData)))),
             y = c(testData, 
                   baseData),
             group = c(rep("Position",length(testData)), 
                       rep("Base", length(baseData))))
}

#' @rdname mergePositionsOfReplicates
#'
#' @description
#' \code{mod_m7G}
#'
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "mergePositionsOfReplicates",
  signature = signature(object = "mod_m7G",
                        gff = "GRanges",
                        seq = "FaFile",
                        data = "list"),
  definition = function(object,
                        gff,
                        seq,
                        data) {
    return(NA)
  }
)
