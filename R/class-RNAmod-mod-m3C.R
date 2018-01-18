#' @include class-RNAmod-mod-type.R
NULL

RNAMOD_M3C_ROLLING_MEAN_WINDOW_WIDTH <- 9
RNAMOD_M3C_ARREST_RATE <- 0.95
RNAMOD_M3C_P_THRESHOLD <- 0.05
RNAMOD_M3C_SIGMA_THRESHOLD <- 3

#' @rdname mod
#'
#' @description 
#' \code{mod_m3C}
#'
#' @return
#' @export
#'
#' @examples
setClass("mod_m3C",
         contains = "mod",
         prototype = list(modType = "m3C")
)

#' @rdname parseMod
#' 
#' @description 
#' \code{mod_m3C}
#' 
#' @return
#' @export
#' 
#' @importFrom stringr str_locate_all
#' @importFrom BiocParallel bplapply
#'
#' @examples
setMethod(
  f = "checkForModification",
  signature = signature(object = "mod_D",
                        data = "DataFrame",
                        locations = "numeric"),
  definition = function(object,
                        data,
                        locations) {
    
    
    
    
    
    
    # short cut if amount of data is not sufficient
    if( is.null(.do_M3C_pretest(location,
                                data))) return(NULL)
    # If potential modification right in front of current location
    if(length(locs[locs == (location-1)]) != 0) {
      locTestPre <- .check_for_M3C((location-1), 
                                   .mask_data(data, location),
                                   locs[locs != location])
      if(!is.null(locTestPre)){
        # udpate data accordingly
        data <- .mask_data(data, (location-1))
      }
    }
    # If potential modification right after of current location
    if(length(locs[locs == (location+1)]) != 0) {
      locTestPost <- .check_for_M3C((location+1), 
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
    locTest <- .calc_M3C_test_values(location,
                                     data,
                                     arrestData)
    # If insufficient data is present
    if(is.null(locTest)) return(NULL)
    # dynamic threshold based on the noise of the signal (high sd)
    if(!.validate_M3C_pos(RNAMOD_M3C_SIGMA_THRESHOLD, 
                          RNAMOD_M3C_P_THRESHOLD, 
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
)


# check if any data is available to proceed with test
.do_M3C_pretest <- function(location,
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

# test for m3C at current location
.calc_M3C_test_values <- function(location,
                                  data,
                                  arrestData){
  # short cut if amount of data is not sufficient
  pretestData <- .do_M3C_pretest(location,
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
  testArrest <- length(testArrestData[testArrestData >= RNAMOD_M3C_ARREST_RATE])
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
.validate_M3C_pos <- function(sig.threshold, 
                              p.threshold, 
                              sig, 
                              p.value){
  ((sig > sig.threshold &&
      p.value <= p.threshold) ||
     (sig > sig.threshold &&
        !.get_use_p()))
}

# create data.frame for sample plotting
.create_M3C_plot_data <- function(testData, baseData, name){
  data.frame(x = c(rep(name,(length(testData) + 
                               length(baseData)))),
             y = c(testData, 
                   baseData),
             group = c(rep("Position",length(testData)), 
                       rep("Baseline", length(baseData))))
}
