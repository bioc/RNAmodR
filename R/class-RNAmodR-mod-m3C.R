#' @include class-RNAmodR-mod-type.R
NULL

RNAMODR_M3C_NUCLEOTIDE <- "C"
RNAMODR_M3C_ARREST_RATE <- 0.65
RNAMODR_M3C_Z_THRESHOLD <- 3
RNAMODR_M3C_SIG_THRESHOLD <- 5


#' @rdname RNAmodR-mod-class
#'
#' @description 
#' \code{mod_m3C}:
#' 
#' @export
setClass("mod_m3C",
         contains = "mod",
         prototype = list(modType = "m3C")
)

#' @rdname preTest
#' 
#' @description 
#' \code{mod_m3C}
#' 
#' @export
setMethod(
  f = "preTest",
  signature = signature(object = "mod_m3C",
                        location = "numeric",
                        locations = "numeric",
                        data = "list"),
  definition = function(object,
                        location,
                        locations,
                        data) {
    # do pretest
    res <- .do_M3C_pretest(location,
                           locations,
                           data)
    return(res)
  }
)

# check if any data is available to proceed with test
# this is in a seperate function since it is also called by checkForModification
.do_M3C_pretest <- function(location,
                            locations,
                            data){
  # if non C position skip position
  if( names(locations[locations == location]) 
      != RNAMODR_M3C_NUCLEOTIDE){
    return(NULL)
  }
  # do not take into account position 1
  if(location == 1) return(NULL)
  # number of replicates
  n <- length(data)
  # merge data for positions
  # data on the N+1 location
  testData <- .aggregate_location_data(data, (location + 1))
  # if number of data points is not high enough
  if(length(testData) < n) return(NULL)
  # if stop coverage is to low
  if(length(testData[testData > RNAMODR_M3C_ARREST_RATE]) < n) return(NULL)
  # base data to compare against
  baseData <- .aggregate_area_data(data, 
                                   (location + 1), 
                                   50)
  baseData <- baseData[baseData > 0]
  # if not enough data is present
  if(length(baseData) < (3 * n)) return(NULL)
  return(list(n = n,
              testData = testData,
              baseData = baseData))
}

#' @rdname checkForModification
#' 
#' @description 
#' \code{mod_m3C}
#' 
#' @export
setMethod(
  f = "checkForModification",
  signature = signature(object = "mod_m3C",
                        location = "numeric",
                        locations = "numeric",
                        data = "list"),
  definition = function(object,
                        location,
                        locations,
                        data) {
    # get test result for the current location
    locTest <- .calc_M3C_test_values(location,
                                     locations,
                                     data)
    # If insufficient data is present
    if(is.null(locTest)) return(NULL)
    # dynamic threshold based on the noise of the signal (high sd)
    if(!.validate_M3C_pos(RNAMODR_M3C_SIG_THRESHOLD, 
                          RNAMODR_M3C_Z_THRESHOLD, 
                          locTest$sig.mean, 
                          locTest$z) ) {
      # debug
      if( getOption("RNAmodR_debug") ){
        .print_location_info(paste(location,"_no"),locations)
      }
      return(NULL)
    }
    # debug
    if( getOption("RNAmodR_debug") ){
      .print_location_info(paste(location,"_yes"), locations)
    }
    # Return data
    return(list(location = location,
                type = getModType(object),
                signal = locTest$sig.mean,
                signal.sd = locTest$sig.sd,
                z = locTest$z,
                nbsamples = locTest$n))
  }
)

# test for m3C at current location
.calc_M3C_test_values <- function(location,
                                  locations,
                                  data){
  # short cut if amount of data is not sufficient
  pretestData <- .do_M3C_pretest(location,
                                 locations,
                                 data)
  if(is.null(pretestData)) return(NULL)
  # data from pretest
  testData <- pretestData$testData
  baseData <- pretestData$baseData
  n <- pretestData$n
  # difference to threshold by relative percent
  # minimal value of 1 since y is one percent lower than arrest rate threshold
  sig <- floor((testData - RNAMODR_M3C_ARREST_RATE) / 
                 (1 - RNAMODR_M3C_ARREST_RATE) * 100)
  sig.mean <- mean(sig)
  sig.sd <- stats::sd(sig)
  # generate z score
  z <- mean((testData - mean(baseData))/sd(baseData))
  return(list(sig = sig,
              sig.mean = sig.mean,
              sig.sd = sig.sd,
              z = z,
              n = n))
}

# call yes or nor position
.validate_M3C_pos <- function(sig.threshold, 
                              p.threshold, 
                              sig, 
                              z){
  ((sig > sig.threshold &&
      z >= p.threshold) ||
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
