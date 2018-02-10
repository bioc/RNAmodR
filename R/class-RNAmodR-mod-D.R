#' @include class-RNAmodR-mod-type.R
NULL

RNAMODR_D_NUCLEOTIDE <- "T"
RNAMODR_D_ARREST_RATE <- 0.6
RNAMODR_D_Z_THRESHOLD <- 3
RNAMODR_D_SIG_THRESHOLD <- 5


#' @rdname RNAmodR-mod-class
#'
#' @description 
#' \code{mod_D}:
#' 
#' @export
setClass("mod_D",
         contains = "mod",
         prototype = list(modType = "D",
                          positionOffset = 1)
)

#' @rdname preTest
#' 
#' @description 
#' \code{mod_D}
#' 
#' @export
setMethod(
  f = "preTest",
  signature = signature(object = "mod_D",
                        location = "numeric",
                        locations = "numeric",
                        data = "list"),
  definition = function(object,
                        location,
                        locations,
                        data) {
    # do pretest
    res <- .do_D_pretest(location,
                         locations,
                         data)
    return(res)
  }
)

# check if any data is available to proceed with test
# this is in a seperate function since it is also called by checkForModification
.do_D_pretest <- function(location,
                          locations,
                          data){
  # if non D position skip position
  if( names(locations[locations == location]) != RNAMODR_D_NUCLEOTIDE){
    return(NULL)
  }
  # do not take into account position 1
  if(location == 1) return(NULL)
  # split into conditions
  res <- mapply(FUN = .get_data_per_condition_D,
                split(data,names(data)),
                MoreArgs = list(location = location),
                SIMPLIFY = FALSE)
  res <- res[!vapply(res,is.null,logical(1))]
  # check if Treated condition has valid data. Otherwise return NULL (Abort)
  if(is.null(res$Treated)) return(NULL)
  return(res)
}

.get_data_per_condition_D <- function(data, location){
  # number of replicates
  n <- length(data)
  # merge data for positions
  # data on the N+1 location
  testData <- .aggregate_location_data(data, 
                                       (location + 1))
  # if number of data points is not high enough
  if(any(is.na(testData))) return(NULL)
  # base data to compare against
  baseData <- .aggregate_area_data(data, 
                                   (location + 1), 
                                   50)
  # if not enough data is present
  if(sum(vapply(lapply(baseData, 
                       function(x){x[!is.na(x)]}),
                length,
                numeric(1))) < n) {
    return(NULL)
  }
  return(list(n = n,
              testData = testData,
              baseData = baseData))
}

#' @rdname checkForModification
#' 
#' @description 
#' \code{mod_D}
#' 
#' @export
setMethod(
  f = "checkForModification",
  signature = signature(object = "mod_D",
                        location = "numeric",
                        locations = "numeric",
                        data = "list"),
  definition = function(object,
                        location,
                        locations,
                        data) {
    # get test result for the current location
    locTest <- .calc_D_test_values(location,
                                   locations,
                                   data)
    # If insufficient data is present
    if(is.null(locTest)) return(NULL)
    # dynamic threshold based on the noise of the signal (high sd)
    if(!.validate_D_pos(RNAMODR_D_SIG_THRESHOLD, 
                        RNAMODR_D_Z_THRESHOLD, 
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

# test for D at current location
.calc_D_test_values <- function(location,
                                locations,
                                data){
  if(location == 20) browser()
  # short cut if amount of data is not sufficient
  pretestData <- .do_D_pretest(location,
                               locations,
                               data)
  if(is.null(pretestData)) return(NULL)
  # If Control sample is available
  if(!is.null(pretestData$Control)){
    # data from pretest
    testData <- pretestData$Treated$testData
    baseData <- pretestData$Treated$baseData
    n <- pretestData$Treated$n
    testDataC <- pretestData$Control$testData
    baseDataC <- pretestData$Control$baseData
    nC <- pretestData$Control$n
    testDataCombined <- unlist(testData) - mean(unlist(testDataC))
    # if stop coverage is to low
    if(any(testDataCombined < (RNAMODR_D_ARREST_RATE/2))) return(NULL)
    # difference to threshold by relative percent
    # minimal value of 1 since y is one percent lower than arrest rate threshold
    sig <- floor((testDataCombined - (RNAMODR_D_ARREST_RATE / 2)) / 
                   (1 - (RNAMODR_D_ARREST_RATE / 2)) * 100)
    sig.mean <- mean(sig)
    # approx. sd since covariance of testData and testDataC cannot not be 
    # assumed in case of unequal observations
    sig.sd <- sqrt(stats::sd(unlist(testData))^2 + 
                     stats::sd(unlist(testDataC))^2)
    #
    baseData <- .merge_base_data_D(baseData,
                                   baseDataC)
    # generate z score
    z <- mean((testDataCombined - mean(baseData)) / sd(baseData))
  } else {
    # data from pretest
    testData <- unlist(pretestData$Treated$testData)
    baseData <- unlist(pretestData$Treated$baseData)
    n <- pretestData$Treated$n
    # if stop coverage is to low
    if(any(testData < RNAMODR_D_ARREST_RATE)) return(NULL)
    browser()
    # difference to threshold by relative percent
    # minimal value of 1 since y is one percent lower than arrest rate threshold
    sig <- floor((testData - RNAMODR_D_ARREST_RATE) / 
                   (1 - RNAMODR_D_ARREST_RATE) * 100)
    sig.mean <- mean(sig)
    sig.sd <- stats::sd(sig)
    # generate z score
    z <- mean((testData - mean(baseData)) / sd(baseData))
  }
  return(list(sig = sig,
              sig.mean = sig.mean,
              sig.sd = sig.sd,
              z = z,
              n = n))
}
# iterates on every position and calculates the difference of the means
.merge_base_data_D <- function(treated,
                                 control){
  df <- data.frame(append(treated,control))
  colnames(df) <- c(names(treated),names(control))
  df <- df[complete.cases(df),]
  treated <- df[,colnames(df) == "Treated"]
  control <- df[,colnames(df) == "Control"]
  res <- rowMeans(treated) - rowMeans(control)
  res[res < 0] <- 0
  return(res)
}

# call yes or nor position
.validate_D_pos <- function(sig.threshold,
                            p.threshold,
                            sig,
                            z){
  ((sig > sig.threshold &&
      z >= p.threshold) ||
     (sig > sig.threshold &&
        !.get_use_p()))
}

# create data.frame for sample plotting
.create_D_plot_data <- function(testData, baseData, name){
  data.frame(x = c(rep(name,(length(testData) + 
                               length(baseData)))),
             y = c(testData, 
                   baseData),
             group = c(rep("Position",length(testData)), 
                       rep("Baseline", length(baseData))))
}
