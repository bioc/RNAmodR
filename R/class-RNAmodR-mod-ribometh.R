#' @include class-RNAmodR-mod-type.R
NULL

RNAMODR_RIBOMETH_SCORE_A_THRESHOLD <- 0.5
RNAMODR_RIBOMETH_SCORE_B_THRESHOLD <- 3.2

#' @rdname RNAmodR-mod-class
#'
#' @description 
#' \code{mod_ribometh}:
#' 
#' @export
setClass("mod_ribometh",
         contains = "mod",
         prototype = list(dataType = "allends",
                          modType = "ribometh",
                          positionOffset = 0)
)

#' @rdname checkForModification
#' 
#' @description 
#' \code{mod_ribometh}
#' 
#' @export
setMethod(
  f = "checkForModification",
  signature = signature(x = "mod_ribometh",
                        location = "numeric",
                        locations = "numeric",
                        data = "list"),
  definition = function(x,
                        location,
                        locations,
                        data) {
    # get test result for the current location
    locTest <- .calc_RiboMeth_test_values(location,
                                          locations,
                                          data)
    # If insufficient data or tow low scores is present return NULL
    if(is.null(locTest)) return(NULL)
    if(!.validate_RiboMeth_pos(locTest) ) {
      return(NULL)
    }
    score_c <- .aggregate_score_C(data,
                                  location,
                                  .get_ribometh_score_weights())
    # Return data
    return(list(location = location,
                type = getModType(x),
                score_a = locTest$score_a,
                score_a.sd = locTest$score_a.sd,
                score_b = locTest$score_b,
                score_b.sd = locTest$score_b.sd,
                score_c = score_c$mean,
                score_c.sd = score_c$sd,
                nbsamples = locTest$n))
  }
)

# test for RiboMeth at current location
.calc_RiboMeth_test_values <- function(location,
                                       locations,
                                       data){
  # split data
  data <- lapply(data,"[[","data")
  # short cut if amount of data is not sufficient
  pretestData <- .do_RiboMeth_pretest(location,
                                      locations,
                                      data)
  if(is.null(pretestData)) {
    return(NULL)
  }
  # Return data
  return(list(score_a = pretestData$score_a$mean,
              score_a.sd = pretestData$score_a$sd,
              score_b = pretestData$score_b$mean,
              score_b.sd = pretestData$score_b$sd,
              nbsamples = pretestData$n))
}

# check if any data is available to proceed with test
# this is in a seperate function since it is also called by checkForModification
.do_RiboMeth_pretest <- function(location,
                                 locations,
                                 data,
                                 coverage){
  # do not take into account position 1 or the other end
  if(location == 1 || location > (max(locations) - 10) ) return(NULL)
  # number of replicates
  n <- length(data)
  a <- .aggregate_score_A(data,
                          location)
  b <- .aggregate_score_B(data,
                          location,
                          .get_ribometh_score_weights())
  if(is.null(a) || is.null(b)) {
    return(NULL)
  }
  return(list(n = n,
              score_a = a,
              score_b = b))
}

# call yes or nor position
.validate_M7G_pos <- function(locData){
  (locData.score_a > RNAMODR_RIBOMETH_SCORE_A_THRESHOLD &&
     locData.score_b > RNAMODR_RIBOMETH_SCORE_B_THRESHOLD)
}

# RiboMeth scores --------------------------------------------------------------

# calculates score A according to Birkedal et al. 2014
.aggregate_score_A <- function(data,
                               pos){
  # subset to position
  locData <- .aggregate_location_data(data,
                                      pos)
  areaL <- .aggregate_area_data_left(data,
                                     pos,
                                     .get_ribometh_score_width())
  areaR <- .aggregate_area_data_left(data,
                                     pos,
                                     .get_ribometh_score_width())
  # check that data
  # all data points populated
  if( any(vapply(areaL, length, numeric(1)) < .get_ribometh_score_width()) ||
      any(vapply(areaR, length, numeric(1)) < .get_ribometh_score_width()) ) {
    return(NULL)
  }
  # average data meets minimum required
  if( any( (vapply(areaL, sum, numeric(1)) / .get_ribometh_score_width()) < RNAMODR_ALLENDS_COVERAGE_MIN ) ||
      any( (vapply(areaR, sum, numeric(1)) / .get_ribometh_score_width()) < RNAMODR_ALLENDS_COVERAGE_MIN ) ) {
    return(NULL)
  }
  # calculate means and sd
  meanLeft <- lapply(areaL, mean)
  meanRight <- lapply(areaR, mean)
  sdLeft <- lapply(areaL, mean)
  sdRight <- lapply(areaR, mean)
  # calc score per replicate
  scores <- mapply(FUN = .calculate_ribometh_score_A,
                   locData,
                   meanLeft,
                   sdLeft,
                   meanRight,
                   sdRight)
  return(scores = scores,
         mean = mean(scores),
         sd = sd(scores))
}

# calculates score B according to Birkedal et al. 2014
.aggregate_score_B <- function(data,
                               pos,
                               weights){
  # subset to position
  locData <- .aggregate_location_data(data,
                                      pos)
  areaL <- .aggregate_area_data_left(data,
                                     pos,
                                     .get_ribometh_score_width())
  areaR <- .aggregate_area_data_right(data,
                                      pos,
                                      .get_ribometh_score_width())
  # check that datas
  # all data points populated
  if( any(vapply(areaL, length, numeric(1)) < .get_ribometh_score_width()) ||
      any(vapply(areaR, length, numeric(1)) < .get_ribometh_score_width())) {
    return(NULL)
  }
  # average data meets minimum required
  if( any( (vapply(areaL, sum, numeric(1)) / .get_ribometh_score_width() ) < RNAMODR_ALLENDS_COVERAGE_MIN ) ||
      any( (vapply(areaR, sum, numeric(1)) / .get_ribometh_score_width() ) < RNAMODR_ALLENDS_COVERAGE_MIN )) {
    return(NULL)
  }
  # calc score per replicate
  scores <- mapply(FUN = .calculate_ribometh_score_B,
                   locData,
                   areaL,
                   areaR,
                   MoreArgs = list(weights))
  return(scores = scores,
         mean = mean(scores),
         sd = sd(scores))
}

# calculates score C according to Birkedal et al. 2014
.aggregate_score_C <- function(data,
                               pos,
                               weights){
  # subset to position
  locData <- .aggregate_location_data(data,
                                      pos)
  areaL <- .aggregate_area_data_left(data,
                                     pos,
                                     .get_ribometh_score_width())
  areaR <- .aggregate_area_data_right(data,
                                      pos,
                                      .get_ribometh_score_width())
  # calc score per replicate
  scores <- mapply(FUN = .calculate_ribometh_score_C,
                   locData,
                   areaL,
                   areaR,
                   MoreArgs = list(weights))
  return(scores = scores,
         mean = mean(scores),
         sd = sd(scores))
}

# calculates score MAX according to Marchand et al. 2016
# .aggregate_score_MAX <- function(data,
#                                  pos){
#   locData <- .aggregate_location_data(data,
#                                       pos)
#   areaData <- .aggregate_area_data2(data,
#                                     pos,
#                                     (.get_ribometh_score_width() + 1))
#   areaData <- .get_relative(areaData,
#                             pos,
#                             .get_ribometh_score_width())
# }
# 
# .get_relative <- function(areaData,
#                           pos,
#                           width){
#   lapply(areaData,
#          function(y){
#            a <- lapply(1:(length(y)-1), function(i){
#              y[i] / y[i+1]
#            })
#            a <- c(a,0)
#            names(a) <- names(y)
#            b <- lapply((length(y)-1):2, function(i){
#              y[i] / y[i-1]
#            })
#            b <- c(0,b)
#            names(b) <- names(y)
#            
#            df <- data.frame(a,
#                             b)
#            rownames(df) <- names(y)
#            df <- df[2:(nrow(df) - 1),]
#            avgA <- mean(df[rownames(df) != pos,1])
#            avgB <- mean(df[rownames(df) != pos,2])
#            df[,1] <- df[,1]/avgA
#            df[,2] <- df[,2]/avgB
#            
#            
#          })
# }
