#' @include RNAmodR.R
#' @include Modifier-class.R
# #' @include RiboMethSeq.R
NULL

#' @name ModRiboMethSeq
#' @aliases RiboMethSeq ModRiboMethSeq
#' 
#' @title ModRiboMethSeq
#' @description 
#' title
#' 

NULL

#' @rdname ModRiboMethSeq
#' @export
setClass("ModRiboMethSeq",
         contains = c("Modifier"),
         prototype = list(mod = c("Am","Cm","Gm","Um"),
                          dataClass = "EndSequenceData"))


setMethod(
  f = "initialize", 
  signature = signature(.Object = "ModRiboMethSeq"),
  definition = function(.Object,
                        bamfiles,
                        fasta,
                        gff) {
    .Object <- callNextMethod(.Object,
                              bamfiles,
                              fasta,
                              gff)
    return(.Object)
  }
)

# constructors -----------------------------------------------------------------

.norm_rms_args <- function(input){
  args <- RNAmodR:::.norm_args(input)
  args
}


setGeneric( 
  name = "ModRiboMethSeq",
  def = function(x,
                 ...) standardGeneric("ModRiboMethSeq")
)
# Create Modifier class from file character, fasta and gff file
#' @rdname ModRiboMethSeq
#' @export
setMethod("ModRiboMethSeq",
          signature = c(x = "character"),
          function(x,
                   fasta,
                   gff,
                   modifications = NULL,
                   ...){
            args <- .norm_rms_args(list(...))
            ans <- RNAmodR:::.ModFromCharacter("ModRiboMethSeq",
                                               x,
                                               fasta,
                                               gff,
                                               args)
            ans <- RNAmodR:::.norm_modifications(ans,
                                                 args)
            ans
          }
)

# Create Modifier class from bamfiles, fasta and gff file
#' @rdname ModRiboMethSeq
#' @export
setMethod("ModRiboMethSeq",
          signature = c(x = "BamFileList"),
          function(x,
                   fasta,
                   gff,
                   modifications = NULL,
                   ...){
            args <- .norm_rms_args(list(...))
            ans <- RNAmodR:::.ModFromCharacter("ModRiboMethSeq",
                                               x,
                                               fasta,
                                               gff,
                                               args)
            ans <- RNAmodR:::.norm_modifications(ans,
                                                 args)
            ans
          }
)

# Create Modifier class from existing SequenceData
#' @rdname ModRiboMethSeq
#' @export
setMethod("ModRiboMethSeq",
          signature = c(x = "SequenceData"),
          function(x,
                   modifications = NULL,
                   ...){
            args <- .norm_rms_args(list(...))
            ans <- RNAmodR:::.ModFromSequenceData("ModRiboMethSeq",
                                                  x,
                                                  args)
            ans <- RNAmodR:::.norm_modifications(ans,
                                                 args)
            ans
          }
)

# functions --------------------------------------------------------------------

.get_rms_scores <- function(data){
  
}

# calculates score A according to Birkedal et al. 2015
.calculate_ribometh_score_A_c <- function(n,
                                          meanL,
                                          sdL,
                                          meanR,
                                          sdR){
  if(is.nan(meanL) || is.na(sdL) || is.na(meanR) || is.na(sdR)){
    return(NA)
  }
  dividend <- (2 * n  + 1)
  divisor <- (0.5 * abs(meanL - sdL) + n + 
                0.5 * abs(meanR - sdR) + 1 ) 
  # return(max(0, (1 - (dividend / divisor))))
  return((1 - (dividend / divisor)))
}
.calculate_ribometh_score_A <- compiler::cmpfun(.calculate_ribometh_score_A_c)

.calculate_ribometh_score_B_c <- function(n,
                                          areaL,
                                          areaR,
                                          weights){
  # split the weights
  weightsL <- weights[names(weights) < 0]
  weightsR <- weights[names(weights) > 0]
  if(length(weightsL) != length(areaL) || length(weightsR) != length(areaR)){
    return(NA)
  }
  # calculates score
  dividend <- abs(n - 0.5 * (stats::weighted.mean(areaL, weightsL) +
                               stats::weighted.mean(areaR, weightsR)))
  divisor <- (n + 1)
  return(dividend / divisor)
}
.calculate_ribometh_score_B <- compiler::cmpfun(.calculate_ribometh_score_B_c)

# calculates score C according to Birkedal et al. 2014
.calculate_ribometh_score_meth_c <- function(n,
                                             areaL,
                                             areaR,
                                             weights){
  # split the weights
  weightsL <- weights[names(weights) < 0]
  weightsR <- weights[names(weights) > 0]
  if(length(weightsL) != length(areaL) || length(weightsR) != length(areaR)){
    return(NA)
  }
  # calculates score
  return((1 - (n / 
                 ( 0.5 * (stats::weighted.mean(areaL, weightsL) +
                            stats::weighted.mean(areaR, weightsR))
                 )
  )
  )
  )
}
.calculate_ribometh_score_meth <- compiler::cmpfun(.calculate_ribometh_score_meth_c)

# calculates score mean according to Marchand et al. 2016
.calculate_ribometh_score_mean_c <- function(locationData5,
                                             locationData3,
                                             mean5,
                                             mean3,
                                             sd5,
                                             sd3){
  if(length(locationData5) == 0 || length(locationData3) == 0){
    return(NA)
  }
  value5 <- (locationData5 - mean5)/sd5
  value3 <- (locationData3 - mean3)/sd3
  valueMean <- mean(c(value5,value3))
  return(max(abs(value5 - valueMean),abs(value3 - valueMean)))
}
.calculate_ribometh_score_mean <- compiler::cmpfun(.calculate_ribometh_score_mean_c)

# calculates score A according to Birkedal et al. 2015
.aggregate_score_A <- function(data,
                               locations,
                               weights){
  min <- min(as.numeric(names(weights)))
  max <- max(as.numeric(names(weights)))
  counts <- data$counts
  # names(counts) <- pos(data)
  names(counts) <- data$position
  # get data
  locationData <- counts[as.numeric(names(counts)) %in% locations]
  area5 <- lapply(locations, function(l){counts[as.numeric(names(counts)) < l & as.numeric(names(counts)) >= (l + min)]})
  area3 <- lapply(locations, function(l){counts[as.numeric(names(counts)) <= (l + max) & as.numeric(names(counts)) > l]})
  # calculate means and sd
  mean5 <- lapply(area5, mean)
  mean3 <- lapply(area3, mean)
  sd5 <- lapply(area5, stats::sd)
  sd3 <- lapply(area3, stats::sd)
  # calc score per replicate
  scores <- mapply(FUN = .calculate_ribometh_score_A,
                   locationData,
                   mean5,
                   sd5,
                   mean3,
                   sd3)
  return(scores)
}

# calculates score B according to Birkedal et al. 2015
.aggregate_score_B <- function(data,
                               locations,
                               weights){
  min <- min(as.numeric(names(weights)))
  max <- max(as.numeric(names(weights)))
  counts <- data$counts
  # names(counts) <- pos(data)
  names(counts) <- data$position
  # get data
  locationData <- counts[as.numeric(names(counts)) %in% locations]
  area5 <- lapply(locations, function(l){counts[as.numeric(names(counts)) < l & as.numeric(names(counts)) >= (l + min)]})
  area3 <- lapply(locations, function(l){counts[as.numeric(names(counts)) <= (l + max) & as.numeric(names(counts)) > l]})
  # calc score per replicate
  scores <- mapply(FUN = .calculate_ribometh_score_B,
                   locationData,
                   area5,
                   area3,
                   MoreArgs = list(weights))
  return(scores)
}

# calculates score C according to Birkedal et al. 2015
.aggregate_score_meth <- function(data,
                                  locations,
                                  weights){
  min <- min(as.numeric(names(weights)))
  max <- max(as.numeric(names(weights)))
  counts <- data$counts
  # names(counts) <- pos(data)
  names(counts) <- data$position
  # get data
  locationData <- counts[as.numeric(names(counts)) %in% locations]
  area5 <- lapply(locations, function(l){counts[as.numeric(names(counts)) < l & as.numeric(names(counts)) >= (l + min)]})
  area3 <- lapply(locations, function(l){counts[as.numeric(names(counts)) <= (l + max) & as.numeric(names(counts)) > l]})
  # calc score per replicate
  scores <- mapply(FUN = .calculate_ribometh_score_meth,
                   locationData,
                   area5,
                   area3,
                   MoreArgs = list(weights))
  return(scores)
}

# calculates score mean according to Marchand et al. 2016
.aggregate_score_mean <- function(data,
                                  locations,
                                  windowSize){
  counts <- data$counts
  # names(counts) <- pos(data)
  names(counts) <- data$position
  #
  locationData5 <- 
    lapply(locations, 
           function(l){
             pos5 <- counts[as.numeric(names(counts)) == l]
             pos <- counts[as.numeric(names(counts)) == (l + 1)]
             return(pos/pos5)
           })
  locationData3 <- 
    lapply(locations, 
           function(l){
             pos3 <- counts[as.numeric(names(counts)) == l]
             pos <- counts[as.numeric(names(counts)) == (l - 1)]
             return(pos/pos3)
           })
  area5 <- 
    lapply(locations, 
           function(l){
             locs <- (l - windowSize):(l + windowSize)
             locs <- locs[locs > 3 & locs != l & locs < max(locations)]
             pos5 <- counts[as.numeric(names(counts)) %in% (locs - 2)]
             pos <- counts[as.numeric(names(counts)) %in% (locs - 1)]
             return(pos/pos5)
           })
  area3 <- 
    lapply(locations, 
           function(l){
             locs <- (l - windowSize):(l + windowSize)
             locs <- locs[locs > 1 & locs != l & locs < max(locations)]
             pos3 <- counts[as.numeric(names(counts)) %in% locs]
             pos <- counts[as.numeric(names(counts)) %in% (locs - 1)]
             return(pos/pos3)
           })
  mean5 <- lapply(area5, mean)
  mean3 <- lapply(area3, mean)
  sd5 <- lapply(area5, stats::sd)
  sd3 <- lapply(area3, stats::sd)
  # calc score per replicate
  scores <- mapply(FUN = .calculate_ribometh_score_mean,
                   locationData5,
                   locationData3,
                   mean5,
                   mean3,
                   sd5,
                   sd3)
  return(scores)
}

.find_rms <- function(x,
                      args){
  message("Searching for 2'-O methylations...")
  # parameter data
  weights <- c(0.5,0.6,0.7,0.8,0.9,1,0,1,0.9,0.8,0.7,0.6,0.5)
  weightPositions <- c(-6L,-5L,-4L,-3L,-2L,-1L,0L,1L,2L,3L,4L,5L,6L)
  # ToOo check for continuity
  if(length(weights) != length(weightPositions)){
    stop("Something went wrong.")
  }
  minWeightPos <- min(weightPositions)
  maxWeightPos <- max(weightPositions)
  #
  browser()
  data <- seqdata(x)
  mod <- aggregate(data)
  means <- SplitDataFrameList(mod@unlistData[grepl("mean",colnames(mod@unlistData))])
  means@partitioning <- mod@partitioning
  ncol <- ncol(means[[1]])
  ncolV <- seq_len(ncol)
  n <- length(mod)
  nV <- seq_len(n)
  seq <- weightPositions - minWeightPos
  pos <- lapply(lengths(mod),seq_len)
  # subset to neightbouring positions
  neighborCounts <- lapply(ncolV,
                           function(i){
                             z <- IntegerList(mod@unlistData[,i])
                             z@partitioning <- mod@partitioning 
                             lapply(nV,
                                    function(j){
                                      IntegerList(lapply(pos[[j]],
                                             function(k){
                                               f <- weightPositions + k - 1L
                                               f <- f[f >= 0]
                                               z[[j]][f]
                                             }))
                                    })
                           })
  # calculate mean for neighbouring position
  neighborAmean <- lapply(ncolV,
                          function(i){
                            z <- neighborCounts[[i]]
                            IntegerList(lapply(nV,
                                               function(j){
                                                 unlist(lapply(z[[j]],mean,na.rm = TRUE))
                                               }))
                          })
  # calculate sd for neighbouring positions
  neighborAsd <- lapply(ncolV,
                        function(i){
                          z <- neighborCounts[[i]]
                          IntegerList(lapply(nV,
                                             function(j){
                                               unlist(lapply(z[[j]],sd,na.rm = TRUE))
                                             }))
                        })
  # calculate weighted mean for neighbouring position
  neighborBmean <- lapply(ncolV,
                      function(i){
                        z <- neighborCounts[[i]]
                        IntegerList(lapply(nV,
                                           function(j){
                                             unlist(lapply(lapply(z[[j]],"*",weights),mean,na.rm = TRUE))
                                           }))
                      })
  
  neighborRMS <- NULL
  neighborMAX <- NULL
  
  
  
  
  # Calculate ScoreMax for each position, make a ratio of number
  # of 5′-reads ends between preceding and following position and
  # calculate ScoreMAX as a ratio of a drop for a given position
  # compared to the average for 12 neighboring positions (−6/+6).
  # 6. Calculate the RiboMethScore, use the following formula:
  #   RiboMethScore = 1 – ni/(0.5*(SUM(nj*Wj)/SUM(Wj)+
#   SUM(nk*Wk)/SUM(Wk)), where ni – 5′-end count for a given
# position, j—varies from i−6 to i−1, k varies from i+1 to i+6,
# weight parameters are defined as following Wj = 1–0.1*(j−1)
# idem for Wk. Thus, Wj or Wk parameter varies from 0.5 for
# positions −6/+6 to 1 for positions −1/+1.
  
  
  NULL
}

#' @rdname ModRiboMethSeq
#' @export
setMethod("modify",
          signature = c(x = "ModRiboMethSeq"),
          function(x,
                   ...){
            x@modifications <- .find_rms(x,
                                         .norm_rms_args(list(...)))
            message("done.")
            x
          }
)
