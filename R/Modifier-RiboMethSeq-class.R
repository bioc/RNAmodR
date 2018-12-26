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
#' score MAX as described by publiccations from the motorin lab are not 
#' implemented since an unambigeous description is not available from the 
#' literature.
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
  minSignal <- 50L
  minScoreA <- 0.6
  minScoreB <- 1
  minScoreRMS <- 0.65
  scoreOperator <- "&"
  if(!is.null(input[["minSignal"]])){
    minSignal <- input[["minSignal"]]
    if(!is.integer(minSignal) | minSignal < 0){
      stop("'minSignal' must be integer with a value higher than 0.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["minScoreA"]])){
    minScoreA <- input[["minScoreA"]]
    if(!is.numeric(minScoreA) | minScoreA < 0 | minScoreA > 1){
      stop("'minScoreA' must be numeric with a value between 0 and 1.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["minScoreB"]])){
    minScoreB <- input[["minScoreB"]]
    if(!is.numeric(minScoreB) | minScoreB < 0){
      stop("'minScoreB' must be numeric and greater then 0.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["minScoreRMS"]])){
    minScoreRMS <- input[["minScoreRMS"]]
    if(!is.numeric(minScoreRMS) | minScoreRMS < 0 | minScoreRMS > 1){
      stop("'minScoreRMS' must be numeric with a value between 0 and 1.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["scoreOperator"]])){
    scoreOperator <- input[["scoreOperator"]]
    if(!(scoreOperator %in% c("|","&"))){
      stop("'scoreOperator' must be either '|' or '&'.",
           call. = FALSE)
    }
  }
  args <- .norm_args(input)
  args <- c(args,
            list(minSignal = minSignal,
                 minScoreA = minScoreA,
                 minScoreB = minScoreB,
                 minScoreRMS = minScoreRMS,
                 scoreOperator = scoreOperator))
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


# calculates score A according to Birkedal et al. 2015
# it is simplified to use the mean/sd for the neighboring positions
# and not the left and right mean/sd seperatly
.calculate_ribometh_score_A_c <- function(n,
                                          mean,
                                          sd){
  dividend <- (2 * n  + 1)
  divisor <- (abs(mean - sd) + n + 1)
  ans <- 1 - (dividend / divisor)
  return(max(0, ans))
}
.calculate_ribometh_score_A <- compiler::cmpfun(.calculate_ribometh_score_A_c)

# calculates score B according to Birkedal et al. 2015
# it is simplified to use the weighted neighboring positions
# and not the left and right area seperatly
.calculate_ribometh_score_B_c <- function(n,
                                          area,
                                          weights){
  dividend <- abs(n - (sum(area) / sum(weights) ) )
  divisor <- (n + 1)
  return(dividend / divisor)
}
.calculate_ribometh_score_B <- compiler::cmpfun(.calculate_ribometh_score_B_c)

# calculates score C according to Birkedal et al. 2015
# it is simplified to use the weighted neighboring positions
# and not the left and right area seperatly
.calculate_ribometh_score_meth_c <- function(n,
                                             area,
                                             weights){
  dividend <- n
  divisor <- sum(area) / sum(weights)
  ans <- 1 - (dividend / divisor)
  return(max(0, ans))
}
.calculate_ribometh_score_meth <- compiler::cmpfun(.calculate_ribometh_score_meth_c)

# calculates score MAX according to Marchand et al. 2016
# not used since no clear description available in the literature
# .calculate_ribometh_score_max_c <- function(n,
#                                              area,
#                                              weights){
# }
# .calculate_ribometh_score_max <- compiler::cmpfun(.calculate_ribometh_score_max_c)


#' @name Modifier
#' @export
setMethod(
  f = "aggregate", 
  signature = signature(x = "ModRiboMethSeq"),
  definition = 
    function(x,
             args){
      if(!hasAggregateData(x)){
        message("Aggregating data and calculating scores...")
        # parameter data
        weights <- c(0.5,0.6,0.7,0.8,0.9,1,0,1,0.9,0.8,0.7,0.6,0.5)
        weightPositions <- c(-6L,-5L,-4L,-3L,-2L,-1L,0L,1L,2L,3L,4L,5L,6L)
        # ToOo check for continuity
        if(length(weights) != length(weightPositions)){
          stop("Something went wrong.")
        }
        # get the means. the sds arecurrently disregarded for this analysis
        mod <- aggregate(seqData(x), condition = "Treated")
        means <- IntegerList(mod@unlistData[,which(grepl("mean",
                                                   colnames(mod@unlistData)))])
        means@partitioning <- mod@partitioning
        # browser()
        # set up variables
        n <- length(mod)
        nV <- seq_len(n)
        lengths <- lengths(mod)
        pos <- lapply(lengths,seq_len)
        # subset to neightbouring positions and set position to zero
        neighborCounts <- lapply(nV,
                                function(j){
                                  IntegerList(lapply(pos[[j]],
                                                     function(k){
                                                       f <- weightPositions + k
                                                       ans <- means[[j]][f[f > 0]]
                                                       ans[weightPositions[f > 0] == 0] <- 0
                                                       ans <- ans[!is.na(ans)]
                                                       ans
                                                     }))
                                })
        # create list of weights vector alongside the neighbor counts
        weightsList <- lapply(nV,
                              function(j){
                                NumericList(mapply(
                                  function(k,l){
                                    f <- weightPositions + k
                                    ans <- weights[f > 0 & f <= l]
                                    ans
                                  },
                                  pos[[j]],
                                  MoreArgs = list(l = lengths[j])))
                              })
        # remove zero from neighborCounts
        neighborCountsWoZero <- lapply(nV,
                                 function(j){
                                   IntegerList(lapply(neighborCounts[[j]],
                                                      function(k){
                                                        k[k > 0]
                                                      }))
                                 })
        # calculate mean and sd for neighbouring position
        neighborMean <- IntegerList(
          lapply(nV,
                 function(j){
                   unlist(lapply(neighborCountsWoZero[[j]],
                                 mean,
                                 na.rm = TRUE))
                 }))
        neighborSd <- IntegerList(
          lapply(nV,
                 function(j){
                   unlist(lapply(neighborCountsWoZero[[j]],
                                 sd,
                                 na.rm = TRUE))
                 }))
        # calculate weighted mean and sd for neighbouring position
        # neighborWeightedMean <- IntegerList(
        #   lapply(nV,
        #          function(j){
        #            unlist(lapply(neighborCounts[[j]] * weights[[j]],
        #                          mean,
        #                          na.rm = TRUE))
        #          }))
        # neighborWeightedSd <- IntegerList(
        #   lapply(nV,
        #          function(j){
        #            unlist(lapply(neighborCounts[[j]] * weights[[j]],
        #                          sd,
        #                          na.rm = TRUE))
        #          }))
        # apply the weights to the counts
        neighborWeightedArea <- lapply(nV,
                                       function(j){
                                         neighborCounts[[j]] * weightsList[[j]]
                                       })
        # calculate the actual scores
        scoreA <- NumericList(mapply(
          function(n,mean,sd){
            unlist(mapply(.calculate_ribometh_score_A,
                          n,
                          mean,
                          sd,
                          SIMPLIFY = FALSE))
          },
          means,
          neighborMean,
          neighborSd))
        scoreB <- NumericList(mapply(
          function(n,area){
            unlist(mapply(.calculate_ribometh_score_B,
                          n,
                          area,
                          MoreArgs = list(weights = weights),
                          SIMPLIFY = FALSE))
          },
          means,
          neighborWeightedArea))
        scoreRMS <- NumericList(mapply(
          function(n,area){
            unlist(mapply(.calculate_ribometh_score_meth,
                          n,
                          area,
                          MoreArgs = list(weights = weights),
                          SIMPLIFY = FALSE))
          },
          means,
          neighborWeightedArea))
        #
        # scoreMax <- NumericList(mapply(
        #   function(n,area){
        #     unlist(mapply(.calculate_ribometh_score_max,
        #                   n,
        #                   area,
        #                   MoreArgs = list(weights = weights),
        #                   SIMPLIFY = FALSE))
        #   },
        #   means,
        #   neighborWeightedArea))
        # construct the result DataFrameList
        ans <- DataFrame(ends = unlist(means),
                         scoreA = unlist(scoreA),
                         scoreB = unlist(scoreB),
                         scoreRMS = unlist(scoreRMS))
        ans <- SplitDataFrameList(ans)
        ans@partitioning <- mod@partitioning
        x@aggregate <- ans
      }
      x
    }
)

.get_rms_scores <- function(data){
  list(score = data$scoreRMS,
       scoreA = data$scoreA,
       scoreB = data$scoreB)
}

.find_rms <- function(x,
                      args){
  message("Searching for 2'-O methylations...")
  #
  data <- seqData(x)
  letters <- CharacterList(strsplit(as.character(sequences(data)),""))
  ranges <- split(.get_parent_annotations(ranges(data)),
                  seq_along(ranges(data)))
  # get the aggregate data
  mod <- aggregateData(x)
  # setup args
  minSignal <- args[["minSignal"]]
  minScoreA <- args[["minScoreA"]]
  minScoreB <- args[["minScoreB"]]
  minScoreRMS <- args[["minScoreRMS"]]
  scoreOperator <- args[["scoreOperator"]]
  # find modifications
  modifications <- mapply(
    function(m,l,r){
      rownames(m) <- seq_len(width(r))
      m <- m[!is.na(m$scoreA) &
               !is.na(m$scoreB) &
               !is.na(m$scoreRMS),]
      if(nrow(m) == 0L) return(NULL)
      m <- m[mapply(Reduce,
                    rep(scoreOperator,nrow(m)),
                    m$scoreA >= minScoreA,
                    m$scoreB >= minScoreB,
                    m$scoreRMS >= minScoreRMS),]
      if(nrow(m) == 0L) return(NULL)
      ans <- .construct_mod_ranges(r,m,modType = "Am",.get_rms_scores,
                            "RNAmodR","RNAMOD")
      ans$mod <- paste0(l[start(ans)],"m")
      ans
    },
    mod,
    letters,
    ranges)
  modifications <- GRangesList(modifications[!vapply(modifications,
                                                     is.null,
                                                     logical(1))])
  unname(unlist(modifications))
}

#' @rdname ModRiboMethSeq
#' @export
setMethod("modify",
          signature = c(x = "ModRiboMethSeq"),
          function(x,
                   ...){
            args <- .norm_rms_args(list(...))
            # get the aggregate data
            x <- aggregate(x, args)
            x@modifications <- .find_rms(x, args)
            message("done.")
            x
          }
)
