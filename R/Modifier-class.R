#' @include RNAmodR.R
#' @include SequenceData-class.R
#' @include Modifier-utils.R
NULL

#' @name Modifier
#' @aliases modify
#' 
#' @title Modifier
#' @description 
#' title
NULL

#' @rdname Modifier
#' @export
setClass("Modifier",
         contains = c("VIRTUAL"),
         slots = c(mod = "character", # this have to be populated by subclass
                   score = "character", # this have to be populated by subclass
                   dataClass = "character", # this have to be populated by subclass
                   bamfiles = "BamFileList",
                   conditions = "factor",
                   fasta = "FaFile",
                   gff = "GFFFile",
                   data = "SequenceData",
                   aggregate = "DataFrameList",
                   modifications = "GRanges",
                   arguments = "list",
                   aggregateValidForCurrentArguments = "logical",
                   modificationsValidForCurrentArguments = "logical"),
         prototype = list(aggregateValidForCurrentArguments = FALSE,
                          modificationsValidForCurrentArguments = FALSE))

setMethod(
  f = "initialize", 
  signature = signature(.Object = "Modifier"),
  definition = function(.Object,
                        bamfiles = NULL,
                        fasta = NULL,
                        gff = NULL) {
    # check modification ident
    .Object@mod <- .norm_mod(.Object@mod,className)
    # short cut for creating an empty object
    if(is.null(bamfiles) || is.null(fasta) || is.null(gff)){
      return(.Object)
    }
    className <- class(.Object)[[1]]
    className <- .norm_modifiertype(className)
    # check bam files
    bamfiles <- .norm_bamfiles(bamfiles,className)
    # check genome sequences
    fasta <- .norm_fasta(fasta,className)
    # check genome annotation
    gff <- .norm_gff(gff,className)
    # set clots
    .Object@bamfiles <- bamfiles
    .Object@conditions <- factor(names(bamfiles))
    .Object@fasta <- fasta
    .Object@gff <- gff
    return(.Object)
  }
)

.valid_Modifier <- function(x){
  NULL
}
S4Vectors::setValidity2(Class = "Modifier",.valid_Modifier)

# show -------------------------------------------------------------------------

.show_settings <- function(settings){
  l <- length(settings)
  nc <- 6
  nr <- ceiling(l / nc)
  settings <- lapply(seq_len(nr),
                     function(i){
                       f <- seq.int(from = (i * nc) - nc + 1L,
                                    to = i * nc)
                       f <- f[f <= l]
                       settings[,f]
                     })
  if(is.null(rownames(settings[[1]]))){
    names <- rep(" ",nrow(settings[[1]]))
  } else {
    names <- rownames(settings[[1]])
  }
  for(i in seq_along(settings)){
    out <-
      as.matrix(format(as.data.frame(
        lapply(settings[[i]],showAsCell),
        optional = TRUE)))
    classinfo <-
      matrix(unlist(lapply(settings[[i]], function(x) {
        paste0("<", classNameForDisplay(x)[1],
               ">")
      }), use.names = FALSE), nrow = 1,
      dimnames = list("", colnames(out)))
    out <- rbind(classinfo, out)
    rownames(out) <- c(" ",names)
    print(out, quote = FALSE, right = TRUE)
  }
}

#' @rdname Modifier
#' @export
setMethod(
  f = "show", 
  signature = signature(object = "Modifier"),
  definition = function(object) {
    cat("A", class(object), "object containing",object@dataClass,
        "with",length(object@data),"elements.\n")
    files <- path(object@bamfiles)
    cat("| Input files:\n",paste0("  - ",names(files),": ",files,"\n"))
    cat("| Sequence file:",path(object@fasta),"\n")
    cat("| Annotation file:",path(object@gff),"\n")
    cat("| Modification type(s): ",paste0(object@mod, collapse = " / "),"\n")
    cat("| Modifications found:",ifelse(length(object@modifications) != 0L,
                                      paste0("yes (",
                                             length(object@modifications),
                                             ")"),
                                      "no"),"\n")
    cat("| Settings:\n")
    settings <- settings(object)
    settings <- lapply(settings,
                       function(s){
                         if(length(s) > 1L){
                           ans <- List(s)
                           return(ans)
                         }
                         s
                       })
    settings <- DataFrame(settings)
    .show_settings(settings)
    valid <- c(object@aggregateValidForCurrentArguments,
               object@modificationsValidForCurrentArguments)
    if(!all(valid)){
      warning("Settings were changed after data aggregation or modification ",
              "search. Rerun with modify(x,force = TRUE) to update with ",
              "current settings.", call. = FALSE)
    }
  }
)

# accessors --------------------------------------------------------------------

#' @name Modifier
#' @export
setMethod(f = "modifierType", 
          signature = signature(x = "Modifier"),
          definition = function(x){class(x)[[1]]})
#' @name Modifier
#' @export
setMethod(f = "modType", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@mod})
#' @name Modifier
#' @export
setMethod(f = "mainScore", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@score})

.norm_args <- function(input){
  minCoverage <- 10L
  minReplicate <- 1L
  findMod <- TRUE
  if(!is.null(input[["minCoverage"]])){
    minCoverage <- input[["minCoverage"]]
    if(!is.integer(minCoverage) || 
       minCoverage < 0L ||
       length(minCoverage) != 1){
      stop("'minCoverage' must be a single positive integer value.")
    }
  }
  if(!is.null(input[["minReplicate"]])){
    minReplicate <- input[["minReplicate"]]
    if(!is.integer(minReplicate) || 
       minReplicate < 0L ||
       length(minReplicate) != 1){
      stop("'minReplicate' must be a single positive integer value.")
    }
  }
  if(!is.null(input[["findMod"]])){
    findMod <- input[["findMod"]]
    if(!assertive::is_a_bool(findMod)){
      stop("'findMod' must be a single logical value.")
    }
  }
  args <- list(minCoverage = minCoverage,
               minReplicate = minReplicate,
               findMod = findMod)
  args
}

#' @name Modifier
#' @export
setMethod(f = "settings", 
          signature = signature(x = "Modifier"),
          definition = function(x,name){
            if(missing(name) || is.null(name)){
              return(x@arguments)
            }
            if(!assertive::is_a_string(name)){
              stop("'name' must be a single character value.")
            }
            x@arguments[[name]]
          }
)
#' @name Modifier
#' @export
setReplaceMethod(f = "settings", 
          signature = signature(x = "Modifier"),
          definition = function(x,value){
            if(is.null(names(value)) && length(value) > 0L){
              stop("'value' has to be a named.")
            }
            if(!is.list(value)){
              value <- as.list(value)
            }
            value <- .norm_args(value)
            x@arguments[names(value)] <- unname(value)
            x@aggregateValidForCurrentArguments <- FALSE
            x@modificationsValidForCurrentArguments <- FALSE
            x
          })

#' @name Modifier
#' @export
setMethod(f = "gff", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@gff})
  
#' @name Modifier
#' @export
setMethod(f = "fasta", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@fasta})
 
#' @name Modifier
#' @export
setMethod(f = "sequences", 
          signature = signature(x = "Modifier"),
          definition = 
            function(x,
                     modified = FALSE,
                     with.qualities = FALSE){
              if(!assertive::is_a_bool(modified)){
                stop("'modified' has to be a single logical value.")
              }
              if(!assertive::is_a_bool(with.qualities)){
                stop("'with.qualities' has to be a single logical value.")
              }
              if(modified == FALSE){
                return(x@data@sequences)
              }
              mod <- .get_modifications_per_transcript(x)
              mod <- split(mod,names(mod))
              ans <- ModRNAStringSet(sequences(seqData(x)))
              modSeqList <- ans[names(ans) %in% names(mod)]
              mod <- mod[match(names(mod),names(modSeqList))]
              ans[names(ans) %in% names(mod)] <- 
                combineIntoModstrings(modSeqList,
                                      mod)
              if(with.qualities == TRUE){
                browser()
              }
              ans
            }
)
  
#' @name Modifier
#' @export
setMethod(f = "ranges", 
          signature = signature(x = "Modifier"),
          definition = function(x){ranges(seqData(x))})

#' @name Modifier
#' @export
setMethod(f = "bamfiles", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@bamfiles})

#' @name Modifier
#' @export
setMethod(f = "seqData", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@data})

.check_score_name <- function(data,score){
  if(is(data,"CompressedSplitDataFrameList")){
    columns <- colnames(data@unlistData)
  } else {
    columns <- colnames(data[[1]])
  }
  if(!(score %in% columns)){
    stop("The default score is not present in the aggregate data. Contact the ",
         "maintainer of the class used.",
         call. = FALSE)
  }
  data
}

#' @name aggregate
#' @export
setMethod(f = "aggregateData", 
          signature = signature(x = "Modifier"),
          definition = function(x){
            x <- aggregate(x)
            .check_score_name(x@aggregate,x@score)
          })

# converts the genomic coordinates to transcript based coordinates
.get_modifications_per_transcript <- function(x){
  ranges <- .get_parent_annotations(ranges(x))
  modifications <- modifications(x)
  modRanges <- 
    ranges[as.character(ranges$ID) %in% as.character(modifications$Parent),]
  modRanges <- modRanges[match(as.character(modRanges$ID),
                               as.character(modifications$Parent))]
  # modify modifcation positions from genome centric to transcript centric
  start(modifications[strand(modifications) == "+"]) <- 
    start(modifications[strand(modifications) == "+"]) - 
    start(modRanges[strand(modRanges) == "+"]) + 1L
  end(modifications[strand(modifications) == "+"]) <- 
    end(modifications[strand(modifications) == "+"]) - 
    start(modRanges[strand(modRanges) == "+"]) + 1L
  end(modifications[strand(modifications) == "-"]) <- 
    end(modRanges[strand(modRanges) == "-"]) - 
    end(modifications[strand(modifications) == "-"]) + 1L
  start(modifications[strand(modifications) == "-"]) <- 
    end(modRanges[strand(modRanges) == "-"]) - 
    start(modifications[strand(modifications) == "-"]) + 1L
  names(modifications) <- as.character(modifications$Parent)
  modifications
}

#' @rdname modify
#' @export
setMethod(f = "modifications", 
          signature = signature(x = "Modifier"),
          definition = 
            function(x,
                     perTranscript = FALSE){
              if(!assertive::is_a_bool(perTranscript)){
                stop("'perTranscript' has to be a single logical value.")
              }
              if(perTranscript){
                return(.get_modifications_per_transcript(x))
              }
              x@modifications
            }
)

# constructors -----------------------------------------------------------------

.ModFromCharacter <- function(class,
                              bamfiles,
                              fasta,
                              gff,
                              args){
  ans <- new(class,
             bamfiles,
             fasta,
             gff)
  settings(ans) <- args
  modName <- fullName(ModRNAString())[
    which(shortName(ModRNAString()) %in% ans@mod)]
  message("Starting to search for '",
          paste(tools::toTitleCase(modName), collapse = "', '"),
          "'...")
  ans@data <- do.call(ans@dataClass,
                      c(list(bamfiles = ans@bamfiles,
                           fasta = ans@fasta,
                           gff = ans@gff),
                        settings(ans)))
  ans@data@sequences <- RNAStringSet(ans@data@sequences)
  if(settings(ans,"findMod")){
    ans <- do.call(modify,
                   list(ans))
  }
  ans
}

.ModFromSequenceData <- function(class,
                                 x,
                                 args){
  ans <- new(class,
             x@bamfiles,
             x@fasta,
             x@gff)
  settings(ans) <- args
  # check data type
  ans@data <- .norm_data_type(ans,x)
  if(settings(ans,"findMod")){
    ans <- do.call(modify,list(ans))
  }
  ans
}

# dummy functions --------------------------------------------------------------
# these need to be implemented by each subclass

#' @name modify
#' 
#' @title modify
#' 
#' @description 
#' title
#' 
#' @export
setMethod(f = "modify", 
          signature = signature(x = "Modifier"),
          definition = 
            function(x){
              stop("This functions needs to be implemented by '",class(x),"'.",
                   call. = FALSE)
            }
)

#' @name aggregate
#' @aliases hasAggregateData aggregateData
#' 
#' @title aggregate
#' 
#' @description 
#' title
#' 
#' @export
setMethod(f = "aggregate", 
          signature = signature(x = "Modifier"),
          definition = 
            function(x){
              stop("This functions needs to be implemented by '",class(x),"'.",
                   call. = FALSE)
            }
)

# check function ---------------------------------------------------------------

#' @rdname aggregate
#' @export
setMethod(f = "hasAggregateData", 
          signature = signature(x = "Modifier"),
          definition = 
            function(x){
              if(is.null(nrow(x@aggregate))){
                return(FALSE)
              }
              return(TRUE)
            }
)
