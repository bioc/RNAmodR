#' @include RNAmodR.R
#' @include SequenceData-class.R
#' @include SequenceDataList-class.R
#' @include Modifier-utils.R
NULL

#' @name Modifier-class
#' @aliases Modifier
#' 
#' @title The Modifier class
#' 
#' @description
#' The \code{Modifier} class is a virtual class, which provides the central 
#' functionality for searching for patterns of post-transcriptional RNA 
#' modifications in high throughput sequencing data.
#' 
#' Each subclass has to implement the following functions:
#' 
#' \itemize{
#' \item{\code{\link{aggregate}}: }{used for specific data aggregation}
#' \item{\code{\link{modify}}: }{used for specific search for modifications}
#' }
#' 
#' Optionally the function \code{\link[=Modifier-functions]{settings<-}} can be
#' implemented to store additional arguments, which the base class does not
#' recognize.
#' 
#' \code{Modifier} objects are constructed centrally by calling 
#' \code{Modifier()} with a \code{className} matching the specific class to be
#' constructed. Thei will trigger the immediate analysis, if \code{findMod} is
#' not set to \code{TRUE}.
#' 
#' 
#' @param className The name of the class which should be constructed.
#' @param x the input which can be of the following types
#' \itemize{
#' \item{\code{SequenceData}:} {a single \code{SequenceData} or a list containg
#' only \code{SequenceData} objects. The input will just be used as elements of
#' the \code{Modifier} and must match the requirements of specific
#' \code{Modifier} class }
#' \item{\code{BamFileList}:} {a named \code{BamFileList}}
#' \item{\code{character}:} {a \code{character} vector, which must be coercible
#' to a named \code{BamFileList} referencing existing bam files. Valid names are
#' \code{control} and \code{treated} to define conditions and replicates}
#' }
#' @param annotation annotation data, which must match the information contained
#' in the BAM files. This is parameter is only required if \code{x} if not a 
#' \code{Modifier} object.
#' @param sequences sequences matching the target sequences the reads were 
#' mapped onto. This must match the information contained in the BAM files. This
#' is parameter is only required if \code{x} if not a \code{Modifier} object.
#' @param seqinfo optional \code{\link[GenomeInfoDb:Seqinfo]{Seqinfo}} to 
#' subset the transcripts analyzed on a chromosome basis.
#' @param ... Additional otpional parameters:
#' \itemize{
#' \item{findMod: }{\code{TRUE} or \code{FALSE}: should the search for for 
#' modifications be triggered upon construction? If not the search can be 
#' started by calling the \code{modify()} function.}
#' }
#' All other arguments will be passed onto the \code{SequenceData} objects.
#' 
#' @slot mod a \code{character} value, which needs to contain one or more 
#' elements from the alphabet of a 
#' \code{\link[Modstrings:ModRNAString]{ModRNAString}} class.
#' @slot score the main score identifier used for visualizations
#' @slot dataType the class name(s) of the \code{SequenceData} class used 
#' @slot bamfiles the input bam files as \code{BamFileList}
#' @slot condition conditions along the \code{BamFileList}: Either 
#' \code{control} or \code{treated}
#' @slot replicate replicate number along the \code{BamFileList} for each of the
#' condition types.
#' @slot data The sequence data object: Either a \code{SequenceData} or a 
#' \code{SequenceDataList} object, if more than one \code{dataType} is used.
#' @slot aggregate the aggregated data as a \code{SplitDataFrameList}
#' @slot modifications the found modifications as a \code{GRanges} object
#' @slot arguments arguments used for the analysis as a \code{list}
#' @slot aggregateValidForCurrentArguments \code{TRUE} or \code{FALSE} whether
#' the aggregate data was constructed with the current arguments
#' @slot modificationsValidForCurrentArguments \code{TRUE} or \code{FALSE} 
#' whether the modifications were found with the current arguments
NULL

#' @name Modifier-functions
#' 
#' @title Modifier/ModifierSet functions
#' 
#' @description
#' For the \code{Modifier} and  \code{ModifierSet} classes a number of functions
#' are implemented to access the data stored by the object.
#' 
#' @param x,object a \code{Modifier} or \code{ModifierSet} class
#' @param name For \code{settings}: name of the setting to be returned or set
#' @param value For \code{settings}: value of the setting to be set
#' @param modified For \code{sequences}: \code{TRUE} or \code{FALSE}: Should
#' the sequences be returned as a \code{ModRNAString} with the found 
#' modifications added on top of the \code{RNAString}? See 
#' \code{\link[Modstrings:separate]{combineIntoModstrings}}.
#' @param ... Additional arguments. 
NULL

setClassUnion("SequenceData_OR_SequenceDataList",
              c("SequenceData", "SequenceDataList"))

#' @rdname Modifier-class
#' @export
setClass("Modifier",
         contains = c("VIRTUAL"),
         slots = c(mod = "character", # this have to be populated by subclass
                   score = "character", # this have to be populated by subclass
                   dataType = "character", # this have to be populated by subclass
                   bamfiles = "BamFileList",
                   condition = "factor",
                   replicate = "factor",
                   data = "SequenceData_OR_SequenceDataList",
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
  definition = function(.Object, bamfiles = NULL) {
    # check modification ident
    .Object@mod <- .norm_mod(.Object@mod,className)
    # short cut for creating an empty object
    if(is.null(bamfiles)){
      return(.Object)
    }
    className <- class(.Object)[[1]]
    className <- .norm_modifiertype(className)
    # check bam files
    bamfiles <- .norm_bamfiles(bamfiles,className)
    # set clots
    .Object@bamfiles <- bamfiles
    .Object@condition <- factor(names(bamfiles))
    .Object@replicate <- .get_replicate_number(bamfiles, .Object@condition)
    return(.Object)
  }
)

# validity ---------------------------------------------------------------------

.norm_SequenceData_elements <- function(x,list){
  if(is(list,"SequenceData")){
    list <- list(list)
  } else {
    if(is(list,"list")){
      browser()
      elementTypeMatch <- !vapply(list,is,logical(1),"SequenceData")
      if(any(elementTypeMatch)){
        stop("Not all elements are 'SequenceData' objects.", call. = FALSE)
      }
    }
  }
  elementTypes <- vapply(list,class,character(1))
  if(length(elementTypes) != length(x@dataType)){
    stop("Number of 'SequenceData' elements does not match the requirements of",
         " ",class(x),". '",paste(x@dataType, collapse = "','"),"' are ",
         "required", call. = FALSE)
  }
  elementTypes <- elementTypes[match(elementTypes,x@dataType)]
  if(any(elementTypes != x@dataType)){
    stop("Type of SequenceData elements does not match the requirements of ",
         class(x),". '",paste(x@dataType, collapse = "','"),"' are ",
         "required", call. = FALSE)
  }
}

.valid_SequenceData <- function(x){
  tmp <- try(.norm_SequenceData_elements(x,x@data))
  if (inherits(tmp, "try-error")){
    return(tmp)
  }
  NULL
}

.valid_Modifier <- function(x){
  c(.valid_SequenceData(x))
}
S4Vectors::setValidity2(Class = "Modifier",.valid_Modifier)

# show -------------------------------------------------------------------------

.show_settings <- function(settings){
  l <- length(settings)
  if(l == 0L){
    return(NULL)
  }
  nc <- 6
  nr <- ceiling(l / nc)
  settings <- lapply(seq_len(nr),
                     function(i){
                       f <- seq.int(from = (i * nc) - nc + 1L,
                                    to = i * nc)
                       f <- f[f <= l]
                       settings[,f,drop = FALSE]
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

#' @rdname Modifier-functions
setMethod(
  f = "show", 
  signature = signature(object = "Modifier"),
  definition = function(object) {
    cat("A", class(object), "object containing",object@dataType,
        "with",length(object@data),"elements.\n")
    files <- path(object@bamfiles)
    cat("| Input files:\n",paste0("  - ",names(files),": ",files,"\n"))
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

#' @rdname Modifier-functions
#' @export
setMethod(f = "modifierType", 
          signature = signature(x = "Modifier"),
          definition = function(x){class(x)[[1]]})
#' @rdname Modifier-functions
#' @export
setMethod(f = "modType", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@mod})
#' @rdname Modifier-functions
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

#' @rdname Modifier-functions
#' @export
setMethod(f = "settings", 
          signature = signature(x = "Modifier"),
          definition = function(x, name){
            if(missing(name) || is.null(name)){
              return(x@arguments)
            }
            if(!assertive::is_a_string(name)){
              stop("'name' must be a single character value.")
            }
            x@arguments[[name]]
          }
)
#' @rdname Modifier-functions
#' @export
setReplaceMethod(f = "settings", 
          signature = signature(x = "Modifier"),
          definition = function(x, value){
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

#' @rdname Modifier-functions
#' @export
setMethod(f = "seqData", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@data})

#' @rdname Modifier-functions
#' @export
setMethod(f = "names", 
          signature = signature(x = "Modifier"),
          definition = function(x){names(seqData(x))})

#' @rdname Modifier-functions
#' @export
setMethod(f = "seqinfo", 
          signature = signature(x = "Modifier"),
          definition = function(x){seqinfo(seqData(x))})

#' @rdname Modifier-functions
#' @export
setMethod(f = "sequences", 
          signature = signature(x = "Modifier"),
          definition = 
            function(x, modified = FALSE){
              if(!assertive::is_a_bool(modified)){
                stop("'modified' has to be a single logical value.")
              }
              if(modified == FALSE){
                return(sequences(seqData(x)))
              }
              mod <- .get_modifications_per_transcript(x)
              mod <- split(mod,names(mod))
              ans <- ModRNAStringSet(sequences(seqData(x)))
              modSeqList <- ans[names(ans) %in% names(mod)]
              mod <- mod[match(names(mod),names(modSeqList))]
              ans[names(ans) %in% names(mod)] <- 
                combineIntoModstrings(modSeqList,
                                      mod)
              ans
            }
)

#' @rdname Modifier-functions
#' @export
setMethod(f = "ranges", 
          signature = signature(x = "Modifier"),
          definition = function(x){ranges(seqData(x))})

#' @rdname Modifier-functions
#' @export
setMethod(f = "bamfiles", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@bamfiles})

# converts the genomic coordinates to transcript based coordinates
.get_modifications_per_transcript <- function(x){
  grl <- ranges(x)
  strand_u <- .get_strand_u_GRangesList(grl)
  seqs <- .seqs_rl_by(grl, 1L)
  seqs[strand_u == "-"] <- .seqs_rl_by(grl[strand_u == "-"], -1L)
  modifications <- modifications(x)
  start_mod <- start(modifications)
  parent_mod <- as.character(modifications$Parent)
  new_start_mod <- which(seqs[parent_mod] == start_mod)
  # reset strand since it is now transcipt centric
  strand(modifications) <- "*"
  ranges(modifications) <- IRanges::IRanges(start = unlist(new_start_mod),
                                            end = unlist(new_start_mod))
  modifications
}

#' @rdname modify
#' @export
setMethod(f = "modifications", 
          signature = signature(x = "Modifier"),
          definition = 
            function(x, perTranscript = FALSE){
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

.new_ModFromCharacter <- function(className, x, annotation, sequences, seqinfo,
                                  ...){
  ans <- new(className, x)
  # settings
  settings(ans) <- list(...)
  #
  annotation <- .norm_annotation(annotation, className)
  sequences <- .norm_sequences(sequences, className)
  seqinfo <- .norm_seqnames(bamfiles(ans), annotation, sequences, seqinfo,
                            className)
  # get SequenceData
  data <- lapply(ans@dataType,
                 function(class){
                   do.call(class, c(list(bamfiles = bamfiles(ans),
                                         annotation = annotation,
                                         sequences = sequences,
                                         seqinfo = seqinfo),
                                    settings(ans)))
                 })
  if(length(data) > 1L){
    ans@data <- as(data,"SequenceDataList")
  } else {
    ans@data <- data[[1]]
  }
  #
  f <- which(shortName(ModRNAString()) %in% ans@mod)
  modName <- fullName(ModRNAString())[f]
  message("Starting to search for '", paste(tools::toTitleCase(modName), 
                                            collapse = "', '"),
          "' ... ", appendLF = FALSE)
  if(settings(ans,"findMod")){
    ans <- do.call(modify,
                   list(ans))
  }
  validObject(ans)
  ans <- .norm_modifications(ans, settings(ans))
  ans
}

.check_list_for_SequenceData_elements <- function(ans, list){
  if(is(ans,"character") && extends(ans,"Modifier")){
    ans <- getClass(ans)@prototype
  } else if(!is(ans,"Modifier")) {
    stop("Something went wrong.")
  }
  if(!is(list,"list")){
    list <- list(list)
  }
  .norm_SequenceData_elements(ans,list)
  if(length(list) == 1L){
    return(list[[1]])
  }
  as(list,"SequenceDataList")
}

.new_ModFromSequenceData <- function(className, x, ...){
  browser()
  ans <- new(className, bamfiles(x))
  # settings
  settings(ans) <- list(...)
  # check data type, length
  ans@data <- .check_list_for_SequenceData_elements(ans,x)
  # validate
  validObject(ans)
  #
  f <- which(shortName(ModRNAString()) %in% ans@mod)
  modName <- fullName(ModRNAString())[f]
  message("Starting to search for '", paste(tools::toTitleCase(modName),
                                            collapse = "', '"),
          "' ... ", appendLF = FALSE)
  if(settings(ans,"findMod")){
    ans <- do.call(modify,list(ans))
  }
  validObject(ans)
  ans <- .norm_modifications(ans, settings(ans))
  ans
}

#' @rdname Modifier-class
#' @export
setMethod("Modifier",
          signature = c(x = "SequenceData"),
          function(className, x, annotation = NULL, sequences = NULL, 
                   seqinfo = NULL, ...){
            .new_ModFromSequenceData(className, x, ...)
          })
#' @rdname Modifier-class
#' @export
setMethod("Modifier",
          signature = c(x = "character"),
          function(className, x, annotation = NULL, sequences = NULL, 
                   seqinfo = NULL, ...){
            .new_ModFromCharacter(className, x, annotation, sequences, seqinfo,
                                  ...)
          })
#' @rdname Modifier-class
#' @export
setMethod("Modifier",
          signature = c(x = "BamFileList"),
          function(className, x, annotation = NULL, sequences = NULL, 
                   seqinfo = NULL, ...){
            .new_ModFromCharacter(className, x, annotation, sequences, seqinfo,
                                  ...)
          })

# modify -----------------------------------------------------------------------

#' @name modify
#' @aliases modifications
#' 
#' @title Searching for modifications in \code{SequenceData}
#' 
#' @description 
#' \code{modify} triggers the search for modifications for a 
#' \code{\link[=Modifier-class]{Modifier}} class. Usually this is done 
#' automatically during construction of a \code{Modifier} object.
#' 
#' \code{modifications} is the accessor for the found modifications.
#' 
#' @param x a \code{Modifier} object.
#' @param force force to run \code{aggregate} again, if data is already stored
#' in \code{x}.
#' @param perTranscript For \code{modifications>} \code{TRUE} or \code{FALSE}:
#'   Should the coordinates be returned as local per transcript coordinates?
#' @param ... additional arguments
#' 
#' @return 
#' \itemize{
#' \item{\code{modify}: }{the updated \code{Modifier} object.}
#' \item{\code{modifications}: }{the modifications found as a \code{GRanges}
#' object.}
#' }
NULL

#' @rdname modify
#' @export
setMethod(f = "modify", 
          signature = signature(x = "Modifier"),
          definition = 
            function(x, force){
              stop("This functions needs to be implemented by '",class(x),"'.",
                   call. = FALSE)
            }
)


# aggregate --------------------------------------------------------------------

#' @name aggregate
#' @aliases hasAggregateData aggregateData
#' 
#' @title Aggreagte data per positions
#' 
#' @description 
#' The aggregate function is defined per
#' \code{\link[=SequenceData-class]{SequenceData}} object and can be triggered
#' by either using the a \code{\link[=SequenceData-class]{SequenceData}} or
#' \code{\link[=Modifier-class]{Modifier}}. For the letter it will just redirect
#' to the \code{\link[=SequenceData-class]{SequenceData}} object, but will store
#' the result. The data is then used for subsequent tasks, such as search for
#' modifications and visualization of the results.
#' 
#' @param x a \code{\link[=SequenceData-class]{SequenceData}} or
#'   \code{\link[=Modifier-class]{Modifier}} object.
#' @param force whether to recreate the aggregated data, if it is already stored
#' inside the \code{Modifier} object.
#' @param condition character value, which selects, for which condition the data
#' should be aggregated. One of the following values: \code{Both}, 
#' \code{Control}, \code{Treated} 
#' @param ... additional arguments
#' 
#' @return 
#' \itemize{
#' \item{\code{aggregate}: }{for \code{SequenceData} the aggregated data is
#' returned as a \code{SplitDataFrameList} with an element per transcript, 
#' whereas for a \code{Modifier} the modified input object is returned, 
#' containing the aggregated data, which can be accessed using 
#' \code{aggregateData}.}
#' \item{\code{aggregateData}: }{only for \code{Modifier}: a 
#' \code{SplitDataFrameList} with an element per transcript is returned. If the 
#' aggregated data is not stored in the object, it is generated on the fly, but 
#' does not persist}
#' \item{\code{hasAggregateData}: }{TRUE or FALSE. Does the \code{Modifier} 
#' object already contain aggregated data?}
#' }
NULL

#' @rdname aggregate
#' @export
setMethod(f = "aggregate", 
          signature = signature(x = "Modifier"),
          definition = 
            function(x, force = FALSE){
              stop("This functions needs to be implemented by '",class(x),"'.",
                   call. = FALSE)
            }
)

.add_positions_as_rownames <- function(x){
  if(is.null(rownames(x@aggregate@unlistData))){
    seqs <- IRanges::CharacterList(.seqs_rl(ranges(x)))
    rownames(x@aggregate@unlistData) <- unlist(seqs)
  }
  x
}
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
#' @rdname aggregate
#' @export
setMethod(f = "aggregateData", 
          signature = signature(x = "Modifier"),
          definition = function(x){
            x <- aggregate(x)
            x <- .add_positions_as_rownames(x)
            .check_score_name(x@aggregate, x@score)
          })

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
