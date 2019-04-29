#' @include RNAmodR.R
#' @include SequenceData-class.R
#' @include SequenceDataSet-class.R
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
#' functionality to search for post-transcriptional RNA modification patterns in
#' high throughput sequencing data.
#'
#' Each subclass has to implement the following functions:
#'
#' \itemize{
#' \item{\code{\link{aggregateData}}: }{used for specific data aggregation}
#' \item{\code{\link{findMod}}: }{used for specific search for modifications}
#' }
#'
#' Optionally the function \code{\link[=Modifier-functions]{settings<-}} can be
#' implemented to store additional arguments, which the base class does not
#' recognize.
#'
#' \code{Modifier} objects are constructed centrally by calling
#' \code{Modifier()} with a \code{className} matching the specific class to be
#' constructed. This will trigger the immediate analysis, if \code{find.mod} is
#' not set to \code{FALSE}.
#'
#'
#' @section Creation:
#' \code{Modifier} objects can be created in two ways, either by providing a
#' list of bamfiles or
#' \code{SequenceData}/\code{SequenceDataSet}/\code{SequenceDataList} objects,
#' which match the structure in \code{dataType()}.
#'
#' \code{dataType()} can be a \code{character} vector or a \code{list} of
#' \code{character} vectors and depending on this the input files have to
#' follow this structure:
#'
#' \itemize{
#' \item{a single \code{character}:} {a \code{SequenceData} is
#' constructed/expected.}
#' \item{a \code{character} vector:} {a \code{SequenceDataSet} is
#' constructed/expected.}
#' \item{a \code{list} of \code{character} vectors:} {a \code{SequenceDataList}
#' is constructed/expected.}
#' }
#'
#' The cases for a \code{SequenceData} or \code{SequenceDataSet} are straight
#' forward, since the input remains the same. The last case is special, since it
#' is a hypothetical option, in which bam files from two or more different
#' methods have to be combined to reliably detect a single modification (The
#' elements of a \code{SequenceDataList} don't have to be created from the
#' bamfiles, whereas from a \code{SequenceDataSet} they have to be).
#'
#' For this example a \code{list} of \code{character} vectors is expected.
#' Each element must be named according to the names of \code{dataType()} and
#' contain a \code{character} vector for creating a \code{SequenceData} object.
#'
#' @param className The name of the class which should be constructed.
#' @param x the input which can be of the following types
#' \itemize{
#' \item{\code{SequenceData}:} {a single \code{SequenceData} or a list
#' containing only \code{SequenceData} objects. The input will just be used to
#' file the \code{data} slot of the \code{Modifier} and must match the
#' requirements of specific \code{Modifier} class.}
#' \item{\code{BamFileList}:} {a named \code{BamFileList}}
#' \item{\code{character}:} {a \code{character} vector, which must be coercible
#' to a named \code{BamFileList} referencing existing bam files. Valid names are
#' \code{control} and \code{treated} to define conditions and replicates}
#' }
#' @param annotation annotation data, which must match the information contained
#' in the BAM files. This parameter is only required if \code{x} is not a
#' \code{SequenceData} object or a list of \code{SequenceData} objects.
#' @param sequences sequences matching the target sequences the reads were
#'   mapped onto. This must match the information contained in the BAM files.
#'   TThis parameter is only required if \code{x} is not a \code{SequenceData}
#'   object or a list of \code{SequenceData} objects.
#' @param seqinfo An optional \code{\link[GenomeInfoDb:Seqinfo-class]{Seqinfo}}
#' argument or character vector, which can be coerced to one, to subset the
#' sequences to be analyzed on a per chromosome basis.
#' @param ... Additional otpional parameters:
#' \itemize{
#' \item{\code{find.mod}:} {\code{TRUE} or \code{FALSE}: should the search for
#' for modifications be triggered upon construction? If not the search can be
#' started by calling the \code{modify()} function.}
#' }
#' All other arguments will be passed onto the \code{SequenceData} objects, if
#' \code{x} is not a \code{SequenceData} object or a list of \code{SequenceData}
#' objects.
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
#' @slot data The sequence data object: Either a \code{SequenceData},
#' \code{SequenceDataSet} or a \code{SequenceDataList} object, if more than one
#' \code{dataType} is used.
#' @slot aggregate the aggregated data as a \code{SplitDataFrameList}
#' @slot modifications the found modifications as a \code{GRanges} object
#' @slot arguments arguments used for the analysis as a \code{list}
#' @slot aggregateValidForCurrentArguments \code{TRUE} or \code{FALSE} whether
#' the aggregate data was constructed with the current arguments
#' @slot modificationsValidForCurrentArguments \code{TRUE} or \code{FALSE}
#' whether the modifications were found with the current arguments
#'
#' @return a \code{Modifier} object of type \code{className}
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
#' @param perTranscript \code{TRUE} or \code{FALSE}: Should the positions shown
#' per transcript? (default: \code{perTranscript = FALSE})
#' @param ... Additional arguments.
#'
#' @return
#' \itemize{
#' \item{\code{modifierType}:} {a character vector with the appropriate class
#' Name of a \code{\link[=Modifier-class]{Modifier}}.}
#' \item{\code{mainScore}:} {a character vector.}
#' \item{\code{settings}:} {a \code{Seqinfo} object.}
#' \item{\code{sequenceData}:} {a \code{SequenceData} object.}
#' \item{\code{modifications}:} {a \code{GRanges} or \code{GRangesList} object
#' describing the found modifications.}
#' \item{\code{seqinfo}:} {a \code{Seqinfo} object.}
#' \item{\code{sequences}:} {a \code{RNAStingSet} object.}
#' \item{\code{ranges}:} {a \code{GRangesList} object with each element per
#' transcript.}
#' \item{\code{bamfiles}:} {a \code{BamFileList} object.}
#' }
#'
#' @examples
#' data(msi,package="RNAmodR")
#' mi <- msi[[1]]
#' modifierType(mi) # The class name of the Modifier object
#' modifierType(msi) #
#' mainScore(mi)
#' settings(mi)
#' sequenceData(mi)
#' modifications(mi)
#' # general accessors
#' seqinfo(mi)
#' sequences(mi)
#' ranges(mi)
#' bamfiles(mi)
NULL

setClassUnion("list_OR_character",
              c("list", "character"))
#' @importClassesFrom Rsamtools BamFileList PileupFiles
setClassUnion("list_OR_BamFileList",
              c("list", "BamFileList"))

#' @rdname Modifier-class
#' @export
setClass("Modifier",
         contains = c("VIRTUAL"),
         slots = c(mod = "character", # this have to be populated by subclass
                   score = "character", # this have to be populated by subclass
                   dataType = "list_OR_character", # this have to be populated by subclass
                   bamfiles = "list_OR_BamFileList",
                   condition = "factor",
                   replicate = "factor",
                   data = "SD_or_SDS_or_SDL",
                   aggregate = "CompressedSplitDataFrameList",
                   modifications = "GRanges",
                   arguments = "list",
                   aggregateValidForCurrentArguments = "logical",
                   modificationsValidForCurrentArguments = "logical"),
         prototype = list(aggregateValidForCurrentArguments = FALSE,
                          modificationsValidForCurrentArguments = FALSE))

# validity ---------------------------------------------------------------------

.check_SequenceData_elements <- function(x, list){
  if(is(list,"SequenceData")){
    list <- list(list)
  } else if(is(list,"list")){
    elementTypeMatch <- !vapply(list,is,logical(1),"SequenceData")
    if(any(elementTypeMatch)){
      stop("Not all elements are 'SequenceData' objects.", call. = FALSE)
    }
  } else if(!is(list,"SequenceDataSet")){
   stop("Something went wrong.")
  }
  elementTypes <- vapply(list,class,character(1))
  if(length(elementTypes) != length(dataType(x))){
    stop("Number of 'SequenceData' elements does not match the requirements of",
         " ",class(x),". '",paste(dataType(x), collapse = "','"),"' are ",
         "required", call. = FALSE)
  }
  elementTypes <- elementTypes[match(elementTypes,dataType(x))]
  if(is.na(elementTypes) || any(elementTypes != dataType(x))){
    stop("Type of SequenceData elements does not match the requirements of ",
         class(x),". '",paste(dataType(x), collapse = "','"),"' are ",
         "required", call. = FALSE)
  }
  NULL
}

.check_SequenceDataList_data_elements <- function(x, list){
  ans <- lapply(list, .check_SequenceData_elements, x)
  if(all(vapply(ans,is.null))) {
    return(NULL)
  }
  ans
}

.check_Modifier_data_elements <- function(x, data){
  if(is(data,"SequenceData") || is(data,"SequenceDataSet")){
    return(.check_SequenceData_elements(x, data))
  }
  .check_SequenceDataList_data_elements(x, data)
}

.valid_SequenceData <- function(x){
  tmp <- try(.check_Modifier_data_elements(x,x@data))
  if (inherits(tmp, "try-error")){
    return(tmp)
  }
  NULL
}

.valid_Modifier <- function(x){
  seqdata <- x@data
  if(is.list(x@bamfiles)){
    test <- !vapply(x@bamfiles,is,logical(1),"BamFileList")
    if(any(test)){
      return("Invalid slot: 'bamfiles'")
    }
  }
  if(is(seqdata,"SequenceData")){
    seqdata <- list(seqdata)
  }
  if(!assertive::is_a_bool(x@aggregateValidForCurrentArguments)){
    return("Invalid slot: 'aggregateValidForCurrentArguments'")
  }
  if(!assertive::is_a_bool(x@aggregateValidForCurrentArguments)){
    return("Invalid slot: 'modificationsValidForCurrentArguments'")
  }
  if(hasAggregateData(x)){
    data <- getAggregateData(x)
    if(is.null(rownames(data@unlistData))){
      return("rownames of aggregate data is not set.")
    } else {
      seqs <- .seqs_rl_strand(ranges(x))
      if(any(any(seqs != IRanges::IntegerList(rownames(data))))){
        return(paste0("Out of range rownames of aggregate data. The rownames ",
                      "do not match the genomic coordinate covered by the ",
                      "annotation data."))
      }
    }
  }
  c(.valid_SequenceData(x), unlist(lapply(seqdata, .valid.SequenceData)))
}
S4Vectors::setValidity2(Class = "Modifier", .valid_Modifier)

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
#' @importFrom BiocGenerics path
setMethod(
  f = "show",
  signature = signature(object = "Modifier"),
  definition = function(object) {
    cat("A", class(object), "object containing",dataType(object),
        "with",length(object@data),"elements.\n")
    files <- BiocGenerics::path(object@bamfiles)
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
    valid <- c(validAggregate(object), validModification(object))
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
setMethod(f = "bamfiles",
          signature = signature(x = "Modifier"),
          definition = function(x){x@bamfiles})
#' @rdname Modifier-functions
#' @export
setMethod(f = "conditions",
          signature = signature(object = "Modifier"),
          definition = function(object){
            conditions(sequenceData(object))
          })
#' @rdname Modifier-functions
#' @export
setMethod(f = "mainScore",
          signature = signature(x = "Modifier"),
          definition = function(x){x@score})

# converts the genomic coordinates to transcript based coordinates
.get_modifications_per_transcript <- function(x){
  grl <- ranges(x)
  strand_u <- .get_strand_u_GRangesList(grl)
  seqs <- .seqs_rl_by(grl, 1L)
  seqs[strand_u == "-"] <- .seqs_rl_by(grl[strand_u == "-"], -1L)
  modifications <- modifications(x)
  if(is(modifications,"GRangesList")){
    modifications <- unlist(modifications)
    modifications <- modifications[!duplicated(modifications)]
  }
  start_mod <- start(modifications)
  parent_mod <- as.character(modifications$Parent)
  new_start_mod <- BiocGenerics::which(seqs[parent_mod] == start_mod)
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
setMethod(f = "dataType",
          signature = signature(x = "Modifier"),
          definition = function(x){x@dataType})
#' @rdname Modifier-functions
#' @export
setMethod(f = "names",
          signature = signature(x = "Modifier"),
          definition = function(x){names(sequenceData(x))})
#' @rdname Modifier-functions
#' @export
setMethod(f = "ranges",
          signature = signature(x = "Modifier"),
          definition = function(x){ranges(sequenceData(x))})
#' @rdname Modifier-functions
#' @export
setMethod(f = "replicates",
          signature = signature(x = "Modifier"),
          definition = function(x){
            replicates(sequenceData(x))
          })
#' @rdname Modifier-functions
#' @export
setMethod(f = "sequenceData",
          signature = signature(x = "Modifier"),
          definition = function(x){x@data})
#' @rdname Modifier-functions
#' @export
setMethod(f = "sequences",
          signature = signature(x = "Modifier"),
          definition =
            function(x, modified = FALSE){
              if(!assertive::is_a_bool(modified)){
                stop("'modified' has to be a single logical value.",
                     call. = FALSE)
              }
              if(modified == FALSE){
                return(sequences(sequenceData(x)))
              }
              mod <- .get_modifications_per_transcript(x)
              mod <- .rebase_seqnames(mod, mod$Parent)
              mod <- split(mod,factor(mod$Parent, levels = mod$Parent))
              ans <- ModRNAStringSet(sequences(sequenceData(x)))
              modSeqList <- ans[names(ans) %in% names(mod)]
              mod <- mod[match(names(mod),names(modSeqList))]
              ans[names(ans) %in% names(mod)] <-
                Modstrings::combineIntoModstrings(modSeqList, mod)
              ans
            }
)
#' @rdname Modifier-functions
#' @export
setMethod(f = "seqinfo",
          signature = signature(x = "Modifier"),
          definition = function(x){seqinfo(sequenceData(x))})

.norm_args <- function(input){
  minCoverage <- 10L
  minReplicate <- 1L
  find.mod <- TRUE
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
  if(!is.null(input[["find.mod"]])){
    find.mod <- input[["find.mod"]]
    if(!assertive::is_a_bool(find.mod)){
      stop("'find.mod' must be a single logical value.")
    }
  }
  args <- list(minCoverage = minCoverage,
               minReplicate = minReplicate,
               find.mod = find.mod)
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
setMethod(f = "validAggregate",
          signature = signature(x = "Modifier"),
          definition = function(x) x@aggregateValidForCurrentArguments
)
#' @rdname Modifier-functions
#' @export
setMethod(f = "validModification",
          signature = signature(x = "Modifier"),
          definition = function(x) x@modificationsValidForCurrentArguments
)

# constructors -----------------------------------------------------------------

.norm_Modifier_input_SequenceData_elements <- function(list, ans){
  if(is(ans,"character") && extends(ans,"Modifier")){
    ans <- getClass(ans)@prototype
  } else if(!is(ans,"Modifier")) {
    stop("Something went wrong.")
  }
  .check_Modifier_data_elements(ans, list)
  if(length(list) == 1L){
    return(list[[1]])
  }
  list
}

.Modifier <- function(className, data){
  proto <- new(className)  # create prototype object for mod normalization only
  data <- .norm_Modifier_input_SequenceData_elements(data, proto)
  bamfiles <- bamfiles(data)
  condition <- factor(names(bamfiles))
  new2(className,
       mod = .norm_mod(proto@mod, className),
       bamfiles = bamfiles,
       condition = condition,
       replicate = .get_replicate_number(bamfiles, condition),
       data = data)
}

.load_SequenceData <- function(classes, bamfiles, annotation, sequences,
                               seqinfo, args){
  if(is.list(classes)){
    if(is.list(bamfiles)){
      if(length(classes) != length(bamfiles)){
        stop("'x' has invalid length. '",paste(classes, collapse = "' and '"),
             "' ",ifelse(length(classes) > 1L,"are","is")," required.",
             call. = FALSE)
      }
      if(is.null(names(bamfiles))){
        stop("If 'x' is a list, it must be named. The names must match the ",
             "required SequenceData class names.",
             call. = FALSE)
      }
      if(all(classes %in% names(bamfiles))){
        stop("If 'x' is named, names must match the required SequenceData ",
             "class names.",
             call. = FALSE)
      }
      bamfiles <- bamfiles[match(class,names(bamfiles))]
      data <- BiocParallel::bpmapply(.load_SequenceData, classes, bamfiles,
                                     MoreArgs = list(annotation, sequences,
                                                     seqinfo, args))
    } else {
      data <- BiocParallel::bplapply(classes, .load_SequenceData, bamfiles,
                                     annotation, sequences, seqinfo, args)
    }
    data <- as(data,"SequenceDataList")
  } else if(is.character(classes)){
    data <- lapply(classes,
                   function(class){
                     do.call(class, c(list(bamfiles = bamfiles,
                                           annotation = annotation,
                                           sequences = sequences,
                                           seqinfo = seqinfo),
                                      args))
                   })
    data <- as(data,"SequenceDataSet")
  } else {
    stop("Something went wrong.")
  }
  data
}

.new_ModFromCharacter <- function(className, x, annotation, sequences, seqinfo,
                                  ...){
  # Check that external classes are implemented correctly
  className <- .norm_modifiertype(className)
  # create prototype object for mod normalization and settings only
  proto <- new(className)
  # short cut for creating an empty object
  if(is.null(x)){
    return(new2(className, mod = .norm_mod(proto@mod, className)))
  }
  bamfiles <- .norm_bamfiles(x, className) # check bam files
  # settings
  settings(proto) <- list(...)
  #
  annotation <- .norm_annotation(annotation, className)
  sequences <- .norm_sequences(sequences, className)
  seqinfo <- .norm_seqnames(bamfiles, annotation, sequences, seqinfo, className)
  annotation <- .load_annotation(annotation)
  annotation <- .subset_by_seqinfo(annotation, seqinfo)
  # get SequenceData
  data <- .load_SequenceData(dataType(proto), bamfiles = bamfiles,
                             annotation = annotation, sequences = sequences,
                             seqinfo = seqinfo, args = settings(proto))
  .new_ModFromSequenceData(className, data, ...)
}

.new_ModFromSequenceData <- function(className, x, ...){
  # create Modifier object
  ans <- .Modifier(className, x)
  # settings
  settings(ans) <- list(...)
  # aggregate data
  message("Aggregating data and calculating scores ... ", appendLF = FALSE)
  ans <- aggregate(ans)
  # search for modifications
  if(settings(ans,"find.mod")){
    f <- which(Modstrings::shortName(Modstrings::ModRNAString()) %in% ans@mod)
    modName <- Modstrings::fullName(Modstrings::ModRNAString())[f]
    message("Starting to search for '", paste(tools::toTitleCase(modName),
                                              collapse = "', '"),
            "' ... ", appendLF = FALSE)
    ans <- modify(ans)
  }
  message("done.")
  validObject(ans)
  ans
}

#' @rdname Modifier-class
#' @export
setGeneric(
  name = "Modifier",
  signature = c("x"),
  def = function(className, x, annotation, sequences, seqinfo, ...)
    standardGeneric("Modifier")
)

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
          signature = c(x = "SequenceDataSet"),
          function(className, x, annotation = NULL, sequences = NULL,
                   seqinfo = NULL, ...){
            .new_ModFromSequenceData(className, x, ...)
          })
#' @rdname Modifier-class
#' @export
setMethod("Modifier",
          signature = c(x = "SequenceDataList"),
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
          signature = c(x = "list"),
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

# aggregate --------------------------------------------------------------------

#' @name aggregate
#' @aliases hasAggregateData aggregateData getAggregateData
#'
#' @title Aggregate data per positions
#'
#' @description
#' The \code{aggregate} function is defined for each
#' \code{\link[=SequenceData-class]{SequenceData}} object and can be used
#' directly on a \code{\link[=SequenceData-class]{SequenceData}} object or
#' indirectly via a \code{\link[=Modifier-class]{Modifier}} object.
#'
#' For the letter the call is redirect to the
#' \code{\link[=SequenceData-class]{SequenceData}} object, the result summarized
#' as defined for the individual \code{Modifier} class and stored in the
#' \code{aggregate} slot of the \code{Modifier} object. The data is then used
#' for subsequent tasks, such as search for modifications and visualization of
#' the results.
#'
#' The summarization is implemented in the \code{aggregateData} for each type of
#' \code{Modifier} class. The stored data from the \code{aggregate} slot can be
#' retrieved using the \code{getAggregateData} function.
#'
#' Whether the aggrgeated data is already present in the \code{aggregate} slot
#' can be checked using the \code{hasAggregateData} function.
#'
#' For \code{SequenceDataSet}, \code{SequenceDataList} and \code{ModfierSet}
#' classes wrapper of the \code{aggregate} function exist as well.
#'
#' @param x a \code{\link[=SequenceData-class]{SequenceData}},
#' \code{SequenceDataSet}, \code{SequenceDataList},
#' \code{\link[=Modifier-class]{Modifier}} or
#' \code{\link[=Modifier-class]{ModfierSet}}  object.
#' @param force whether to recreate the aggregated data, if it is already stored
#' inside the \code{Modifier} object.
#' @param condition character value, which selects, for which condition the data
#' should be aggregated. One of the following values: \code{Both},
#' \code{Control}, \code{Treated}
#' @param ... additional arguments
#'
#' @return
#' \itemize{
#' \item{\code{aggregate}: }{for \code{SequenceData} object the aggregated data
#' is returned as a \code{SplitDataFrameList} with an element per transcript,
#' whereas for a \code{Modifier} the modified input object is returned,
#' containing the aggregated data, which can be accessed using
#' \code{getAggregateData}.}
#' \item{\code{getAggregateData}: }{only for \code{Modifier}: a
#' \code{SplitDataFrameList} with an element per transcript is returned. If the
#' aggregated data is not stored in the object, it is generated on the fly, but
#' does not persist.}
#' \item{\code{hasAggregateData}: }{TRUE or FALSE. Does the \code{Modifier}
#' object already contain aggregated data?}
#' }
#'
#' @return
#' If 'x' is a
#' \itemize{
#' \item{\code{\link[=SequenceData-class]{SequenceData}}} {a
#' \code{SplitDataFrameList} with elments per transcript.}
#' \item{\code{\link[=SequenceDataSet-class]{SequenceDataSet}} or
#' \code{\link[=SequenceDataList-class]{SequenceDataList}}} {a \code{SimpleList}
#' with \code{SplitDataFrameList} as elements.}
#' \item{\code{\link[=Modifier-class]{Modifier}} or
#' \code{\link[=ModifierSet-class]{ModifierSet}}} {an updated \code{Modifier}
#' object. The data can be accessed by using the \code{aggregateData} function.}
#' }
#'
#' @examples
#' data(e5sd,package="RNAmodR")
#' data(msi,package="RNAmodR")
#' # modify() triggers the search for modifications in the data contained in
#' # the Modifier or ModifierSet object
#' sdfl <- aggregate(e5sd)
#' mi <- aggregate(msi[[1]])
NULL

.check_aggregate_modifier <- function(data, x){
  score <- x@score
  if(is(data,"CompressedSplitDataFrameList")){
    columns <- colnames(data@unlistData)
  } else {
    stop("aggregate data has to be a 'CompressedSplitDataFrameList' object. ",
         "Contact the maintainer of the class used.",
         call. = FALSE)
  }
  if(!(score %in% columns)){
    stop("The default score is not present in the aggregate data. Contact the ",
         "maintainer of the class used.",
         call. = FALSE)
  }
  seqs <- .seqs_rl_strand(ranges(x))
  if(any(any(seqs != IRanges::IntegerList(rownames(data))))){
    stop("rownames() of aggregate data is not named or not named along the ",
         "genomic coordinates. Contact the maintainer of the class used.",
         call. = FALSE)
  }
  data
}

#' @rdname aggregate
#' @export
setMethod(f = "aggregate",
          signature = signature(x = "Modifier"),
          definition =
            function(x, force = FALSE){
              if(missing(force)){
                force <- FALSE
              }
              assertive::assert_is_a_bool(force)
              if(!hasAggregateData(x) || force){
                x@aggregate <- .check_aggregate_modifier(aggregateData(x), x)
                x@aggregateValidForCurrentArguments <- TRUE
              }
              x
            }
)

#' @rdname aggregate
#' @export
setMethod(f = "aggregateData",
          signature = signature(x = "Modifier"),
          definition =
            function(x){
              stop("The 'aggregateData' functions needs to be implemented by
                   '",class(x),"'.",call. = FALSE)
            }
)
#' @rdname aggregate
#' @export
setMethod(f = "getAggregateData",
          signature = signature(x = "Modifier"),
          definition = function(x){
            x <- aggregate(x)
            x@aggregate
          })

#' @rdname aggregate
#' @export
setMethod(f = "hasAggregateData",
          signature = signature(x = "Modifier"),
          definition =
            function(x){
              if(is.null(x@aggregate) || nrow(x@aggregate@unlistData) == 0L){
                return(FALSE)
              }
              return(TRUE)
            }
)

# modify -----------------------------------------------------------------------

#' @name modify
#' @aliases modifications
#'
#' @title Searching for modifications in \code{SequenceData}
#'
#' @description
#' The \code{modify} function executes the search for modifications for a
#' \code{\link[=Modifier-class]{Modifier}} class. Usually this is done
#' automatically during construction of a \code{Modifier} object.
#'
#' When the \code{modify} functions is called, the aggregated data is checked
#' for validity for the current settings and the search for modifications is
#' performed using the \code{findMod}. The results are stored in the
#' \code{modification} slot of the \code{Modifier} object, which is returned by
#' \code{modify}. The results can be accessed via the \code{modifications()}
#' function.
#'
#' \code{findMod} returns the found modifications as a \code{GRanges}
#' object and has to be implemented for each individual \code{Modifier} class.
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
#'
#' @examples
#' data(msi,package="RNAmodR")
#' # modify() triggers the search for modifications in the data contained in
#' # the Modifier or ModifierSet object
#' mi <- modify(msi[[1]])
NULL

#' @rdname modify
#' @export
setMethod(f = "modify",
          signature = signature(x = "Modifier"),
          definition =
            function(x, force = FALSE){
              if(missing(force)){
                force <- FALSE
              }
              assertive::assert_is_a_bool(force)
              if(!validAggregate(x) | force){
                x <- aggregate(x, force = TRUE)
              }
              x@modifications <- findMod(x)
              x@modificationsValidForCurrentArguments <- TRUE
              x
            }
)

#' @rdname modify
#' @export
setMethod(f = "findMod",
          signature = signature(x = "Modifier"),
          definition =
            function(x){
              stop("The 'findMod' functions needs to be implemented by
                   '",class(x),"'.",call. = FALSE)
            }
)
