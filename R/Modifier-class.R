#' @include RNAmodR.R
#' @include SequenceData-class.R
#' @include SequenceDataList-class.R
#' @include Modifier-utils.R
NULL

#' @name Modifier
#' @aliases modify
#' 
#' @title Modifier
#' @description 
#' title
NULL

setClassUnion("SequenceData_OR_SequenceDataList",
              c("SequenceData", "SequenceDataList"))

#' @rdname Modifier
#' @export
setClass("Modifier",
         contains = c("VIRTUAL"),
         slots = c(mod = "character", # this have to be populated by subclass
                   score = "character", # this have to be populated by subclass
                   dataType = "character", # this have to be populated by subclass
                   bamfiles = "BamFileList",
                   conditions = "factor",
                   fasta = "FaFile",
                   gff = "GFFFile",
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
    if(is.null(bamfiles) || is.null(fasta) || is.null(gff)){
      return(.Object)
    }
    className <- class(.Object)[[1]]
    className <- .norm_modifiertype(className)
    # check bam files
    bamfiles <- .norm_bamfiles(bamfiles,className)
    # set clots
    .Object@bamfiles <- bamfiles
    .Object@conditions <- factor(names(bamfiles))
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
    cat("A", class(object), "object containing",object@dataType,
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
#' @name Modifier
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

#' @name Modifier
#' @export
setMethod(f = "seqinfo", 
          signature = signature(x = "Modifier"),
          definition = function(x){x@seqinfo})
 
#' @name Modifier
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
  browser()
  ans <- new(className, x)
  # settings
  settings(ans) <- list(...)
  #
  f <- which(shortName(ModRNAString()) %in% ans@mod)
  modName <- fullName(ModRNAString())[f]
  message("Starting to search for '", paste(tools::toTitleCase(modName), 
                                            collapse = "', '"),
          "' ... ", appendLF = FALSE)
  #
  annotation <- .norm_annotation(annotation, className)
  sequences <- .norm_sequences(sequences, className)
  seqinfo <- .norm_seqnames(bamfiles(ans), annotation, sequences, seqinfo,
                            className)
  # get SequenceData
  data <- lapply(ans@dataType,
                 function(class){
                   do.call(class, c(list(files = bamfiles(ans),
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
  if(settings(ans,"findMod")){
    ans <- do.call(modify,
                   list(ans))
  }
  message("OK")
  validObject(ans)
  ans <- .norm_modifications(ans, settings(ans))
  ans
}

.check_list_for_SequenceData_elements <- function(ans, list){
  if(is(ans,"character") && extends(ans,"Modifier")){
    ans <- getclass(ans)@prototype
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
  ans <- new(className, x)
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
  message("OK")
  validObject(ans)
  ans <- .norm_modifications(ans, settings(ans))
  ans
}

setMethod("Modifier",
          signature = c(x = "SequenceData"),
          function(className, x, annotation = NULL, sequences = NULL, 
                   seqinfo = NULL, ...){
            .new_ModFromSequenceData(className, x, ...)
          })
setMethod("Modifier",
          signature = c(x = "character"),
          function(className, x, annotation = NULL, sequences = NULL, 
                   seqinfo = NULL, ...){
            .new_ModFromCharacter(className, x, annotation, sequences, seqinfo,
                                  ...)
          })
setMethod("Modifier",
          signature = c(x = "BamFileList"),
          function(className, x, annotation = NULL, sequences = NULL, 
                   seqinfo = NULL, ...){
            .new_ModFromCharacter(className, x, annotation, sequences, seqinfo,
                                  ...)
          })

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
