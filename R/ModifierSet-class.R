#' @include RNAmodR.R
#' @include Modifier-class.R
NULL

#' @name ModifierSet-class
#' @aliases ModifierSet
#' 
#' @title The ModifierSet class
#' 
#' @description 
#' The \code{ModifierSet} class allows multiple
#' \code{\link[=Modifier-class]{Modifier}} objects to be created from the same
#' annotation and sequence data varying only the bam input files.
#' 
#' In addition the comparison of samples is also done via calling functions on 
#' the \code{ModifierSet} objects.
#' 
#' The \code{ModifierSet} is a virtual class, which derives from the 
#' \code{SimpleList} with the slot \code{elementType = "Modifier"}. As 
#' the \code{Modifier} tt must be implemented for the specific analysis.
#' 
#' @param className The name of the class which should be constructed.
#' @param x the input which can be of the following types
#' \itemize{
#' \item{\code{Modifier}:} {a single \code{Modifier} or a list containg only 
#' \code{Modifier} objects. The input will just be used as elements of the
#' \code{ModifierSet}}
#' \item{\code{BamFileList}:} {a named \code{BamFileList} or a list of 
#' named \code{BamFileList}}
#' \item{\code{list}:} {a list of one or more types of elements: 
#' \code{BamFileList}, a named \code{list} or named \code{character} vector. All
#' elements must be or be coercible to a named \code{BamFileList} referencing 
#' existing bam files. Valid names are \code{control} and \code{treated}}
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
#' \item{internalBP} {\code{TRUE} or \code{FALSE}: should 
#' parallilazation used internally during creation of each \code{Modifier} or
#' should the creation of the \code{Modifier} objects be parallalized? (default:
#' \code{internalBP = FALSE})}
#' }
#' All other arguments will be passed onto the \code{Modifier} objects.
#' 
#' @return a \code{ModifierSet} object of type 'className'
NULL

#' @rdname ModifierSet-class
#' @export
setClass("ModifierSet",
         contains = c("VIRTUAL",
                      "SimpleList"),
         prototype = list(elementType = "Modifier"))

setMethod("pcompareRecursively", "ModifierSet", function(x) FALSE)

# validity ---------------------------------------------------------------------

.valid_ModifierSet <- function(x){
  elementTypeX <- modifierType(x)
  if (!all(vapply(as.list(x),
                  function(xi) extends(class(xi), elementTypeX),
                  logical(1)))){
    return(paste("All 'Modifier' in '",class(x),"' must be of ",
                 elementTypeX, "objects"))
  }
  valid_Modifier <- lapply(x@listData, .valid_Modifier)
  valid_Modifier <- valid_Modifier[!vapply(valid_Modifier, is.null, logical(1))]
  if(length(valid_Modifier) != 0L){
    return(paste(paste0(seq_along(valid_Modifier), ". :", valid_Modifier),
                 collapse = "\n"))
  }
  NULL
}
S4Vectors::setValidity2(Class = "ModifierSet", .valid_ModifierSet)

# not supported functions ------------------------------------------------------

setMethod(f = "relistToClass",
          signature = c(x = "ModifierSet"),
          function(x) {
            stop("Relisting not supported.")
          })

# contructor -------------------------------------------------------------------

#' @importFrom stringr str_detect
.norm_classname_ModifierSet <- function(classname){
  if(stringr::str_detect(classname,"ModSet")){
    ans <- classname
  } else {
    ans <- gsub("Mod","ModSet",classname)
  }
  ans
}

.get_classname_for_ModifierSet_from_modifier_type <- function(modifiertype){
  ans <- .norm_classname_ModifierSet(modifiertype)
  class <- try(getClass(ans), silent = TRUE)
  if (is(class, "try-error")){
    stop("Class '",ans,"' is not implemented.",
         call. = FALSE)
  }
  if(isVirtualClass(class)){
    stop("Class '",ans,"' is virtual.")
  }
  if(!("ModifierSet" %in% extends(class))){
    stop("Class '",ans,"' does not extend the 'ModifierSet' class.")
  }
  ans
}

.ModifierSet <- function(className, x){
  new2(.get_classname_for_ModifierSet_from_modifier_type(className),
       listData = x)
}

.norm_ModifierSet_args <- function(input){
  internalBP <- FALSE
  if(!is.null(input[["internalBP"]])){
    internalBP <- input[["internalBP"]]
    if(!assertive::is_a_bool(internalBP)){
      stop("'internalBP' must be TRUE or FALSE.",
           call. = FALSE)
    }
  }
  args <- list(internalBP = internalBP)
  args
}

.contains_only_Modifier <- function(x){
  classNames <- unique(vapply(x, function(z){class(z)[[1]]},character(1)))
  if(length(classNames) != 1L){
    return(FALSE)
  }
  allSameClass <- vapply(x,
                         function(z, c){
                           class(z)[[1]] == c
                         },
                         logical(1),
                         classNames)
  if(!all(allSameClass)){
    return(FALSE)
  }
  x <- try(.norm_modifiertype(classNames), silent = TRUE)
  if (is(x, "try-error")){
    return(FALSE)
  }
  TRUE
}

#' @importFrom Rsamtools BamFileList 
#' @importFrom BiocGenerics path
.contains_only_bamfiles <- function(x){
  x <- unname(x)
  classNames <- vapply(x, function(z){class(z)[[1]]},character(1))
  if(!all(unique(classNames) %in% c("BamFileList","character","list"))){
    return(FALSE)
  }
  namedRequired <- x[classNames %in% c("character","list")]
  names <- unique(names(unlist(namedRequired)))
  if(is.null(names) || 
     !all(tolower(names) %in% c("treated","control"))){
    return(FALSE)
  }
  x <- lapply(x,.norm_bamfiles)
  x <- unlist(x) # a list of BamFileList cannot be unlisted. However this 
                 # normalizes x
  if(!all(vapply(x,is,logical(1),"BamFileList"))){
    return(FALSE)
  }
  if(!all(unlist(lapply(lapply(x,BiocGenerics::path),file.exists)))){
    return(FALSE)
  }
  TRUE
}

#' @importFrom BiocParallel SerialParam register bpmapply bplapply
.bamfiles_to_ModifierSet <- function(className, x, annotation, sequences,
                                     seqinfo, ...){
  # check and normalize input
  args <- .norm_ModifierSet_args(list(...))
  className <- .norm_modifiertype(className)
  if(!is.list(x)){
    x <- list(x)
  }
  names <- as.list(names(x))
  if(length(names) == 0L){
    names <- vector(mode = "list", length = length(x))
  }
  x <- lapply(x, .norm_bamfiles, className)
  annotation <- .norm_annotation(annotation, className)
  annotation <- .load_annotation(annotation)
  sequences <- .norm_sequences(sequences, className)
  if(missing(seqinfo)){
    seqinfo <- NA
  }
  ni <- seq_along(x)
  # choose were to use parallelization
  if(args[["internalBP"]] == TRUE){
    BiocParallel::register(BiocParallel::SerialParam())
  }
  # do analysis by calling the Modifier classes
  FUN <- function(i, z, n, args, className, PACKAGE, CLASSFUN, annotation,
                  sequences, seqinfo, ...){
    suppressPackageStartupMessages({
      requireNamespace(PACKAGE)
    })
    if(!is.null(n)){
      message(i,". ",className," analysis '",n,"':")
    } else {
      message(i,". ",className," analysis:")
    }
    # choose were to use parallelization
    if(args[["internalBP"]] == FALSE){
      BiocParallel::register(BiocParallel::SerialParam())
    }
    # do not pass this argument along to objects
    args[["internalBP"]] <- NULL
    do.call(CLASSFUN,
            c(list(z,
                   annotation = annotation,
                   sequences = sequences,
                   seqinfo = seqinfo),
              list(...)))
  }
  PACKAGE <- getClass(className)@package
  CLASSFUN <- get(className)
  x <- BiocParallel::bpmapply(FUN,
                              ni, x, names,
                              MoreArgs = list(args = args, 
                                              className = className, 
                                              PACKAGE = PACKAGE,
                                              CLASSFUN = CLASSFUN,
                                              annotation = annotation,
                                              sequences = sequences,
                                              seqinfo = seqinfo,
                                              ...),
                              SIMPLIFY = FALSE)
  f <- vapply(names,is.null,logical(1))
  names[f] <- as.list(as.character(seq_along(x))[f])
  names(x) <- unlist(names)
  # pass results to ModifierSet object
  .ModifierSet(className, x)
}

.Modifer_to_ModifierSet <- function(className, x, ...){
  if(!is.list(x)){
    x <- list(x)
  }
  elementType <- modifierType(x[[1]])
  className <- .get_classname_for_ModifierSet_from_modifier_type(className)
  if(className != .norm_classname_ModifierSet(elementType)){
    stop("Something went wrong.")
  }
  if (!all(vapply(x,
                  function(xi) extends(class(xi), elementType),
                  logical(1)))){
    return(paste("All 'Modifier' in '",className,"' must be of ",
                 elementType, " objects"))
  }
  .ModifierSet(className, x)
}

#' @rdname ModifierSet-class
#' @export
setMethod(f = "ModifierSet",
          signature = c(x = "list"),
          function(className, x, annotation = NULL, sequences = NULL, 
                   seqinfo = NULL, ...) {
            if(.contains_only_Modifier(x)){
              return(.Modifer_to_ModifierSet(className, x, ...))
            }
            if(.contains_only_bamfiles(x)){
              return(.bamfiles_to_ModifierSet(className, x, annotation, 
                                              sequences, seqinfo, ...))
            }
            stop("'x' must be a list containing only elements of the same ",
                 "type\nof 'Modifer' or elements of type ('BamFileList', ",
                 "'character', 'list') which are coercible\nto a named ",
                 "BamFileList. In the latter case, elements must contain named",
                 " vectors or lists('treated' or 'control')\nand the files ",
                 "referenced must exist. Please note, that the list a",
                 call. = FALSE)
          })
#' @rdname ModifierSet-class
#' @export
setMethod(f = "ModifierSet",
          signature = c(x = "character"),
          function(className, x, annotation = NULL, sequences = NULL, 
                   seqinfo = NULL, ...) {
            .bamfiles_to_ModifierSet(className, x, annotation, sequences,
                                     seqinfo, ...)
          })
#' @rdname ModifierSet-class
#' @export
setMethod(f = "ModifierSet",
          signature = c(x = "BamFileList"),
          function(className, x, annotation = NULL, sequences = NULL,
                   seqinfo = NULL, ...) {
            .bamfiles_to_ModifierSet(className, x, annotation, sequences,
                                     seqinfo, ...)
          })
#' @rdname ModifierSet-class
#' @export
setMethod(f = "ModifierSet",
          signature = c(x = "Modifier"),
          function(className, x, ...) {
            .Modifer_to_ModifierSet(className, x, ...)
          })

# show -------------------------------------------------------------------------

#' @rdname Modifier-functions
setMethod(
  f = "show", 
  signature = signature(object = "ModifierSet"),
  definition = function(object) {
    callNextMethod()
    cat("| Modification type(s): ",paste0(object[[1]]@mod, collapse = " / "))
    mf <- lapply(seq_along(object),
                 function(i){
                   o <- object[[i]]
                   c(names(object[i]),
                     ifelse(length(o@modifications) != 0L,
                            paste0("yes (",
                                   length(o@modifications),
                                   ")"),
                            "no"))
                 })
    mf <- DataFrame(mf)
    out <- as.matrix(format(as.data.frame(lapply(mf,showAsCell),
                                          optional = TRUE)))
    colnames(out) <- rep(" ",ncol(mf))
    if(is.null(names(object))){
      rownames(out) <- c("| Modifications found:")
    } else {
      rownames(out) <- c("                      ",
                         "| Modifications found:")
    }
    print(out, quote = FALSE, right = TRUE)
    cat("| Settings:\n")
    settings <- lapply(object,
                       function(o){
                         set <- settings(o)
                         set <- lapply(set,
                                       function(s){
                                         if(length(s) > 1L){
                                           ans <- List(s)
                                           return(ans)
                                         }
                                         s
                                       })
                         DataFrame(set)
                       })
    settings <- do.call(rbind,settings)
    rownames(settings) <- names(object)
    .show_settings(settings)
    valid <- unlist(lapply(object,
                           function(o){
                             c(validAggregate(o), validModification(o))
                           }))
    if(!all(valid)){
      warning("Settings were changed after data aggregation or modification ",
              "search. Rerun with modify(x,force = TRUE) to update with ",
              "current settings.", call. = FALSE)
    }
  }
)

# accessors and accessor-like functions ----------------------------------------

#' @rdname Modifier-functions
#' @export
setMethod(f = "bamfiles", 
          signature = signature(x = "ModifierSet"),
          definition = function(x){
            S4Vectors::SimpleList(lapply(x, bamfiles))
          }
)
#' @rdname Modifier-functions
#' @export
setMethod(f = "conditions", 
          signature = signature(object = "ModifierSet"),
          definition = function(object){
            ans <- S4Vectors::SimpleList(lapply(object,conditions))
            ans
          })
#' @rdname Modifier-functions
#' @export
setMethod(f = "mainScore", 
          signature = signature(x = "ModifierSet"),
          definition = function(x) mainScore(new(elementType(x)))
)
#' @rdname Modifier-functions
#' @export
setMethod(f = "modifications", 
          signature = signature(x = "ModifierSet"),
          definition = function(x, perTranscript = FALSE) {
            GenomicRanges::GRangesList(lapply(x,modifications,perTranscript))
          }
)
#' @rdname Modifier-functions
#' @export
setMethod(f = "modifierType", 
          signature = signature(x = "ModifierSet"),
          definition = function(x) modifierType(new(elementType(x)))
)
#' @rdname Modifier-functions
#' @export
setMethod(f = "modType", 
          signature = signature(x = "ModifierSet"),
          definition = function(x) modType(new(elementType(x)))
)
#' @rdname Modifier-functions
#' @export
setMethod(f = "ranges", 
          signature = signature(x = "ModifierSet"),
          definition = function(x){
            ranges(x[[1]])
          }
)
#' @rdname Modifier-functions
#' @export
setMethod(f = "replicates", 
          signature = signature(x = "ModifierSet"),
          definition = function(x){
            ans <- S4Vectors::SimpleList(lapply(x,replicates))
            ans
          })
#' @rdname Modifier-functions
#' @export
setMethod(f = "seqinfo", 
          signature = signature(x = "ModifierSet"),
          definition = function(x){
            S4Vectors::SimpleList(lapply(x, seqinfo))
          }
)
#' @rdname Modifier-functions
#' @export
setMethod(f = "settings", 
          signature = signature(x = "ModifierSet"),
          definition = function(x, name){
            ans <- lapply(x,settings,name)
            names(ans) <- names(x)
            ans
          }
)
#' @rdname Modifier-functions
#' @export
setReplaceMethod(f = "settings", 
                 signature = signature(x = "ModifierSet"),
                 definition = function(x, value){
                   for(i in seq_along(x)){
                     settings(x[[i]]) <- value
                   }
                   x
                 })
#' @rdname Modifier-functions
#' @export
setMethod(f = "sequences", 
          signature = signature(x = "ModifierSet"),
          definition = 
            function(x, modified = FALSE){
              if(!assertive::is_a_bool(modified)){
                stop("'modified' has to be a single logical value.",
                     call. = FALSE)
              }
              if(modified == FALSE){
                return(sequences(sequenceData(x[[1]])))
              }
              mod <- .get_modifications_per_transcript(x)
              mod <- .rebase_seqnames(mod, mod$Parent)
              mod <- split(mod,factor(mod$Parent, levels = mod$Parent))
              ans <- ModRNAStringSet(sequences(sequenceData(x[[1]])))
              modSeqList <- ans[names(ans) %in% names(mod)]
              mod <- mod[match(names(mod),names(modSeqList))]
              ans[names(ans) %in% names(mod)] <- 
                Modstrings::combineIntoModstrings(modSeqList, mod)
              ans
            }
)

# aggregate/modify -------------------------------------------------------------

#' @rdname aggregate
#' @export
setMethod(f = "aggregate", 
          signature = signature(x = "ModifierSet"),
          definition = 
            function(x, force = FALSE){
              ans <- BiocParallel::bplapply(x,aggregate,force)
              ans <- RNAmodR::ModifierSet(class(x),ans)
              ans
            }
)
#' @rdname modify
#' @export
setMethod(f = "modify", 
          signature = signature(x = "ModifierSet"),
          definition = 
            function(x, force = FALSE){
              ans <- BiocParallel::bplapply(x,modify,force)
              ans <- RNAmodR::ModifierSet(class(x),ans)
              ans
            }
)
