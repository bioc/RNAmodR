#' @include RNAmodR.R
#' @include Modifier-class.R
NULL

#' @name ModifierSet
#' 
#' @title ModifierSet
#' @description 
#' title
#' 
#' @param x the input which can be of the following types
#' \itemize{
#' \item{\code{Modifier}}{a single \code{Modifier} or a list containg only 
#' \code{Modifier} objects. The input will just be used as elements of the
#' \code{ModifierSet}}
#' \item{\code{BamFileList}}{a named \code{BamFileList} or a list of 
#' named \code{BamFileList}}
#' \item{\code{list}}{a list of one or more types of elements: 
#' \code{BamFileList}, a named \code{list} or named \code{character} vector. All
#' elements must be or be coercible to a named \code{BamFileList} referencing 
#' existing bam files. In case of a \code{character} vector, it is assumed that
#' each element should be used for creation of a seperate \code{Modifier}
#' object.}
#' }
#' @param fasta sequences matching the target sequences the reads were mapped 
#' onto. This must match the information contained in the BAM files. This is 
#' parameter is only required if \code{x} if not a \code{Modifier} object.
#' @param gff annotation data, which must match the information contained in the
#' BAM files. This is parameter is only required if \code{x} if not a 
#' \code{Modifier} object.
#' @param ... Additional otpional parameters:
#' \itemize{
#' \item{internalBP}{\code{TRUE} or \code{FALSE}: should 
#' parallilazation used internally for creation of each \code{Modifier} or
#' should the creation of the \code{Modifier} objects be parallalized? (default:
#' \code{internalBP = FALSE})}
#' }
#' All other arguments will be passed onto the \code{Modifier} objects.
#' 
NULL

#' @rdname ModifierSet
#' @export
setClass("ModifierSet",
         contains = c("VIRTUAL",
                      "SimpleList"),
         prototype = list(elementType = "Modifier"))


setMethod("pcompareRecursively", "ModifierSet", function(x) FALSE)

setMethod(
  f = "initialize", 
  signature = signature(.Object = "ModifierSet"),
  definition = function(.Object,
                        ...) {
    callNextMethod(.Object,
                   ...)
  }
)

.valid_ModifierSet <- function(x){
  NULL
}
S4Vectors::setValidity2(Class = "ModifierSet",.valid_ModifierSet)

# not supported functions ------------------------------------------------------

setMethod(f = "relistToClass",
          signature = c(x = "ModifierSet"),
          function(x) {
            browser()
            NULL
          })

# contructor -------------------------------------------------------------------

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
  x <- try(.norm_modifiertype(classNames), silent = TRUE)
  if (is(x, "try-error")){
    return(FALSE)
  }
  TRUE
}

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
  x <- unlist(x)
  x[classNames %in% c("BamFileList","BamFile")] <- 
    lapply(x[classNames %in% c("BamFileList","BamFile")],
           path)
  if(!all(file.exists(unlist(x)))){
    return(FALSE)
  }
  TRUE
}

.bamfiles_to_ModifierSet <- function(modifiertype,
                                     x,
                                     fasta,
                                     gff,
                                     modifications,
                                     ...){
  # check and normalize input
  args <- .norm_ModifierSet_args(list(...))
  fasta <- .norm_fasta(fasta, modifiertype)
  gff <- .norm_gff(gff, modifiertype)
  modifiertype <- .norm_modifiertype(modifiertype)
  if(!is.list(x)){
    x <- list(x)
  }
  names <- names(x)
  if(is.null(names)){
    names <- vector(mode = "list", length = length(x))
  }
  x <- lapply(x,
              .norm_bamfiles,
              modifiertype)
  ni <- seq_along(x)
  # choose were to use parallelization
  if(args[["internalBP"]] == TRUE){
    BiocParallel::register(BiocParallel::SerialParam())
  }
  # do analysis by calling the Modifier classes
  FUN <- function(i,
                  z,
                  n,
                  args,
                  modifiertype,
                  PACKAGE){
    suppressPackageStartupMessages({
      requireNamespace(PACKAGE)
    })
    if(!is.null(n)){
      message(i,". ",modifiertype," analysis '",n,"':")
    } else {
      message(i,". ",modifiertype," analysis:")
    }
    # choose were to use parallelization
    if(args[["internalBP"]] == FALSE){
      BiocParallel::register(BiocParallel::SerialParam())
    }
    # do not pass this argument along to objects
    args[["internalBP"]] <- NULL
    #
    if(is(x,"BamFileList")){
      METHOD <- selectMethod(modifiertype,"BamFileList")
    } else {
      METHOD <- selectMethod(modifiertype,"character")
    }
    do.call(METHOD,
            c(list(z,
                   fasta = fasta,
                   gff = gff,
                   modifications = modifications),
              args))
  }
  PACKAGE <- getClass(modifiertype)@package
  x <- BiocParallel::bpmapply(
    FUN,
    ni,
    x,
    names,
    MoreArgs = list(args = args,
                    modifiertype = modifiertype,
                    PACKAGE = PACKAGE),
    SIMPLIFY = FALSE)
  if(!all(vapply(names,is.null,logical(1)))){
    names(x) <- names
  }
  # pass results to ModifierSet object
  new2(.get_class_name_for_set_from_modifier_type(modifiertype),
       listData = x)
}

.modifer_to_ModifierSet <- function(x, ...){
  modifiertype <- modifierType(x[[1]])
  browser()
  new2(.get_class_name_for_set_from_modifier_type(modifiertype),
       listData = x)
}

#' @rdname ModifierSet
#' @export
setMethod(f = "ModifierSet",
          signature = c(x = "list"),
          function(modifiertype, x, ...) {
            if(.contains_only_Modifier(x)){
              return(.modifer_to_ModifierSet(x, ...))
            }
            if(.contains_only_bamfiles(x)){
              args <- list(...)
              return(.bamfiles_to_ModifierSet(modifiertype, 
                                              x, 
                                              args[["fasta"]], 
                                              args[["gff"]], 
                                              ...))
            }
            stop("'x' must be a list containining only elements of the same ",
                 "type\nof 'Modifer' or elements of type ('BamFileList', ",
                 "'character', 'list') which are coercible\nto a named ",
                 "BamFileList. In the latter case, elements must contain named",
                 " vectors or lists('treated' or 'control')\nand the files ",
                 "referenced must exist. Please note, that the list a",
                 call. = FALSE)
          })
#' @rdname ModifierSet
#' @export
setMethod(f = "ModifierSet",
          signature = c(x = "character"),
          function(modifiertype, x, fasta, gff, ...) {
            browser()
            .bamfiles_to_ModifierSet(modifiertype, x, fasta, gff, ...)
          })
#' @rdname ModifierSet
#' @export
setMethod(f = "ModifierSet",
          signature = c(x = "BamFileList"),
          function(modifiertype, x, fasta, gff, ...) {
            browser()
            .bamfiles_to_ModifierSet(modifiertype, x, fasta, gff, ...)
          })
#' @rdname ModifierSet
#' @export
setMethod(f = "ModifierSet",
          signature = c(x = "Modifier"),
          function(modifiertype, x, ...) {
            browser()
            .modifer_to_ModifierSet(x, ...)
          })

# show -------------------------------------------------------------------------

setMethod(
  f = "show", 
  signature = signature(object = "ModifierSet"),
  definition = function(object) {
    callNextMethod()
    cat("| Modification type(s): ",paste0(object[[1]]@mod, collapse = " / "),
        "\n")
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
    out <-
      as.matrix(format(as.data.frame(
        lapply(mf,showAsCell),
        optional = TRUE)))
    colnames(out) <- rep(" ",ncol(mf))
    rownames(out) <- c("| Modifications found:",
                       "                      ")
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
                             c(o@aggregateValidForCurrentArguments,
                               o@modificationsValidForCurrentArguments)
                           }))
    if(!all(valid)){
      warning("Settings were changed after data aggregation or modification ",
              "search. Rerun with modify(x,force = TRUE) to update with ",
              "current settings.", call. = FALSE)
    }
  }
)

# accessors and accessor-like functions ----------------------------------------

#' @name ModifierSet
#' @export
setMethod(f = "modifierType", 
          signature = signature(x = "ModifierSet"),
          definition = function(x) 
            modifierType(new(elementType(x),NULL,NULL,NULL))
)
#' @name ModifierSet
#' @export
setMethod(f = "modType", 
          signature = signature(x = "ModifierSet"),
          definition = function(x) modType(new(elementType(x),NULL,NULL,NULL))
)
#' @name ModifierSet
#' @export
setMethod(f = "mainScore", 
          signature = signature(x = "ModifierSet"),
          definition = function(x) mainScore(new(elementType(x),NULL,NULL,NULL))
)
#' @name ModifierSet
#' @export
setMethod(f = "sequences", 
          signature = signature(x = "ModifierSet"),
          definition = 
            function(x,
                     modified = FALSE,
                     with.qualities = FALSE){
              sequences(x[[1]],
                        modified = modified,
                        with.qualities = with.qualities)
            }
)
#' @name ModifierSet
#' @export
setMethod(f = "ranges", 
          signature = signature(x = "ModifierSet"),
          definition = function(x){
            ranges(seqData(x[[1]]))
          }
)
#' @name ModifierSet
#' @export
setMethod(f = "modifications", 
          signature = signature(x = "ModifierSet"),
          definition = function(x) {
            GRangesList(lapply(x,modifications))
          }
)
