#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name SequenceDataSet-class
#' @aliases SequenceDataSet
#' 
#' @title The SequenceDataSet class
#' 
#' @description 
#' The \code{SequenceDataSet} class is used to hold \code{SequenceData} objects
#' as its elements. It is derived from the 
#' \code{\link[S4Vectors:List-class]{List}} class.
#' 
#' The \code{SequenceDataSet} is used to hold different data types from the of
#' same aligned reads. The same dataset can be used to generate multiple sets of
#' data types. Bam files, annotation and sequence data must be the same for all 
#' elements. 
#' 
#' @param ... The elements to be included in the \code{SequenceDataSet}.
#' 
#' @return a \code{SequenceDataSet}
#' 
#' @examples
#' data(psd,package="RNAmodR")
#' data(e5sd,package="RNAmodR")
#' sdl <- SequenceDataSet(psd,e5sd)
NULL

#' @rdname SequenceDataSet-class
#' @export
setClass("SequenceDataSet",
         contains = "List",
         slots = c(listData = "list"),
         prototype = list(elementType = "SequenceData"))

setClassUnion("SD_or_SDS",
              c("SequenceData", "SequenceDataSet"))

# show method ------------------------------------------------------------------
#' @rdname SequenceData-functions
setMethod("show", "SequenceDataSet",
          function(object)
          {
            lo <- length(object)
            cat(classNameForDisplay(object), " of length ", lo,
                "\n", sep = "")
            if (!is.null(names(object)))
              cat(S4Vectors:::labeledLine("names", names(object)))
            unlisted_ranges <- unlist(ranges(object[[1]]),use.names = FALSE)
            ranges_mcols <- mcols(unlisted_ranges, use.names = FALSE)
            nhead <- S4Vectors::get_showHeadLines()
            ntail <- S4Vectors::get_showTailLines()
            nc <- if (is.null(ranges_mcols)) 0L else ncol(ranges_mcols)
            nr <- if (is.null(ranges_mcols)) 0L else nrow(ranges_mcols)
            nms <- rownames(ranges_mcols)
            if (nr <= (nhead + ntail + 1L)) {
              out <-
                as.matrix(format(as.data.frame(
                  lapply(ranges_mcols, showAsCell),
                  optional = TRUE)))
            } else {
              out <-
                rbind(as.matrix(format(as.data.frame(
                  lapply(ranges_mcols, function(x)
                    showAsCell(head(x, nhead))),
                  optional = TRUE))),
                  rbind(rep.int("...", nc)),
                  as.matrix(format(as.data.frame(
                    lapply(ranges_mcols, function(x) 
                      showAsCell(tail(x, ntail))),
                    optional = TRUE))))
              rownames(out) <- make_rownames_for_RectangularData_display(
                                                 nms, nr,
                                                 nhead, ntail)
            }
            classinfo <-
              matrix(unlist(lapply(ranges_mcols, function(x) {
                paste0("<", classNameForDisplay(x)[1],
                       ">")
              }), use.names = FALSE), nrow = 1,
              dimnames = list("", colnames(out)))
            out <- rbind(classinfo, out)
            cat("- Ranges metadata columns:\n")
            print(out, quote = FALSE, right = TRUE)
          })

# parallel_slot_names ----------------------------------------------------------
#' @rdname RNAmodR-internals
setMethod("parallel_slot_names", "SequenceDataSet",
          function(x) c("listData", callNextMethod())
)

# constructors -----------------------------------------------------------------

.SequenceDataSet <- function(Class, listData, ..., check = FALSE){
  new2(Class, listData = listData, ..., check = check)
}

.compare_element_metadata <- function(input, FUN){
  if(length(input) == 1L){
    return(TRUE)
  }
  input <- lapply(input,FUN)
  first_input <- input[[1]]
  if(is(first_input,"BamFileList")){
    ans <- vapply(input[seq.int(2,length(input))],
                  function(i){
                    all(path(first_input) == path(i))
                  },
                  logical(1))
  } else {
    first_input <- unlist(first_input)
    ans <- vapply(input[seq.int(2,length(input))],
                  function(i){
                    all(all(first_input == unlist(i)))
                  },
                  logical(1))
  }
  all(ans)
}

# not exported. Only used internally
new_SequenceDataSet_from_list <- function(Class, x, ..., mcols){
  if (!extends(Class, "SequenceDataSet")){
    stop("class ", Class, " must extend SequenceDataSet")
  }
  if (!is.list(x)){
    stop("'x' must be a list")
  }
  proto <- new(Class)
  ans_elementType <- elementType(proto)
  if (is(S4Vectors::mcols(proto, use.names = FALSE), "DataFrame")){
    mcols <- S4Vectors:::make_zero_col_DataFrame(length(x))
  }
  extends_elementType <- vapply(x,
                                function(xi){
                                  extends(class(xi), ans_elementType)
                                },
                                logical(1))
  if (!all(extends_elementType)){
    stop("All elements in 'x' must be ", ans_elementType, " objects")
  }
  # check that all sequences and annotation information are the same
  if(!.compare_element_metadata(x,"bamfiles")){
    stop("Annotation data of all SequenceData elements are not equal.",
         call. = FALSE)
  }
  if(!.compare_element_metadata(x,"ranges")){
    stop("Annotation data of all SequenceData elements are not equal.",
         call. = FALSE)
  }
  if(!.compare_element_metadata(x,"sequences")){
    stop("Sequence data of all SequenceData elements are not equal.",
         call. = FALSE)
  }
  # class name as default names
  if(is.null(names(x))){
    names(x) <- vapply(x,class,character(1))
  }
  #
  if (missing(mcols)){
    return(.SequenceDataSet(Class, listData = x, ..., check = FALSE))
  }
  .SequenceDataSet(Class, listData = x, ..., elementMetadata = mcols,
                    check = FALSE)
}

#' @rdname SequenceDataSet-class
#' @export
SequenceDataSet <- function(...){
  new_SequenceDataSet_from_list("SequenceDataSet", list(...))
}

# Validity ---------------------------------------------------------------------
.valid.SequenceDataSet.listData <- function(x){
  elementTypeX <- elementType(x)
  if (!all(vapply(as.list(x),
                  function(xi) extends(class(xi), elementTypeX),
                  logical(1)))){
    return(paste("the 'listData' slot must be a list containing ",
                 elementTypeX, " objects"))
  }
  if(!.compare_element_metadata(x,"bamfiles")){
    return("Bam input files are not equal.")
  }
  if(!.compare_element_metadata(x,"ranges")){
    return("Annotation data is not equal.")
  }
  if(!.compare_element_metadata(x,"sequences")){
    return("Sequence data is not equal.")
  }
  NULL
}
.valid.SequenceDataSet <- function(x){
  c(.valid.SequenceDataSet.listData(x),
    unlist(lapply(x,validObject)))
}
S4Vectors::setValidity2("SequenceDataSet", .valid.SequenceDataSet)


# classNameForDisplay ----------------------------------------------------------
setMethod("classNameForDisplay", "SequenceDataSet",
          function(x) "SequenceDataSet")


# Subsetting -------------------------------------------------------------------
#' @rdname RNAmodR-internals
setMethod("getListElement", "SequenceDataSet",
          function(x, i, exact = TRUE)
            getListElement(x@listData, i, exact = exact)
)

# looping ----------------------------------------------------------------------
setMethod("lapply", "SequenceDataSet",
          function(X, FUN, ...) lapply(as.list(X), match.fun(FUN), ...)
)

# coercion ---------------------------------------------------------------------
# forth and ...
setAs("SequenceDataSet", "list", function(from) as.list(from))
.as.list.SequenceDataSet <- function(x, use.names = TRUE)
{
  if (!isTRUEorFALSE(use.names)){
    stop("'use.names' must be TRUE or FALSE")
  }
  ans <- x@listData
  if (!use.names){
    names(ans) <- NULL
  }
  ans
}
setMethod("as.list", "SequenceDataSet", .as.list.SequenceDataSet)

# ... back
setAs("list", "SequenceDataSet",
      function(from){
        new_SequenceDataSet_from_list("SequenceDataSet", from)
      })
setAs("ANY", "SequenceDataSet",
      function(from) {
        as(as.list(from),"SequenceDataSet")
      })

# additional accessors ---------------------------------------------------------

#' @rdname SequenceData-functions
#' @export
setMethod(f = "bamfiles", 
          signature = signature(x = "SequenceDataSet"),
          definition = function(x){
            bamfiles(x[[1L]])
          })

#' @rdname SequenceData-functions
#' @export
setMethod(f = "conditions", 
          signature = signature(object = "SequenceDataSet"),
          definition = function(object){
            ans <- IRanges::FactorList(
              lapply(object[1L],
                     function(o){
                       ia <- as.integer(interaction(conditions(o),
                                                    replicates(o)))
                       m <- match(unique(ia),ia)
                       conditions(o)[m]
                     }))
            ans[[1L]]
          })

#' @rdname SequenceData-functions
#' @export
setMethod(f = "names", 
          signature = signature(x = "SequenceDataSet"),
          definition = function(x){
            names(x[[1L]])
          })

#' @rdname SequenceData-functions
#' @export
setMethod(f = "ranges", 
          signature = signature(x = "SequenceDataSet"),
          definition = function(x){
            ranges(x[[1L]])
          })

#' @rdname SequenceData-functions
#' @export
setMethod(f = "replicates", 
          signature = signature(x = "SequenceDataSet"),
          definition = function(x){
            ans <- IRanges::FactorList(
              lapply(x[1L],
                     function(o){
                       ia <- as.integer(interaction(conditions(o),
                                                    replicates(o)))
                       m <- match(unique(ia),ia)
                       replicates(o)[m]
                     }))
            ans[[1L]]
          })

#' @rdname SequenceData-functions
#' @export
setMethod(f = "seqinfo", 
          signature = signature(x = "SequenceDataSet"),
          definition = function(x){
            seqinfo(x[[1L]])
          })

#' @rdname SequenceData-functions
#' @export
setMethod(f = "seqtype", 
          signature = signature(x = "SequenceDataSet"),
          definition = function(x){seqtype(x[[1L]])})

#' @rdname SequenceData-functions
#' @export
setReplaceMethod(f = "seqtype", 
                 signature = signature(x = "SequenceDataSet"),
                 definition = function(x, value){
                   as(lapply(x,`seqtype<-`,value),"SequenceDataSet")
                 })

#' @rdname SequenceData-functions
#' @export
setMethod(f = "sequences", 
          signature = signature(x = "SequenceDataSet"),
          definition = function(x){
            sequences(x[[1L]])
          })

# aggregate --------------------------------------------------------------------

#' @rdname aggregate
#' @export
setMethod("aggregate",
          signature = c(x = "SequenceDataSet"),
          function(x, condition = "Treated"){
            ans <- do.call(S4Vectors::SimpleList,
                           lapply(x, aggregate, condition = condition))
            names(ans) <- names(x@listData)
            ans
          })
