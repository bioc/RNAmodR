#' @include RNAmodR.R
#' @include SequenceData-class.R
#' @include SequenceDataSet-class.R
NULL

#' @name SequenceDataList-class
#' @aliases SequenceDataList
#' 
#' @title The SequenceDataList class
#' 
#' @description 
#' The \code{SequenceDataList} class is used to hold \code{SequenceData} or 
#' \code{SequenceDataSet} objects as its elements. It is derived from the 
#' \code{\link[S4Vectors:List-class]{List}} class.
#' 
#' The \code{SequenceDataList} is used to hold data from different sets of
#' aligned reads. This allows multiple methods to be aggregated into one
#' modification detection strategy. Annotation and sequence data must be the 
#' same for all elements, however the bam files can be different. 
#' 
#' @param ... The elements to be included in the \code{SequenceDataList}.
#' 
#' @return a \code{SequenceDataList}
#' 
#' @examples
#' data(psd,package="RNAmodR")
#' data(e5sd,package="RNAmodR")
#' sdl <- SequenceDataList(SequenceDataSet(psd,e5sd),e5sd)
NULL

#' @rdname SequenceDataList-class
#' @export
setClass("SequenceDataList",
         contains = "List",
         slots = c(listData = "list"),
         prototype = list(elementType = "SD_or_SDS"))

setClassUnion("SD_or_SDS_or_SDL",
              c("SequenceData", "SequenceDataSet", "SequenceDataList"))

# show method ------------------------------------------------------------------
#' @rdname SequenceData-functions
setMethod("show", "SequenceDataList",
          function(object)
          {
            lo <- length(object)
            cat(classNameForDisplay(object), " of length ", lo,
                "\n", sep = "")
          })

# parallelSlotNames ------------------------------------------------------------
#' @rdname RNAmodR-internals
setMethod("parallelSlotNames", "SequenceDataList",
          function(x) c("listData", callNextMethod())
)

# constructors -----------------------------------------------------------------

.SequenceDataList <- function(Class, listData, ..., check = FALSE){
  new2(Class, listData = listData, ..., check = check)
}

# not exported. Only used internally
new_SequenceDataList_from_list <- function(Class, x, ..., mcols){
  if (!extends(Class, "SequenceDataList")){
    stop("class ", Class, " must extend SequenceDataList")
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
  if(!.compare_element_metadata(x,"ranges")){
    stop("Annotation data of all SequenceDataList elements are not equal.",
         call. = FALSE)
  }
  if(!.compare_element_metadata(x,"sequences")){
    stop("Sequence data of all SequenceDataList elements are not equal.",
         call. = FALSE)
  }
  # class name as default names
  if(is.null(names(x))){
    names(x) <- vapply(x,class,character(1))
    f <- vapply(x,is,logical(1),"SequenceDataSet")
    if(any(f)){
      names(x)[f] <- vapply(x[f],
                            function(xi){
                              paste0(vapply(xi,class,character(1)),collapse="_")
                            },
                            character(1))
    }
  }
  #
  if (missing(mcols)){
    return(.SequenceDataList(Class, listData = x, ..., check = FALSE))
  }
  .SequenceDataList(Class, listData = x, ..., elementMetadata = mcols,
                   check = FALSE)
}

#' @rdname SequenceDataList-class
#' @export
SequenceDataList <- function(...){
  new_SequenceDataList_from_list("SequenceDataList", list(...))
}

# Validity ---------------------------------------------------------------------
.valid.SequenceDataList.listData <- function(x){
  elementTypeX <- elementType(x)
  if (!all(vapply(as.list(x),
                  function(xi) extends(class(xi), elementTypeX),
                  logical(1)))){
    classes <- getClass("SD_or_SDS")
    if(isClassUnion(classes)){
      classes <- paste(names(classes@subclasses), collapse = " or ")
    } else {
      classes <- classes@className
    }
    return(paste("the 'listData' slot must be a list containing ",
                 classes, " objects"))
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

.valid.SequenceDataList <- function(x){
  c(.valid.SequenceDataList.listData(x),
    unlist(lapply(x,validObject)))
}
S4Vectors::setValidity2("SequenceDataList", .valid.SequenceDataList)

# classNameForDisplay ----------------------------------------------------------
setMethod("classNameForDisplay", "SequenceDataList",
          function(x) "SequenceDataList")


# Subsetting -------------------------------------------------------------------
#' @rdname RNAmodR-internals
setMethod("getListElement", "SequenceDataList",
          function(x, i, exact = TRUE)
            getListElement(x@listData, i, exact = exact)
)

# looping ----------------------------------------------------------------------
setMethod("lapply", "SequenceDataList",
          function(X, FUN, ...) lapply(as.list(X), match.fun(FUN), ...)
)

# coercion ---------------------------------------------------------------------
# forth and ...
setAs("SequenceDataList", "list", function(from) as.list(from))
.as.list.SequenceDataList <- function(x, use.names = TRUE)
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
setMethod("as.list", "SequenceDataList", .as.list.SequenceDataList)

# ... back
setAs("list", "SequenceDataList",
      function(from){
        new_SequenceDataList_from_list("SequenceDataList", from)
      })
setAs("ANY", "SequenceDataList",
      function(from) {
        as(as.list(from),"SequenceDataList")
      })

# additional accessors ---------------------------------------------------------

#' @rdname SequenceData-functions
#' @export
setMethod(f = "bamfiles", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            ans <- do.call(S4Vectors::SimpleList, lapply(x, bamfiles))
            names(ans) <- names(x@listData)
            ans
          })
#' @rdname SequenceData-functions
#' @export
setMethod(f = "conditions", 
          signature = signature(object = "SequenceDataList"),
          definition = function(object){
            ans <- S4Vectors::SimpleList(lapply(object,conditions))
            ans
          })
#' @rdname SequenceData-functions
#' @export
setMethod(f = "names", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            names(x[[1]])
          })
#' @rdname SequenceData-functions
#' @export
setMethod(f = "ranges", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            ranges(x[[1]])
          })
#' @rdname SequenceData-functions
#' @export
setMethod(f = "replicates", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            ans <- S4Vectors::SimpleList(lapply(x,replicates))
            ans
          })
#' @rdname SequenceData-functions
#' @export
setMethod(f = "seqinfo", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            seqinfo(x[[1]])
          })
#' @rdname SequenceData-functions
#' @export
setMethod(f = "sequences", 
          signature = signature(x = "SequenceDataList"),
          definition = function(x){
            sequences(x[[1]])
          })

# aggregate --------------------------------------------------------------------

#' @rdname aggregate
#' @export
setMethod("aggregate",
          signature = c(x = "SequenceDataList"),
          function(x, condition = "Treated"){
            ans <- do.call(S4Vectors::SimpleList,
                           lapply(x, aggregate, condition = condition))
            names(ans) <- names(x@listData)
            ans
          })
