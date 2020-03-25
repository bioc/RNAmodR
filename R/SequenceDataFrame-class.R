#' @include RNAmodR.R
NULL

#' @name SequenceData-functions
#' @aliases show,SequenceDataFrame-method
#' 
#' @title SequenceData/SequenceDataSet/SequenceDataList/SequenceDataFrame 
#' functions
#' 
#' @description 
#' The \code{SequenceData}, \code{SequenceDataSet}, \code{SequenceDataList} and
#' \code{SequenceDataFrame} classes share functionality. Have a look at the 
#' elements listed directly below.
#' 
#' @param x,object a \code{SequenceData}, \code{SequenceDataSet}, 
#' \code{SequenceDataList} or a \code{SequenceDataFrame} object.
#' @param bamfiles a \code{BamFileList}.
#' @param grl a \code{GRangesList} from \code{exonsBy(..., by = "tx")}
#' @param sequences a \code{XStringSet} of type \code{RNAStringSet}, 
#' \code{ModRNAStringSet}, \code{DNAStringSet} or 
#' \code{ModDNAStringSet}
#' @param param a \code{\link[Rsamtools:ScanBamParam-class]{ScanBamParam}} 
#' object
#' @param value a new \code{seqtype}, either "RNA" or "DNA"
#' @param args a list of addition arguments
#' 
#' @return 
#' \itemize{
#' \item{\code{seqinfo}:} {a \code{Seqinfo} object ().}
#' \item{\code{sequences}:} {a \code{RNAStingSet} object or a \code{RNAString} 
#' object for a \code{SequenceDataFrame}.}
#' \item{\code{ranges}:} {a \code{GRangesList} object with each element per 
#' transcript or a \code{GRanges} object for a \code{SequenceDataFrame}.}
#' \item{\code{bamfiles}:} {a \code{BamFileList} object or a SimpleList of 
#' \code{BamFileList} objects for a \code{SequenceDataList}.}
#' }
#' 
#' @examples 
#' data(e5sd,package="RNAmodR")
#' # general accessors
#' seqinfo(e5sd)
#' sequences(e5sd)
#' ranges(e5sd)
#' bamfiles(e5sd)
NULL

#' @name SequenceDataFrame-class
#' @aliases SequenceDataFrame
#' 
#' @title The SequenceDataFrame class
#' 
#' @description 
#' The \code{SequenceDataFrame} class is a virtual class and  contains data for
#' positions along a single transcript. In addition to being used for returning
#' elements from a \code{SequenceData} object, the SequenceDataFrame class is
#' used to store the unlisted data within a
#' \code{\link[=SequenceData-class]{SequenceData}} object. Therefore, a matching
#' \code{SequenceData} and \code{SequenceDataFrame} class must be implemented.
#' 
#' The \code{SequenceDataFrame} class is derived from the
#' \code{\link[S4Vectors:DataFrame-class]{DataFrame}} class. To follow the 
#' functionallity in the \code{S4Vectors} package, \code{SequenceDataFrame} 
#' implements the concept, whereas \code{SequenceDFrame} is the implementation
#' for in-memory data representation from which some specific 
#' \code{*SequenceDataFrame} class derive from, e.g. 
#' \code{\link[=CoverageSequenceData-class]{CoverageSequenceData}}.
#' 
#' Subsetting of a \code{SequenceDataFrame} returns a \code{SequenceDataFrame} or 
#' \code{DataFrame}, if it is subset by a column or row, respectively. The 
#' \code{drop} argument is ignored for column subsetting.
#'
#' @param x,i,j,...,drop,deparse.level arguments used for 
#' \code{\link[S4Vectors:DataFrame-class]{subsetting}} or 
#' \code{\link[base:cbind]{base::cbind}}.
#' 
#' @seealso for an example see
#' \code{\link[=ProtectedEndSequenceData-class]{ProtectedEndSequenceData}}
#' and for more information see \code{\link[=SequenceData-class]{SequenceData}}
#' 
#' @slot ranges a \code{\link[GenomicRanges:GRanges-class]{GRanges}} 
#' object each element describing a transcript including its element. The 
#' \code{GRanges} is constructed from the unlisted results of the
#' \code{\link[GenomicFeatures:transcriptsBy]{exonsBy(x, by="tx")}} function.
#' If during construction a \code{GRangesList} is provided instead of a 
#' character value pointing to a gff3 file or a \code{TxDb} object, it must have
#' a comparable structure. 
#' @slot sequence a \code{\link[Biostrings:XString-class]{XString}} of 
#' type \code{sequencesType} from the parent 
#' \code{\link[=SequenceData-class]{SequenceData}} object.
#' @slot condition conditions along the 
#' \code{\link[Rsamtools:BamFile-class]{BamFileList}}: Either \code{control}
#' or \code{treated}
#' @slot replicate replicate number along the \code{BamFileList} for each of the
#' condition types.
#' @slot bamfiles the input bam files as 
#' \code{\link[Rsamtools:BamFile-class]{BamFileList}}
#' @slot seqinfo a \code{\link[GenomeInfoDb:Seqinfo-class]{Seqinfo}} describing
#' the avialable/used chromosomes.
#' 
#' 
#' @return A \code{SequenceDataFrame} object or if subset to row a 
#' \code{DataFrame}
#'
#' @examples 
#' data(e5sd,package="RNAmodR")
#' # A SequenceDataFrame can is usually constructed by subsetting from 
#' # a SequenceData object
#' sdf <- e5sd[[1]]
#' # Its also used to store the unlisted data in a SequenceData object
#' sdf <- unlist(e5sd) # should probably only used internally
#' e5sd <- relist(sdf,e5sd)
NULL

################################################################################
# SequenceDataFrame (virtual, concept)
################################################################################

#' @rdname SequenceDataFrame-class
#' @export
setClass(Class = "SequenceDataFrame",
         contains = c("VIRTUAL","DataFrame"),
         slots = c(ranges = "GRanges",
                   sequence = "XString",
                   condition = "factor",
                   replicate = "factor",
                   bamfiles = "BamFileList",
                   seqinfo = "Seqinfo"),
         prototype = list(ranges = GRanges(),
                          sequence = RNAString(),
                          condition = factor(),
                          replicate = factor(),
                          bamfiles = Rsamtools::BamFileList(),
                          seqinfo = GenomeInfoDb::Seqinfo()))

setMethod("relistToClass", "SequenceDataFrame",
          function(x) gsub("DataFrame","Data",class(x))
)

# show -------------------------------------------------------------------------

#' @rdname SequenceData-functions
setMethod("show", "SequenceDataFrame",
          function(object){
            callNextMethod(object)
            cat("\ncontaining a ")
            show(object@ranges)
            cat("\nand a ")
            show(object@sequence)
          })

# accessors --------------------------------------------------------------------

#' @rdname SequenceData-functions
#' @export
setMethod(
  f = "conditions", 
  signature = signature(object = "SequenceDataFrame"),
  definition = function(object){object@condition})

#' @rdname SequenceData-functions
#' @export
setMethod(
  f = "bamfiles", 
  signature = signature(x = "SequenceDataFrame"),
  definition = function(x){x@bamfiles})

#' @rdname SequenceData-functions
#' @export
setMethod(f = "dataType",
          signature = signature(x = "SequenceDataFrame"),
          definition = function(x){gsub("SequenceDataFrame","",class(x))})

#' @rdname SequenceData-functions
#' @export
setMethod(
  f = "ranges", 
  signature = signature(x = "SequenceDataFrame"),
  definition = function(x){x@ranges})

#' @rdname SequenceData-functions
#' @export
setMethod(
  f = "replicates", 
  signature = signature(x = "SequenceDataFrame"),
  definition = function(x){x@replicate})

#' @rdname SequenceData-functions
#' @export
setMethod(
  f = "seqinfo", 
  signature = signature(x = "SequenceDataFrame"),
  definition = function(x){x@seqinfo})

#' @rdname SequenceData-functions
#' @export
setMethod(
  f = "seqinfo", 
  signature = signature(x = "SequenceDataFrame"),
  definition = function(x){x@seqinfo})

#' @rdname SequenceData-functions
#' @export
setMethod(
  f = "seqtype", 
  signature = signature(x = "SequenceDataFrame"),
  definition = function(x){seqtype(sequences(x))})

#' @rdname SequenceData-functions
#' @export
setReplaceMethod(
  f = "seqtype", 
  signature = signature(x = "SequenceDataFrame"),
  definition = function(x, value){
    if(!(value %in% c(seqtype(DNAString()),seqtype(RNAString())))){
      stop("Invalid new seqtype.")
    }
    seqtype(x@sequence) <- value
    x
  }
)

#' @rdname SequenceData-functions
#' @export
setMethod(
  f = "sequences", 
  signature = signature(x = "SequenceDataFrame"),
  definition = function(x){x@sequence})

# internals for SequenceDataFrame ----------------------------------------------

#' @importClassesFrom IRanges PartitioningByEnd
#' @importFrom IRanges PartitioningByEnd
setMethod(
  "extractROWS", "SequenceDataFrame",
  function(x, i){
    i <- normalizeSingleBracketSubscript(i, x, exact = FALSE, 
                                         allow.NAs = TRUE, as.NSBS = TRUE)
    start <- which(start(PartitioningByWidth(ranges(x))) == i@subscript[[1L]])
    end <- which(end(PartitioningByWidth(ranges(x))) == i@subscript[[2L]])
    x_ranges <- extractROWS(ranges(x), seq.int(start,end))
    x_sequences <- extractROWS(sequences(x), i)
    # save the other slots, in case they are deleted from the result by calling
    # callNextMethod()
    cl <- class(x)
    x_condition <- conditions(x)
    x_replicate <- replicates(x)
    x_bamfiles <- bamfiles(x)
    x_seqinfo <- seqinfo(x)
    x <- callNextMethod()
    if(!is(x,"SequenceDataFrame")){
      x <- new(cl,
               x,
               ranges = x_ranges,
               sequence = x_sequences,
               condition = x_condition,
               replicate = x_replicate,
               bamfiles = x_bamfiles,
               seqinfo = x_seqinfo)
    } else {
      slot(x, "ranges", check = FALSE) <- x_ranges
      slot(x, "sequence", check = FALSE) <- x_sequences
      validObject(x)
    }
    x
  }
)

setMethod(
  "bindROWS", "SequenceDataFrame",
  function (x, objects = list(), use.names = TRUE, ignore.mcols = FALSE, 
            check = TRUE) 
  {
    objects <- S4Vectors:::prepare_objects_to_bind(x, objects)
    all_objects <- c(list(x), objects)

    ## Call bindROWS() method for DataFrame objects.
    ## Note that the resulting 'ans' might temporarily be an invalid
    ## SequenceDataFrame object (data length and ranges width won't match)
    ## until we update its 'ranges' and 'sequence' slots below.
    ## So we must use 'check=FALSE' to skip validation.
    ans <- callNextMethod(x, objects, use.names = use.names,
                                      ignore.mcols = ignore.mcols,
                                      check = FALSE)

    ## Take care of the 'ranges' and 'sequence' slots.
    ans_ranges <- unlist(GenomicRanges::GRangesList(lapply(all_objects,ranges)))
    ans_sequence <- do.call(xscat,lapply(all_objects,sequences))
    BiocGenerics:::replaceSlots(ans, ranges = ans_ranges,
                                     sequence = ans_sequence,
                                     check = check)
  }
)

#' @rdname SequenceDataFrame-class
#' @export
setMethod(
  "cbind", "SequenceDataFrame",
  function(...){
    args <- list(...)
    if(length(args) == 1L){
      return(args[[1L]])
    }
    # input checks
    classes <- lapply(args,class)
    if(length(unique(classes)) != 1L){
      stop("Inputs must be of the same SequenceDataFrame type.")
    }
    className <- unique(classes)
    lengths <- vapply(args,function(a){sum(lengths(a))},integer(1))
    if(length(unique(lengths)) != 1L){
      stop("Inputs must have the same length.")
    }
    .check_ranges(args)
    .check_sequences(args)
    #
    data <- do.call(cbind,
                    lapply(args,function(a){
                      as(a, .get_first_class_extends(a))
                    }))
    ranges <- ranges(args[[1L]])
    sequences <- sequences(args[[1L]])
    colnames <- IRanges::CharacterList(strsplit(colnames(data),"\\."))
    colnames_conditions <- colnames %in% c("treated","control")
    colnames_replicates <- 
      !is.na(suppressWarnings(IRanges::IntegerList(colnames)))
    colnames_f <- !(colnames_conditions | colnames_replicates)
    conditionsFmultiplier <- length(unique(vapply(colnames[colnames_f],
                                                  paste,character(1),
                                                  collapse=".")))
    condition <- unlist(lapply(args,conditions))
    condition_steps <- seq.int(1,length(condition),by=conditionsFmultiplier)
    replicate <- .get_replicate_number(condition[condition_steps])
    replicate <- rep(replicate, each = conditionsFmultiplier)
    colnames[colnames_conditions] <- IRanges::CharacterList(condition)
    colnames[colnames_replicates] <- IRanges::CharacterList(replicate)
    colnames(data) <- vapply(colnames,paste,character(1),collapse = ".")
    bamfiles <- do.call(c,lapply(args,bamfiles))
    seqinfo <- seqinfo(args[[1L]])
    .SequenceDataFrame(class = gsub("SequenceDataFrame","",className),
                       df = data,
                       ranges = ranges,
                       sequence = sequences,
                       replicate = replicate,
                       condition = condition,
                       bamfiles = bamfiles,
                       seqinfo = seqinfo)
  }
)

#' @importFrom stats setNames
#' @rdname SequenceDataFrame-class
#' @export
setMethod(
  "[", "SequenceDataFrame",
  function(x, i, j, ..., drop = TRUE){
    if (!isTRUEorFALSE(drop)){
      stop("'drop' must be TRUE or FALSE")
    }
    if (length(list(...)) > 0L){
      warning("parameters in '...' not supported")
    }
    classDG <- .get_first_class_extends(x)
    ## We do list-style subsetting when [ was called with no ','.
    ## NOTE: matrix-style subsetting by logical matrix not supported.
    list_style_subsetting <- (nargs() - !missing(drop)) < 3L
    if (list_style_subsetting || !missing(j)) {
      if (list_style_subsetting) {
        if (!missing(drop))
          warning("'drop' argument ignored by list-style subsetting")
        if (missing(i))
          return(x)
        j <- i
      }
      xstub <- stats::setNames(seq_along(x), names(x))
      ia <- interaction(conditions(x), replicates(x))
      if(is.character(j)){
        j <- normalizeSingleBracketSubscript(j, xstub)
        j <- unique(as.integer(ia)[j])
      } else {
        conditionsFmultiplier <- length(ia) / length(unique(ia))
        j <- normalizeSingleBracketSubscript(j, xstub[seq_len(length(ia)/conditionsFmultiplier)])
      }
      j2 <- which(!is.na(match(as.integer(ia), j)))
      x <- initialize(x,
                      as(x,classDG)[, j2, drop = FALSE],
                      ranges = x@ranges,
                      sequence = x@sequence,
                      replicate = factor(x@replicate[j2]),
                      condition = factor(x@condition[j2]),
                      bamfiles = x@bamfiles[j],
                      seqinfo = x@seqinfo)
      if (anyDuplicated(names(x))){
        names(x) <- make.unique(names(x))
      }
      if (list_style_subsetting){
        return(x)
      }
    }
    if (!missing(i)){
      x <- extractROWS(as(x,classDG), i)
    } else {
      return(x) # early exit if subset is column-only
    }
    if (missing(drop)){
      drop <- TRUE
    }  
    if (drop) {
      ## one row left
      if (nrow(x) == 1L){
        return(as(x, "list"))
      }
    }
    x
  }
)

# constructor ------------------------------------------------------------------

# class names must be compatible with this class name generation function
sequenceDataFrameClass <- function(dataType){
  ans <- paste0(dataType,"SequenceDataFrame")
  tmp <- try(getClass(ans))
  if(is(tmp,"try-error")){
    stop("Class '",ans,"' not found: ",tmp)
  }
  ans
}

.SequenceDataFrame <- function(class, df, ranges, sequence, replicate,
                               condition, bamfiles, seqinfo){
  # defaults from function are strangly not set
  if(missing(df)){
    df <- DataFrame()
  }
  if(missing(ranges)){
    ranges <- GRanges()
  }
  if(missing(sequence)){
    sequence <- RNAString()
  }
  if(missing(replicate)){
    replicate <- factor()
  }
  if(missing(condition)){
    condition <- factor()
  }
  if(missing(bamfiles)){
    bamfiles <- Rsamtools::BamFileList()
  }
  if(missing(seqinfo)){
    seqinfo <- GenomeInfoDb::Seqinfo()
  }
  # check inputs
  if(!is(df,"DataFrame")){
    stop("Invalid data object: ", class(df), " found, DataFrame expected.")
  }
  if(ncol(df) != length(replicate) ||
     ncol(df) != length(condition)){
    stop("Replicate and Conditions information must match the DataFrame ",
         "dimensions.")
  }
  if(!is(ranges,"GRanges")){
    stop("Invalid data object: ", class(ranges), " found, GRanges expected.")
  }
  if(!is(sequence,"XString")){
    stop("Invalid data object: ", class(sequence), " found, XString expected.")
  }
  new(paste0(class,"SequenceDataFrame"),
      ranges = ranges,
      sequence = sequence,
      condition = condition,
      replicate = replicate,
      bamfiles = bamfiles,
      seqinfo = seqinfo,
      df)
}

.valid_SequenceDataFrame <-  function(x){
  if(nrow(x) != sum(width(ranges(x)))){
    return("data length and ranges width do not match.")
  }
  if(nrow(x) != length(sequences(x))){
    return("data length and sequence length do not match.")
  }
  if(!is(sequences(x),"RNAString") && !is(sequences(x),"DNAString")){
    stop("")
  }
  S4Vectors:::.valid.DataFrame(x)
  NULL
}

S4Vectors::setValidity2(Class = "SequenceDataFrame", .valid_SequenceDataFrame)

################################################################################
# SequenceDFrame (Virtual, implementation, type independent)
################################################################################

#' @rdname SequenceDataFrame-class
#' @export
setClass(Class = "SequenceDFrame",
         contains = c("VIRTUAL","SequenceDataFrame","DFrame"))
