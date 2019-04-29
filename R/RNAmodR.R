#' @title RNAmodR
#' 
#' @author Felix G M Ernst [aut], Denis L.J. Lafontaine [ctb]
#' 
#' @description
#' Post-transcriptional modifications can be found abundantly in rRNA and tRNA
#' and can be detected classically via several strategies. However, difficulties
#' arise if the identity and the position of the modified nucleotides is to be
#' determined at the same time. Classically, a primer extension, a form of
#' reverse transcription (RT), would allow certain modifications to be accessed
#' by blocks during the RT changes or changes in the cDNA sequences. Other
#' modification would need to be selectively treated by chemical reactions to
#' influence the outcome of the reverse transcription.
#' 
#' With the increased availability of high throughput sequencing, these
#' classical methods were adapted to high throughput methods allowing more RNA
#' molecules to be accessed at the same time. With these advances
#' post-transcriptional modifications were also detected on mRNA. Among these
#' high throughput techniques are for example Pseudo-Seq (Carlile et al. 2014),
#' RiboMethSeq (Birkedal et al. 2015) and AlkAnilineSeq (Marchand et al. 2018)
#' each able to detect a specific type of modification from footprints in
#' RNA-Seq data prepared with the selected methods.
#'     
#' Since similar pattern can be observed from some of these techniques, overlaps
#' of the bioinformatical pipeline already are and will become more frequent
#' with new emerging sequencing techniques.
#' 
#' \code{RNAmodR} implements classes and a workflow to detect
#' post-transcriptional RNA modifications in high throughput sequencing data. It
#' is easily adaptable to new methods and can help during the phase of initial
#' method development as well as more complex screenings.
#' 
#' Briefly, from the \code{SequenceData}, specific subclasses are derived for
#' accessing specific aspects of aligned reads, e.g. 5’-end positions or pileup
#' data. With this a \code{Modifier} class can be used to detect specific
#' patterns for individual types of modifications. The \code{SequenceData}
#' classes can be shared by different \code{Modifier} classes allowing easy
#' adaptation to new methods.
#' 
#' @seealso The \code{RNAmodR.RiboMethSeq} and \code{RNAmodR.AlkAnilineSeq}
#' package.
#'
#' @references 
#' - Carlile TM, Rojas-Duran MF, Zinshteyn B, Shin H, Bartoli KM, Gilbert WV 
#' (2014): "Pseudouridine profiling reveals regulated mRNA pseudouridylation in 
#' yeast and human cells." Nature 515 (7525), P. 143–146. DOI:
#' \href{https://doi.org/10.1038/nature13802}{10.1038/nature13802}.
#' 
#' - Birkedal U, Christensen-Dalsgaard M, Krogh N, Sabarinathan R, Gorodkin J, 
#' Nielsen H (2015): "Profiling of ribose methylations in RNA by high-throughput
#' sequencing." Angewandte Chemie (International ed. in English) 54 (2), 
#' P. 451–455. DOI: 
#' \href{https://doi.org/10.1002/anie.201408362}{10.1002/anie.201408362}.
#' 
#' - Marchand V, Ayadi L, __Ernst FGM__, Hertler J, Bourguignon-Igel V,
#' Galvanin A, Kotter A, Helm M, __Lafontaine DLJ__, Motorin Y (2018): 
#' "AlkAniline-Seq: Profiling of m7 G and m3 C RNA Modifications at Single 
#' Nucleotide Resolution." Angewandte Chemie (International ed. in English) 57 
#' (51), P. 16785–16790. DOI: 
#' \href{https://doi.org/10.1002/anie.201810946}{10.1002/anie.201810946}.
#'
#' @docType package
#' @name RNAmodR
NULL

#' @import methods
#' @import assertive
#' @import XVector
#' @import S4Vectors
#' @import GenomicRanges
#' @importFrom BiocGenerics which conditions
#' @importFrom Biostrings DNAString RNAString DNAStringSet RNAStringSet getSeq
#' @importFrom Modstrings ModRNAString ModRNAStringSet combineIntoModstrings
#' shortName fullName
#' @importClassesFrom IRanges IntegerList CharacterList LogicalList IRanges
#' SplitDataFrameList
NULL

# constants for annotation -----------------------------------------------------

#' @name RNAmodR-internals
#' @aliases .getData
#' 
#' @title RNAmodR internal functions
#' 
#' @description
#' These functions are used internally.
#' 
#' @param object,range,data,modType,scoreFun,source,type,x,i,j,...,exact,value
#' Internal arguments
#' 
#' @keywords internal
#' 
#' @return internally used values
NULL

#' @name RNAmodR-datasets
#' @title Example data in the RNAmodR package
#' @description 
#' The following datasets are contained in the RNAmodR package. They are used
#' in the man page examples.
#' @docType data
#' @usage msi
#' @format 
#' \itemize{
#' \item{msi} {a \code{ModSetInosine} instance}
#' \item{sds} {a \code{SequenceDataSet} instance}
#' \item{sdl} {a \code{SequenceDataList} instance}
#' \item{psd} {a \code{PileupSequenceData} instance}
#' \item{e5sd} {a \code{End5SequenceData} instance}
#' \item{e3sd} {a \code{End3SequenceData} instance}
#' \item{esd} {a \code{EndSequenceData} instance}
#' \item{csd} {a \code{CoverageSequenceData} instance}
#' \item{ne3sd} {a \code{NormEnd3SequenceData} instance}
#' \item{ne5sd} {a \code{NormEnd5SequenceData} instance}
#' \item{pesd} {a \code{ProtectedEndSequenceData} instance}
#' }
#' @keywords datasets
"msi"
#' @rdname RNAmodR-datasets
"sds"
#' @rdname RNAmodR-datasets
"sdl"
#' @rdname RNAmodR-datasets
"psd"
#' @rdname RNAmodR-datasets
"e5sd"
#' @rdname RNAmodR-datasets
"e3sd"
#' @rdname RNAmodR-datasets
"esd"
#' @rdname RNAmodR-datasets
"csd"
#' @rdname RNAmodR-datasets
"ne3sd"
#' @rdname RNAmodR-datasets
"ne5sd"
#' @rdname RNAmodR-datasets
"pesd"

#' @name RNAmodR-development
#' @aliases getData
#' 
#' @title RNAmodR developments functions
#' 
#' @description
#' These functions are not intended for general use, but are used for 
#' additional package development. 
#' 
#' \code{getData} is used to load data into a
#' \code{\link[RNAmodR:SequenceData-class]{SequenceData}} object and must be
#' implented for all \code{SequenceData} classes. The results must match the
#' requirements outlined in the value section.
#' 
#' In addition the following functions should be implemented for complete
#' functionality:
#' 
#' \code{aggregateData} for each \code{SequenceData} and \code{Modifier} class.
#' See also \code{\link[=aggregate]{aggregateData}}
#' 
#' \code{findMod} for each \code{Modifier} class. See also 
#' \code{\link[=modify]{findMod}}.
#' 
#' \code{visualizeData}/\code{visualizeDataByCoord} for each \code{Modifier} 
#' and \code{ModifierSet} class. See also 
#' \code{\link[=visualizeData]{visualizeData}}.
#' 
#' The following helper function can be called from within \code{findMod} to
#' construct a coordinate for each modification found:
#' 
#' \code{constructModRanges} constructs a \code{GRanges} object describing the
#' location, type and associated scores of a modification.
#' \code{constructModRanges} is typically called from the \code{modify}
#' function, which must be implemented for all
#' \code{\link[RNAmodR:Modifier-class]{Modifier}} classes.
#' 
#' @param x for \code{getData}:a \code{SequenceData} object.
#' @param bamfiles for \code{getData}:a \code{BamFileList} object.
#' @param grl for \code{getData}:a \code{GRangesList} object.
#' @param sequences for \code{getData}:a \code{XStringSet} object.
#' @param param for \code{getData}:a \code{ScanBamParam} object.
#' @param args for \code{getData}: a \code{list} with optional arguments.
#' @param range for \code{constructModRanges}: a \code{GRanges} object
#' @param data for \code{constructModRanges}: a \code{DataFrame} object
#' @param modType for \code{constructModRanges}: a valid shortName for the 
#' modification found. Must be present in \code{shortName(ModRNAString())}.
#' @param scoreFun for \code{constructModRanges}: a custom function for 
#' extracting scores from \code{data}. The result must be a \code{list}.
#' @param source for \code{constructModRanges}: a single character vector for
#' populating the source column of the result.
#' @param type for \code{constructModRanges}: a single character vector for
#' populating the source column of the result.
#' 
#' @return 
#' \itemize{
#' \item{\code{getData}:} {returns a list with elements per BamFile in 
#' \code{bamfiles}. Elements can be 
#' \code{\link[IRanges:AtomicList-class]{IntegerList}}, 
#' \code{\link[IRanges:AtomicList-class]{NumericList}} or a
#' \code{\link[IRanges:DataFrameList-class]{CompressedSplitDataFrameList}}. The
#' data in the elements must be order by increasing positions numbers. However, 
#' names and rownames will be discarded.}
#' \item{\code{constructModRanges}:} {returns a \code{GRanges} object with
#' genomic coordinates of modified nucleotides in the associated transcripts.}
#' }
#' 
#' @examples 
#' # new SequenceData class
#' setClass(Class = "ExampleSequenceData",
#'          contains = "SequenceData",
#'          prototype = list(minQuality = 5L))
#' ExampleSequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
#'   RNAmodR:::SequenceData("Example", bamfiles = bamfiles, 
#'                          annotation = annotation, sequences = sequences,
#'                          seqinfo = seqinfo, ...)
#' }
#' setMethod("getData",
#'           signature = c(x = "ExampleSequenceData",
#'                         bamfiles = "BamFileList",
#'                         grl = "GRangesList",
#'                         sequences = "XStringSet",
#'                         param = "ScanBamParam"),
#'           definition = function(x, bamfiles, grl, sequences, param, args){
#'             ###
#'           }
#' )
#' setMethod("aggregateData",
#'           signature = c(x = "ExampleSequenceData"),
#'           function(x, condition = c("Both","Treated","Control")){
#'             ###
#'           }
#' )
#' setMethod(
#'   f = "getDataTrack",
#'   signature = signature(x = "ExampleSequenceData"),
#'   definition = function(x, name, ...) {
#'     ###
#'   }
#' )
#' 
#' # new Modifier class
#' setClass("ModExample",
#'          contains = c("Modifier"),
#'          prototype = list(mod = "X",
#'                           score = "score",
#'                           dataType = "ExampleSequenceData"))
#' ModExample <- function(x, annotation, sequences, seqinfo, ...){
#'   RNAmodR:::Modifier("ModExample", x = x, annotation = annotation,
#'                      sequences = sequences, seqinfo = seqinfo, ...)
#' }
#' 
#' setMethod(f = "aggregateData",
#'           signature = signature(x = "ModExample"),
#'           definition =
#'             function(x, force = FALSE){
#'               # Some data with element per transcript
#'             }
#' )
#' 
#' setMethod("findMod",
#'           signature = c(x = "ModExample"),
#'           function(x){
#'             # an element per modification found.
#'           }
#' )
#' setMethod(
#'   f = "getDataTrack",
#'   signature = signature(x = "ModExample"),
#'   definition = function(x, name, type, ...) {
#'   }
#' )
#' setMethod(
#'   f = "visualizeDataByCoord",
#'   signature = signature(x = "ModExample", coord = "GRanges"),
#'   definition = function(x, coord, type = "score", window.size = 15L, ...) {
#'   }
#' )
#' setMethod(
#'   f = "visualizeData",
#'   signature = signature(x = "ModExample"),
#'   definition = function(x, name, from, to, type = "score", ...) {
#'   }
#' )
#'  
#' # new ModifierSet class
#' setClass("ModSetExample",
#'          contains = "ModifierSet",
#'          prototype = list(elementType = "ModExample"))
#' ModSetExample <- function(x, annotation, sequences, seqinfo, ...){
#'   RNAmodR:::ModifierSet("ModExample", x = x, annotation = annotation,
#'                         sequences = sequences, seqinfo = seqinfo, ...)
#' }
#' 
#' setMethod(
#'   f = "visualizeDataByCoord",
#'   signature = signature(x = "ModSetExample", coord = "GRanges"),
#'   definition = function(x, coord, type = "score", window.size = 15L, ...) {
#'   }
#' )
#' setMethod(
#'   f = "visualizeData",
#'   signature = signature(x = "ModSetExample"),
#'   definition = function(x, name, from, to, type = "score", ...) {
#'   }
#' )
NULL
