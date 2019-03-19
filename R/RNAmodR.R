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
#' These functions are not intended for general use, but are used for 
#' additional package development.
#' 
#' @param object,range,data,modType,scoreFun,source,type,x,i,j,...,exact,value
#' Internal arguments
#' 
#' @return internally used values
NULL

#' @name RNAmodR-datasets
#' @title Example data in the RNAmodR package
#' @description 
#' This contains an example ModifierSet object
#' @docType data
#' @usage msi
#' @usage psd
#' @usage e5sd
#' @format a \code{ModSetInosine} instance
#' @keywords datasets
"msi"
#' @rdname RNAmodR-datasets
#' @format a \code{SequenceDataSet} instance
"sds"
#' @rdname RNAmodR-datasets
#' @format a \code{SequenceDataList} instance
"sdl"
#' @rdname RNAmodR-datasets
#' @format a \code{PileupSequenceData} instance
"psd"
#' @rdname RNAmodR-datasets
#' @format a \code{End5SequenceData} instance
"e5sd"
#' @rdname RNAmodR-datasets
#' @format a \code{End3SequenceData} instance
"e3sd"
#' @rdname RNAmodR-datasets
#' @format a \code{EndSequenceData} instance
"esd"
#' @rdname RNAmodR-datasets
#' @format a \code{CoverageSequenceData} instance
"csd"
#' @rdname RNAmodR-datasets
#' @format a \code{NormEnd3SequenceData} instance
"ne3sd"
#' @rdname RNAmodR-datasets
#' @format a \code{NormEnd5SequenceData} instance
"ne5sd"
#' @rdname RNAmodR-datasets
#' @format a \code{ProtectedEndSequenceData} instance
"pesd"
