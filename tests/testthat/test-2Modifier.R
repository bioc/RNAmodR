
context("Modifier/ModifierSet")
test_that("Modifier/ModifierSet:",{
  data(msi,package = "RNAmodR")
  data(psd,package = "RNAmodR")
  data(e5sd,package = "RNAmodR")
  # arguments
  expect_error(RNAmodR:::.norm_args(),
               'argument "input" is missing, with no default')
  actual <- RNAmodR:::.norm_args(list())
  expect_type(actual,"list")
  expect_named(actual,c("minCoverage","minReplicate","findMod"))
  expect_error(RNAmodR:::.norm_args(list(minCoverage = 1)),
               "'minCoverage' must be a single positive integer value")
  expect_error(RNAmodR:::.norm_args(list(minReplicate = 1)),
               "'minReplicate' must be a single positive integer value")
  expect_error(RNAmodR:::.norm_args(list(minReplicate = -1L)),
               "'minReplicate' must be a single positive integer value")
  expect_error(RNAmodR:::.norm_args(list(findMod = 1)),
               "'findMod' must be a single logical value")
  # .norm_SequenceData_elements
  expect_error(RNAmodR:::.check_SequenceData_elements(),
               'argument "list" is missing, with no default')
  expect_error(RNAmodR:::.check_SequenceData_elements(msi[[1]],character()),
               "Something went wrong")
  expect_error(RNAmodR:::.check_SequenceData_elements(msi[[1]],list()),
               "Number of 'SequenceData' elements does not match")
  expect_error(RNAmodR:::.check_SequenceData_elements(msi[[1]],e5sd),
               "Type of SequenceData elements does not match")
  expect_null(RNAmodR:::.check_SequenceData_elements(msi[[1]],psd))
  # settings
  expect_error(settings(msi[[1]]) <- list(minCoverage = 1),
               "'minCoverage' must be a single positive integer value")
  expect_error(settings(msi[[1]]) <- list(minReplicate = 1),
               "'minReplicate' must be a single positive integer value")
  expect_error(settings(msi[[1]]) <- list(minReplicate = -1L),
               "'minReplicate' must be a single positive integer value")
  expect_error(settings(msi[[1]]) <- list(findMod = 1),
               "'findMod' must be a single logical value")
  actual <- settings(msi[[1]])
  expect_type(actual,"list")
  expect_named(actual,c("minCoverage","minReplicate","findMod","minScore"))
  # Modifier accessors
  expect_type(names(msi[[1]]),"character")
  expect_type(modifierType(msi[[1]]),"character")
  expect_type(modType(msi[[1]]),"character")
  expect_type(mainScore(msi[[1]]),"character")
  expect_s4_class(seqinfo(msi[[1]]),"Seqinfo")
  expect_s4_class(sequences(msi[[1]]),"RNAStringSet")
  actual <- sequences(msi[[1]], modified = TRUE)
  expect_s4_class(actual,"ModRNAStringSet")
  expect_equivalent(unique(unlist(strsplit(as.character(actual),""))),
                    c("U","A","C","G","I"))
  expect_error(sequences(msi[[1]], modified = 1),
               "'modified' has to be a single logical value")
  expect_s4_class(ranges(msi[[1]]),"GRangesList")
  expect_s4_class(modifications(msi[[1]]),"GRanges")
  expect_error(modifications(msi[[1]], perTranscript = 1),
               "'perTranscript' has to be a single logical value")
  expect_equal(unique(unlist(strand(modifications(msi[[1]],
                                                       perTranscript = TRUE)))),
               factor("*", levels = c("+","-","*")))
  expect_true(is.factor(conditions(msi[[1]])))
  expect_equal(conditions(msi[[1]]),factor(rep("treated",ncol(sequenceData(msi[[1]]))[1])))
  # Modifier creation
  library(rtracklayer)
  library(RNAmodR.Data)
  library(rtracklayer)
  sequences <- RNAmodR.Data.example.AAS.fasta()
  annotation <- GFF3File(RNAmodR.Data.example.AAS.gff3())
  files <- c(treated = RNAmodR.Data.example.wt.1(),
             treated = RNAmodR.Data.example.wt.2(),
             treated = RNAmodR.Data.example.wt.3())
  mi <- ModInosine(files, annotation = annotation, sequences = sequences)
  expect_s4_class(mi,"ModInosine")
  mi2 <- ModInosine(sequenceData(mi))
  expect_equal(mi,mi2)
  mi2 <- ModInosine(sequenceData(mi), findMod = FALSE)
  expect_false(validModification(mi2))
  expect_true(validAggregate(mi2))
  mix <- new("ModInosine")
  expect_false(hasAggregateData(mix))
  # ModifierSet accessors
  expect_type(names(msi),"character")
  expect_type(modifierType(msi),"character")
  expect_type(modType(msi),"character")
  expect_type(mainScore(msi),"character")
  expect_s4_class(seqinfo(msi),"SimpleList")
  expect_s4_class(sequences(msi),"RNAStringSet")
  actual <- sequences(msi, modified = TRUE)
  expect_s4_class(actual,"ModRNAStringSet")
  expect_equivalent(unique(unlist(strsplit(as.character(actual),""))),
                    c("U","A","C","G","I"))
  expect_error(sequences(msi[[1]], modified = 1),
               "'modified' has to be a single logical value")
  expect_s4_class(ranges(msi[[1]]),"GRangesList")
  expect_s4_class(modifications(msi[[1]]),"GRanges")
  expect_error(modifications(msi[[1]], perTranscript = 1),
               "'perTranscript' has to be a single logical value")
  expect_equal(unique(unlist(strand(modifications(msi[[1]],
                                                  perTranscript = TRUE)))),
               factor("*", levels = c("+","-","*")))
  actual <- conditions(msi)
  expect_s4_class(actual,"SimpleList")
  expect_true(is.factor(actual[[1]]))
  actual <- replicates(msi)
  expect_s4_class(actual,"SimpleList")
  expect_true(is.factor(actual[[1]]))
  # .get_class_name_for_set_from_modifier_type
  expect_error(
    RNAmodR:::.get_classname_for_ModifierSet_from_modifier_type(),
    'argument "modifiertype" is missing, with no default')
  expect_error(
    RNAmodR:::.get_classname_for_ModifierSet_from_modifier_type("abc"),
    "Class 'abc' is not implemented")
  expect_error(
    RNAmodR:::.get_classname_for_ModifierSet_from_modifier_type("DataFrame"),
    "Class 'DataFrame' does not extend the 'ModifierSet' class")
  expect_equal(
    RNAmodR:::.get_classname_for_ModifierSet_from_modifier_type("ModInosine"),
    "ModSetInosine")
  expect_equal(
    RNAmodR:::.get_classname_for_ModifierSet_from_modifier_type("ModSetInosine"),
    "ModSetInosine")
  # ModifierSet creation
  files <- list("SampleSet1" = c(treated = RNAmodR.Data.example.wt.1(),
                                 treated = RNAmodR.Data.example.wt.2(),
                                 treated = RNAmodR.Data.example.wt.3()),
                "SampleSet2" = c(treated = RNAmodR.Data.example.bud23.1(),
                                 treated = RNAmodR.Data.example.bud23.2()),
                "SampleSet3" = c(treated = RNAmodR.Data.example.trm8.1(),
                                 treated = RNAmodR.Data.example.trm8.2()))
  msi <- ModSetInosine(files, annotation = annotation, sequences = sequences)
  expect_s4_class(msi,"ModSetInosine")
  msi2 <- ModSetInosine(msi[[1]])
  expect_s4_class(msi,"ModSetInosine")
  msi2 <- ModSetInosine(as.list(msi))
  expect_s4_class(msi,"ModSetInosine")
  input <- c(files[[1]],msi[[1]])
  expect_error(ModSetInosine(input),
               "'x' must be a list containing only elements")
  expect_equal(msi,msi2)
  expect_equivalent(msi,aggregate(msi))
  expect_equivalent(msi,modify(msi))
})
