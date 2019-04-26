
context("NormEndSequenceData")
test_that("NormEndSequenceData:",{
  # SequenceData using CoverageSequenceData as test case
  library(Rsamtools)
  library(RNAmodR.Data)
  annotation <- GFF3File(RNAmodR.Data.example.man.gff3())
  sequences <- RNAmodR.Data.example.man.fasta()
  files <- c(control = RNAmodR.Data.example.wt.1(),
             treated = RNAmodR.Data.example.wt.2())
  #
  ne5sd <- NormEnd5SequenceData(files, annotation = annotation,
                                sequences = sequences)
  expect_false(any(lengths(rownames(ne5sd)) == 0L))
  expect_s4_class(ne5sd,"NormEnd5SequenceData")
  expect_named(ne5sd,c("1","2"))
  expect_s4_class(colnames(ne5sd),"CharacterList")
  expect_length(colnames(ne5sd),2)
  expect_equal(lengths(colnames(ne5sd)),c(6,6))
  expect_equal(colnames(ne5sd)[[1]],colnames(ne5sd)[[2]])
  expect_equal(colnames(ne5sd)[[1]],c("normend5.control.1.ends",
                                      "normend5.control.1.norm.tx",
                                      "normend5.control.1.norm.ol",
                                      "normend5.treated.1.ends",
                                      "normend5.treated.1.norm.tx",
                                      "normend5.treated.1.norm.ol"))
  actual <- aggregate(ne5sd)
  expect_false(any(lengths(rownames(actual)) == 0L))
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(12,12))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control.ends","means.control.tx",
                                       "means.control.ol","means.treated.ends",
                                       "means.treated.tx","means.treated.ol",
                                       "sds.control.ends","sds.control.tx",
                                       "sds.control.ol","sds.treated.ends",
                                       "sds.treated.tx","sds.treated.ol"))
  expect_s4_class(seqinfo(ne5sd),"Seqinfo")
  expect_equal(length(seqinfo(ne5sd)),11)
  actual <- aggregate(ne5sd, condition = "Control")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(6,6))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control.ends","means.control.tx",
                                       "means.control.ol","sds.control.ends",
                                       "sds.control.tx","sds.control.ol"))
  expect_s4_class(seqinfo(ne5sd),"Seqinfo")
  expect_equal(length(seqinfo(ne5sd)),11)
  actual <- aggregate(ne5sd, condition = "Treated")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(6,6))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.treated.ends","means.treated.tx",
                                       "means.treated.ol","sds.treated.ends",
                                       "sds.treated.tx","sds.treated.ol"))
  expect_s4_class(seqinfo(ne5sd),"Seqinfo")
  expect_equal(length(seqinfo(ne5sd)),11)
  #
  ne3sd <- NormEnd3SequenceData(files, annotation = annotation,
                                sequences = sequences)
  expect_false(any(lengths(rownames(ne3sd)) == 0L))
  expect_s4_class(ne3sd,"NormEnd3SequenceData")
  expect_named(ne3sd,c("1","2"))
  expect_s4_class(colnames(ne3sd),"CharacterList")
  expect_length(colnames(ne3sd),2)
  expect_equal(lengths(colnames(ne3sd)),c(6,6))
  expect_equal(colnames(ne3sd)[[1]],colnames(ne3sd)[[2]])
  expect_equal(colnames(ne3sd)[[1]],c("normend3.control.1.ends",
                                      "normend3.control.1.norm.tx",
                                      "normend3.control.1.norm.ol",
                                      "normend3.treated.1.ends",
                                      "normend3.treated.1.norm.tx",
                                      "normend3.treated.1.norm.ol"))
  actual <- aggregate(ne3sd)
  expect_false(any(lengths(rownames(actual)) == 0L))
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(12,12))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control.ends","means.control.tx",
                                       "means.control.ol","means.treated.ends",
                                       "means.treated.tx","means.treated.ol",
                                       "sds.control.ends","sds.control.tx",
                                       "sds.control.ol","sds.treated.ends",
                                       "sds.treated.tx","sds.treated.ol"))
  expect_s4_class(seqinfo(ne3sd),"Seqinfo")
  expect_equal(length(seqinfo(ne3sd)),11)
  actual <- aggregate(ne3sd, condition = "Control")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(6,6))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control.ends","means.control.tx",
                                       "means.control.ol","sds.control.ends",
                                       "sds.control.tx","sds.control.ol"))
  expect_s4_class(seqinfo(ne3sd),"Seqinfo")
  expect_equal(length(seqinfo(ne3sd)),11)
  actual <- aggregate(ne3sd, condition = "Treated")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(6,6))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.treated.ends","means.treated.tx",
                                       "means.treated.ol","sds.treated.ends",
                                       "sds.treated.tx","sds.treated.ol"))
  expect_s4_class(seqinfo(ne3sd),"Seqinfo")
  expect_equal(length(seqinfo(ne3sd)),11)
})
