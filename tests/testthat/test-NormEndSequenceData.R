
context("NormEndSequenceData")
test_that("NormEndSequenceData:",{
  # SequenceData using CoverageSequenceData as test case
  annotation <- system.file("extdata","example1.gff3",package = "RNAmodR.Data")
  sequences <- system.file("extdata","example1.fasta",package = "RNAmodR.Data")
  files <- c(control = system.file("extdata","example_wt_1.bam",
                                   package = "RNAmodR.Data"),
             treated = system.file("extdata","example_wt_2.bam",
                                   package = "RNAmodR.Data"))
  #
  ne5sd <- NormEnd5SequenceData(files, annotation = annotation,
                                sequences = sequences)
  expect_s4_class(ne5sd,"NormEnd5SequenceData")
  expect_named(ne5sd,c("1","2"))
  expect_s4_class(colnames(ne5sd),"CharacterList")
  expect_length(colnames(ne5sd),2)
  expect_equal(lengths(colnames(ne5sd)),c(6,6))
  expect_equal(colnames(ne5sd)[[1]],colnames(ne5sd)[[2]])
  expect_equal(colnames(ne5sd)[[1]],c("norm.end5.control.1.ends",
                                      "norm.end5.control.1.norm.tx",
                                      "norm.end5.control.1.norm.ol",
                                      "norm.end5.treated.1.ends",
                                      "norm.end5.treated.1.norm.tx",
                                      "norm.end5.treated.1.norm.ol"))
  actual <- aggregate(ne5sd)
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
  expect_s4_class(ne3sd,"NormEnd3SequenceData")
  expect_named(ne3sd,c("1","2"))
  expect_s4_class(colnames(ne3sd),"CharacterList")
  expect_length(colnames(ne3sd),2)
  expect_equal(lengths(colnames(ne3sd)),c(6,6))
  expect_equal(colnames(ne3sd)[[1]],colnames(ne3sd)[[2]])
  expect_equal(colnames(ne3sd)[[1]],c("norm.end3.control.1.ends",
                                      "norm.end3.control.1.norm.tx",
                                      "norm.end3.control.1.norm.ol",
                                      "norm.end3.treated.1.ends",
                                      "norm.end3.treated.1.norm.tx",
                                      "norm.end3.treated.1.norm.ol"))
  actual <- aggregate(ne3sd)
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
