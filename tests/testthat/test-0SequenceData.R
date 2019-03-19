
context("CoverageSequenceData")
test_that("CoverageSequenceData:",{
  # SequenceData using CoverageSequenceData as test case
  annotation <- system.file("extdata","example1.gff3",package = "RNAmodR.Data")
  sequences <- system.file("extdata","example1.fasta",package = "RNAmodR.Data")
  files <- c(control = system.file("extdata","example_wt_1.bam",
                                   package = "RNAmodR.Data"),
             treated = system.file("extdata","example_wt_2.bam",
                                   package = "RNAmodR.Data"))
  csd <- CoverageSequenceData(files, annotation = annotation,
                              sequences = sequences)
  expect_false(any(lengths(rownames(csd)) == 0L))
  expect_s4_class(csd,"CoverageSequenceData")
  expect_named(csd,c("1","2"))
  expect_s4_class(colnames(csd),"CharacterList")
  expect_length(colnames(csd),2)
  expect_equal(lengths(colnames(csd)),c(2,2))
  expect_equal(colnames(csd)[[1]],colnames(csd)[[2]])
  expect_equal(colnames(csd)[[1]],c("coverage.control.1","coverage.treated.1"))
  actual <- aggregate(csd)
  expect_false(any(lengths(rownames(actual)) == 0L))
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(4,4))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control","means.treated",
                                       "sds.control","sds.treated"))
  expect_s4_class(seqinfo(csd),"Seqinfo")
  expect_equal(length(seqinfo(csd)),11)
  actual <- aggregate(csd, condition = "Control")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(2,2))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control","sds.control"))
  expect_s4_class(seqinfo(csd),"Seqinfo")
  expect_equal(length(seqinfo(csd)),11)
  actual <- aggregate(csd, condition = "Treated")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(2,2))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.treated","sds.treated"))
  expect_s4_class(seqinfo(csd),"Seqinfo")
  expect_equal(length(seqinfo(csd)),11)
})
context("PileupSequenceData")
test_that("PileupSequenceData:",{
  annotation <- system.file("extdata","example1.gff3",package = "RNAmodR.Data")
  sequences <- system.file("extdata","example1.fasta",package = "RNAmodR.Data")
  files <- c(control = system.file("extdata","example_wt_1.bam",
                                   package = "RNAmodR.Data"),
             treated = system.file("extdata","example_wt_2.bam",
                                   package = "RNAmodR.Data"))
  #
  psd <- PileupSequenceData(files, annotation = annotation,
                            sequences = sequences)
  expect_false(any(lengths(rownames(psd)) == 0L))
  expect_s4_class(psd,"PileupSequenceData")
  expect_named(psd,c("1","2"))
  expect_s4_class(colnames(psd),"CharacterList")
  expect_length(colnames(psd),2)
  expect_equal(lengths(colnames(psd)),c(10,10))
  expect_equal(colnames(psd)[[1]],colnames(psd)[[2]])
  expect_equal(colnames(psd)[[1]],c("pileup.control.1.-","pileup.control.1.G",
                                    "pileup.control.1.A","pileup.control.1.T",
                                    "pileup.control.1.C","pileup.treated.1.-",
                                    "pileup.treated.1.G","pileup.treated.1.A",
                                    "pileup.treated.1.T","pileup.treated.1.C"))
  actual <- aggregate(psd)
  expect_false(any(lengths(rownames(actual)) == 0L))
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(20,20))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control..","means.control.G",
                                       "means.control.A","means.control.T",
                                       "means.control.C","means.treated..",
                                       "means.treated.G","means.treated.A",
                                       "means.treated.T","means.treated.C",
                                       "sds.control..","sds.control.G",
                                       "sds.control.A","sds.control.T",
                                       "sds.control.C","sds.treated..",
                                       "sds.treated.G","sds.treated.A",
                                       "sds.treated.T","sds.treated.C"))
  expect_s4_class(seqinfo(psd),"Seqinfo")
  expect_equal(length(seqinfo(psd)),11)
  actual <- aggregate(psd, condition = "Control")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(10,10))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control..","means.control.G",
                                       "means.control.A","means.control.T",
                                       "means.control.C","sds.control..",
                                       "sds.control.G","sds.control.A",
                                       "sds.control.T","sds.control.C"))
  expect_s4_class(seqinfo(psd),"Seqinfo")
  expect_equal(length(seqinfo(psd)),11)
  actual <- aggregate(psd, condition = "Treated")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(10,10))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.treated..","means.treated.G",
                                       "means.treated.A","means.treated.T",
                                       "means.treated.C","sds.treated..",
                                       "sds.treated.G","sds.treated.A",
                                       "sds.treated.T","sds.treated.C"))
  expect_s4_class(seqinfo(psd),"Seqinfo")
  expect_equal(length(seqinfo(psd)),11)
})

context("ProtectedEndSequenceData")
test_that("ProtectedEndSequenceData:",{
  annotation <- system.file("extdata","example1.gff3",package = "RNAmodR.Data")
  sequences <- system.file("extdata","example1.fasta",package = "RNAmodR.Data")
  files <- c(control = system.file("extdata","example_wt_1.bam",
                                   package = "RNAmodR.Data"),
             treated = system.file("extdata","example_wt_2.bam",
                                   package = "RNAmodR.Data"))
  #
  pesd <- ProtectedEndSequenceData(files, annotation = annotation,
                                   sequences = sequences)
  expect_false(any(lengths(rownames(pesd)) == 0L))
  expect_s4_class(pesd,"ProtectedEndSequenceData")
  expect_named(pesd,c("1","2"))
  expect_s4_class(colnames(pesd),"CharacterList")
  expect_length(colnames(pesd),2)
  expect_equal(lengths(colnames(pesd)),c(2,2))
  expect_equal(colnames(pesd)[[1]],colnames(pesd)[[2]])
  expect_equal(colnames(pesd)[[1]],c("protectedend.control.1",
                                     "protectedend.treated.1"))
  actual <- aggregate(pesd)
  expect_false(any(lengths(rownames(actual)) == 0L))
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(4,4))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control","means.treated",
                                       "sds.control","sds.treated"))
  expect_s4_class(seqinfo(pesd),"Seqinfo")
  expect_equal(length(seqinfo(pesd)),11)
  actual <- aggregate(pesd, condition = "Control")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(2,2))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control","sds.control"))
  expect_s4_class(seqinfo(pesd),"Seqinfo")
  expect_equal(length(seqinfo(pesd)),11)
  actual <- aggregate(pesd, condition = "Treated")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(2,2))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.treated","sds.treated"))
  expect_s4_class(seqinfo(pesd),"Seqinfo")
  expect_equal(length(seqinfo(pesd)),11)
})
