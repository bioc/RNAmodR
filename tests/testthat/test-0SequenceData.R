
context("SequenceData")
test_that("SequenceData:",{
  # SequenceData using CoverageSequenceData as test case
  library(RNAmodR.Data)
  library(rtracklayer)
  annotation <- GFF3File(RNAmodR.Data.example.man.gff3())
  sequences <- RNAmodR.Data.example.man.fasta()
  files <- c(treated = RNAmodR.Data.example.wt.2())
  #
  e5sd <- End5SequenceData(files, annotation = annotation,
                              sequences = sequences)
  expect_false(any(lengths(rownames(e5sd)) == 0L))
  expect_s4_class(e5sd,"End5SequenceData")
  expect_named(e5sd,c("1","2"))
  expect_s4_class(colnames(e5sd),"CharacterList")
  expect_length(colnames(e5sd),2)
  expect_equal(lengths(colnames(e5sd)),c(1,1))
  expect_equal(colnames(e5sd)[[1]],colnames(e5sd)[[2]])
  expect_equal(colnames(e5sd)[[1]],c("end5.treated.1"))
  ##############################################################################
  skip_on_bioc()
  annotation <- GFF3File(RNAmodR.Data.example.man.gff3())
  sequences <- RNAmodR.Data.example.man.fasta()
  files <- c(control = RNAmodR.Data.example.wt.1(),
             treated = RNAmodR.Data.example.wt.2())
  #
  e5sd <- End5SequenceData(files, annotation = annotation,
                           sequences = sequences)
  expect_false(any(lengths(rownames(e5sd)) == 0L))
  expect_s4_class(e5sd,"End5SequenceData")
  expect_named(e5sd,c("1","2"))
  expect_s4_class(colnames(e5sd),"CharacterList")
  expect_length(colnames(e5sd),2)
  expect_equal(lengths(colnames(e5sd)),c(2,2))
  expect_equal(colnames(e5sd)[[1]],colnames(e5sd)[[2]])
  expect_equal(colnames(e5sd)[[1]],c("end5.control.1","end5.treated.1"))
  actual <- aggregate(e5sd)
  expect_false(any(lengths(rownames(actual)) == 0L))
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(4,4))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control","means.treated",
                                       "sds.control","sds.treated"))
  expect_s4_class(seqinfo(e5sd),"Seqinfo")
  expect_equal(length(seqinfo(e5sd)),11)
  actual <- aggregate(e5sd, condition = "Control")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(2,2))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control","sds.control"))
  expect_s4_class(seqinfo(e5sd),"Seqinfo")
  expect_equal(length(seqinfo(e5sd)),11)
  actual <- aggregate(e5sd, condition = "Treated")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(2,2))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.treated","sds.treated"))
  expect_s4_class(seqinfo(e5sd),"Seqinfo")
  expect_equal(length(seqinfo(e5sd)),11)
  # General accessors
  expect_type(names(e5sd),"character")
  expect_s4_class(seqinfo(e5sd),"Seqinfo")
  expect_s4_class(sequences(e5sd),"RNAStringSet")
  expect_s4_class(ranges(e5sd),"GRangesList")
  expect_true(is.factor(conditions(e5sd)))
  expect_equal(conditions(e5sd),factor(c("control","treated")))
  expect_true(is.factor(replicates(e5sd)))
  expect_equal(replicates(e5sd),factor(c(1,1)))
  #
  e3sd <- End3SequenceData(files, annotation = annotation,
                          sequences = sequences)
  expect_false(any(lengths(rownames(e3sd)) == 0L))
  expect_s4_class(e3sd,"End3SequenceData")
  expect_named(e3sd,c("1","2"))
  expect_s4_class(colnames(e3sd),"CharacterList")
  expect_length(colnames(e3sd),2)
  expect_equal(lengths(colnames(e3sd)),c(2,2))
  expect_equal(colnames(e3sd)[[1]],colnames(e3sd)[[2]])
  expect_equal(colnames(e3sd)[[1]],c("end3.control.1","end3.treated.1"))
  #############################################################################
  actual <- aggregate(e3sd)
  expect_false(any(lengths(rownames(actual)) == 0L))
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(4,4))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control","means.treated",
                                       "sds.control","sds.treated"))
  expect_s4_class(seqinfo(e3sd),"Seqinfo")
  expect_equal(length(seqinfo(e3sd)),11)
  actual <- aggregate(e3sd, condition = "Control")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(2,2))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control","sds.control"))
  expect_s4_class(seqinfo(e3sd),"Seqinfo")
  expect_equal(length(seqinfo(e3sd)),11)
  actual <- aggregate(e3sd, condition = "Treated")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(2,2))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.treated","sds.treated"))
  expect_s4_class(seqinfo(e3sd),"Seqinfo")
  expect_equal(length(seqinfo(e3sd)),11)
  #
  esd <- EndSequenceData(files, annotation = annotation,
                          sequences = sequences)
  expect_false(any(lengths(rownames(esd)) == 0L))
  expect_s4_class(esd,"EndSequenceData")
  expect_named(esd,c("1","2"))
  expect_s4_class(colnames(esd),"CharacterList")
  expect_length(colnames(esd),2)
  expect_equal(lengths(colnames(esd)),c(2,2))
  expect_equal(colnames(esd)[[1]],colnames(esd)[[2]])
  expect_equal(colnames(esd)[[1]],c("end.control.1","end.treated.1"))
  ##############################################################################
  actual <- aggregate(esd)
  expect_false(any(lengths(rownames(actual)) == 0L))
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(4,4))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control","means.treated",
                                       "sds.control","sds.treated"))
  expect_s4_class(seqinfo(esd),"Seqinfo")
  expect_equal(length(seqinfo(esd)),11)
  actual <- aggregate(esd, condition = "Control")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(2,2))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.control","sds.control"))
  expect_s4_class(seqinfo(esd),"Seqinfo")
  expect_equal(length(seqinfo(esd)),11)
  actual <- aggregate(esd, condition = "Treated")
  expect_s4_class(actual,"CompressedSplitDataFrameList")
  expect_s4_class(actual,"SplitDataFrameList")
  expect_equal(length(actual),2)
  expect_length(colnames(actual),2)
  expect_equal(lengths(colnames(actual)),c(2,2))
  expect_equal(colnames(actual)[[1]],colnames(actual)[[2]])
  expect_equal(colnames(actual)[[1]],c("means.treated","sds.treated"))
  expect_s4_class(seqinfo(esd),"Seqinfo")
  expect_equal(length(seqinfo(esd)),11)
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
  ##############################################################################
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
  ##############################################################################
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
  # SequenceData using CoverageSequenceData as test case
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
  ##############################################################################
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
  ##############################################################################
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
  ##############################################################################
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
