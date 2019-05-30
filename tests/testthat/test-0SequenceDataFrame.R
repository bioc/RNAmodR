context("SequenceDataFrame")
test_that("SequenceDataFrame:",{
  data(psd,package = "RNAmodR")
  sdf <- psd[[1]]
  #Accessors
  expect_type(names(sdf),"character")
  expect_s4_class(sequences(sdf),"RNAString")
  expect_s4_class(ranges(sdf),"GRanges")
  expect_true(is.factor(conditions(sdf)))
  expect_equal(replicates(sdf),factor(c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3)))
  expect_equal(conditions(sdf),factor(rep("treated",ncol(sdf))))
  #
  sdf2 <- PileupSequenceDataFrame(as(sdf,"DataFrame"),
                                  ranges(sdf),
                                  sequences(sdf),
                                  replicates(sdf),
                                  conditions(sdf))
  expect_equal(sdf,sdf2)
  ##############################################################################
  # errors
  skip_on_bioc()
  expect_error(PileupSequenceDataFrame(as(sdf,"DataFrame"),
                                       BString(),
                                       sequences(sdf),
                                       replicates(sdf),
                                       conditions(sdf)),
               "Invalid data object: BString found, GRanges expected")
  expect_error(PileupSequenceDataFrame(as(sdf,"DataFrame"),
                                       ranges(sdf),
                                       GRanges(),
                                       replicates(sdf),
                                       conditions(sdf)),
               "Invalid data object: GRanges found, XString expected")
  expect_error(PileupSequenceDataFrame(BString(),
                                       ranges(sdf),
                                       sequences(sdf),
                                       replicates(sdf),
                                       conditions(sdf)),
               "Invalid data object: BString found, DataFrame expected")
  expect_error(PileupSequenceDataFrame(as(sdf,"DataFrame"),
                                       ranges(sdf),
                                       sequences(sdf),
                                       replicates(sdf)[1],
                                       conditions(sdf)[1]),
               "Replicate and Conditions information must match the DataFrame")
  # subsetting
  expect_type(sdf[1,],"list")
  expect_equal(length(sdf[1,]),ncol(sdf))
  expect_s4_class(sdf[1,,drop = FALSE],"DataFrame")
  expect_equal(ncol(sdf[1,,drop = FALSE]),ncol(sdf))
  expect_s4_class(sdf[,1],"PileupSequenceDataFrame")
  expect_equal(ncol(sdf[,1]),1L)
  expect_type(sdf["1",],"list")
  expect_equal(length(sdf["1",]),ncol(sdf))
  expect_s4_class(sdf["1",,drop = FALSE],"DataFrame")
  expect_equal(ncol(sdf["1",,drop = FALSE]),ncol(sdf))
  expect_s4_class(sdf[,"pileup.treated.1.G"],"PileupSequenceDataFrame")
  expect_equal(ncol(sdf[,"pileup.treated.1.G"]),1L)
})
