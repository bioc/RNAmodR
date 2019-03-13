
context("SequenceDataSet")
test_that("SequenceDataSet:",{
  data(psd,package = "RNAmodR")
  data(e5sd,package = "RNAmodR")
  sds <- SequenceDataSet(e5sd,psd)
  sds2 <- as(list(e5sd,psd),"SequenceDataSet")
  expect_equal(sds,sds2)
  sds2 <- as(SimpleList(e5sd,psd),"SequenceDataSet")
  expect_equal(sds,sds2)
  expect_equal(seqinfo(sds),seqinfo(sds[[1]]))
  expect_equal(names(sds),names(sds[[1]]))
  expect_equal(ranges(sds),ranges(sds[[1]]))
  expect_equal(sequences(sds),sequences(sds[[1]]))
  expect_equal(bamfiles(sds),bamfiles(sds[[1]]))
  expect_equal(aggregate(sds),
               SimpleList(End5SequenceData = aggregate(sds[[1]]),
                          PileupSequenceData = aggregate(sds[[2]])))
  actual <- as.list(sds)
  expect_type(actual,"list")
  expect_named(actual,c(class(e5sd),class(psd)))
  expect_equal(actual[[1]],e5sd)
  expect_equal(actual[[2]],psd)
  actual <- as.list(sds, use.names = FALSE)
  expect_null(names(actual))
  expect_true(validObject(actual))
  # error
  expect_error(SequenceDataSet(e5sd,c(1,2,3)),
               "All elements in 'x' must be SequenceData objects")
  expect_error(as.list(sds, use.names = 1))
})

context("SequenceDataList")
test_that("SequenceDataList:",{
  data(psd,package = "RNAmodR")
  data(e5sd,package = "RNAmodR")
  sds <- SequenceDataSet(e5sd,psd)
  sdl <- SequenceDataList(sds,e5sd,psd)
  sdl2 <- as(list(sds,e5sd,psd),"SequenceDataList")
  expect_equal(sdl,sdl2)
  sdl2 <- as(SimpleList(sds,e5sd,psd),"SequenceDataList")
  expect_equal(sdl,sdl2)
  # expect_equal(seqinfo(sdl),seqinfo(sdl[[1]]))
  # expect_equal(names(sdl),names(sdl[[1]]))
  # expect_equal(ranges(sdl),ranges(sdl[[1]]))
  # expect_equal(sequences(sdl),sequences(sdl[[1]]))
  # expect_equal(bamfiles(sdl),bamfiles(sdl[[1]]))
  # expect_equal(aggregate(sdl),
  #              SimpleList(End5SequenceData = aggregate(sdl[[1]]),
  #                         PileupSequenceData = aggregate(sdl[[2]])))
  # actual <- as.list(sdl)
  # expect_type(actual,"list")
  # expect_named(actual,c(class(e5sd),class(psd)))
  # expect_equal(actual[[1]],e5sd)
  # expect_equal(actual[[2]],psd)
  # actual <- as.list(sdl, use.names = FALSE)
  # expect_null(names(actual))
  # expect_true(validObject(actual))
  # # error
  # expect_error(SequenceDataSet(e5sd,c(1,2,3)),
  #              "All elements in 'x' must be SequenceData objects")
  # expect_error(as.list(sdl, use.names = 1))
})
