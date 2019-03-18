
context("Export")
test_that("Export:",{
  data(msi, package = "RNAmodR")
  data(sds, package = "RNAmodR")
  data(psd, package = "RNAmodR")
  file <- tempfile()
  # SequenceData
  actual <- export.wig(psd, file)
  expect_null(actual)
  actual <- export.bw(psd, file)
  expect_s4_class(actual,"BigWigFile")
  # SequenceDataSet
  actual <- export.wig(sds, file)
  expect_null(actual)
  actual <- export.bw(sds, file)
  expect_s4_class(actual,"BigWigFile")
  # Modifier
  actual <- export.wig(msi[[1]], file)
  expect_null(actual)
  actual <- export.bw(msi[[1]], file)
  expect_s4_class(actual,"BigWigFile")
})
