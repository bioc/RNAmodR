
context("ROC")
test_that("ROC:",{
  data(msi, package = "RNAmodR")
  # arguments
  actual <- RNAmodR:::.norm_prediction_args()
  pred <- actual
  expect_type(actual,"list")
  expect_equal(actual,list())
  expect_warning(RNAmodR:::.norm_prediction_args(list("abc")),
                 "Unnamed list for 'prediction.args'")
  #
  expect_error(RNAmodR:::.norm_performance_args(),
               'argument "x" is missing, with no default')
  expect_error(RNAmodR:::.norm_performance_args(list()),
               'argument "x" is missing, with no default')
  actual <- RNAmodR:::.norm_performance_args(list(), msi)
  perf <- actual
  expect_type(actual,"list")
  expect_named(actual,c("measure","x.measure","score"))
  expect_warning(RNAmodR:::.norm_performance_args(list("abc"), msi),
                 "Unnamed list for 'performance.args'")
  expect_error(RNAmodR:::.norm_performance_args(list(measure = 1), msi),
               "'measure' must a single character compatible with")
  expect_error(RNAmodR:::.norm_performance_args(list(x.measure = 1), msi),
               "'x.measure' must a single character compatible with")
  expect_error(RNAmodR:::.norm_performance_args(list(score = 1), msi),
               "'score' must a single character and a valid column")
  #
  actual <- RNAmodR:::.norm_plot_args()
  plot <- actual
  expect_type(actual,"list")
  expect_named(actual,c("colorize","lwd","colorize.palette","abline","AUC"))
  expect_error(RNAmodR:::.norm_plot_args(list(colorize.palette = 1)),
               "'colorize.palette' must a single character compatible")
  expect_error(RNAmodR:::.norm_plot_args(list(lwd = "abc")),
               "'lwd' must be a single numeric value")
  expect_error(RNAmodR:::.norm_plot_args(list(colorize = 1)),
               "'colorize' must a single logical value")
  expect_error(RNAmodR:::.norm_plot_args(list(AUC = 1)),
               "'AUC' must a single logical value")
  expect_error(RNAmodR:::.norm_plot_args(list(abline = 1)),
               "'abline' must a single logical value")
  #
  expect_error(RNAmodR:::.readjust_plot_args(),
               'argument "performance.args" is missing')
  actual <- RNAmodR:::.readjust_plot_args(plot, perf)
  expect_type(actual,"list")
  #
  coord <- modifications(msi[[1]])
  expect_error(RNAmodR:::.get_prediction_data_Modifier(),
               'argument "x" is missing, with no default')
  actual <- RNAmodR:::.get_prediction_data_Modifier(msi[[1]],coord,"score")
  expect_type(actual,"list")
  expect_named(actual,c("score"))
  expect_s4_class(actual[[1]],"DataFrame")
  expect_equal(colnames(actual[[1]]),c("labels","predictions"))
  expect_error(RNAmodR:::.get_prediction_data_Modifier(msi[[1]],coord,"score1"),
               "Score identifier 'score1' not found in the data")
  #
  expect_error(RNAmodR:::.get_prediction_data_Modifier(),
               'argument "x" is missing, with no default')
  expect_error(RNAmodR:::.get_prediction_data_ModifierSet(msi,coord,"score1"),
               "Score identifier 'score1' not found in the data")
  actual <- RNAmodR:::.get_prediction_data_ModifierSet(msi,coord,"score")
  expect_type(actual,"list")
  expect_named(actual,c("score"))
  expect_type(actual[[1]],"list")
  expect_named(actual[[1]],c("predictions","labels"))
  expect_s3_class(actual[[1]][[1]],"data.frame")
  expect_s3_class(actual[[1]][[2]],"data.frame")
  #
  expect_null(plotROC(msi[[1]], coord))
  expect_null(plotROC(msi, coord))
})
