
context("Settings")
test_that("Settings:",{
  variable1 <- 1
  variable2 <- 2
  variable3 <- 3
  expect_equal(RNAmodR:::.get_name_in_parent_list(),character(0))
  expect_equal(RNAmodR:::.get_name_in_parent_list(variable1),c("variable1"))
  expect_equal(RNAmodR:::.get_name_in_parent_list(variable1,variable2),
               c("variable1","variable2"))
  .test_settings_df <- data.frame(
    variable = c("variable1","variable2"),
    stringsAsFactors = FALSE)
  expect_error(RNAmodR:::.test_setting("variable1",.test_settings_df,list(variable1 = 2),
                                        list(variable1 = 3)),
               "invalid first argument")
  expect_error(RNAmodR:::.norm_settings(list(variable2 = 4),.test_settings_df,
                                        variable1),
               "Invalid columns in settings test definition")
  .test_settings_df <- data.frame(
    variable = c("variable1","variable1"),
    testFUN = c(".test_test_TRUE",".test_test_FALSE"),
    errorValue = c(FALSE,TRUE),
    errorMessage = c("1","2"),
    stringsAsFactors = FALSE)
  expect_error(RNAmodR:::.norm_settings(list(variable2 = 4),.test_settings_df,
                                        variable1),
               "Duplicated variable names in settings test definition")
  .test_settings_df <- data.frame(
    variable = c("variable1","variable2"),
    testFUN = c(".test_test_TRUE",".test_test_FALSE"),
    errorValue = c(FALSE,TRUE),
    errorMessage = c("1","2"),
    stringsAsFactors = FALSE)
  expect_equal(RNAmodR:::.test_setting("variable1",.test_settings_df,list(variable1 = 2),
                                       list(variable1 = 3)),
               3)
  expect_equal(RNAmodR:::.test_setting("variable1",.test_settings_df,list(variable1 = 2),
                                       list()),
               2)
  .test_settings_df <- data.frame(
    variable = c("variable1","variable2"),
    testFUN = c(".test_test_TRUE",".test_test_FALSE"),
    errorValue = c(TRUE,TRUE),
    errorMessage = c("1","2"),
    stringsAsFactors = FALSE)
  expect_error(RNAmodR:::.test_setting("variable1",.test_settings_df,list(variable1 = 2),
                                       list(variable1 = 3)),
               "1")
  expect_equal(RNAmodR:::.test_setting("variable2",.test_settings_df,list(variable2 = 2),
                                       list(variable2 = 3)),
               3)
  expect_equal(RNAmodR:::.norm_settings(list(variable2 = 4),.test_settings_df,
                                        variable1),
               list(variable1 = 1))
  expect_equal(RNAmodR:::.norm_settings(list(variable2 = 4),.test_settings_df,
                                        variable2),
               list(variable2 = 4))
  expect_error(RNAmodR:::.norm_settings(list(variable1 = 4),.test_settings_df,
                                        variable1),
               "1")
  expect_error(RNAmodR:::.norm_settings(list(variable3 = 4),.test_settings_df,
                                        variable3),
               "Test for variables 'variable3' not found.")
})
