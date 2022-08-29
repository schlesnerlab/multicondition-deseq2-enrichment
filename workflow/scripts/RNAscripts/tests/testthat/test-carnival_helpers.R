test_that("ConvertHumanGEneHomologs works", {
  testthat::expect_equal(
    RNAscripts::convertHumanGeneHomologs("APLN"),
    "Apln"
  )
})
