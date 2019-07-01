context("ReactomeAnalysisRequest functions")
library(ReactomeGSA)

test_that("Request method is set by the constructor", {
  request_obj <- new("ReactomeAnalysisRequest", method = "Camera")
  expect_equal(request_obj@method, "Camera")

  request_obj <- new("ReactomeAnalysisRequest", method = "ssGSEA")
  expect_equal(request_obj@method, "ssGSEA")

  # invalid methods possible as well
  request_obj <- new("ReactomeAnalysisRequest", method = "asda")
  expect_equal(request_obj@method, "asda")

  # make sure the method is changed
  request_obj <- set_method(request_obj, method = "osda")
  expect_equal(request_obj@method, "osda")
})

test_that("Set parameters changes parameters", {
  request_obj <- new("ReactomeAnalysisRequest", method = "Camera")

  # parameters are missing
  expect_equal("parameters" %in% names(request_obj@request_object), FALSE)

  # setting a parameter creates the parameters object
  request_obj <- set_parameters(request_obj, min_missing_values = 0.01)

  expect_equal("parameters" %in% names(request_obj@request_object), TRUE)
  expect_equal(nrow(request_obj@request_object[["parameters"]]), 1)

  request_obj <- set_parameters(request_obj, create_reactome_visualization = FALSE)
  expect_equal(nrow(request_obj@request_object[["parameters"]]), 2)
  # bool converted to string
  expect_equal(request_obj@request_object[["parameters"]][2, "value"], "FALSE")
})

test_that("add_dataset adds different datasets", {
  request_obj <- new("ReactomeAnalysisRequest", method = "Camera")
  expect_equal("data" %in% names(request_obj@request_object), FALSE)
})
