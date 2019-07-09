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

test_that("remove_dataset removes datasets", {
  request_obj <- new("ReactomeAnalysisRequest", method = "Camera")

  # create test data
  test_expr <- data.frame(
    sample.1 = c(1,2,3,4),
    sample.2 = c(1,2,3,4),
    sample.3 = c(1,2,3,4),
    row.names = c("CD19", "MS4A1", "MITF", "SDC1")
  )

  sample_data <- data.frame(
    treatment = c("Control", "Treatment", "Control"),
    lab = c("L1", "L1", "L2"),
    row.names = paste0("sample.", 1:3)
  )

  request_obj <- add_dataset(request_obj, expression_values = test_expr, name = "Test 1", type = "rnaseq_counts",
                             comparison_factor = "treatment", comparison_group_1 = "Control", comparison_group_2 = "Treatment",
                             sample_data = sample_data)

  expect_equal("datasets" %in% names(request_obj@request_object), TRUE)
  expect_equal(nrow(request_obj@request_object[["datasets"]]), 1)

  # remove the dataset again
  request_obj <- remove_dataset(request_obj, "Test 1")

  expect_equal(nrow(request_obj@request_object[["datasets"]]), 0)
})

test_that("add_dataset adds different datasets", {
  request_obj <- new("ReactomeAnalysisRequest", method = "Camera")
  expect_equal("datasets" %in% names(request_obj@request_object), FALSE)

  # create test data
  test_expr <- data.frame(
    sample.1 = c(1,2,3,4),
    sample.2 = c(1,2,3,4),
    sample.3 = c(1,2,3,4),
    row.names = c("CD19", "MS4A1", "MITF", "SDC1")
  )

  sample_data <- data.frame(
    treatment = c("Control", "Treatment", "Control"),
    lab = c("L1", "L1", "L2"),
    row.names = paste0("sample.", 1:3)
  )

  request_obj <- add_dataset(request_obj, expression_values = test_expr, name = "Test 1", type = "rnaseq_counts",
                             comparison_factor = "treatment", comparison_group_1 = "Control", comparison_group_2 = "Treatment",
                             sample_data = sample_data)

  expect_equal("datasets" %in% names(request_obj@request_object), TRUE)
  expect_equal(nrow(request_obj@request_object[["datasets"]]), 1)

  # add a second dataset
  request_obj <- add_dataset(request_obj, expression_values = test_expr, name = "Test 2", type = "rnaseq_counts",
                             comparison_factor = "treatment", comparison_group_1 = "Control", comparison_group_2 = "Treatment",
                             sample_data = sample_data)

  expect_equal(nrow(request_obj@request_object[["datasets"]]), 2)

  # overwrite the second dataset
  request_obj <- add_dataset(request_obj, expression_values = test_expr, name = "Test 2", type = "rnaseq_counts",
                             comparison_factor = "treatment", comparison_group_1 = "Control", comparison_group_2 = "Treatment",
                             sample_data = sample_data, overwrite = TRUE)

  expect_equal(nrow(request_obj@request_object[["datasets"]]), 2)
})
