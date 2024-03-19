#Loading of packages
library(testthat)

#Source code of the function to test
source("R/gprofiler2_analysis.R")

de_res <- get_demo_de_res()
gp <- gprofiler2_analysis(de_res)

is_valid_gp <- function(gp){

  #Test expected objects in gp to be lists
  expect_true(is.list(gp))
  expect_true(is.list(gp$result$parents))
  expect_true(is.list(gp$meta))
  expect_true(!any(is.na(gp)))

  # Test expected columns in the list result
  expect_true(("query" %in% colnames(gp$result)))
  expect_true(("significant" %in% colnames(gp$result)))
  expect_true(("p_value" %in% colnames(gp$result)))
  expect_true(("term_size" %in% colnames(gp$result)))
  expect_true(("query_size" %in% colnames(gp$result)))
  expect_true(("intersection_size" %in% colnames(gp$result)))
  expect_true(("precision" %in% colnames(gp$result)))
  expect_true(("recall" %in% colnames(gp$result)))
  expect_true(("term_id" %in% colnames(gp$result)))
  expect_true(("source" %in% colnames(gp$result)))
  expect_true(("term_name" %in% colnames(gp$result)))
  expect_true(("effective_domain_size" %in% colnames(gp$result)))
  expect_true(("source_order" %in% colnames(gp$result)))
  expect_true(("parents" %in% colnames(gp$result)))
  expect_true(("evcodes" %in% colnames(gp$result)))
  expect_true(("intersection" %in% colnames(gp$result)))

  #test expected values to be characters in the list result
  expect_true(is.character(gp$result$query))
  expect_true(is.character(gp$result$term_id))
  expect_true(is.character(gp$result$source))
  expect_true(is.character(gp$result$term_name))
  expect_true(is.logical(gp$result$significant))

  #test expected values to be double in the list result
  expect_true(is.double(gp$result$p_value))
  expect_true(is.double(gp$result$precision))
  expect_true(is.double(gp$result$recall))
  #test expected values to be integer in the list result
  int_vector <- c(gp$result$term_size, gp$result$query_size, gp$result$intersection_size, gp$result$effective_domain_size, gp$result$source_order)
  expect_true(is.integer(int_vector))

}

test_that("gprofiler2_analysis works with valid data", {
  is_valid_gp(gp)
})
