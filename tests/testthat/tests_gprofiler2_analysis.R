#Loading of packages
library(rnaseq)

de_res <- get_demo_de_res()
enr <- gprofiler2_analysis(de_res)

is_valid_gp <- function(enr){

  #Test expected objects in gp to be lists
  expect_true(is.list(enr))
  expect_true(is.list(enr$result))
  expect_true(is.list(enr$result$parents))
  expect_true(is.list(enr$meta))
  expect_true(!any(is.na(enr)))

  # Test expected columns in the list result
  expect_true(("query" %in% colnames(enr$result)))
  expect_true(("significant" %in% colnames(enr$result)))
  expect_true(("p_value" %in% colnames(enr$result)))
  expect_true(("term_size" %in% colnames(enr$result)))
  expect_true(("query_size" %in% colnames(enr$result)))
  expect_true(("intersection_size" %in% colnames(enr$result)))
  expect_true(("precision" %in% colnames(enr$result)))
  expect_true(("recall" %in% colnames(enr$result)))
  expect_true(("term_id" %in% colnames(enr$result)))
  expect_true(("source" %in% colnames(enr$result)))
  expect_true(("term_name" %in% colnames(enr$result)))
  expect_true(("effective_domain_size" %in% colnames(enr$result)))
  expect_true(("source_order" %in% colnames(enr$result)))
  expect_true(("parents" %in% colnames(enr$result)))
  expect_true(("evidence_codes" %in% colnames(enr$result)))
  expect_true(("intersection" %in% colnames(enr$result)))
  expect_true(("GeneRatio" %in% colnames(enr$result)))

  #test expected values to be characters in the list result
  expect_true(is.character(enr$result$query))
  expect_true(is.character(enr$result$term_id))
  expect_true(is.character(enr$result$source))
  expect_true(is.character(enr$result$term_name))
  expect_true(is.character(enr$result$evidence_codes))
  expect_true(is.character(enr$result$intersection))
  expect_true(is.character(enr$result$GeneRatio))
  expect_true(is.logical(enr$result$significant))

  #test expected values to be double in the list result
  expect_true(is.double(enr$result$p_value))
  expect_true(is.double(enr$result$precision))
  expect_true(is.double(enr$result$recall))
  #test expected values to be integer in the list result
  int_vector <- c(enr$result$term_size, enr$result$query_size, enr$result$intersection_size, enr$result$effective_domain_size, enr$result$source_order)
  expect_true(is.integer(int_vector))

}

test_that("gprofiler2_analysis works with valid data", {
  is_valid_gp(enr$gp$de_res)
})

is_valid_cp <- function(enr){

  #Test expected type of objects in cp
  expect_true(class(enr$enrichResult) == "enrichResult")
  expect_true(is.list(enr$df))
  expect_true(is.list(enr$df$parents))

  #Test columns in cp$x$result
  expect_true(("ID" %in% colnames(enr$df)))
  expect_true(("query" %in% colnames(enr$df)))
  expect_true(("significant" %in% colnames(enr$df)))
  expect_true(("pvalue" %in% colnames(enr$df)))
  expect_true(("term_size" %in% colnames(enr$df)))
  expect_true(("query_size" %in% colnames(enr$df)))
  expect_true(("intersection_size" %in% colnames(enr$df)))
  expect_true(("precision" %in% colnames(enr$df)))
  expect_true(("recall" %in% colnames(enr$df)))
  expect_true(("source" %in% colnames(enr$df)))
  expect_true(("Description" %in% colnames(enr$df)))
  expect_true(("effective_domain_size" %in% colnames(enr$df)))
  expect_true(("source_order" %in% colnames(enr$df)))
  expect_true(("parents" %in% colnames(enr$df)))
  expect_true(("evidence_codes" %in% colnames(enr$df)))
  expect_true(("geneID" %in% colnames(enr$df)))
  expect_true(("GeneRatio" %in% colnames(enr$df)))
  expect_true(("p.adjust" %in% colnames(enr$df)))
  expect_true(("Count" %in% colnames(enr$df)))
  expect_true(("setSize" %in% colnames(enr$df)))

  #test expected values to be characters in cp$x$df
  expect_true(is.character(enr$df$query))
  expect_true(is.character(enr$df$ID))
  expect_true(is.character(enr$df$source))
  expect_true(is.character(enr$df$Description))
  expect_true(is.character(enr$df$evidence_codes))
  expect_true(is.character(enr$df$geneID))
  expect_true(is.character(enr$df$GeneRatio))

  #test expected values to be double in cp$x$df
  expect_true(is.double(enr$df$pvalue))
  expect_true(is.double(enr$df$precision))
  expect_true(is.double(enr$df$recall))
  expect_true(is.double(enr$df$p.adjust))
  expect_true(is.double(enr$df$setSize))

  #test expected values to be logical in cp$x$df
  expect_true(is.logical(enr$df$significant))

  #test expected values to be integer in cp$x$df
  expect_true(is.integer(enr$df$term_size))
  expect_true(is.integer(enr$df$query_size))
  expect_true(is.integer(enr$df$intersection_size))
  expect_true(is.integer(enr$df$effective_domain_size))
  expect_true(is.integer(enr$df$source_order))
  expect_true(is.integer(enr$df$Count))
}

test_that("gprofiler2_analysis works with valid data", {
  is_valid_cp(enr$cp$de_res)
})
