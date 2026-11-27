fixture_cache <- new.env(parent = emptyenv())

read_input_fixture <- function(input_name) {
  if (exists(input_name, fixture_cache, inherits = FALSE)) {
    return(get(input_name, fixture_cache, inherits = FALSE))
  }

  input <- utils::read.csv(
    system.file("extdata", input_name, package = "cancerRC"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  assign(input_name, input, fixture_cache)
  input
}

read_expected_fixture <- function(expected_name) {
  utils::read.csv(
    system.file("extdata", expected_name, package = "cancerRC"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

test_that("internal bundle loaders validate the frozen contracts", {
  main_loader <- getFromNamespace(".load_main_bundle", "cancerRC")
  panel_loader <- getFromNamespace(".load_panel_bundle", "cancerRC")
  pathway_loader <- getFromNamespace(".load_pathway_bundle", "cancerRC")
  score_matrix <- getFromNamespace(".score_matrix", "cancerRC")

  main <- main_loader()
  panel <- panel_loader()
  pathway <- pathway_loader()

  expect_true(is.list(main))
  expect_true(is.list(panel))
  expect_true(is.list(pathway))
  pathway_input <- read_input_fixture("pathway_example_input.csv")
  expect_equal(
    dim(score_matrix(pathway_input, pathway)),
    c(nrow(pathway_input), length(pathway$features))
  )
})

test_that("Full-Matrix Multinomial Classifier returns the frozen expected classes", {
  input <- read_input_fixture("main_example_input.csv")
  expected <- read_expected_fixture("main_expected_output.csv")
  result <- predict_main(input)

  expect_type(result, "character")
  expect_length(result, nrow(input))
  expect_identical(result, as.character(expected$predicted_class))
})

test_that("15-Feature Panel Classifier returns the frozen expected classes", {
  input <- read_input_fixture("panel_example_input.csv")
  expected <- read_expected_fixture("panel_expected_output.csv")
  result <- suppressWarnings(predict_panel(input))

  expect_type(result, "character")
  expect_length(result, nrow(input))
  expect_identical(result, as.character(expected$predicted_class))
})

test_that("50-Pathway-Score Classifier returns the frozen expected classes", {
  input <- read_input_fixture("pathway_example_input.csv")
  expected <- read_expected_fixture("pathway_expected_output.csv")
  result <- suppressWarnings(predict_pathway(input))

  expect_type(result, "character")
  expect_length(result, nrow(input))
  expect_identical(result, as.character(expected$predicted_class))
})

test_that("metadata columns are never used as predictors", {
  fixtures <- list(
    list("main_example_input.csv", predict_main),
    list("panel_example_input.csv", predict_panel),
    list("pathway_example_input.csv", predict_pathway)
  )

  for (item in fixtures) {
    input <- read_input_fixture(item[[1]])
    with_metadata <- input
    with_metadata$class <- rep("not_a_predictor", nrow(with_metadata))
    with_metadata$free_text <- rep("ignored", nrow(with_metadata))
    expect_identical(
      suppressWarnings(item[[2]](input)),
      suppressWarnings(item[[2]](with_metadata))
    )
  }
})

test_that("input row order is preserved for all models", {
  fixtures <- list(
    list("main_example_input.csv", predict_main),
    list("panel_example_input.csv", predict_panel),
    list("pathway_example_input.csv", predict_pathway)
  )

  for (item in fixtures) {
    input <- read_input_fixture(item[[1]])
    reversed <- input[nrow(input):1, , drop = FALSE]
    expect_identical(
      suppressWarnings(item[[2]](reversed)),
      rev(suppressWarnings(item[[2]](input)))
    )
  }
})

test_that("Full-Matrix Multinomial Classifier rejects a partial matrix", {
  input <- read_input_fixture("main_example_input.csv")
  selected <- utils::read.csv(
    system.file("extdata", "features.csv", package = "cancerRC"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )$feature_name
  partial <- input[, c("sample_id", selected), drop = FALSE]

  expect_error(predict_main(partial), "Missing")
})

test_that("invalid feature values are rejected", {
  input <- read_input_fixture("panel_example_input.csv")
  feature_names <- setdiff(names(input), "sample_id")

  non_numeric <- input
  non_numeric[[feature_names[1]]] <- as.character(non_numeric[[feature_names[1]]])
  expect_error(predict_panel(non_numeric), "Non-numeric feature")

  non_finite <- input
  non_finite[[feature_names[1]]][1] <- NA_real_
  expect_error(predict_panel(non_finite), "finite")

  duplicated_names <- input
  names(duplicated_names)[2] <- names(duplicated_names)[3]
  expect_error(predict_panel(duplicated_names), "duplicated column")

  expect_error(predict_pathway(list()), "data.frame")
})

test_that("empty inputs return empty character vectors", {
  fixtures <- list(
    list("main_example_input.csv", predict_main),
    list("panel_example_input.csv", predict_panel),
    list("pathway_example_input.csv", predict_pathway)
  )

  for (item in fixtures) {
    input <- read_input_fixture(item[[1]])
    empty <- input[FALSE, , drop = FALSE]
    expect_identical(suppressWarnings(item[[2]](empty)), character(0))
  }
})
