# Shared prediction API for the three frozen BIO215 classifiers.
# The helpers below are package-private and are not exported.
# Their structured comments use @noRd and do not create user help pages.

.bio_classes <- c("BRCA", "KIRC", "COAD", "LUAD", "PRAD")
.bio_cache <- new.env(parent = emptyenv())

#' Resolve a package resource in an installed package or source tree.
#'
#' This internal helper avoids machine-specific paths.
#' Tests: `tests/testthat/test-models.R`.
#' @param filename A file name under `inst/extdata`.
#' @returns An existing local path.
#' @keywords internal
#' @noRd
.bio_resource <- function(filename) {
  installed_path <- system.file(
    "extdata", filename, package = "cancerRC", mustWork = FALSE
  )
  if (nzchar(installed_path) && file.exists(installed_path)) {
    return(installed_path)
  }

  source_path <- file.path("inst", "extdata", filename)
  if (file.exists(source_path)) {
    return(source_path)
  }

  stop(
    "Package resource not found: ", filename,
    ". Install the package or run from its source directory.",
    call. = FALSE
  )
}

#' Validate a model input table.
#'
#' This internal helper checks feature presence, uniqueness, numeric type, and
#' finite values. Non-feature columns are ignored.
#' Tests: `tests/testthat/test-models.R`.
#' @param newdata A candidate input data frame.
#' @param required_features Ordered feature names required by a model.
#' @returns `TRUE` invisibly, or an error describing the failed check.
#' @keywords internal
#' @noRd
.validate_input <- function(newdata, required_features) {
  if (!is.data.frame(newdata)) {
    stop("newdata must be a data.frame.", call. = FALSE)
  }
  if (anyDuplicated(names(newdata))) {
    stop("newdata has duplicated column names.", call. = FALSE)
  }

  missing_features <- setdiff(required_features, names(newdata))
  if (length(missing_features)) {
    stop(
      "Missing ", length(missing_features), " required feature column(s): ",
      paste(utils::head(missing_features, 5), collapse = ", "),
      if (length(missing_features) > 5L) ", ..." else "", ".",
      call. = FALSE
    )
  }

  feature_data <- newdata[, required_features, drop = FALSE]
  non_numeric <- required_features[
    !vapply(feature_data, is.numeric, logical(1))
  ]
  if (length(non_numeric)) {
    stop(
      "Non-numeric feature column(s): ",
      paste(utils::head(non_numeric, 5), collapse = ", "), ".",
      call. = FALSE
    )
  }

  values <- as.matrix(feature_data)
  storage.mode(values) <- "double"
  if (any(!is.finite(values))) {
    stop("Feature values must be finite (no NA, NaN or Inf).", call. = FALSE)
  }

  invisible(TRUE)
}

#' Warn about values outside a fitted reference range.
#'
#' This check is advisory. It never blocks a valid numeric prediction.
#' Tests: `tests/testthat/test-models.R`.
#' @param newdata A validated input data frame.
#' @param feature_names Ordered model feature names.
#' @param reference_range A table containing feature names, `p01`, and `p99`.
#' @param label A model label used in the warning.
#' @returns `NULL` invisibly.
#' @keywords internal
#' @noRd
.warn_reference_range <- function(newdata, feature_names, reference_range, label) {
  if (!is.data.frame(reference_range) ||
      !all(c("p01", "p99") %in% names(reference_range))) {
    return(invisible(NULL))
  }

  ref_name <- if ("feature" %in% names(reference_range)) {
    "feature"
  } else if ("pathway" %in% names(reference_range)) {
    "pathway"
  } else {
    return(invisible(NULL))
  }

  index <- match(feature_names, as.character(reference_range[[ref_name]]))
  if (anyNA(index)) {
    return(invisible(NULL))
  }

  values <- as.matrix(newdata[, feature_names, drop = FALSE])
  lower <- as.numeric(reference_range$p01[index])
  upper <- as.numeric(reference_range$p99[index])
  outside <- sweep(values, 2L, lower, "<") |
    sweep(values, 2L, upper, ">")
  count <- sum(outside)

  if (count > 0L) {
    warning(
      label, ": ", count,
      " value(s) fall outside the fitted P1-P99 reference bands. ",
      "They are retained and scored.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

#' Load and validate the Full-Matrix Multinomial Classifier bundle.
#'
#' The loader checks the selected feature order, full matrix schema, frozen
#' preprocessing values, class order, and the serialized model class.
#' Tests: `tests/testthat/test-models.R`.
#' @returns A validated internal model object.
#' @keywords internal
#' @noRd
.load_main_bundle <- function() {
  if (exists("m0", .bio_cache, inherits = FALSE)) {
    return(get("m0", .bio_cache, inherits = FALSE))
  }
  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("The nnet package is required for the Full-Matrix Multinomial Classifier model.", call. = FALSE)
  }

  bundle <- readRDS(.bio_resource("main_model_v01.rds"))
  selected_schema <- utils::read.csv(
    .bio_resource("features.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  full_schema <- utils::read.csv(
    .bio_resource("full_matrix_features.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required_bundle_fields <- c("model", "preprocess", "classes", "n_training")
  missing_bundle_fields <- setdiff(required_bundle_fields, names(bundle))
  if (length(missing_bundle_fields)) {
    stop(
      "The Full-Matrix Multinomial Classifier bundle is missing field(s): ",
      paste(missing_bundle_fields, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!inherits(bundle$model, "multinom")) {
    stop("The Full-Matrix Multinomial Classifier bundle is not an nnet::multinom object.", call. = FALSE)
  }

  model_features <- as.character(bundle$preprocess$genes)
  center <- as.numeric(bundle$preprocess$mean)
  scale <- as.numeric(bundle$preprocess$sd)
  full_features <- as.character(full_schema$feature_name)

  if (length(full_features) != 20531L || anyDuplicated(full_features)) {
    stop("The full-matrix feature schema is invalid.", call. = FALSE)
  }
  if (!all(model_features %in% full_features)) {
    stop("The Full-Matrix Multinomial Classifier selected features are absent from the full-matrix schema.", call. = FALSE)
  }
  if (length(center) != length(model_features) ||
      length(scale) != length(model_features) ||
      any(!is.finite(center)) || any(!is.finite(scale)) || any(scale <= 0)) {
    stop("The Full-Matrix Multinomial Classifier preprocessing values are invalid.", call. = FALSE)
  }
  if (!is.data.frame(selected_schema) ||
      !identical(as.character(selected_schema$feature_name), model_features) ||
      max(abs(as.numeric(selected_schema$center) - center)) > 1e-12 ||
      max(abs(as.numeric(selected_schema$scale) - scale)) > 1e-12) {
    stop("features.csv does not match the Full-Matrix Multinomial Classifier bundle.", call. = FALSE)
  }
  if (!identical(as.character(bundle$classes), c("BRCA", "COAD", "KIRC", "LUAD", "PRAD"))) {
    stop("The Full-Matrix Multinomial Classifier class order is invalid.", call. = FALSE)
  }

  model <- list(
    bundle = bundle,
    model_features = model_features,
    full_features = full_features,
    selected_schema = selected_schema
  )
  assign("m0", model, .bio_cache)
  model
}

#' Predict cancer type with the Full-Matrix Multinomial Classifier
#'
#' Apply the BIO215 `main_v01` multinomial model. The input must contain the
#' complete 20,531-column gene-expression matrix. The model uses the 100
#' selected features stored in `features.csv`.
#'
#' @param newdata A data frame containing `sample_id` optionally and all 20,531
#'   `gene_XX` columns. The input values must remain on the original supplied
#'   `log2(norm_count + 1)` scale.
#'
#' @returns A character vector with one predicted tumor class per input row.
#'   The input row order is preserved.
#'
#' @export
#'
#' @examples
#' example_file <- system.file(
#'   "extdata", "main_example_input.csv", package = "cancerRC"
#' )
#' example_input <- utils::read.csv(
#'   example_file, check.names = FALSE, stringsAsFactors = FALSE
#' )
#' predict_main(example_input)
predict_main <- function(newdata) {
  model <- .load_main_bundle()
  .validate_input(newdata, model$full_features)

  if (nrow(newdata) == 0L) {
    return(character(0))
  }

  x <- as.matrix(newdata[, model$model_features, drop = FALSE])
  storage.mode(x) <- "double"
  z <- sweep(
    sweep(x, 2L, as.numeric(model$bundle$preprocess$mean), "-"),
    2L,
    as.numeric(model$bundle$preprocess$sd),
    "/"
  )

  as.character(
    stats::predict(
      model$bundle$model,
      newdata = as.data.frame(z, check.names = FALSE),
      type = "class"
    )
  )
}

#' Load and validate a coefficient-based classifier bundle.
#'
#' This helper checks schema order, centering, scaling, coefficients, and class
#' order before prediction.
#' Tests: `tests/testthat/test-models.R`.
#' @param bundle_file Bundle file name.
#' @param schema_file Feature schema file name.
#' @param feature_field Bundle field containing ordered features.
#' @param center_field Bundle field containing frozen centers.
#' @param scale_field Bundle field containing frozen scales.
#' @param expected_features Expected feature count.
#' @returns A validated internal model object.
#' @keywords internal
#' @noRd
.load_linear_bundle <- function(
    bundle_file,
    schema_file,
    feature_field,
    center_field,
    scale_field,
    expected_features
) {
  bundle <- readRDS(.bio_resource(bundle_file))
  required <- c(
    "model_id", "bundle_version", "classes", "coefficients",
    feature_field, center_field, scale_field, "reference_range"
  )
  missing <- setdiff(required, names(bundle))
  if (length(missing)) {
    stop(
      "The classifier bundle is missing field(s): ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }

  features <- as.character(bundle[[feature_field]])
  center <- as.numeric(bundle[[center_field]])
  scale <- as.numeric(bundle[[scale_field]])
  coefficients <- bundle$coefficients

  if (length(features) != expected_features || anyDuplicated(features)) {
    stop("The classifier feature schema is invalid.", call. = FALSE)
  }
  if (length(center) != length(features) || length(scale) != length(features) ||
      any(!is.finite(center)) || any(!is.finite(scale)) || any(scale <= 0)) {
    stop("The classifier centering and scaling values are invalid.", call. = FALSE)
  }
  if (!identical(as.character(bundle$classes), .bio_classes)) {
    stop("The classifier class order is invalid.", call. = FALSE)
  }
  if (!is.matrix(coefficients) || !is.numeric(coefficients) ||
      !identical(rownames(coefficients), c("(Intercept)", features)) ||
      !identical(colnames(coefficients), .bio_classes) ||
      any(!is.finite(coefficients))) {
    stop("The classifier coefficient matrix is invalid.", call. = FALSE)
  }

  schema <- utils::read.csv(
    .bio_resource(schema_file),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  schema_tolerance <- 1e-8
  if (!all(c("feature_name", "center", "scale") %in% names(schema)) ||
      nrow(schema) != length(features) ||
      !identical(as.character(schema$feature_name), features) ||
      max(abs(as.numeric(schema$center) - center)) > schema_tolerance ||
      max(abs(as.numeric(schema$scale) - scale)) > schema_tolerance) {
    stop("The feature schema does not match the classifier bundle.", call. = FALSE)
  }

  list(
    bundle = bundle,
    features = features,
    center = center,
    scale = scale
  )
}

#' Load the 15-Feature Panel Classifier bundle.
#'
#' This is an internal wrapper around `.load_linear_bundle()`.
#' Tests: `tests/testthat/test-models.R`.
#' @returns A validated internal model object.
#' @keywords internal
#' @noRd
.load_panel_bundle <- function() {
  if (exists("m1", .bio_cache, inherits = FALSE)) {
    return(get("m1", .bio_cache, inherits = FALSE))
  }
  model <- .load_linear_bundle(
    bundle_file = "panel_bundle_v01.rds",
    schema_file = "panel_features.csv",
    feature_field = "feature_names",
    center_field = "center",
    scale_field = "scale",
    expected_features = 15L
  )
  assign("m1", model, .bio_cache)
  model
}

#' Load the 50-Pathway-Score Classifier bundle.
#'
#' This is an internal wrapper around `.load_linear_bundle()`.
#' Tests: `tests/testthat/test-models.R`.
#' @returns A validated internal model object.
#' @keywords internal
#' @noRd
.load_pathway_bundle <- function() {
  if (exists("m2", .bio_cache, inherits = FALSE)) {
    return(get("m2", .bio_cache, inherits = FALSE))
  }
  model <- .load_linear_bundle(
    bundle_file = "pathway_bundle_v01.rds",
    schema_file = "pathway_features.csv",
    feature_field = "pathway_names",
    center_field = "score_center",
    scale_field = "score_scale",
    expected_features = 50L
  )
  assign("m2", model, .bio_cache)
  model
}

#' Score a validated coefficient-based classifier.
#'
#' This helper standardizes the selected inputs, applies the frozen coefficient
#' matrix, and returns class labels. Reference-range warnings are advisory.
#' Tests: `tests/testthat/test-models.R`.
#' @param newdata A candidate input data frame.
#' @param model A validated internal model object.
#' @param label A model label used in range warnings.
#' @returns A character vector of class labels.
#' @keywords internal
#' @noRd
.predict_linear_bundle <- function(newdata, model, label) {
  .validate_input(newdata, model$features)
  .warn_reference_range(
    newdata,
    model$features,
    model$bundle$reference_range,
    label
  )

  if (nrow(newdata) == 0L) {
    return(character(0))
  }

  x <- as.matrix(newdata[, model$features, drop = FALSE])
  storage.mode(x) <- "double"
  z <- sweep(sweep(x, 2L, model$center, "-"), 2L, model$scale, "/")
  coefficients <- model$bundle$coefficients
  eta <- z %*% coefficients[model$features, , drop = FALSE]
  eta <- sweep(eta, 2L, coefficients["(Intercept)", ], "+")
  eta <- sweep(eta, 1L, apply(eta, 1L, max), "-")
  probabilities <- exp(eta)
  probabilities <- probabilities / rowSums(probabilities)

  as.character(.bio_classes[max.col(probabilities, ties.method = "first")])
}

#' Predict cancer type with the 15-Feature Panel Classifier
#'
#' Apply the BIO215 `panel_v01` classifier to the 15 selected gene-expression
#' features used by the Shiny panel routes.
#'
#' @param newdata A data frame containing the 15 required panel features. A
#'   `sample_id` column and other non-feature columns are allowed and ignored.
#'
#' @returns A character vector with one predicted tumor class per input row.
#'   Values outside the fitted P1-P99 bands generate a warning and are retained.
#'
#' @export
#'
#' @examples
#' example_file <- system.file(
#'   "extdata", "panel_example_input.csv", package = "cancerRC"
#' )
#' example_input <- utils::read.csv(
#'   example_file, check.names = FALSE, stringsAsFactors = FALSE
#' )
#' predict_panel(example_input)
predict_panel <- function(newdata) {
  model <- .load_panel_bundle()
  .predict_linear_bundle(newdata, model, "15-Feature Panel Classifier panel input")
}

#' Predict cancer type with the 50-Pathway-Score Classifier
#'
#' Apply the BIO215 `pathway_v01` classifier to 50 precomputed MSigDB Hallmark
#' ssGSEA scores. The function is the lightweight 50-Pathway-Score Classifier entry point.
#'
#' @param newdata A data frame containing the 50 required `HALLMARK_*` score
#'   columns. A `sample_id` column and other non-feature columns are allowed and
#'   ignored.
#'
#' @returns A character vector with one predicted tumor class per input row.
#'   Values outside the fitted P1-P99 bands generate a warning and are retained.
#'
#' @export
#'
#' @examples
#' example_file <- system.file(
#'   "extdata", "pathway_example_input.csv", package = "cancerRC"
#' )
#' example_input <- utils::read.csv(
#'   example_file, check.names = FALSE, stringsAsFactors = FALSE
#' )
#' predict_pathway(example_input)
predict_pathway <- function(newdata) {
  model <- .load_pathway_bundle()
  .predict_linear_bundle(newdata, model, "50-Pathway-Score Classifier pathway-score input")
}

#' Standardize pathway scores for the pathway classifier.
#'
#' This internal helper mirrors the score-matrix preprocessing step. The public
#' pathway API applies the same operation inside `.predict_linear_bundle()`.
#' Tests: `tests/testthat/test-models.R`.
#' @param newdata A data frame containing pathway score columns.
#' @param model A validated pathway model object.
#' @returns A standardized numeric matrix.
#' @keywords internal
#' @noRd
.score_matrix <- function(newdata, model) {
  x <- as.matrix(newdata[, model$features, drop = FALSE])
  storage.mode(x) <- "double"
  sweep(sweep(x, 2L, model$center, "-"), 2L, model$scale, "/")
}
