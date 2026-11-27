# Cancer RNA-seq Classifier model contract

Package name: `cancerRC`

The package carries the three frozen classifiers used by the BIO215 Shiny
application. Each model has a separate schema and a separate public function.

## Full-Matrix Multinomial Classifier

- Model ID: `main_v01`.
- Engine: `nnet::multinom` multinomial logistic regression.
- Training rows: 801.
- Input: all 20,531 `gene_XX` columns.
- Model features: 100 columns listed in `features.csv`.
- Input scale: original supplied `log2(norm_count + 1)` scale.
- Function: `predict_main()`.
- Independent test: not available for this full-801 refit.

The target label is not required. It is never used as a predictor.

## 15-Feature Panel Classifier

- Model ID: `panel_v01`.
- Engine: frozen multinomial `glmnet` coefficients.
- Training split: 638 training rows and 163 holdout rows.
- Input: 15 selected `gene_XX` columns.
- Function: `predict_panel()`.

The bundle stores the feature order, centering, scaling, coefficients, and
reference ranges. Values outside the fitted P1-P99 bands produce warnings.
They are not blocked.

## 50-Pathway-Score Classifier

- Model ID: `pathway_v01`.
- Engine: frozen multinomial `glmnet` coefficients.
- Training split: 638 training rows and 163 holdout rows.
- Input: 50 precomputed MSigDB Hallmark ssGSEA scores.
- Function: `predict_pathway()`.
- ssGSEA training setting: `normalize = FALSE`.

The package exports the lightweight score-matrix entry point. It does not
export the full-profile pathway-scoring workflow. The Shiny application owns
that route outside the package contract.

## Shared input rules

All public functions require a `data.frame`. Required features must be present,
numeric, non-duplicated, and finite. Missing columns, `NA`, `NaN`, and `Inf`
values produce errors. Input row order is preserved.

`sample_id` and other non-feature columns are ignored. A `class` column is never
used as a predictor. All package toy inputs omit the target label.

## Example fixtures

- `main_example_input.csv` demonstrates Full-Matrix Multinomial Classifier.
- `panel_example_input.csv` demonstrates 15-Feature Panel Classifier.
- `pathway_example_input.csv` demonstrates 50-Pathway-Score Classifier.
- The matching `*_expected_output.csv` files record expected class labels.

The examples are software fixtures. They are not new independent evaluation
data.

## Performance scope

Full-Matrix Multinomial Classifier has no independent test. Its independent-test metrics are unavailable.
15-Feature Panel Classifier and 50-Pathway-Score Classifier have frozen 163-row holdout metrics. Metrics remain model-specific.
15-Feature Panel Classifier or 50-Pathway-Score Classifier values must never be assigned to Full-Matrix Multinomial Classifier.

## Safety scope

This package is a BIO215 teaching and research-software demonstration. It is
not a diagnostic tool. It is not a clinical risk score. It is not validated for
unknown tumor classes or clinical deployment.
