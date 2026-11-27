---

editor_options: 
  markdown: 
    wrap: 72
---

# Cancer RNA-seq Classifier {#cancer-rna-seq-classifier}

- [Cancer RNA-seq Classifier](#cancer-rna-seq-classifier)
  - [Acknowledgments](#acknowledgments)
  - [Installation](#installation)
  - [Exported functions](#exported-functions)
  - [Full-Matrix Multinomial Classifier example](#full-matrix-multinomial-classifier-example)
  - [15-Feature Panel Classifier example](#15-feature-panel-classifier-example)
  - [50-Pathway-Score Classifier example](#50-pathway-score-classifier-example)
  - [Model contracts](#model-contracts)
  - [Model Performance Summary](#model-performance-summary)
  - [Reference Evaluation of the Original Full-Matrix Model](#reference-evaluation-of-the-original-full-matrix-model)
  - [Package resources](#package-resources)
  - [Scope and limitations](#scope-and-limitations)
  - [Data source](#data-source)
  - [License](#license)

# Cancer RNA-seq Classifier

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/license/MIT)

**Package name:** `cancerRC`

[Cancer RNA-seq datasets](https://www.kaggle.com/datasets/waalbannyantudre/gene-expression-cancer-rna-seq-donated-on-682016) contain thousands of gene-expression variables, which makes tumor-class prediction and model deployment difficult. This package addresses the problem by providing prediction functions for three pre-trained classifiers. Its inputs include 20,531 `gene_XX` expression columns, 15 selected gene features, or 50 Hallmark ssGSEA pathway scores (depending on the classifier). The cancer classes are `BRCA`, `KIRC`, `COAD`, `LUAD`, and `PRAD`.

The package covers the three classifiers used in the prediction module of our [web app](https://bio215-leyang.shinyapps.io/Group28_Cancer-rnaseq-db/):

- the Full-Matrix Multinomial Classifier;
- the 15-Feature Panel Classifier;
- the 50-Pathway-Score Classifier.

The package provides these pre-trained models. It does not train any additional model at runtime. It returns one predicted tumor class per sample row.

Note: This package was developed for the group coursework of the BIO215 project.

## Acknowledgments {#acknowledgments}

We acknowledge the frozen full-matrix model and its associated original 80/20 evaluation materials. The reference ROC and precision-recall curves, AUROC, AUPRC, and accuracy values describe that earlier full-matrix model evaluation. They do not describe the full-data refit shipped in `predict_main()`.

## Installation {#installation}

Install the package from GitHub with `remotes`:

Note: Install the package `remotes` if you haven’t.

``` r
install.packages("remotes")
remotes::install_github(
  "rakuyonn28-ux/cancerRC",
  upgrade = "never"
)
```

Alternatively, install the package from GitHub with `devtools`:

``` r
devtools::install_github("rakuyonn28-ux/cancerRC")
```

## Exported functions {#exported-functions}

| Function | Model | Input | Output |
|:-----------------|:-----------------|:-----------------|:-----------------|
| `predict_main()` | Full-Matrix Multinomial Classifier `main_v01` | 20,531 gene-expression columns | Character class labels |
| `predict_panel()` | 15-Feature Panel Classifier `panel_v01` | 15 selected gene columns | Character class labels |
| `predict_pathway()` | 50-Pathway-Score Classifier `pathway_v01` | 50 Hallmark ssGSEA scores | Character class labels |

All functions preserve input row order. The `sample_id` column is optional. It is not used as a predictor. A `class` column is also ignored.

## Full-Matrix Multinomial Classifier example {#full-matrix-multinomial-classifier-example}

This model’s input is all 20,531 `gene_XX` columns. It uses the 100 columns listed in `features.csv`.

``` r
library(cancerRC)

full_matrix_file <- system.file(
  "extdata", "main_example_input.csv", package = "cancerRC"
)
full_matrix_input <- utils::read.csv(
  full_matrix_file, check.names = FALSE, stringsAsFactors = FALSE
)

full_matrix_prediction <- predict_main(full_matrix_input)
data.frame(
  sample_id = full_matrix_input$sample_id,
  predicted_class = full_matrix_prediction
)
```

The Full-Matrix Multinomial Classifier input remains on the original supplied `log2(norm_count + 1)` scale. The function applies the frozen means and standard deviations internally. Please do not log-transform or standardize the data a second time.

## 15-Feature Panel Classifier example

15-Feature Panel Classifier uses 15 selected `gene_XX` features.

``` r
panel_file <- system.file(
  "extdata", "panel_example_input.csv", package = "cancerRC"
)
panel_input <- utils::read.csv(
  panel_file, check.names = FALSE, stringsAsFactors = FALSE
)

panel_prediction <- predict_panel(panel_input)
data.frame(
  sample_id = panel_input$sample_id,
  predicted_class = panel_prediction
)
```

The model bundle stores the feature order, centers, scales, coefficients, and reference ranges. Values outside the fitted P1-P99 bands will generate a warning but are retained for scoring. Missing, duplicated, non-numeric, and non-finite values will stop the prediction with an error message.

## 50-Pathway-Score Classifier example

50-Pathway-Score Classifier uses 50 precomputed MSigDB Hallmark ssGSEA scores. The package exposes this lightweight score-matrix entry point.

``` r
pathway_file <- system.file(
  "extdata", "pathway_example_input.csv", package = "cancerRC"
)
pathway_input <- utils::read.csv(
  pathway_file, check.names = FALSE, stringsAsFactors = FALSE
)

pathway_prediction <- predict_pathway(pathway_input)
data.frame(
  sample_id = pathway_input$sample_id,
  predicted_class = pathway_prediction
)
```

The package does not export the full-profile pathway-scoring route. It accepts precomputed pathway scores only. The Shiny application maintains its own full-profile transformation outside this package contract. The package keeps the classifier entry point small and independent of application runtime code.

## Model contracts {#model-contracts}

| Model | Algorithm | Training or holdout scope | Required input |
|:-----------------|:-----------------|:-----------------|:-----------------|
| Full-Matrix Multinomial Classifier `main_v01` | `nnet::multinom` | 801-row full-data refit | All 20,531 gene columns |
| 15-Feature Panel Classifier `panel_v01` | Frozen multinomial `glmnet` coefficients | 638 train, 163 holdout | 15 selected gene columns |
| 50-Pathway-Score Classifier `pathway_v01` | Frozen multinomial `glmnet` coefficients | 638 train, 163 holdout | 50 Hallmark score columns |

The five fixed classes are `BRCA`, `KIRC`, `COAD`, `LUAD`, and `PRAD`.

## Model Performance Summary {#model-performance-summary}

The Full-Matrix Multinomial Classifier has no independent test set. The other two classifiers use the fixed-seed stratified 80/20 fallback split.

Macro-F1 is the primary metric because it weights the five classes equally when their sample counts differ. The test set contains 163 samples. All reported values are point estimates from a single internal holdout split.

| Model | Macro-F1 (primary) | Accuracy | Balanced accuracy | Macro AUROC | Macro AUPRC (AP) | Test size |
|:----------|----------:|----------:|----------:|----------:|----------:|----------:|
| Full-Matrix Multinomial Classifier | Not available | Not available | Not available | Not available | Not available | No independent test |
| 15-Feature Panel Classifier | 0.9854 | 0.9877 | 0.9898 | 0.9996 | 0.9972 | 163 |
| 50-Pathway-Score Classifier | 0.9853 | 0.9877 | 0.9895 | 0.9996 | 0.9974 | 163 |

AUROC uses one-vs-rest evaluation. Macro AUROC is the unweighted mean across the five classes. AUPRC uses average precision for each one-vs-rest problem. Macro AUPRC is the unweighted mean across the five classes. The random precision-recall baseline is the class prevalence in the holdout set.

[Read the detailed performance report](PERFORMANCE_REPORT.md), including the metric definitions, aggregate model metrics, and evaluation scope.

The Full-Matrix Multinomial Classifier used all 801 rows for the full-data refit. Historical metrics from another version must not be assigned to this object.

The 15-Feature Panel Classifier and 50-Pathway-Score Classifier use a fixed-seed stratified fallback split. This split is distinct from the original full-matrix model evaluation split. Their metrics should not be compared directly with the full-data refit.

## Reference Evaluation of the Original Full-Matrix Model {#reference-evaluation-of-the-original-full-matrix-model}

This section is reserved for the earlier full-matrix model evaluation. It will include aggregate ROC and precision-recall figures, AUROC, AUPRC, accuracy, and balanced accuracy after the source evaluation record is supplied. No per-class performance is required for this reference section.

| Metric | Value |
|:-----------------------------------|:-----------------------------------|
| AUROC | To be added from the original full-matrix model evaluation |
| AUPRC | To be added from the original full-matrix model evaluation |
| Accuracy | To be added from the original full-matrix model evaluation |
| Balanced accuracy | To be added from the original full-matrix model evaluation |
| Test size | To be confirmed from the evaluation record |
| Aggregate ROC figure | To be added |
| Aggregate precision-recall figure | To be added |

The reference results will describe the original full-matrix model evaluation. They will not describe the full-data refit shipped in `predict_main()`.

## Package resources {#package-resources}

The installed package includes:

- Full-Matrix Multinomial Classifier model, selected-feature schema, and full-matrix feature schema;
- 15-Feature Panel Classifier panel bundle and 15-feature schema;
- 50-Pathway-Score Classifier pathway bundle and 50-score schema;
- label-free toy inputs and expected class outputs;
- model cards, metric records, and the model contract.

Use `system.file("extdata", ..., package = "cancerRC")` to locate installed resources. The detailed performance report and its figures are available in the source repository; their presence in an installed package is not required for prediction.

## Scope and limitations {#scope-and-limitations}

This package was developed for the BIO215 group project.

It is not a diagnostic tool. It is not a clinical risk score. It is not validated for unknown tumor classes or clinical deployment.

## Data source {#data-source}

The project uses the UCI Gene Expression Cancer RNA-seq dataset distributed through Kaggle:

<https://www.kaggle.com/datasets/waalbannyantudre/gene-expression-cancer-rna-seq-donated-on-682016>

## License {#license}

MIT. See `LICENSE.md` for the full terms.
