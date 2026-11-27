# 15-Feature Panel Classifier (panel_v01)

Trained by the app developer with `scripts/20_train_panel.R` on 2026-11-09 (created 2026-11-09T14:14:17Z).
This card records the frozen training run used by the package.

## What it is

A multinomial glmnet classifier over **K = 15** features (`gene_XX` dummy keys), one of the five classes BRCA, KIRC, COAD, LUAD, PRAD.
The app never trains: it reads the frozen bundle from `app/models/` and scores
with the stored center, scale and coefficient matrix.

## Features (bundle order)

1. `gene_9175` - stability 100%, training P1-P99 0.00-17.72
2. `gene_9176` - stability 100%, training P1-P99 0.00-19.45
3. `gene_220` - stability 100%, training P1-P99 0.00-13.08
4. `gene_219` - stability 100%, training P1-P99 0.00-13.29
5. `gene_18135` - stability 100%, training P1-P99 0.00-13.93
6. `gene_13818` - stability 100%, training P1-P99 0.00-9.52
7. `gene_14114` - stability 100%, training P1-P99 0.00-10.03
8. `gene_3737` - stability 100%, training P1-P99 0.00-12.81
9. `gene_450` - stability 100%, training P1-P99 0.00-11.48
10. `gene_7964` - stability 100%, training P1-P99 0.00-12.48
11. `gene_15895` - stability 100%, training P1-P99 0.00-12.68
12. `gene_6733` - stability 100%, training P1-P99 0.00-14.60
13. `gene_3439` - stability 100%, training P1-P99 0.00-13.21
14. `gene_12995` - stability 100%, training P1-P99 0.00-14.68
15. `gene_16169` - stability 100%, training P1-P99 0.00-13.57

## How it was selected

* Split: fixed-seed stratified 80/20 fallback (seed 20261115) - 638 train / 163 test rows.
* Variance filter on the training rows removed 278 constant features; a one-way F test ranked the rest.
* Stability selection: 20 subsamples of 80% of the training rows, top 50 features per repeat, frequency kept for the ranking.
* Panel sizes 10, 15, 20 were compared by 5-fold CV balanced accuracy inside the training split;
  the smallest K within 0.01 of the best was kept (K = 15).
* Lasso penalty chosen by the same CV: lambda = 0.002.

## Metrics (test split, scored once)

* Balanced accuracy: **0.9898**
* Accuracy: 0.9877 over 163 test samples
* 5-fold CV balanced accuracy inside training: 0.9845 (sd 0.0108)


## Limits

* The split is the fixed-seed fallback (20261115), distinct from the original full-matrix model evaluation split; comparison with the Full-Matrix Multinomial Classifier is indicative only.
* The classes are imbalanced (300/146/78/141/136); accuracy and balanced accuracy differ.
* `gene_XX` is a placeholder key, not a gene symbol; the panel is a teaching artifact.
* Probabilities are model output, not diagnoses and not clinical risk scores.
* The runtime scores with frozen coefficients; parity with glmnet on the training rows is 2.22e-16.
