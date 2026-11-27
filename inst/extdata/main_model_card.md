# Full-Matrix Multinomial Classifier (main_v01)

Derived from the frozen full-matrix model object and its training record.
No number on this card is an accuracy estimate.

## What it is

- Engine: `nnet::multinom` (multinomial logistic regression, softmax).
- Input: 100 fixed `gene_XX` features, listed in `features.csv` in model order.
- Regularisation: `decay = 0.1`.
- Fitted on: all 801 samples of the local dataset (full-data refit).
- Preprocessing: per-feature mean/SD refitted on those 801 samples; the app
  freezes them in `features.csv` and verifies them against the RDS at start-up.
- Convergence: 0. Seed: 20261111.
- Class order inside the RDS: BRCA, COAD, KIRC, LUAD, PRAD. The app maps
  probabilities by name to BRCA, KIRC, COAD, LUAD, PRAD.

## Evaluation status

**No independent test exists for this object.** A used all 801 samples to fit,
and the 100 genes had already been chosen on the training split (641 rows), so
no sample was ever held out. The earlier 160-sample evaluation belongs to the
previous 641-sample model.
Therefore this card quotes no accuracy, balanced accuracy or confusion matrix,
and the Models page shows `metrics: pending (no independent test)`.

The shipped example rows are real samples, but they are part of the 801 rows
used to fit. Their near-certain probabilities are a property of the training
rows, not evidence of generalization.

## Input scale

The original supplied `log2(norm_count+1)` scale; do not log-transform again.
Do not standardize before upload: the app applies the frozen mean/SD.

## Provenance

- Model RDS MD5: 715f38209e4729ba5cdad0dfb714bd09
- Model RDS SHA-256: f10bad87384f71ad7db235afbbadbcac6e36a37c88a7a5b8bad5c6358de6ebc6
- Training script SHA-256: d22296fd3c4d25fbb5eca318f4e9300d8b7b34d30978e01eedb2ff435c43f8f3
- Source delivery: handoff/A/full801_model_delivery (2026-11-11)
- Derived: 2026-11-13T19:30:11+0800 by scripts/11_derive_handoff.R

## Scope

Teaching classification into five known tumor classes. Not a diagnosis, not a
clinical risk score, and not validated for unknown classes.
