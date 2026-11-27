# 50-Pathway-Score Classifier (pathway_v01)

Trained by the app developer with `scripts/30_train_pathway.R` on 2026-11-10 (created 2026-11-10T12:06:47Z).
This card records the frozen training run used by the package.

## What it is

A multinomial glmnet classifier over the **50 MSigDB Hallmark (H) ssGSEA scores** of one sample, predicting one of BRCA, KIRC, COAD, LUAD, PRAD.

Pipeline: expression profile -> `gene_XX` to Entrez (verified map) -> ssGSEA
per sample -> 50 scores -> linear classifier. Training happens in the script;
the app reads the frozen bundle and never fits anything.

## Gene sets (locked)

* Collection: MSigDB Hallmark (H), release **2024.1.Hs**
* Source: https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.entrez.gmt
* SHA-256 of the GMT: `fcf3a5471f973874f3c37996dc93a45a9c84edf7d7828bba54a7eb95af53bc45`
* Gene sets: 50; genes per set in this universe: 32-200 (median 179.5)

## ssGSEA parameters

* `normalize = FALSE` (single-sample scores must not depend on the batch), alpha = 0.25, minSize = 10, maxSize = 500.
* Single sample vs batch, 5 probe samples: max |diff| = 0e+00.

## Mapping

* Map `uci_synapse_gene_20261108` (verified): 20531 gene_XX columns -> 20531 distinct Entrez ids, 0 duplicates.
* Hallmark genes present in this universe: 7317 of 4384 set members.

## Metrics (test split, scored once)

* Balanced accuracy: **0.9895**
* Accuracy: 0.9877 over 163 test samples
* 5-fold CV balanced accuracy inside training: 0.9852 (sd 0.0112)
* Lasso lambda chosen by the same CV: 0.00317


## Limits

* The split is the fixed-seed fallback (20261115), distinct from the original full-matrix model evaluation split; comparison with the Full-Matrix Multinomial Classifier is indicative only.
* ssGSEA scores are relative ranks within one sample; they are not expression
  levels and cannot be compared across cohorts without care.
* The classifier is linear: a loading is an association, not a pathway that
  causes a cancer type.
* Probabilities are model output, not diagnoses and not clinical risk scores.
* The runtime scores with frozen coefficients; parity with glmnet on the training rows is 3.33e-16.
