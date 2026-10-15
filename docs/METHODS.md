# What the supplied code implements

This document describes the supplied implementation. Source data and participant-level results are not distributed.

## Compilation

The compilation script uses enrolment clinical records, normalises study identifiers and joins outcomes and workbook modalities. It builds the complete intersection across the study-ID reference, clinical/outcome data and all required molecular/pathology sources. Numeric replicates are collapsed using the existing helper. Features with more than 50% missingness across this compiled cohort are removed before modelling. Both filtered and unfiltered cohort files are retained locally.

The binary `mace` field is constructed from the logical OR of the two exported MACE composite fields. The endpoint and any prevalent-event exclusions must be verified against the study records; they are not inferred from previous descriptions of the project.

## Elastic-net comparisons

The scripts create five stratified global validation folds. For each held-out fold, the remaining development data are divided into five stratified development folds. Numeric medians, categorical encoding, zero-variance removal, centring, scaling and outcome-informed molecular AUC screening are fitted using the corresponding development training subset.

Clinical predictors are retained outside the molecular AUC screen. Eligible molecular variables are ranked by `max(AUC, 1 - AUC)` and the selected K applies across eligible molecular predictors in that model specification. The original alpha grid is 0.25, 0.5, 0.75 and 1; the original K grid is 20, 30, 50, 75, 100, 150 and 200.

For each alpha/K fit, `cv.glmnet()` performs its own five-fold CV to choose `lambda.min` using AUC. This additional lambda-selection step receives an already preprocessed and screened training matrix. Each development fold retains its winning alpha/K pair. The final pair is chosen by aggregating those winning rows, rather than averaging every grid combination over every development fold. These behaviours have been preserved.

The selected strategy is refitted on the full development subset and evaluated on the untouched global validation fold. Positive observations receive a weight equal to the training negative/positive class ratio; negative observations receive weight one. Saved out-of-fold predictions supply pooled performance summaries and downstream comparisons.

## Downstream analyses

Candidate characterisation combines out-of-fold performance and feature-selection summaries with an interpretation-only full-cohort refit. Exploratory classification uses the pooled Youden threshold. The three-class framework converts model scores to cohort-relative percentiles and adds a distributional lower boundary. MASLD subgroup evaluation reuses the candidate predictions; optional MASLD-only development is a separate nested analysis.

The optional random forest uses its own tuning for K, mtry and minimum node size and retains the original tree counts. It is not required by the main elastic-net workflow.

## Implementation qualifications

- The compilation uses `any(..., na.rm = TRUE)` to collapse exported outcomes. Confirm all-missing composite records are handled as intended; this construction can otherwise produce zero.
- Clinical predictor lists are intersected with the available columns. Confirm retained counts and coding against the study dictionary.
- Several ROC calls use automatic direction selection. This is retained from the original analysis and should be reported consistently when interpreting probability predictions.
- The internal model search, pooled Youden threshold and percentile classes are exploratory. Validation of prospective clinical use requires a separate study.

These details are preserved for fidelity to the original implementation. Changes to statistical calculations require a separate, re-evaluated analysis version.
