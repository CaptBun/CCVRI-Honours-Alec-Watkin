# CCRVI research code

R code for *Development of a Multiomics Machine Learning Framework for Cardiovascular Risk Stratification in Liver Cirrhosis*.

The workflow compiles clinical and molecular data, compares weighted elastic-net models for recorded MACE status, and characterises the selected candidate. Supporting analyses cover medications, MASLD subgroups, biomarker recurrence and an optional random forest. The candidate combines expanded clinical predictors excluding aspirin, lipoproteins, lipids and small molecules.

## Setup

Use **R 4.5.3** and open `CCRVI.Rproj`, or run commands from the repository root. Install renv if necessary and restore the supplied package versions:

```r
install.packages("renv")
renv::restore(lockfile = "renv.lock")
```

The lockfile records the package environment; R itself and required system libraries must be installed separately. It includes 210 package records for the registered analyses and their dependencies. A clean restore on another computer remains to be checked.

## Restricted inputs

Place authorised source files in `data/raw/`:

- `clindata.csv`
- `outcomedata.csv`
- `Alec Honours.xlsx`

See [data/README.md](data/README.md) for sheet names, identifiers and expected fields. For alternative locations or resource settings, copy `config/local.example.R` to ignored `config/local.R`.

Study inputs, compiled cohorts, participant predictions, fitted models and generated outputs are excluded from Git. The repository provides no participant data or assumed public data-access permission.

## Check and run

From the project terminal:

```bash
Rscript tools/validate_project.R
Rscript tools/smoke_test.R
Rscript run_pipeline.R --list
Rscript run_pipeline.R main
```

The last command prints a plan. Add `--execute` only to start the selected analysis:

```bash
Rscript run_pipeline.R main --execute
```

The default exhaustive main workflow trains **191 original plus 64 aspirin-excluded model specifications**, each with nested fitting. It can take hours and is unnecessary for the structural/synthetic checks above.

| Profile | Analyses |
|---|---|
| `main` | Compilation, both elastic-net comparisons, selection, candidate, risk classes, MASLD evaluation, biomarkers and final synthesis |
| `descriptive` | Cohort EDA and endpoint characterisation |
| `sensitivity` | Clinical interrogation, medication investigation and medication sensitivity fits |
| `rf` | Complete-modality random forest |
| `masld-exploratory` | MASLD subgroup evaluation plus MASLD-only model development |
| `results` | Downstream characterisation and synthesis using completed comparison results |
| `figures` | Final synthesis and presentation export |
| `all` | All registered analyses, including random forest and MASLD-only development |

Supporting profiles require the relevant compiled cohort/results. `results` includes the candidate's interpretation refit. `figures` fits no new prediction model. MASLD-only retraining and SHAP are disabled by default; the named MASLD exploratory profile enables retraining. Optional SHAP requires separately installing fastshap and recording that addition in the lockfile.

## Resume an existing private run

With matching inputs, software, configuration and checkpoint folders, skip the completed original comparison:

```bash
Rscript run_pipeline.R main --from no-aspirin --execute
```

After both comparison result files are complete:

```bash
Rscript run_pipeline.R results --execute
```

Individual stages can be selected with `--stage NAME`. Each stage runs in a fresh R process; private logs and captured warning conditions are written under `outputs/logs/`. Changing tracked code, settings, inputs, software or the lockfile changes checkpoint provenance and may prevent reuse. Keep existing private working runs separate from the publication copy.

## Methods and validation

See [docs/METHODS.md](docs/METHODS.md) for preprocessing, validation layers, feature screening, tuning and implementation qualifications. The 28-variable expanded clinical set becomes 27 after aspirin exclusion, before dummy encoding and fold-specific preprocessing. The exploratory 0–100 score is a cohort-relative percentile, not an absolute event probability.

The tests check R syntax, workflow dependencies, warning logging and synthetic helper behaviour. Study-result reproduction requires authorised inputs and the corresponding saved results. CI runs structural and source-file checks without study data.

## Submission and citation

Review the intended Git index with `Rscript tools/check_submission.R`. Use `--all-files` only on a clean publication copy, and add `--release` to require the actual licence and lockfile. The licence remains to be agreed and added as `LICENSE`. Citation metadata are supplied in `CITATION.cff`.
