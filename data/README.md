# Restricted data inputs

No study data are distributed with the code. Paths are set in `config/default.R` and can be overridden in ignored `config/local.R`.

| Default input | Purpose |
|---|---|
| `raw/clindata.csv` | Raw clinical export, including enrolment records and event-level fields |
| `raw/outcomedata.csv` | Clinical-team composite MACE and CVD outcome export |
| `raw/Alec Honours.xlsx` | Study-ID reference and all molecular/pathology sheets |

The clinical compiler accepts study-ID columns named `Study.ID` or `Study.ID.` and selects `Event.Name == "Enrolment"` with no repeat instrument. It requires the clinical fields explicitly requested by `get_clin()` in `scripts/01_compile_cohort.R`. Some legacy column aliases are accepted by that helper. The endpoint-characterisation script additionally searches event-level clinical fields and writes missing-column audits.

The outcome compiler expects the two exported MACE fields corresponding to `MACE.5..HF.` and `MACE.5..revasc.`, plus `Any.CVD`, with the aliases listed in the script. These are study exports; their names are not universal definitions of MACE.

| Workbook sheet | Raw patient identifier column | Prefix in compiled predictors |
|---|---|---|
| Study IDs for Alec Project | First column, read without a header | Reference list only |
| Lipids | sourceID | `lipid_` |
| Small Molecules | Study ID | `sm_` |
| Lipoproteins | sourceID | `lipo_` |
| Bile Acids | sourceID | `ba_` |
| Pathology | Study ID | `path_` |
| Cytokines | sourceID | `cyto_` |

The compiler explicitly excludes technical metadata such as sample, plate and tube identifiers. It normalises patient IDs and collapses duplicate observations using its existing helpers. Verify that these rules are appropriate for the actual source exports.

Generated cohort files contain one row per retained participant, `Study.ID`, binary `mace` and `cvd`, clinical variables and prefixed molecular/pathology predictors. These files remain restricted. A complete validated column dictionary, units and assay metadata should be supplied by the study team; the prefix table alone is not a substitute for that dictionary.

Raw CSV encoding is handled by the original scripts: the compiler uses latin1, while endpoint characterisation uses Windows-1252. Confirm that the actual exports are compatible.
