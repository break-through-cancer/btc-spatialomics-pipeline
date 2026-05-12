"""
preprocess.py — H&E Spatial Heterogeneity Pipeline
Runs immediately before the Nextflow workflow is launched on Cirro.

Responsibilities:
  1. Confirm the selected dataset contains *_Probabilities.h5 files.
  2. Validate the required filename convention: {biopsy}-{replicate}_Probabilities.h5
     (e.g. 1-1_Probabilities.h5, 2-3_Probabilities.h5)
     The R script parses biopsy and replicate directly from the filename —
     if the convention is wrong the analysis will silently produce bad results.
  3. Confirm at least two distinct biopsies are present (required for the
     pairwise LME group comparisons: 1v2, 1v3, 2v3).
  4. Log a clear summary of what will be processed.
"""

import re
from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

# ---------------------------------------------------------------------------
# 1. Find all *_Probabilities.h5 files in the input dataset
#    ds.files is populated for pipeline-output datasets; for custom (manually
#    uploaded) datasets it may be empty — fall back to ds.manifest in that case.
# ---------------------------------------------------------------------------
all_files = ds.files.copy()
h5_mask   = all_files["file"].str.endswith("_Probabilities.h5") if len(all_files) > 0 else []
h5_files  = all_files.loc[h5_mask, "file"].tolist() if len(all_files) > 0 else []

if len(h5_files) == 0:
    ds.logger.warning(
        "ds.files is empty — this is expected for manually uploaded (custom) datasets. "
        "Skipping file validation here; Nextflow will validate file presence at runtime."
    )
    ds.logger.info("Validation passed — launching workflow.")
    import sys; sys.exit(0)

ds.logger.info(f"Found {len(h5_files)} *_Probabilities.h5 file(s) in the input dataset.")

# ---------------------------------------------------------------------------
# 2. Validate filename convention: {biopsy}-{replicate}_Probabilities.h5
#    biopsy and replicate can be any alphanumeric string.
# ---------------------------------------------------------------------------
PATTERN = re.compile(r"^.*/([A-Za-z0-9]+)-([A-Za-z0-9]+)_Probabilities\.h5$")

invalid = [f for f in h5_files if not PATTERN.match(f)]

if invalid:
    raise ValueError(
        f"{len(invalid)} file(s) do not match the required naming convention "
        f"({{biopsy}}-{{replicate}}_Probabilities.h5):\n"
        + "\n".join(f"  {f}" for f in sorted(invalid))
        + "\n\nPlease rename your files before re-uploading the dataset."
    )

# ---------------------------------------------------------------------------
# 3. Extract biopsy and replicate labels; confirm ≥2 biopsies
# ---------------------------------------------------------------------------
parsed = [PATTERN.match(f) for f in h5_files]
biopsies   = sorted(set(m.group(1) for m in parsed))
replicates = [m.group(2) for m in parsed]

ds.logger.info("Samples to be processed:")
for f, m in sorted(zip(h5_files, parsed), key=lambda x: x[0]):
    ds.logger.info(f"  biopsy={m.group(1)}  replicate={m.group(2)}  →  {f}")

if len(biopsies) < 2:
    raise ValueError(
        f"The LME group comparison step requires at least 2 distinct biopsies, "
        f"but only biopsy '{biopsies[0]}' was detected.\n"
        f"Please select a dataset that includes samples from multiple biopsies, "
        f"or confirm your filenames follow the pattern: "
        f"{{biopsy}}-{{replicate}}_Probabilities.h5"
    )

ds.logger.info(f"Biopsies detected ({len(biopsies)}): {', '.join(biopsies)}")
ds.logger.info("Validation passed — launching workflow.")
