# Reproducibility checklist

## Static validation

From the repository root, run:

```bash
python -m json.tool notebooks/16s_colab_pipeline_tested.ipynb >/dev/null
python -m py_compile third_party/diet-microbiome-depression-ml/*.py
```

The notebook contains Colab shell and magic cells, so ordinary Python compilation is not a complete notebook execution test. Open and run it from a fresh Colab session for an end-to-end check.

## Record for each run

Save a small run manifest containing:

- Git commit SHA;
- notebook filename and modification date;
- QIIME 2 and Nextflow versions;
- classifier release and SHA-256 checksum;
- input metadata checksum;
- branch and all key quality/feature thresholds;
- random states;
- `pipeline_info/trace.txt`, `timeline.html`, and `report.html`;
- output checksum manifest.

## Clean-room test

1. Clone the repository into a new directory.
2. Copy the private preferences file into the Colab runtime.
3. Provide private input data and the classifier through the configured paths.
4. Run cells from top to bottom without relying on an old notebook kernel.
5. Confirm that the DADA2 artifacts validate before using archive restoration.
6. Confirm that `vsearch`, `mafft`, and `FastTree` are visible before the corresponding steps.
7. Compare the generated trace, key table dimensions, and model-output schemas with the expected run manifest.

## What is not validated locally

This repository assembly does not execute the biological pipeline because it requires the private reads, metadata, classifier, substantial compute, and external QIIME/Nextflow environments. Structural checks confirm that the notebook and copied Python files are present and parseable; they do not replace a full biological rerun.
