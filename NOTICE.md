# Source and attribution notice

This repository combines an integration notebook with code originating from two public repositories. The source repositories did not contain license files when audited for this assembly. Their authors and contributors retain their respective rights.

## QIIME 2 / Nextflow workflow

Source repository: [Luna-Maan/CSthesis](https://github.com/Luna-Maan/CSthesis)

Included or adapted files:

- `workflows/CSthesis/main.nf`
- `workflows/CSthesis/nextflow.config`
- `third_party/CSthesis/README.upstream.md`
- `third_party/CSthesis/.gitignore`

The workflow defines the QIIME 2 stages, including read import, Cutadapt, DADA2, taxonomy assignment, phylogeny, diversity analysis, and table normalisation.

## Diet/microbiome machine-learning and statistical analysis

Source repository: [R-joumaa/diet-microbiome-depression-ml](https://github.com/R-joumaa/diet-microbiome-depression-ml)

Included files are retained under `third_party/diet-microbiome-depression-ml/`, including the Python analysis scripts, R analysis scripts, and upstream dependency lists.

## Integration notebook

`notebooks/16s_colab_pipeline_tested.ipynb` is the supplied integration notebook. It orchestrates environment setup, input validation, archive restoration, Nextflow execution, ASV/OTU processing, downstream ML, provenance capture, and result archiving. Its redistribution and licensing status should be confirmed with the person or institution that owns the notebook and data-analysis work.

## Data and dependencies

The repository does not redistribute the biological data, participant metadata, QIIME classifier binaries, or external software environments. Each dependency remains subject to its own license. QIIME 2, Nextflow, SILVA, Python packages, R packages, and container images should be checked against their respective license and citation requirements before publication.
