# Pipeline overview

## Stages

The workflow has two connected layers.

### Biology layer

The Nextflow/QIIME 2 layer imports paired-end reads, removes primers, denoises with DADA2, assigns taxonomy, filters features, builds a phylogeny, calculates diversity, and exports tables. The ASV result is the primary output. The optional OTU branch clusters the ASVs at the configured identity threshold for comparison.

### Statistical and ML layer

The downstream layer joins the taxonomic table with phenotype metadata, applies the configured filtering and transformation rules, performs feature selection, and evaluates models using repeated random states. The notebook keeps data-dependent operations inside each split where required so that test-set information does not influence training preprocessing.

## Important file relationships

- `notebooks/16s_colab_pipeline_tested.ipynb` orchestrates the complete run.
- `workflows/CSthesis/main.nf` contains the biology processes.
- `workflows/CSthesis/nextflow.config` defines the upstream container and reporting defaults.
- `third_party/diet-microbiome-depression-ml/` contains the downstream analysis scripts used by the notebook.
- `config/preferences_colombian.example.yaml` is the safe template for runtime parameters.

## ASV and OTU branches

ASVs are inferred sequence variants produced by DADA2. OTUs are clusters of sequences grouped at a similarity threshold. The two branches share the upstream processing and downstream analysis structure, but they do not represent identical features; branch choice must therefore be recorded with the results.

## Leakage control

Feature filtering, transformations, and selection must be fitted using the training portion of a split and then applied to the held-out portion. This prevents information from the test samples from influencing the learned feature set. The notebook also records random states and output files to support an audit of the modelling process.
