# Integrated 16S Microbiome and Machine-Learning Pipeline

An end-to-end, reproducible workflow for processing paired-end 16S rRNA amplicon data and evaluating microbiome-based classification models. The project combines a QIIME 2/Nextflow biology workflow with downstream taxonomic preprocessing, statistical analysis, feature selection, machine learning, visualisation, and result archiving.

The project is designed for two audiences:

- **Researchers and students** who want to run the complete analysis in Google Colab with minimal setup.
- **Developers and reviewers** who want to inspect the workflow, understand each decision, reproduce a result, or adapt one stage without reading a long notebook cell by cell.

## What the workflow does

At a high level, the pipeline:

1. Installs and checks the required computational environments.
2. Mounts Google Drive when running in Colab and creates isolated project, result, cache, and log directories.
3. Obtains or validates paired-end FASTQ reads and the associated metadata.
4. Imports the reads into QIIME 2 and visualises demultiplexing and quality information.
5. Removes primers with Cutadapt.
6. Denoises the reads with DADA2 and produces an amplicon sequence variant (ASV) feature table, representative sequences, and denoising statistics.
7. Optionally clusters the ASVs into operational taxonomic units (OTUs) at the configured similarity threshold.
8. Assigns taxonomy with a SILVA classifier, filters unwanted features, and exports analysis-ready tables.
9. Calculates alpha- and beta-diversity metrics and creates diversity visualisations.
10. Prevents information leakage during machine-learning preprocessing by fitting transformations only on the training data within each split.
11. Applies the configured statistical/feature-selection workflow and evaluates classifiers over repeated random splits.
12. Writes plots, Excel workbooks, checksums, provenance files, and an archive of the generated results.

The current tested configuration uses the `asv` branch. Set `run.branch` to `otu` in the preferences file when the OTU comparison is required.

## Repository contents

```text
.
├── README.md
├── CITATION.cff
├── CONTRIBUTING.md
├── LICENSING.md
├── NOTICE.md
├── requirements.txt
├── notebooks/
│   └── 16s_colab_pipeline_tested.ipynb
├── config/
│   └── preferences_colombian.example.yaml
├── workflows/
│   └── CSthesis/
│       ├── main.nf
│       └── nextflow.config
├── third_party/
│   └── diet-microbiome-depression-ml/
│       ├── dataset1_ancombc2_analysis.R
│       ├── dataset1_diet_correlation.py
│       ├── dataset1_filtered_diet_ml_pipeline.py
│       ├── dataset1_posthoc_diet_and_boxplots.py
│       ├── dataset1_svm__plot.R
│       ├── dataset1_taxa_preprocessing.py
│       ├── requirements.upstream.txt
│       └── R_packages.upstream.txt
├── docs/
│   ├── pipeline_overview.md
│   └── reproducibility.md
└── examples/
    └── metadata.template.tsv
```

The files under `workflows/CSthesis` and `third_party/diet-microbiome-depression-ml` are retained so that the code used by the two source repositories is present and inspectable. Their provenance and current licensing status are documented in [NOTICE.md](NOTICE.md) and [LICENSING.md](LICENSING.md).

## Quick start: Google Colab

The supported entry point is `notebooks/16s_colab_pipeline_tested.ipynb`.

### 1. Prepare the input files

Before starting a run, prepare:

- paired-end FASTQ files, following the configured `_1`/`_2` naming convention;
- a QIIME 2 metadata file containing a `sample-id` column;
- phenotype metadata containing the configured participant identifier and group columns;
- the participant mapping workbook if the selected data source needs it;
- a SILVA 138 classifier, or permission for the notebook to download it;
- optionally, a previously created ASV-results archive containing the DADA2 artifacts.

Do not commit participant-level metadata, sequencing reads, QIIME artifacts, archives, or generated results to a public repository. Keep them in private storage and provide their paths through the configuration file.

### 2. Configure the run

Copy `config/preferences_colombian.example.yaml` to `/content/preferences_colombian.yaml` in Colab, then edit only the values appropriate for your data. The notebook reads this file near the beginning of the run.

The most important settings are:

- `run.branch`: `asv` or `otu`;
- `run.qiime2_release`: the QIIME 2 release to install;
- `input.qiime_metadata_source`: the relative location of QIIME metadata;
- `input.phenotype_metadata`: the phenotype metadata filename;
- `reads.patterns`: the forward and reverse FASTQ filename patterns;
- `primers`: primer sequences and whether primer removal is enabled;
- `quality`: truncation lengths and the sampling depth;
- `taxonomy.classifier`: the classifier path;
- `ml`: feature-selection and repeated-split settings.

The example configuration uses `/content` paths because that is the layout expected by the Colab notebook. These are runtime paths, not repository paths.

### 3. Run the notebook in order

Run the cells from top to bottom. The notebook contains explicit checks for missing tools and invalid QIIME artifacts. If an earlier run produced valid DADA2 artifacts and an archive was configured, the archive restoration step can reuse them; otherwise, the Nextflow biology workflow is run.

The notebook may ask for access to Google Drive. This is required only for the configured private inputs, caches, archives, and result backup locations.

### 4. Choose the biology branch

The `asv` branch uses the DADA2 feature table directly. The `otu` branch clusters the ASVs with VSEARCH at the configured similarity, currently `0.97`. The user-facing branch name remains `otu`; compatibility translation is handled inside the notebook where necessary.

### 5. Inspect the outputs

The main result directory is configured by `run.results_dir`. Typical outputs include:

- QIIME 2 artifacts for imported reads, denoising, taxonomy, diversity, and phylogeny;
- exported feature tables and taxonomy tables;
- filtered taxonomic matrices used by the ML stage;
- alpha- and beta-diversity plots;
- model metrics and repeated-split results;
- feature-selection and coefficient plots;
- Excel workbooks for tables and summaries;
- Nextflow `trace.txt`, `timeline.html`, and `report.html` when the workflow is run;
- checksums and a compressed results archive when reporting is enabled.

## Running the Nextflow workflow separately

The biology workflow can also be run outside the notebook with Nextflow. It expects the input layout and metadata described in the upstream workflow README and uses the container configuration in `workflows/CSthesis/nextflow.config`.

Install Nextflow and either Apptainer/Singularity or Docker, then adapt the paths for the local machine. A representative command is:

```bash
nextflow run workflows/CSthesis/main.nf \
  --metadata_url /path/to/metadata.tsv \
  --classifier_url /path/to/silva-138-99-nb-classifier.qza \
  --outdir /path/to/results \
  -work-dir /path/to/nextflow-work \
  -resume
```

The exact parameters available are defined in `main.nf`. For a full Colab execution, use the notebook because it also performs environment setup, input validation, archive restoration, downstream ML, and result packaging.

## Computational resources

The notebook's Colab Nextflow configuration requests a local executor with:

| Setting | Value |
|---|---:|
| Executor | `local` |
| CPUs requested | 8 |
| Memory requested | 36 GB |
| Maximum parallel processes | 1 |
| GPU | 0 |
| Work directory | `/content/nxf_work` |

The standalone upstream `nextflow.config` contains a separate default allocation for its original execution environment. Do not confuse a requested allocation with measured use: actual CPU, memory, duration, and RSS values come from the Nextflow trace generated by the run.

## Reproducibility and interpretation

Reproducibility depends on more than the Python package versions. Record all of the following with every analysis:

- repository commit and notebook version;
- QIIME 2 release and external executable versions;
- SILVA classifier release and file checksum;
- primer and quality-trimming settings;
- sampling depth and feature-filtering thresholds;
- branch (`asv` or `otu`) and OTU similarity, if applicable;
- metadata version and the mapping between sample and participant identifiers;
- random states used for repeated model splits;
- Nextflow trace and report files;
- the checksum manifest for important inputs and outputs.

The ML stage is exploratory research code, not a clinical diagnostic system. Performance metrics should be interpreted with regard to sample size, class balance, repeated-split variability, preprocessing choices, and the possibility of cohort-specific confounding. A high cross-validation score does not establish clinical utility or causation.

## Data protection and repository hygiene

This repository intentionally does not include:

- raw FASTQ files;
- participant-level metadata or mapping workbooks;
- QIIME `.qza`/`.qzv` artifacts;
- Google Drive archives;
- trained models, generated plots, or result workbooks;
- private credentials, access tokens, or classifier binaries.

The `.gitignore` file excludes common versions of these files. Check `git status` and inspect the complete diff before every push. If a sensitive file has ever been committed, removing it from the working tree is not sufficient; rotate exposed credentials and rewrite the repository history as appropriate.

## Licensing summary

There is currently no blanket open-source license in this assembled repository. This is deliberate: the two upstream repositories currently expose no `LICENSE` file, and a public GitHub repository without a license does not grant permission to copy, modify, or redistribute their code.

My recommendation is:

1. Obtain written permission from the authors of both upstream repositories, or confirm that their code is already covered by another agreement.
2. Keep the upstream files clearly attributed and unchanged where possible.
3. License only the integration code and original documentation that you own, preferably under **Apache License 2.0** if you want explicit patent protection, or **MIT** if maximum simplicity is more important.
4. Add the chosen license only after the ownership and permission question is resolved.

See [LICENSING.md](LICENSING.md) for the detailed decision and a release checklist. This documentation is practical project guidance, not legal advice.

## Citation

If you use this repository in academic work, cite the repository and the upstream projects listed in [NOTICE.md](NOTICE.md). The machine-readable citation information is in [CITATION.cff](CITATION.cff).

## Troubleshooting

### `vsearch`, `mafft`, or `FastTree` cannot be found

Run the environment setup cells again and confirm that the QIIME 2 environment's `bin` directory is on `PATH`. The notebook performs a preflight check before OTU clustering and phylogeny.

### DADA2 artifacts are missing or invalid

Check that the archive contains valid `dada2/table.qza`, `dada2/rep-seqs.qza`, and `dada2/denoising-stats.qza` files. If any artifact fails `qiime tools peek` or QIIME validation, rerun the DADA2 stage rather than using a partial archive.

### The ML stage cannot find a module

Confirm that the ML repository snapshot is present and that the Python dependencies in `requirements.txt` are installed. The notebook uses dynamic imports because the upstream filenames are analysis scripts rather than an installed package.

### The run stops during metadata validation

Compare the metadata column names and sample identifiers with `config/preferences_colombian.example.yaml`. The QIIME `sample-id` values and the phenotype participant mapping must be consistent before modelling.

## Project status

The attached notebook is the tested integration entry point at the time this repository was assembled. The source snapshots are included for transparency and reproducibility. Before publishing the repository publicly, resolve the upstream licensing question, review the notebook's embedded outputs for sensitive information, and perform a clean test run from a fresh environment.
