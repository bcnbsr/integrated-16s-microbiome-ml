# Running the Pipeline

## Prerequisites

Install Nextflow and Apptainer (a drop-in replacement for Singularity):

```bash
# Nextflow
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
nextflow -version

# Apptainer
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
apptainer --version
```

## Setup

In the same directory as `main.nf`, place the following:

- A folder named `fastq_files/` containing your `.fastq` or `.fastq.gz` files
- The corresponding metadata file from the `metadata_files/` directory, renamed to `metadata.tsv`

## Run

```bash
nextflow run main.nf --reads_dir fastq_files/
```

## Using a Different Dataset

Replace the contents of `fastq_files/` with your new dataset (Either `.fastq` or `.fastq.gz` files) and provide a matching `metadata.tsv`.

The metadata file must have a sample identifier column (`sample-id` or any column name QIIME2 recognises as a sample identifier), where each value matches the sample names of the `.fastq` files in `fastq_files/`. You can include additional columns with sample information, these will be used in downstream QIIME2 analyses such as taxa bar plots and alpha/beta diversity comparisons.
