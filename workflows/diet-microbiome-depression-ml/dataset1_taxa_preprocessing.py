#!/usr/bin/env python3
"""Dataset 1 taxonomy filtering, hierarchical roll-up, and CLR transformation."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

HEADER_ROWS_TO_REMOVE = 14
MIN_READS = 5_000
MIN_PREVALENCE_FRACTION = 0.10
MIN_TOTAL_COUNT = 10
PSEUDOCOUNT = 0.5
SAMPLE_PREFIX = "c16S"

DEFAULT_INPUT_FILENAME = "taxonomy completo paper.xlsx"
DEFAULT_PREPROCESSED_DIR = "Preprocessed Files"

def running_in_colab() -> bool:
    return importlib.util.find_spec("google.colab") is not None


def resolve_project_dir(project_dir: Path | None) -> Path:
    if project_dir is not None:
        return project_dir.expanduser().resolve()

    if running_in_colab():
        from google.colab import drive
        drive.mount("/content/drive")
        return Path("/content/drive/MyDrive/project")

    return Path.cwd().resolve()


def validate_required_columns(df: pd.DataFrame) -> None:
    required = {"taxlevel", "rankID", "taxon", "daughterlevels", "total"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")


def identify_sample_columns(df: pd.DataFrame) -> list[str]:
    sample_cols = [col for col in df.columns if str(col).startswith(SAMPLE_PREFIX)]
    if not sample_cols:
        raise ValueError(f"No sample columns beginning with '{SAMPLE_PREFIX}' were found.")
    return sample_cols


def load_raw_table(input_file: Path):
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    df_raw = pd.read_excel(input_file, dtype={"rankID": str})
    validate_required_columns(df_raw)
    sample_cols = identify_sample_columns(df_raw)

    bacteria_row = df_raw.loc[df_raw["taxon"].eq("Bacteria")].copy()
    if len(bacteria_row) != 1:
        raise ValueError("Exactly one row with taxon == 'Bacteria' is required.")

    df_taxa = df_raw.iloc[HEADER_ROWS_TO_REMOVE:].reset_index(drop=True).copy()
    df_taxa[sample_cols] = df_taxa[sample_cols].apply(pd.to_numeric, errors="raise")
    bacteria_row[sample_cols] = bacteria_row[sample_cols].apply(
        pd.to_numeric, errors="raise"
    )

    print(f"Raw table: {df_raw.shape[0]} rows x {df_raw.shape[1]} columns")
    print(
        f"After removing {HEADER_ROWS_TO_REMOVE} rows: "
        f"{df_taxa.shape[0]} rows x {df_taxa.shape[1]} columns"
    )
    print(f"Samples detected: {len(sample_cols)}")
    return df_taxa, bacteria_row, sample_cols


def filter_low_depth_samples(df_taxa, bacteria_row, sample_cols):
    sample_depths = bacteria_row.loc[:, sample_cols].iloc[0]
    low_depth = sample_depths.index[sample_depths < MIN_READS].tolist()

    print(f"\nMinimum read threshold: {MIN_READS}")
    print(f"Low-depth samples: {len(low_depth)}")

    if low_depth:
        print("Removing:", ", ".join(low_depth))
        df_taxa = df_taxa.drop(columns=low_depth).copy()
        sample_cols = [sample for sample in sample_cols if sample not in low_depth]
    else:
        print("No samples were removed.")

    df_taxa["total"] = df_taxa.loc[:, sample_cols].sum(axis=1)
    print(f"Samples retained: {len(sample_cols)}")
    return df_taxa, sample_cols


def filter_taxa(df_taxa, sample_cols, output_dir):
    min_prevalence_samples = int(
        np.ceil(MIN_PREVALENCE_FRACTION * len(sample_cols))
    )

    prevalence = (df_taxa.loc[:, sample_cols] > 0).sum(axis=1)
    keep_prevalence = prevalence >= min_prevalence_samples

    after_prevalence = df_taxa.loc[keep_prevalence].reset_index(drop=True)
    removed_prevalence = df_taxa.loc[~keep_prevalence].reset_index(drop=True)
    removed_prevalence.to_excel(
        output_dir / "filtered_out_prevalence.xlsx", index=False
    )

    print(
        f"\nPrevalence threshold: at least {min_prevalence_samples} "
        f"of {len(sample_cols)} samples"
    )
    print(f"Taxa retained after prevalence filtering: {len(after_prevalence)}")
    print(f"Taxa removed by prevalence filtering: {len(removed_prevalence)}")

    keep_count = after_prevalence["total"] >= MIN_TOTAL_COUNT
    filtered = after_prevalence.loc[keep_count].reset_index(drop=True)
    removed_low_count = after_prevalence.loc[~keep_count].reset_index(drop=True)
    removed_low_count.to_excel(
        output_dir / "filtered_out_lowcount.xlsx", index=False
    )

    print(f"\nMinimum total-count threshold: {MIN_TOTAL_COUNT}")
    print(f"Taxa retained after count filtering: {len(filtered)}")
    print(f"Taxa removed by count filtering: {len(removed_low_count)}")
    return filtered


def recalculate_internal_nodes(df_filtered, sample_cols):
    df_rollup = df_filtered.copy()
    df_rollup["daughterlevels"] = (
        df_rollup["daughterlevels"].fillna(0).astype(int)
    )

    leaf_mask = df_rollup["daughterlevels"].eq(0)
    leaf_df = df_rollup.loc[leaf_mask].copy()
    internal_indices = df_rollup.index[~leaf_mask]
    leaf_rank_ids = leaf_df["rankID"].astype(str).str.strip()

    print(f"\nLeaf nodes: {leaf_mask.sum()}")
    print(f"Internal nodes: {(~leaf_mask).sum()}")

    for index in internal_indices:
        parent_rank_id = str(df_rollup.at[index, "rankID"]).strip()
        descendant_mask = leaf_rank_ids.str.startswith(
            f"{parent_rank_id}.", na=False
        )
        descendants = leaf_df.loc[descendant_mask, sample_cols]

        if descendants.empty:
            aggregated = pd.Series(0.0, index=sample_cols)
        else:
            aggregated = descendants.sum(axis=0)

        df_rollup.loc[index, sample_cols] = aggregated.to_numpy()
        df_rollup.at[index, "total"] = float(aggregated.sum())

    bacteria = df_rollup.loc[df_rollup["taxon"].eq("Bacteria")]
    if not bacteria.empty:
        bacteria_total = float(bacteria.loc[:, sample_cols].iloc[0].sum())
        print(f"Bacteria-row total after roll-up: {bacteria_total:,.0f}")
    else:
        print("Bacteria-row sanity check skipped because the row was not retained.")

    return df_rollup


def apply_clr_and_select_leaves(df_rollup, sample_cols):
    counts = df_rollup.loc[:, sample_cols].to_numpy(dtype=float)
    if np.any(counts < 0):
        raise ValueError("Negative counts were found.")

    log_counts = np.log(counts + PSEUDOCOUNT)
    clr_counts = log_counts - log_counts.mean(axis=0, keepdims=True)

    clr_values = pd.DataFrame(clr_counts, columns=sample_cols)
    metadata_cols = [col for col in df_rollup.columns if col not in sample_cols]

    clr_df = pd.concat(
        [
            df_rollup.loc[:, metadata_cols].reset_index(drop=True),
            clr_values.reset_index(drop=True),
        ],
        axis=1,
    )
    clr_df["normalization"] = f"CLR (pseudocount={PSEUDOCOUNT}, no TSS)"

    if not np.allclose(
        clr_df.loc[:, sample_cols].mean(axis=0).to_numpy(),
        0.0,
        atol=1e-10,
    ):
        raise RuntimeError("CLR validation failed: sample means are not zero.")

    return clr_df.loc[clr_df["daughterlevels"].eq(0)].reset_index(drop=True)


def run_pipeline(project_dir: Path) -> None:
    input_dir = project_dir / "input files"
    output_dir = project_dir / DEFAULT_PREPROCESSED_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    input_file = input_dir / DEFAULT_INPUT_FILENAME
    rollup_output = output_dir / "taxonomy_completo_paper_v5.xlsx"
    clr_output = output_dir / "taxonomy_leaf_CLR_v2.xlsx"

    print(f"Project directory: {project_dir}")
    print(f"Input file: {input_file}")
    print(f"Output directory: {output_dir}")

    df_taxa, bacteria_row, sample_cols = load_raw_table(input_file)
    df_taxa, sample_cols = filter_low_depth_samples(
        df_taxa, bacteria_row, sample_cols
    )
    df_filtered = filter_taxa(df_taxa, sample_cols, output_dir)
    df_rollup = recalculate_internal_nodes(df_filtered, sample_cols)
    df_rollup.to_excel(rollup_output, index=False)

    clr_leaf = apply_clr_and_select_leaves(df_rollup, sample_cols)
    clr_leaf.to_excel(clr_output, index=False)

    values = clr_leaf.loc[:, sample_cols].to_numpy(dtype=float)
    print("\nPipeline completed successfully.")
    print(f"Filtered hierarchy: {rollup_output}")
    print(f"Leaf CLR table: {clr_output}")
    print(f"Final leaf taxa: {len(clr_leaf)}")
    print(f"Final CLR range: {values.min():.4f} to {values.max():.4f}")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Preprocess Dataset 1 taxonomy counts and apply CLR."
    )
    parser.add_argument(
        "--project-dir",
        type=Path,
        default=None,
        help=(
            "Project directory containing 'input files'. Defaults to "
            "/content/drive/MyDrive/project in Colab and the current "
            "directory elsewhere."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    run_pipeline(resolve_project_dir(args.project_dir))

