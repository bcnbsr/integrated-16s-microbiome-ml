#!/usr/bin/env python3
"""Dataset 1 post-hoc dietary comparisons and selected-feature boxplots.

The script performs two analyses used for biological interpretation:

1. Mann–Whitney U tests for all nine Dataset 1 dietary variables, with
   Benjamini–Hochberg correction across all nine tests.
2. Group-wise boxplots for the SVM-prioritized taxa and dietary variables.

Microbial features are displayed as log10(TSS relative abundance + 1e-6).
Dietary variables are displayed in their original units.

Example
-------
python src/dataset1_posthoc_diet_and_boxplots.py \
    --metadata-file "prepared_data/participants_v3_corrected_depflag.tsv" \
    --tss-file "prepared_data/taxonomy_leaf_TSS_v1.xlsx" \
    --output-dir "results/dataset1/posthoc"
"""

from __future__ import annotations

import argparse
import re
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


# =============================================================================
# Analysis settings
# =============================================================================

SAMPLE_ID_COLUMN = "participant_id"
OUTCOME_COLUMN = "dep_flag_corrected"
PSEUDOCOUNT = 1e-6
JITTER_RANDOM_STATE = 42

ALL_DIET_VARIABLES = {
    "Energy intake": "Kcal",
    "Protein": "Prot",
    "Total fat": "TF",
    "Saturated fat": "SatFat",
    "Monounsaturated fat": "MUF",
    "Polyunsaturated fat": "Poly_InsFat",
    "Cholesterol": "Chol",
    "Carbohydrate": "Carb",
    "Fibre": "Fiber",
}

SELECTED_DIET_VARIABLES = (
    "Cholesterol",
    "Polyunsaturated fat",
    "Saturated fat",
    "Fibre",
    "Carbohydrate",
    "Protein",
)

SELECTED_TAXA = (
    "Dielma",
    "Lachnospiraceae_ND3007_group",
    "Lachnospiraceae_FCS020_group",
    "Coprobacter",
    "Intestinibacter",
    "Oscillibacter",
    "Rhizobiales_unclassified",
    "Catabacter",
    "Anaerofustis",
    "Eisenbergiella",
)

DIET_Y_AXIS_LABELS = {
    variable: f"{variable}\n(raw value)"
    for variable in SELECTED_DIET_VARIABLES
}


# =============================================================================
# General helpers
# =============================================================================

def canonical_sample_key(value: object) -> str:
    """Convert Dataset 1 sample-ID formats to a common comparison key."""
    text = str(value).strip()

    match = re.search(
        r"MIC[^0-9]*(\d+)",
        text,
        flags=re.IGNORECASE,
    )
    if match:
        return f"MIC{int(match.group(1))}"

    text = re.sub(
        r"^(sub-|c16s)",
        "",
        text,
        flags=re.IGNORECASE,
    )
    return text.strip().upper()


def normalise_taxon_name(value: object) -> str:
    """Normalise taxon labels without changing their biological meaning."""
    text = str(value).strip().lower()
    text = text.replace("\\_", "_")
    text = re.sub(
        r"(?:^|[;| ])(?:k|p|c|o|f|g|s)__",
        "_",
        text,
    )
    text = re.sub(r"[^a-z0-9]+", "_", text)
    text = re.sub(r"_+", "_", text)
    return text.strip("_")


def benjamini_hochberg(p_values: pd.Series) -> np.ndarray:
    """Return Benjamini–Hochberg adjusted p-values."""
    values = np.asarray(p_values, dtype=float)

    if (
        np.any(~np.isfinite(values))
        or np.any(values < 0)
        or np.any(values > 1)
    ):
        raise ValueError(
            "P-values must be finite values between zero and one."
        )

    n_tests = len(values)
    order = np.argsort(values)
    ranked = values[order]

    adjusted_ranked = (
        ranked * n_tests / np.arange(1, n_tests + 1)
    )
    adjusted_ranked = np.minimum.accumulate(
        adjusted_ranked[::-1]
    )[::-1]
    adjusted_ranked = np.minimum(adjusted_ranked, 1.0)

    adjusted = np.empty(n_tests, dtype=float)
    adjusted[order] = adjusted_ranked
    return adjusted


# =============================================================================
# Metadata loading
# =============================================================================

def load_metadata(metadata_file: Path) -> pd.DataFrame:
    """Load participant metadata and create the binary group label."""
    if not metadata_file.exists():
        raise FileNotFoundError(
            f"Participant metadata file not found: {metadata_file}"
        )

    metadata = pd.read_csv(metadata_file, sep="\t")
    metadata.columns = metadata.columns.astype(str).str.strip()

    required_columns = {
        SAMPLE_ID_COLUMN,
        OUTCOME_COLUMN,
        *ALL_DIET_VARIABLES.values(),
    }
    missing_columns = required_columns.difference(metadata.columns)

    if missing_columns:
        raise ValueError(
            "Participant metadata is missing required columns: "
            f"{sorted(missing_columns)}"
        )

    metadata = metadata.copy()
    metadata[SAMPLE_ID_COLUMN] = (
        metadata[SAMPLE_ID_COLUMN].astype(str).str.strip()
    )

    if metadata[SAMPLE_ID_COLUMN].duplicated().any():
        duplicated_ids = metadata.loc[
            metadata[SAMPLE_ID_COLUMN].duplicated(keep=False),
            SAMPLE_ID_COLUMN,
        ].tolist()
        raise ValueError(
            "Duplicated participant IDs were found: "
            f"{duplicated_ids}"
        )

    outcome = pd.to_numeric(
        metadata[OUTCOME_COLUMN],
        errors="coerce",
    )

    if outcome.isna().any():
        raise ValueError(
            "Missing or non-numeric depression labels were found."
        )

    unique_outcomes = set(outcome.astype(int).unique())
    if not unique_outcomes.issubset({0, 1}):
        raise ValueError(
            f"{OUTCOME_COLUMN} must contain only binary values 0 and 1."
        )

    metadata[OUTCOME_COLUMN] = outcome.astype(int)
    metadata["Group"] = np.where(
        metadata[OUTCOME_COLUMN].eq(1),
        "Depressive",
        "Control",
    )

    for column in ALL_DIET_VARIABLES.values():
        metadata[column] = pd.to_numeric(
            metadata[column],
            errors="coerce",
        )

    metadata["_sample_key"] = metadata[
        SAMPLE_ID_COLUMN
    ].map(canonical_sample_key)

    if metadata["_sample_key"].duplicated().any():
        duplicated_keys = metadata.loc[
            metadata["_sample_key"].duplicated(keep=False),
            [SAMPLE_ID_COLUMN, "_sample_key"],
        ]
        raise ValueError(
            "Participant IDs produced duplicated matching keys:\n"
            f"{duplicated_keys.to_string(index=False)}"
        )

    return metadata


# =============================================================================
# Mann–Whitney dietary comparisons
# =============================================================================

def run_dietary_comparisons(
    metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare all nine dietary variables and return full and selected tables."""
    rows: list[dict[str, object]] = []

    for variable_name, column_name in ALL_DIET_VARIABLES.items():
        control = metadata.loc[
            metadata["Group"].eq("Control"),
            column_name,
        ].dropna()

        depressive = metadata.loc[
            metadata["Group"].eq("Depressive"),
            column_name,
        ].dropna()

        if control.empty or depressive.empty:
            raise ValueError(
                "At least one group has no valid values for "
                f"{variable_name}."
            )

        statistic, p_value = mannwhitneyu(
            control,
            depressive,
            alternative="two-sided",
            method="auto",
        )

        control_q1 = control.quantile(0.25)
        control_q3 = control.quantile(0.75)
        depressive_q1 = depressive.quantile(0.25)
        depressive_q3 = depressive.quantile(0.75)

        control_median = control.median()
        depressive_median = depressive.median()

        if depressive_median > control_median:
            higher_group = "Depressive"
        elif control_median > depressive_median:
            higher_group = "Control"
        else:
            higher_group = "Equal median"

        rows.append(
            {
                "Diet variable": variable_name,
                "Column name": column_name,
                "n Control": len(control),
                "n Depressive": len(depressive),
                "Control median": control_median,
                "Control IQR": control_q3 - control_q1,
                "Depressive median": depressive_median,
                "Depressive IQR": depressive_q3 - depressive_q1,
                "Higher median in": higher_group,
                "Mann-Whitney U": statistic,
                "p-value": p_value,
            }
        )

    full_results = pd.DataFrame(rows)
    full_results["BH-FDR q-value"] = benjamini_hochberg(
        full_results["p-value"]
    )
    full_results["FDR significant at 0.05"] = (
        full_results["BH-FDR q-value"] < 0.05
    )

    selected_results = (
        full_results.loc[
            full_results["Diet variable"].isin(
                SELECTED_DIET_VARIABLES
            )
        ]
        .assign(
            _order=lambda frame: frame["Diet variable"].map(
                {
                    variable: index
                    for index, variable in enumerate(
                        SELECTED_DIET_VARIABLES
                    )
                }
            )
        )
        .sort_values("_order")
        .drop(columns="_order")
        .reset_index(drop=True)
    )

    return full_results, selected_results


# =============================================================================
# TSS loading and sample alignment
# =============================================================================

def load_and_align_tss(
    tss_file: Path,
    metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    """Load the TSS table and align its sample columns with metadata rows."""
    if not tss_file.exists():
        raise FileNotFoundError(
            f"TSS taxonomy file not found: {tss_file}"
        )

    tss = pd.read_excel(tss_file)
    tss.columns = tss.columns.astype(str).str.strip()

    if "taxon" not in tss.columns:
        raise ValueError(
            "The TSS table must contain a column named 'taxon'."
        )

    metadata_keys = set(metadata["_sample_key"])
    sample_column_map: dict[str, str] = {}

    for column in tss.columns:
        key = canonical_sample_key(column)

        if key not in metadata_keys:
            continue

        if key in sample_column_map:
            raise ValueError(
                "Multiple TSS columns produced the same sample key "
                f"{key}: {sample_column_map[key]} and {column}"
            )

        sample_column_map[key] = column

    if not sample_column_map:
        raise ValueError(
            "No sample columns in the TSS table matched participant IDs."
        )

    common_keys = [
        key
        for key in metadata["_sample_key"]
        if key in sample_column_map
    ]

    metadata_aligned = (
        metadata
        .set_index("_sample_key")
        .loc[common_keys]
        .copy()
    )

    ordered_sample_columns = [
        sample_column_map[key]
        for key in common_keys
    ]

    tss[ordered_sample_columns] = tss[
        ordered_sample_columns
    ].apply(
        pd.to_numeric,
        errors="coerce",
    )

    if tss[ordered_sample_columns].isna().any().any():
        raise ValueError(
            "The matched TSS sample columns contain missing or "
            "non-numeric values."
        )

    tss["_normalised_taxon"] = tss[
        "taxon"
    ].map(normalise_taxon_name)

    return (
        tss,
        metadata_aligned,
        common_keys,
        ordered_sample_columns,
    )


# =============================================================================
# Selected-feature plotting data
# =============================================================================

def find_taxon_row(
    tss: pd.DataFrame,
    selected_taxon: str,
) -> tuple[pd.Series, str]:
    """Find one taxon row using exact normalised matching and a strict fallback."""
    target = normalise_taxon_name(selected_taxon)

    exact_matches = tss.loc[
        tss["_normalised_taxon"].eq(target)
    ]

    if len(exact_matches) == 1:
        return exact_matches.iloc[0], "exact_normalised_match"

    if len(exact_matches) > 1:
        raise ValueError(
            f"Multiple exact TSS rows matched {selected_taxon}."
        )

    target_tokens = {
        token
        for token in target.split("_")
        if token
    }

    fallback_mask = tss["_normalised_taxon"].map(
        lambda value: target_tokens.issubset(
            set(str(value).split("_"))
        )
    )
    fallback_matches = tss.loc[fallback_mask]

    if len(fallback_matches) == 1:
        return fallback_matches.iloc[0], "unique_token_match"

    if len(fallback_matches) > 1:
        candidates = fallback_matches["taxon"].astype(str).tolist()
        raise ValueError(
            f"Multiple fallback TSS rows matched {selected_taxon}: "
            f"{candidates}"
        )

    raise ValueError(
        f"No TSS row matched selected taxon: {selected_taxon}"
    )


def build_plotting_table(
    tss: pd.DataFrame,
    metadata_aligned: pd.DataFrame,
    common_keys: list[str],
    ordered_sample_columns: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build participant-level long data for the selected taxa and diet variables."""
    plot_frames: list[pd.DataFrame] = []
    match_rows: list[dict[str, object]] = []

    for selected_taxon in SELECTED_TAXA:
        taxon_row, match_method = find_taxon_row(
            tss,
            selected_taxon,
        )

        raw_values = pd.Series(
            taxon_row[ordered_sample_columns].to_numpy(
                dtype=float
            ),
            index=common_keys,
        )

        if (raw_values < 0).any():
            raise ValueError(
                f"Negative TSS values were found for {selected_taxon}."
            )

        plot_values = np.log10(
            raw_values + PSEUDOCOUNT
        )

        plot_frames.append(
            pd.DataFrame(
                {
                    "Participant ID": metadata_aligned[
                        SAMPLE_ID_COLUMN
                    ].to_numpy(),
                    "Feature": selected_taxon,
                    "Feature group": "Taxa",
                    "Group": metadata_aligned[
                        "Group"
                    ].to_numpy(),
                    "Raw value": raw_values.to_numpy(),
                    "Plot value": plot_values.to_numpy(),
                    "Plot scale": (
                        "log10(TSS relative abundance + 1e-6)"
                    ),
                    "Y-axis label": (
                        "log10(TSS relative abundance + 1e-6)"
                    ),
                }
            )
        )

        match_rows.append(
            {
                "Selected taxon": selected_taxon,
                "Matched TSS taxon": taxon_row["taxon"],
                "Match method": match_method,
            }
        )

    for variable_name in SELECTED_DIET_VARIABLES:
        column_name = ALL_DIET_VARIABLES[variable_name]
        raw_values = metadata_aligned[
            column_name
        ].astype(float)

        plot_frames.append(
            pd.DataFrame(
                {
                    "Participant ID": metadata_aligned[
                        SAMPLE_ID_COLUMN
                    ].to_numpy(),
                    "Feature": variable_name,
                    "Feature group": "Diet",
                    "Group": metadata_aligned[
                        "Group"
                    ].to_numpy(),
                    "Raw value": raw_values.to_numpy(),
                    "Plot value": raw_values.to_numpy(),
                    "Plot scale": "Original dietary value",
                    "Y-axis label": DIET_Y_AXIS_LABELS[
                        variable_name
                    ],
                }
            )
        )

    plotting_data = pd.concat(
        plot_frames,
        ignore_index=True,
    )
    match_report = pd.DataFrame(match_rows)

    feature_order = list(SELECTED_TAXA) + list(
        SELECTED_DIET_VARIABLES
    )

    plotting_data["Feature"] = pd.Categorical(
        plotting_data["Feature"],
        categories=feature_order,
        ordered=True,
    )

    plotting_data = (
        plotting_data
        .sort_values(["Feature", "Participant ID"])
        .reset_index(drop=True)
    )

    return plotting_data, match_report


# =============================================================================
# Boxplot
# =============================================================================

def create_combined_boxplot(
    plotting_data: pd.DataFrame,
    output_dir: Path,
) -> None:
    """Create and save the combined selected-feature boxplot figure."""
    available_features = [
        feature
        for feature in (
            list(SELECTED_TAXA)
            + list(SELECTED_DIET_VARIABLES)
        )
        if feature in plotting_data["Feature"].astype(str).unique()
    ]

    n_columns = 4
    n_rows = int(
        np.ceil(len(available_features) / n_columns)
    )

    figure, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(n_columns * 4.5, n_rows * 4.4),
        squeeze=False,
    )
    flat_axes = axes.flatten()

    group_order = ("Control", "Depressive")
    box_colours = ("#F8766D", "#00BFC4")
    random_generator = np.random.default_rng(
        JITTER_RANDOM_STATE
    )

    for index, feature in enumerate(available_features):
        axis = flat_axes[index]
        subset = plotting_data.loc[
            plotting_data["Feature"].astype(str).eq(feature)
        ].copy()

        feature_group = subset[
            "Feature group"
        ].iloc[0]
        y_axis_label = subset[
            "Y-axis label"
        ].iloc[0]

        values_by_group = [
            subset.loc[
                subset["Group"].eq(group),
                "Plot value",
            ]
            .dropna()
            .to_numpy()
            for group in group_order
        ]

        boxplot = axis.boxplot(
            values_by_group,
            tick_labels=group_order,
            patch_artist=True,
            widths=0.6,
            showfliers=False,
        )

        for patch, colour in zip(
            boxplot["boxes"],
            box_colours,
        ):
            patch.set_facecolor(colour)
            patch.set_alpha(0.65)

        for median in boxplot["medians"]:
            median.set_color("black")
            median.set_linewidth(1.3)

        for x_position, values in enumerate(
            values_by_group,
            start=1,
        ):
            jitter = random_generator.normal(
                loc=0,
                scale=0.045,
                size=len(values),
            )
            axis.scatter(
                np.full(len(values), x_position) + jitter,
                values,
                s=14,
                alpha=0.75,
                color="black",
                linewidths=0,
            )

        title = feature.replace("_", " ")
        title = "\n".join(
            textwrap.wrap(title, width=28)
        )

        axis.set_title(
            title,
            fontweight="bold",
            pad=8,
        )
        axis.set_ylabel(y_axis_label)
        axis.set_xlabel("")
        axis.grid(True, axis="y", alpha=0.25)
        axis.grid(False, axis="x")

        scale_note = (
            "Taxa: log10(TSS + 1e-6)"
            if feature_group == "Taxa"
            else "Diet: original value"
        )
        axis.text(
            0.5,
            0.98,
            scale_note,
            transform=axis.transAxes,
            ha="center",
            va="top",
            fontsize=8,
        )

    for index in range(
        len(available_features),
        len(flat_axes),
    ):
        figure.delaxes(flat_axes[index])

    figure.suptitle(
        "Dataset 1 Selected Machine-Learning Features by Group",
        y=1.005,
        fontweight="bold",
        fontsize=14,
    )
    figure.tight_layout()

    figure.savefig(
        output_dir
        / "dataset1_selected_feature_boxplots.png",
        dpi=300,
        bbox_inches="tight",
    )
    figure.savefig(
        output_dir
        / "dataset1_selected_feature_boxplots.pdf",
        bbox_inches="tight",
    )
    plt.close(figure)


# =============================================================================
# Main
# =============================================================================

def run_analysis(
    metadata_file: Path,
    tss_file: Path,
    output_dir: Path,
) -> None:
    """Run the complete Dataset 1 post-hoc dietary and plotting workflow."""
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    metadata = load_metadata(metadata_file)

    print("Group counts:")
    print(
        metadata["Group"]
        .value_counts()
        .to_string()
    )

    full_diet_results, selected_diet_results = (
        run_dietary_comparisons(metadata)
    )

    full_diet_results.to_csv(
        output_dir
        / "dataset1_all_diet_mannwhitney_results.csv",
        index=False,
    )
    selected_diet_results.to_csv(
        output_dir
        / "dataset1_selected_diet_mannwhitney_results.csv",
        index=False,
    )

    (
        tss,
        metadata_aligned,
        common_keys,
        ordered_sample_columns,
    ) = load_and_align_tss(
        tss_file,
        metadata,
    )

    plotting_data, taxon_match_report = (
        build_plotting_table(
            tss,
            metadata_aligned,
            common_keys,
            ordered_sample_columns,
        )
    )

    plotting_data.to_csv(
        output_dir
        / "dataset1_selected_feature_boxplot_values.csv",
        index=False,
    )
    taxon_match_report.to_csv(
        output_dir
        / "dataset1_selected_taxa_match_report.csv",
        index=False,
    )

    create_combined_boxplot(
        plotting_data,
        output_dir,
    )

    print("\nSelected dietary comparison results:")
    print(
        selected_diet_results[
            [
                "Diet variable",
                "Control median",
                "Control IQR",
                "Depressive median",
                "Depressive IQR",
                "Higher median in",
                "p-value",
                "BH-FDR q-value",
            ]
        ]
        .round(3)
        .to_string(index=False)
    )

    print(
        "\nAnalysis completed. Outputs saved to:"
    )
    print(output_dir)


def parse_arguments() -> argparse.Namespace:
    """Parse input and output paths."""
    parser = argparse.ArgumentParser(
        description=(
            "Run Dataset 1 post-hoc dietary comparisons and "
            "selected-feature boxplots."
        )
    )
    parser.add_argument(
        "--metadata-file",
        type=Path,
        required=True,
        help=(
            "Path to participants_v3_corrected_depflag.tsv"
        ),
    )
    parser.add_argument(
        "--tss-file",
        type=Path,
        required=True,
        help="Path to taxonomy_leaf_TSS_v1.xlsx",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for generated tables and figures",
    )
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_arguments()
    run_analysis(
        metadata_file=arguments.metadata_file,
        tss_file=arguments.tss_file,
        output_dir=arguments.output_dir,
    )

