#!/usr/bin/env python3
"""Dataset 1 feature selection and classification using the filtered diet set.

The script expects the outputs of the preceding preprocessing steps:
1. participant metadata;
2. the CLR-transformed leaf-taxa table.

Within each of five repeated stratified splits, 40% of participants are used
for taxa feature selection, 40% for model training, and 20% for testing.
Only the filtered dietary set is analysed.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from imblearn.ensemble import BalancedRandomForestClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import f_classif, mutual_info_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC


# ---------------------------------------------------------------------------
# Final analysis settings
# ---------------------------------------------------------------------------

RANDOM_STATES = (42, 43, 44, 45, 46)
CLASSIFIER_RANDOM_STATE = 42

ANOVA_P_THRESHOLD = 0.05
MI_PERCENTILE = 95
RF_TREES = 500
TOP_N_FEATURES = 20

SAMPLE_PREFIX = "c16S"
PARTICIPANT_PREFIX = "sub-"

FILTERED_DIET_COLUMNS = (
    "Prot",
    "SatFat",
    "Poly_InsFat",
    "Chol",
    "Carb",
    "Fiber",
)

DIET_NAME_MAP = {
    "Prot": "Protein",
    "SatFat": "Saturated fat",
    "Poly_InsFat": "Polyunsaturated fat",
    "Chol": "Cholesterol",
    "Carb": "Carbohydrate",
    "Fiber": "Fibre",
}


# ---------------------------------------------------------------------------
# Data loading and alignment
# ---------------------------------------------------------------------------

def load_and_align_data(
    metadata_file: Path,
    clr_file: Path,
) -> tuple[
    pd.DataFrame,
    np.ndarray,
    np.ndarray,
    list[str],
    dict[str, str],
]:
    """Load metadata and CLR taxa, then align all rows by participant ID."""
    if not metadata_file.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_file}")
    if not clr_file.exists():
        raise FileNotFoundError(f"CLR taxa file not found: {clr_file}")

    metadata = pd.read_csv(metadata_file, sep="\t")
    clr_raw = pd.read_excel(clr_file, dtype={"rankID": str})

    required_metadata = {
        "participant_id",
        "dep_flag_corrected",
        *FILTERED_DIET_COLUMNS,
    }
    missing_metadata = required_metadata.difference(metadata.columns)
    if missing_metadata:
        raise ValueError(
            f"Metadata is missing required columns: {sorted(missing_metadata)}"
        )

    required_taxa = {"rankID", "taxon"}
    missing_taxa = required_taxa.difference(clr_raw.columns)
    if missing_taxa:
        raise ValueError(
            f"CLR table is missing required columns: {sorted(missing_taxa)}"
        )

    sample_columns = [
        column
        for column in clr_raw.columns
        if str(column).startswith(SAMPLE_PREFIX)
    ]
    if not sample_columns:
        raise ValueError(
            f"No CLR sample columns beginning with '{SAMPLE_PREFIX}' were found."
        )

    taxa_table = clr_raw.set_index("rankID")
    taxa_by_participant = taxa_table.loc[:, sample_columns].T
    taxa_by_participant.index = (
        taxa_by_participant.index.astype(str).str.replace(
            SAMPLE_PREFIX,
            PARTICIPANT_PREFIX,
            regex=False,
        )
    )

    metadata = metadata.copy()
    metadata["participant_id"] = metadata["participant_id"].astype(str)
    metadata = metadata.set_index("participant_id", drop=False)

    common_ids = taxa_by_participant.index[
        taxa_by_participant.index.isin(metadata.index)
    ]
    if common_ids.empty:
        raise ValueError("No participant IDs matched between metadata and CLR data.")

    # Explicit .loc alignment prevents labels and diet rows from becoming
    # misaligned with the taxa matrix.
    taxa_aligned = taxa_by_participant.loc[common_ids]
    metadata_aligned = metadata.loc[common_ids].reset_index(drop=True)

    if metadata_aligned["dep_flag_corrected"].isna().any():
        raise ValueError("Missing depression labels were found.")
    if metadata_aligned.loc[:, FILTERED_DIET_COLUMNS].isna().any().any():
        raise ValueError("Missing values were found in the filtered diet variables.")

    y = metadata_aligned["dep_flag_corrected"].astype(int).to_numpy()
    if set(np.unique(y)) != {0, 1}:
        raise ValueError("The depression outcome must contain binary values 0 and 1.")

    x_taxa = taxa_aligned.to_numpy(dtype=float)
    taxa_ids = taxa_aligned.columns.astype(str).tolist()

    taxon_lookup = (
        clr_raw.loc[:, ["rankID", "taxon"]]
        .drop_duplicates(subset="rankID")
        .set_index("rankID")["taxon"]
        .astype(str)
        .to_dict()
    )

    return metadata_aligned, x_taxa, y, taxa_ids, taxon_lookup


# ---------------------------------------------------------------------------
# Feature selection
# ---------------------------------------------------------------------------

def select_taxa(
    x_fs: np.ndarray,
    y_fs: np.ndarray,
    random_state: int,
) -> tuple[list[int], dict[str, set[int]], float, float]:
    """Apply ANOVA, mutual information, and RF importance on the FS subset."""
    n_taxa = x_fs.shape[1]
    if n_taxa == 0:
        raise ValueError("The taxa matrix contains no features.")

    _, p_values = f_classif(x_fs, y_fs)
    p_values = np.nan_to_num(p_values, nan=1.0)
    anova_selected = set(np.flatnonzero(p_values < ANOVA_P_THRESHOLD))

    mi_values = mutual_info_classif(
        x_fs,
        y_fs,
        random_state=random_state,
    )
    mi_threshold = float(np.percentile(mi_values, MI_PERCENTILE))
    mi_selected = set(np.flatnonzero(mi_values >= mi_threshold))

    rf_selector = RandomForestClassifier(
        n_estimators=RF_TREES,
        class_weight="balanced",
        random_state=random_state,
        n_jobs=-1,
    )
    rf_selector.fit(x_fs, y_fs)

    rf_threshold = 1.0 / n_taxa
    rf_selected = set(
        np.flatnonzero(rf_selector.feature_importances_ > rf_threshold)
    )

    method_sets = {
        "ANOVA": anova_selected,
        "MI": mi_selected,
        "RF": rf_selected,
    }

    consensus = [
        index
        for index in range(n_taxa)
        if sum(index in selected for selected in method_sets.values()) >= 2
    ]

    if not consensus:
        raise RuntimeError(
            f"No taxa passed the two-of-three consensus rule for "
            f"random_state={random_state}."
        )

    return consensus, method_sets, mi_threshold, rf_threshold


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def build_classifiers() -> dict[str, object]:
    """Return fresh classifier instances with the final thesis settings."""
    return {
        "LR": LogisticRegression(
            C=1.0,
            max_iter=1_000,
            class_weight="balanced",
            random_state=CLASSIFIER_RANDOM_STATE,
        ),
        "SVM": LinearSVC(
            C=1.0,
            max_iter=5_000,
            class_weight="balanced",
            random_state=CLASSIFIER_RANDOM_STATE,
        ),
        "RF": RandomForestClassifier(
            n_estimators=RF_TREES,
            class_weight="balanced",
            random_state=CLASSIFIER_RANDOM_STATE,
            n_jobs=-1,
        ),
        "BRF": BalancedRandomForestClassifier(
            n_estimators=RF_TREES,
            random_state=CLASSIFIER_RANDOM_STATE,
            n_jobs=-1,
        ),
    }


def continuous_scores(model: object, x_test: np.ndarray) -> np.ndarray:
    """Return continuous model scores for ROC-AUC calculation."""
    if hasattr(model, "decision_function"):
        return np.asarray(model.decision_function(x_test), dtype=float)

    if hasattr(model, "predict_proba"):
        probabilities = np.asarray(model.predict_proba(x_test), dtype=float)
        return probabilities[:, 1]

    return np.asarray(model.predict(x_test), dtype=float)


def calculate_metrics(
    y_true: np.ndarray,
    predictions: np.ndarray,
    scores: np.ndarray,
) -> dict[str, float]:
    """Calculate the classification metrics reported in the thesis."""
    return {
        "ROC-AUC": roc_auc_score(y_true, scores),
        "ACC": accuracy_score(y_true, predictions),
        "BAL-ACC": balanced_accuracy_score(y_true, predictions),
        "SEN": recall_score(
            y_true,
            predictions,
            pos_label=1,
            zero_division=0,
        ),
        "SPE": recall_score(
            y_true,
            predictions,
            pos_label=0,
            zero_division=0,
        ),
        "PRE": precision_score(
            y_true,
            predictions,
            pos_label=1,
            zero_division=0,
        ),
        "F1": f1_score(
            y_true,
            predictions,
            pos_label=1,
            zero_division=0,
        ),
    }


# ---------------------------------------------------------------------------
# Main repeated-split pipeline
# ---------------------------------------------------------------------------

def run_analysis(
    metadata_file: Path,
    clr_file: Path,
    output_dir: Path,
) -> None:
    """Run filtered-diet, taxa-only, and combined classification analyses."""
    output_dir.mkdir(parents=True, exist_ok=True)

    (
        metadata,
        x_all_taxa,
        y,
        taxa_ids,
        taxon_lookup,
    ) = load_and_align_data(metadata_file, clr_file)

    x_diet = metadata.loc[:, FILTERED_DIET_COLUMNS].to_numpy(dtype=float)
    n_taxa = x_all_taxa.shape[1]

    print(
        f"Participants: {len(y)} | "
        f"Depressive: {(y == 1).sum()} | "
        f"Control: {(y == 0).sum()}"
    )
    print(f"Candidate taxa: {n_taxa}")
    print(f"Filtered diet variables: {list(FILTERED_DIET_COLUMNS)}")

    result_rows: list[dict[str, object]] = []
    selection_rows: list[dict[str, object]] = []
    coefficient_rows: list[dict[str, object]] = []
    svm_prediction_rows: list[dict[str, object]] = []

    for split_number, random_state in enumerate(RANDOM_STATES, start=1):
        first_split = StratifiedShuffleSplit(
            n_splits=1,
            test_size=0.40,
            random_state=random_state,
        )
        train_test_indices, fs_indices = next(
            first_split.split(x_all_taxa, y)
        )

        second_split = StratifiedShuffleSplit(
            n_splits=1,
            test_size=1.0 / 3.0,
            random_state=random_state,
        )
        train_relative, test_relative = next(
            second_split.split(
                x_all_taxa[train_test_indices],
                y[train_test_indices],
            )
        )

        train_indices = train_test_indices[train_relative]
        test_indices = train_test_indices[test_relative]

        consensus_indices, method_sets, mi_threshold, rf_threshold = select_taxa(
            x_all_taxa[fs_indices],
            y[fs_indices],
            random_state,
        )

        selected_taxa_ids = [taxa_ids[index] for index in consensus_indices]
        selected_taxa_names = [
            taxon_lookup.get(taxon_id, taxon_id)
            for taxon_id in selected_taxa_ids
        ]

        for taxon_index, taxon_id in enumerate(taxa_ids):
            selection_rows.append(
                {
                    "Split": split_number,
                    "random_state": random_state,
                    "rankID": taxon_id,
                    "Taxon": taxon_lookup.get(taxon_id, taxon_id),
                    "ANOVA": int(taxon_index in method_sets["ANOVA"]),
                    "MI": int(taxon_index in method_sets["MI"]),
                    "RF": int(taxon_index in method_sets["RF"]),
                    "Consensus": int(taxon_index in consensus_indices),
                    "MI threshold": mi_threshold,
                    "RF threshold": rf_threshold,
                }
            )

        x_selected_taxa = x_all_taxa[:, consensus_indices]

        feature_sets = {
            "Diet": {
                "matrix": x_diet,
                "ids": list(FILTERED_DIET_COLUMNS),
                "names": [
                    DIET_NAME_MAP[column]
                    for column in FILTERED_DIET_COLUMNS
                ],
                "groups": ["Diet"] * len(FILTERED_DIET_COLUMNS),
            },
            "Taxa": {
                "matrix": x_selected_taxa,
                "ids": selected_taxa_ids,
                "names": selected_taxa_names,
                "groups": ["Taxa"] * len(selected_taxa_ids),
            },
            "Diet + taxa": {
                "matrix": np.hstack((x_diet, x_selected_taxa)),
                "ids": list(FILTERED_DIET_COLUMNS) + selected_taxa_ids,
                "names": (
                    [DIET_NAME_MAP[column] for column in FILTERED_DIET_COLUMNS]
                    + selected_taxa_names
                ),
                "groups": (
                    ["Diet"] * len(FILTERED_DIET_COLUMNS)
                    + ["Taxa"] * len(selected_taxa_ids)
                ),
            },
        }

        print(
            f"Split {split_number}: "
            f"FS={len(fs_indices)}, Train={len(train_indices)}, "
            f"Test={len(test_indices)}, "
            f"Consensus taxa={len(consensus_indices)}"
        )

        for feature_set_name, feature_data in feature_sets.items():
            matrix = feature_data["matrix"]

            scaler = StandardScaler()
            x_train = scaler.fit_transform(matrix[train_indices])
            x_test = scaler.transform(matrix[test_indices])
            y_train = y[train_indices]
            y_test = y[test_indices]

            for model_name, model in build_classifiers().items():
                model.fit(x_train, y_train)
                predictions = model.predict(x_test)
                scores = continuous_scores(model, x_test)
                metrics = calculate_metrics(y_test, predictions, scores)

                result_rows.append(
                    {
                        "Split": split_number,
                        "random_state": random_state,
                        "Input features": feature_set_name,
                        "Model": model_name,
                        "n selected taxa": len(consensus_indices),
                        "n diet features": len(FILTERED_DIET_COLUMNS),
                        "n total features": matrix.shape[1],
                        **metrics,
                    }
                )

                if model_name == "SVM" and feature_set_name == "Diet + taxa":
                    for y_value, prediction, score in zip(
                        y_test,
                        predictions,
                        scores,
                    ):
                        svm_prediction_rows.append(
                            {
                                "Split": split_number,
                                "random_state": random_state,
                                "y_true": int(y_value),
                                "y_pred": int(prediction),
                                "score": float(score),
                            }
                        )

                    coefficients = np.asarray(model.coef_[0], dtype=float)
                    absolute_coefficients = np.abs(coefficients)
                    coefficient_order = np.argsort(-absolute_coefficients)
                    coefficient_ranks = np.empty_like(coefficient_order)
                    coefficient_ranks[coefficient_order] = (
                        np.arange(len(coefficient_order)) + 1
                    )

                    for index, (
                        feature_id,
                        feature_name,
                        feature_group,
                    ) in enumerate(
                        zip(
                            feature_data["ids"],
                            feature_data["names"],
                            feature_data["groups"],
                        )
                    ):
                        coefficient_rows.append(
                            {
                                "Split": split_number,
                                "random_state": random_state,
                                "Feature ID": feature_id,
                                "Feature name": feature_name,
                                "Feature group": feature_group,
                                "Coefficient": coefficients[index],
                                "Absolute coefficient": absolute_coefficients[index],
                                "Rank by absolute coefficient": int(
                                    coefficient_ranks[index]
                                ),
                                "Ranked in top 20": bool(
                                    coefficient_ranks[index] <= TOP_N_FEATURES
                                ),
                            }
                        )

    results = pd.DataFrame(result_rows)
    results.to_csv(
        output_dir / "classification_results_all_splits.csv",
        index=False,
    )

    metric_columns = [
        "ROC-AUC",
        "ACC",
        "BAL-ACC",
        "SEN",
        "SPE",
        "PRE",
        "F1",
    ]
    performance_summary = (
        results
        .groupby(["Input features", "Model"], as_index=False)[metric_columns]
        .mean()
    )
    performance_summary.loc[:, metric_columns] = (
        performance_summary.loc[:, metric_columns].round(3)
    )
    performance_summary.to_csv(
        output_dir / "performance_summary.csv",
        index=False,
    )

    selection_details = pd.DataFrame(selection_rows)
    selection_details.to_csv(
        output_dir / "taxa_feature_selection_all_splits.csv",
        index=False,
    )

    taxa_frequency = (
        selection_details
        .groupby(["rankID", "Taxon"], as_index=False)
        .agg(
            Splits_selected_by_consensus=("Consensus", "sum")
        )
        .sort_values(
            ["Splits_selected_by_consensus", "Taxon"],
            ascending=[False, True],
        )
        .reset_index(drop=True)
    )
    taxa_frequency.to_csv(
        output_dir / "taxa_selection_frequency.csv",
        index=False,
    )

    svm_coefficients = pd.DataFrame(coefficient_rows)
    svm_coefficients.to_csv(
        output_dir / "svm_coefficients_all_splits.csv",
        index=False,
    )

    coefficient_summary = (
        svm_coefficients
        .groupby(
            ["Feature ID", "Feature name", "Feature group"],
            as_index=False,
        )
        .agg(
            Mean_coefficient=("Coefficient", "mean"),
            Mean_absolute_coefficient=("Absolute coefficient", "mean"),
            Splits_ranked_in_top_20=("Ranked in top 20", "sum"),
            Splits_present_in_model=("Split", "nunique"),
        )
        .sort_values(
            "Mean_absolute_coefficient",
            ascending=False,
        )
        .reset_index(drop=True)
    )

    taxa_selection_counts = taxa_frequency.set_index(
        "rankID"
    )["Splits_selected_by_consensus"]

    coefficient_summary["Splits_selected_by_FS_consensus"] = (
        coefficient_summary["Feature ID"].map(taxa_selection_counts)
    )
    coefficient_summary.loc[
        coefficient_summary["Feature group"].eq("Diet"),
        "Splits_selected_by_FS_consensus",
    ] = np.nan

    coefficient_summary.to_csv(
        output_dir / "svm_coefficient_summary.csv",
        index=False,
    )

    pd.DataFrame(svm_prediction_rows).to_csv(
        output_dir / "svm_predictions_all_splits.csv",
        index=False,
    )

    print(f"\nAnalysis completed. Outputs saved to: {output_dir}")


def parse_arguments() -> argparse.Namespace:
    """Parse command-line paths."""
    parser = argparse.ArgumentParser(
        description=(
            "Run Dataset 1 taxa feature selection and classification using "
            "the filtered dietary-variable set."
        )
    )
    parser.add_argument(
        "--metadata-file",
        type=Path,
        required=True,
        help="Path to participants_v3_corrected_depflag.tsv",
    )
    parser.add_argument(
        "--clr-file",
        type=Path,
        required=True,
        help="Path to taxonomy_leaf_CLR_v2.xlsx",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for generated CSV outputs",
    )
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_arguments()
    run_analysis(
        metadata_file=arguments.metadata_file,
        clr_file=arguments.clr_file,
        output_dir=arguments.output_dir,
    )

