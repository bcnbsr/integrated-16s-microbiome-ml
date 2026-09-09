"""Dataset 1 dietary-variable correlation analysis."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def analyse_diet_correlations(
    meta_aligned: pd.DataFrame,
    diet_cols: list[str],
    output_dir: str | Path,
    threshold: float = 0.90,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate Spearman correlations among dietary variables, identify pairs
    with an absolute correlation at or above the specified threshold, and
    save a correlation heatmap.

    Parameters
    ----------
    meta_aligned
        Participant-level metadata containing the dietary variables.
    diet_cols
        Names of the dietary-variable columns.
    output_dir
        Directory in which the heatmap will be saved.
    threshold
        Absolute Spearman-correlation threshold used to flag variable pairs.

    Returns
    -------
    correlation_matrix
        Full Spearman correlation matrix.
    highly_correlated_pairs
        Unique variable pairs with absolute correlation at or above threshold.
    """
    missing_columns = [
        column for column in diet_cols
        if column not in meta_aligned.columns
    ]
    if missing_columns:
        raise ValueError(f"Missing dietary variables: {missing_columns}")

    diet_df = meta_aligned.loc[:, diet_cols].copy()
    correlation_matrix = diet_df.corr(method="spearman")

    # Use the upper triangle only, excluding duplicate pairs and the diagonal.
    upper_triangle = np.triu(
        np.ones(correlation_matrix.shape, dtype=bool),
        k=1,
    )

    correlation_pairs = (
        correlation_matrix
        .where(upper_triangle)
        .stack()
        .rename("spearman_rho")
        .reset_index()
        .rename(
            columns={
                "level_0": "diet_variable_1",
                "level_1": "diet_variable_2",
            }
        )
    )

    correlation_pairs["abs_rho"] = correlation_pairs[
        "spearman_rho"
    ].abs()

    highly_correlated_pairs = (
        correlation_pairs.loc[
            correlation_pairs["abs_rho"] >= threshold
        ]
        .sort_values("abs_rho", ascending=False)
        .reset_index(drop=True)
    )

    print("\nSpearman correlation matrix:")
    print(correlation_matrix.round(3).to_string())

    print(f"\nDietary-variable pairs with |rho| >= {threshold:.2f}:")
    if highly_correlated_pairs.empty:
        print("No pairs met the threshold.")
    else:
        print(highly_correlated_pairs.round(3).to_string(index=False))

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(9, 7))
    sns.heatmap(
        correlation_matrix,
        annot=True,
        fmt=".2f",
        cmap="coolwarm",
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
        linewidths=0.5,
        cbar_kws={"label": "Spearman correlation coefficient"},
    )
    plt.title("Spearman Correlations Among Dietary Variables")
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(
        output_path / "dataset1_diet_correlation_heatmap.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()

    return correlation_matrix, highly_correlated_pairs

