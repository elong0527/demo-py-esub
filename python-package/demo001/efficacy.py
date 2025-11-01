"""
Efficacy analysis functions for clinical trials.

Functions for performing ANCOVA, summary statistics, and efficacy endpoint analysis.
"""

from __future__ import annotations

import polars as pl
import numpy as np
import statsmodels.api as sm
from patsy import dmatrices, dmatrix
from scipy import stats as scipy_stats


def prepare_locf_data(
    adlbc: pl.DataFrame,
    adsl_eff: pl.DataFrame,
    param: str = "GLUC",
    endpoint_week: int = 24,
) -> pl.DataFrame:
    """
    Prepare efficacy data with Last Observation Carried Forward (LOCF) imputation.

    Parameters:
    -----------
    adlbc : pl.DataFrame
        Laboratory data (ADLBC)
    adsl_eff : pl.DataFrame
        Efficacy population from ADSL
    param : str
        Parameter code (e.g., "GLUC" for glucose)
    endpoint_week : int
        Primary endpoint week

    Returns:
    --------
    pl.DataFrame
        Analysis-ready dataset with LOCF values
    """
    # Filter to efficacy population and parameter
    data = (
        adlbc.join(adsl_eff, on="USUBJID", how="inner")
        .filter((pl.col("PARAMCD") == param) & (pl.col("AVISITN") <= endpoint_week))
        .sort(["USUBJID", "AVISITN"])
    )

    # Apply LOCF strategy
    locf_data = (
        data.group_by("USUBJID")
        .agg(
            [
                pl.col("TRTP").first(),
                pl.col("BASE").first(),
                pl.col("AVAL").filter(pl.col("AVISITN") == 0).first().alias("Baseline"),
                pl.col("AVAL").last().alias(f"Week_{endpoint_week}_LOCF"),
                pl.col("AVISITN").max().alias("Last_Visit"),
            ]
        )
        .filter(
            pl.col("Baseline").is_not_null()
            & pl.col(f"Week_{endpoint_week}_LOCF").is_not_null()
        )
        .with_columns(
            (pl.col(f"Week_{endpoint_week}_LOCF") - pl.col("Baseline")).alias("CHG")
        )
    )

    return locf_data


def calculate_descriptive_stats(
    data: pl.DataFrame, treatments: list[str]
) -> list[dict[str, any]]:
    """
    Calculate descriptive statistics by treatment group.

    Parameters:
    -----------
    data : pl.DataFrame
        Analysis dataset with Baseline, Week_X_LOCF, and CHG columns
    treatments : list[str]
        List of treatment group names

    Returns:
    --------
    list[dict[str, any]]
        Descriptive statistics for each treatment group
    """
    desc_stats = []

    for trt in treatments:
        trt_data = data.filter(pl.col("TRTP") == trt)

        if trt_data.height > 0:
            stats_dict = {
                "Treatment": trt,
                "N_Analysis": trt_data.height,
                "Baseline_Mean": float(trt_data["Baseline"].mean()),
                "Baseline_SD": float(trt_data["Baseline"].std()),
                "Endpoint_Mean": float(trt_data["Week_24_LOCF"].mean()),
                "Endpoint_SD": float(trt_data["Week_24_LOCF"].std()),
                "Change_Mean": float(trt_data["CHG"].mean()),
                "Change_SD": float(trt_data["CHG"].std()),
            }
        else:
            stats_dict = {
                "Treatment": trt,
                "N_Analysis": 0,
                "Baseline_Mean": np.nan,
                "Baseline_SD": np.nan,
                "Endpoint_Mean": np.nan,
                "Endpoint_SD": np.nan,
                "Change_Mean": np.nan,
                "Change_SD": np.nan,
            }

        desc_stats.append(stats_dict)

    return desc_stats


def perform_ancova(data: pl.DataFrame, treatments: list[str]) -> dict[str, any]:
    """
    Perform ANCOVA analysis for efficacy endpoint.

    Parameters:
    -----------
    data : pl.DataFrame
        Analysis dataset with CHG, TRTP, and BASE columns
    treatments : list[str]
        List of treatment group names

    Returns:
    --------
    dict[str, any]
        ANCOVA results including model, LS means, and comparisons
    """
    # Prepare data for patsy
    data_dict = {
        "CHG": data["CHG"].to_numpy(),
        "BASE": data["BASE"].to_numpy(),
        "TRTP": data["TRTP"].to_numpy(),
    }

    # Set treatment reference level
    formula = f"CHG ~ C(TRTP, Treatment(reference='{treatments[0]}')) + BASE"

    # Create design matrices
    y, X = dmatrices(formula, data=data_dict, return_type="matrix")

    # Fit ANCOVA model
    model = sm.OLS(y, X).fit()

    # Get coefficient names from design info
    coef_names = X.design_info.column_names

    # Calculate LS means
    base_mean = data["BASE"].mean()
    ls_means = []

    for trt in treatments:
        # Create a new data point for prediction
        new_data = {"TRTP": [trt], "BASE": [base_mean]}
        # Create the design matrix for the new data point
        X_pred_patsy = dmatrix(X.design_info, data=new_data, return_type="matrix")

        # Get the prediction
        prediction = model.get_prediction(X_pred_patsy)
        ls_mean = prediction.predicted_mean[0]
        se_pred = prediction.se_mean[0]

        ls_means.append(
            {
                "Treatment": trt,
                "LS_Mean": ls_mean,
                "SE": se_pred,
                "CI_Lower": ls_mean - 1.96 * se_pred,
                "CI_Upper": ls_mean + 1.96 * se_pred,
            }
        )

    # Pairwise comparisons vs placebo
    comparisons = []
    if len(treatments) > 1:
        placebo = treatments[0]

        for trt in treatments[1:]:
            coef_name = f"C(TRTP, Treatment(reference='{placebo}'))[T.{trt}]"
            if coef_name in coef_names:
                idx = coef_names.index(coef_name)
                coef = model.params[idx]
                se = model.bse[idx]
                t_stat = model.tvalues[idx]
                p_value = model.pvalues[idx]
                ci = model.conf_int()[idx]

                comparisons.append(
                    {
                        "Comparison": f"{trt} vs. {placebo}",
                        "Estimate": coef,
                        "SE": se,
                        "CI_Lower": ci[0],
                        "CI_Upper": ci[1],
                        "t_stat": t_stat,
                        "p_value": p_value,
                    }
                )

    return {
        "model": model,
        "ls_means": ls_means,
        "comparisons": comparisons,
        "baseline_mean": base_mean,
    }


def format_efficacy_table(
    desc_stats: list[dict[str, any]], ls_means: list[dict[str, any]]
) -> pl.DataFrame:
    """
    Format efficacy results into a publication-ready table.

    Parameters:
    -----------
    desc_stats : list[dict[str, any]]
        Descriptive statistics from calculate_descriptive_stats()
    ls_means : list[dict[str, any]]
        LS means from perform_ancova()

    Returns:
    --------
    pl.DataFrame
        Formatted efficacy table
    """
    table_data = []

    for i, (desc, ls) in enumerate(zip(desc_stats, ls_means)):
        row = [
            desc["Treatment"],
            str(desc["N_Analysis"]),
            f"{desc['Baseline_Mean']:.1f} ({desc['Baseline_SD']:.2f})"
            if not np.isnan(desc["Baseline_Mean"])
            else "",
            str(desc["N_Analysis"]),
            f"{desc['Endpoint_Mean']:.1f} ({desc['Endpoint_SD']:.2f})"
            if not np.isnan(desc["Endpoint_Mean"])
            else "",
            str(desc["N_Analysis"]),
            f"{desc['Change_Mean']:.1f} ({desc['Change_SD']:.2f})"
            if not np.isnan(desc["Change_Mean"])
            else "",
            f"{ls['LS_Mean']:.2f} ({ls['CI_Lower']:.2f}, {ls['CI_Upper']:.2f})",
        ]
        table_data.append(row)

    return pl.DataFrame(
        table_data,
        schema=[
            "Treatment",
            "N_Base",
            "Mean_SD_Base",
            "N_End",
            "Mean_SD_End",
            "N_Chg",
            "Mean_SD_Chg",
            "LS_Mean_CI",
        ],
        orient="row",
    )


def format_comparison_table(comparisons: list[dict[str, any]]) -> pl.DataFrame:
    """
    Format pairwise comparison results.

    Parameters:
    -----------
    comparisons : list[dict[str, any]]
        Comparison results from perform_ancova()

    Returns:
    --------
    pl.DataFrame
        Formatted comparison table
    """
    table_data = []

    for comp in comparisons:
        row = [
            comp["Comparison"],
            f"{comp['Estimate']:.2f} ({comp['CI_Lower']:.2f}, {comp['CI_Upper']:.2f})",
            f"{comp['p_value']:.4f}" if comp["p_value"] >= 0.0001 else "<0.0001",
        ]
        table_data.append(row)

    return pl.DataFrame(
        table_data, schema=["Comparison", "Diff_CI", "P_Value"], orient="row"
    )
