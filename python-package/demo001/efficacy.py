"""
Efficacy analysis functions for clinical trials.

Functions for performing ANCOVA, summary statistics, and efficacy endpoint analysis.
"""

from __future__ import annotations

import polars as pl
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
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
                "Endpoint_Mean": float(
                    trt_data.select(pl.col("*").exclude("USUBJID", "TRTP")).columns[2]
                ).mean()
                if len(trt_data.columns) > 3
                else float(trt_data["Week_24_LOCF"].mean()),
                "Endpoint_SD": float(
                    trt_data.select(pl.col("*").exclude("USUBJID", "TRTP")).columns[2]
                ).std()
                if len(trt_data.columns) > 3
                else float(trt_data["Week_24_LOCF"].std()),
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
    # Convert to pandas for statsmodels
    ancova_df = data.to_pandas()
    ancova_df["TRTP"] = pd.Categorical(ancova_df["TRTP"], categories=treatments)

    # Fit ANCOVA model
    model = smf.ols("CHG ~ TRTP + BASE", data=ancova_df).fit()

    # Calculate LS means
    base_mean = ancova_df["BASE"].mean()
    var_cov = model.cov_params()
    ls_means = []

    for i, trt in enumerate(treatments):
        # Create prediction vector
        x_pred = np.array([1, int(i == 1), int(i == 2), base_mean])

        # Calculate LS mean
        ls_mean = model.predict(pd.DataFrame({"TRTP": [trt], "BASE": [base_mean]}))[0]

        # Calculate standard error
        se_pred = np.sqrt(x_pred @ var_cov @ x_pred.T)

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
        placebo = treatments[0]  # Assume first treatment is placebo

        for i, trt in enumerate(treatments[1:], 1):
            coef_name = f"TRTP[T.{trt}]"
            if coef_name in model.params:
                coef = model.params[coef_name]
                se = model.bse[coef_name]
                t_stat = coef / se
                df = model.df_resid
                p_value = 2 * (1 - scipy_stats.t.cdf(abs(t_stat), df))

                ci_lower = coef - scipy_stats.t.ppf(0.975, df) * se
                ci_upper = coef + scipy_stats.t.ppf(0.975, df) * se

                comparisons.append(
                    {
                        "Comparison": f"{trt} vs. {placebo}",
                        "Estimate": coef,
                        "SE": se,
                        "CI_Lower": ci_lower,
                        "CI_Upper": ci_upper,
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
