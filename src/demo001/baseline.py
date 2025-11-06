"""
Baseline characteristics analysis functions.

Functions for summarizing continuous and categorical variables in clinical trials.
"""

from __future__ import annotations

import polars as pl


def summarize_continuous(df: pl.DataFrame, var: str) -> pl.DataFrame:
    """
    Calculate summary statistics for continuous variables.

    Provides standard summary statistics used in baseline characteristics tables:
    mean, standard deviation, median, minimum, and maximum values.

    Parameters:
    -----------
    df : pl.DataFrame
        DataFrame containing the data
    var : str
        Name of the continuous variable to summarize

    Returns:
    --------
    pl.DataFrame
        DataFrame with summary statistics by treatment group

    Examples:
    ---------
    >>> # Summarize age by treatment
    >>> age_stats = summarize_continuous(adsl, "AGE")
    """
    return df.group_by("TRT01P").agg(
        [
            pl.col(var).mean().round(1).alias("mean"),
            pl.col(var).std().round(2).alias("sd"),
            pl.col(var).median().alias("median"),
            pl.col(var).min().alias("min"),
            pl.col(var).max().alias("max"),
            pl.len().alias("n"),
        ]
    )


def summarize_categorical(df: pl.DataFrame, var: str) -> pl.DataFrame:
    """
    Calculate counts and percentages for categorical variables.

    Provides standard summary for categorical variables used in baseline
    characteristics tables.

    Parameters:
    -----------
    df : pl.DataFrame
        DataFrame containing the data
    var : str
        Name of the categorical variable to summarize

    Returns:
    --------
    pl.DataFrame
        DataFrame with counts and percentages by treatment and category

    Examples:
    ---------
    >>> # Summarize sex by treatment
    >>> sex_stats = summarize_categorical(adsl, "SEX")
    """
    # Get counts by treatment and category
    counts = df.group_by(["TRT01P", var]).len()

    # Get treatment totals for percentage calculations
    totals = df.group_by("TRT01P").len().rename({"len": "total"})

    # Calculate percentages
    result = counts.join(totals, on="TRT01P").with_columns(
        [(100.0 * pl.col("len") / pl.col("total")).round(1).alias("pct")]
    )

    return result


def get_value(df: pl.DataFrame, treatment: str) -> str:
    """
    Get value for a specific treatment group or return default.

    Helper function to extract values from summary statistics for table building.

    Parameters:
    -----------
    df : pl.DataFrame
        DataFrame with treatment groups and values
    treatment : str
        Treatment group name to extract

    Returns:
    --------
    str
        Value for the treatment group or default "0 (0.0%)"

    Examples:
    ---------
    >>> # Get mean age for placebo group
    >>> placebo_age = get_value(age_formatted, "Placebo")
    """
    result = df.filter(pl.col("TRT01P") == treatment)
    return result[result.columns[-1]][0] if result.height > 0 else "0 (0.0%)"


def format_continuous_stats(stats: pl.DataFrame) -> pl.DataFrame:
    """
    Format continuous variable statistics for display.

    Parameters:
    -----------
    stats : pl.DataFrame
        Output from summarize_continuous()

    Returns:
    --------
    pl.DataFrame
        Formatted statistics with mean_sd and median_range columns
    """
    return stats.with_columns(
        [
            pl.format("{} ({})", pl.col("mean"), pl.col("sd")).alias("mean_sd"),
            pl.format(
                "{} [{}, {}]", pl.col("median"), pl.col("min"), pl.col("max")
            ).alias("median_range"),
        ]
    ).select(["TRT01P", "mean_sd", "median_range"])


def format_categorical_stats(stats: pl.DataFrame, var_col: str) -> pl.DataFrame:
    """
    Format categorical variable statistics for display.

    Parameters:
    -----------
    stats : pl.DataFrame
        Output from summarize_categorical()
    var_col : str
        Name of the categorical variable column

    Returns:
    --------
    pl.DataFrame
        Formatted statistics with n_pct column
    """
    return stats.with_columns(
        pl.format("{} ({}%)", pl.col("len"), pl.col("pct")).alias("n_pct")
    ).select(["TRT01P", var_col, "n_pct"])


def create_baseline_table(
    adsl: pl.DataFrame,
    continuous_vars: list[str],
    categorical_vars: list[str],
    treatments: list[str],
) -> pl.DataFrame:
    """
    Create a complete baseline characteristics table.

    Parameters:
    -----------
    adsl : pl.DataFrame
        Subject-level analysis dataset
    continuous_vars : list[str]
        List of continuous variable names
    categorical_vars : list[str]
        List of categorical variable names
    treatments : list[str]
        List of treatment group names

    Returns:
    --------
    pl.DataFrame
        Complete baseline characteristics table
    """
    table_rows = []

    # Process continuous variables
    for var in continuous_vars:
        var_stats = summarize_continuous(adsl, var)
        var_formatted = format_continuous_stats(var_stats)

        # Add variable header
        label = f"{var.title()} (years)" if var == "AGE" else var.title()
        table_rows.append([label, "", "", ""])

        # Add mean (SD) row
        mean_row = ["  Mean (SD)"] + [
            get_value(var_formatted.select(["TRT01P", "mean_sd"]), trt).replace(
                "0 (0.0%)", ""
            )
            for trt in treatments
        ]
        table_rows.append(mean_row)

        # Add median [min, max] row
        median_row = ["  Median [Min, Max]"] + [
            get_value(var_formatted.select(["TRT01P", "median_range"]), trt).replace(
                "0 (0.0%)", ""
            )
            for trt in treatments
        ]
        table_rows.append(median_row)

    # Process categorical variables
    for var in categorical_vars:
        var_stats = summarize_categorical(adsl, var)
        var_formatted = format_categorical_stats(var_stats, var)

        # Add variable header
        table_rows.append([var.title(), "", "", ""])

        # Add category rows
        # Get categories from the original adsl data to preserve order
        categories = adsl[var].unique(maintain_order=True).drop_nulls().to_list()

        for category in categories:
            if category not in var_formatted[var].unique().to_list():
                continue
            cat_data = var_formatted.filter(pl.col(var) == category)
            cat_row = [f"  {category}"] + [
                get_value(cat_data, trt) for trt in treatments
            ]
            table_rows.append(cat_row)

    # Create DataFrame from table rows
    return pl.DataFrame(
        table_rows, schema=["Characteristic"] + treatments, orient="row"
    )
