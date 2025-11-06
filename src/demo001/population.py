"""
Population analysis functions for clinical trials.

Functions for analyzing study populations, enrollment, and disposition.
"""

from __future__ import annotations

import polars as pl


def count_by_treatment(data: pl.DataFrame, population_name: str) -> pl.DataFrame:
    """
    Count participants by treatment group and add population label.

    This function is commonly used in population analysis tables to count
    participants in different analysis populations (ITT, Safety, Efficacy, etc.).

    Parameters:
    -----------
    data : pl.DataFrame
        DataFrame containing participant data with treatment information
    population_name : str
        Name/label for the population being counted

    Returns:
    --------
    pl.DataFrame
        DataFrame with columns: TRT01P, n, population

    Examples:
    ---------
    >>> # Count all randomized participants
    >>> count_by_treatment(adsl, "Participants in population")

    >>> # Count ITT population
    >>> adsl_itt = adsl.filter(pl.col("ITTFL") == "Y")
    >>> count_by_treatment(adsl_itt, "Intent-to-treat population")
    """
    return (
        data.group_by("TRT01P")
        .agg(n=pl.len())
        .with_columns(population=pl.lit(population_name))
    )


def create_population_summary(adsl: pl.DataFrame) -> pl.DataFrame:
    """
    Create a comprehensive population summary table.

    Parameters:
    -----------
    adsl : pl.DataFrame
        Subject-level analysis dataset (ADSL)

    Returns:
    --------
    pl.DataFrame
        Population summary with all standard populations
    """
    populations = []

    # All randomized participants
    pop_all = count_by_treatment(adsl, "Participants in population")
    populations.append(pop_all)

    # Intent-to-treat population
    if "ITTFL" in adsl.columns:
        adsl_itt = adsl.filter(pl.col("ITTFL") == "Y")
        pop_itt = count_by_treatment(adsl_itt, "Participants included in ITT population")
        populations.append(pop_itt)

    # Efficacy population
    if "EFFFL" in adsl.columns:
        adsl_eff = adsl.filter(pl.col("EFFFL") == "Y")
        pop_eff = count_by_treatment(adsl_eff, "Participants included in efficacy population")
        populations.append(pop_eff)

    # Safety population
    if "SAFFL" in adsl.columns:
        adsl_saf = adsl.filter(pl.col("SAFFL") == "Y")
        pop_saf = count_by_treatment(adsl_saf, "Participants included in safety population")
        populations.append(pop_saf)

    return pl.concat(populations, how="diagonal")


def format_population_table(
    pop_summary: pl.DataFrame, totals: pl.DataFrame
) -> pl.DataFrame:
    """
    Format population summary with percentages for display.

    Parameters:
    -----------
    pop_summary : pl.DataFrame
        Population summary from create_population_summary()
    totals : pl.DataFrame
        Total counts by treatment group

    Returns:
    --------
    pl.DataFrame
        Formatted table with n (%) display values
    """
    # Calculate percentages
    stats_with_pct = pop_summary.join(totals, on="TRT01P").with_columns(
        pct=(100.0 * pl.col("n") / pl.col("total")).round(1)
    )

    # Format display values
    formatted_stats = stats_with_pct.with_columns(
        display=pl.when(pl.col("population") == "Participants in population")
        .then(pl.col("n").cast(str))
        .otherwise(
            pl.concat_str(
                [
                    pl.col("n").cast(str),
                    pl.lit(" ("),
                    pl.col("pct").round(1).cast(str),
                    pl.lit(")"),
                ]
            )
        )
    )

    # Pivot to wide format
    return formatted_stats.pivot(
        values="display", index="population", on="TRT01P", maintain_order=True
    )
