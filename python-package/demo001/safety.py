"""
Safety analysis functions for adverse events.

Functions for analyzing adverse events, safety populations, and related safety endpoints.
"""

from __future__ import annotations

import polars as pl


def count_participants(
    df: pl.DataFrame, condition: pl.Expr | None = None
) -> pl.DataFrame:
    """
    Count unique participants meeting a condition.

    This function is the foundation for adverse event summary tables, counting
    unique participants rather than events to avoid double-counting.

    Parameters:
    -----------
    df : pl.DataFrame
        DataFrame with adverse events data
    condition : pl.Expr, optional
        Polars expression for filtering (None = count all participants)

    Returns:
    --------
    pl.DataFrame
        DataFrame with counts by treatment group

    Examples:
    ---------
    >>> # Count participants with any AE
    >>> count_participants(adae_safety)

    >>> # Count participants with serious AEs
    >>> count_participants(adae_safety, pl.col("AESER") == "Y")

    >>> # Count participants with drug-related AEs
    >>> count_participants(adae_safety,
    ...     pl.col("AEREL").is_in(["POSSIBLE", "PROBABLE", "DEFINITE", "RELATED"]))
    """
    if condition is not None:
        df = df.filter(condition)

    return df.group_by("TRT01A").agg(n=pl.col("USUBJID").n_unique())


def create_ae_summary(
    adae_safety: pl.DataFrame, pop_counts: pl.DataFrame, treatment_levels: pl.DataFrame
) -> pl.DataFrame:
    """
    Create a comprehensive adverse events summary table.

    Generates standard AE summary categories used in regulatory submissions:
    - Participants with any adverse event
    - Participants with drug-related adverse events
    - Participants with serious adverse events
    - Participants with serious drug-related adverse events
    - Participants who died
    - Participants who discontinued due to AE

    Parameters:
    -----------
    adae_safety : pl.DataFrame
        Adverse events data for safety population
    pop_counts : pl.DataFrame
        Population counts by treatment group
    treatment_levels : pl.DataFrame
        Treatment group reference for consistent ordering

    Returns:
    --------
    pl.DataFrame
        Complete AE summary table with counts and percentages
    """
    categories = []

    # 1. Participants in population
    pop_row = pop_counts.with_columns(
        category=pl.lit("Participants in population")
    ).rename({"N": "n"})
    categories.append(pop_row)

    # 2. With any adverse event
    any_ae = count_participants(adae_safety).with_columns(
        category=pl.lit("With any adverse event")
    )
    categories.append(any_ae)

    # 3. With drug-related adverse event
    drug_related = count_participants(
        adae_safety,
        pl.col("AEREL").is_in(["POSSIBLE", "PROBABLE", "DEFINITE", "RELATED"]),
    ).with_columns(category=pl.lit("With drug-related adverse event"))
    categories.append(drug_related)

    # 4. With serious adverse event
    serious = count_participants(adae_safety, pl.col("AESER") == "Y").with_columns(
        category=pl.lit("With serious adverse event")
    )
    categories.append(serious)

    # 5. With serious drug-related adverse event
    serious_drug_related = count_participants(
        adae_safety,
        (pl.col("AESER") == "Y")
        & pl.col("AEREL").is_in(["POSSIBLE", "PROBABLE", "DEFINITE", "RELATED"]),
    ).with_columns(category=pl.lit("With serious drug-related adverse event"))
    categories.append(serious_drug_related)

    # 6. Who died
    deaths = count_participants(adae_safety, pl.col("AEOUT") == "FATAL").with_columns(
        category=pl.lit("Who died")
    )
    categories.append(deaths)

    # 7. Discontinued due to adverse event
    discontinued = count_participants(
        adae_safety, pl.col("AEACN") == "DRUG WITHDRAWN"
    ).with_columns(category=pl.lit("Discontinued due to adverse event"))
    categories.append(discontinued)

    return pl.concat(categories, how="diagonal")


def format_ae_summary(
    ae_summary: pl.DataFrame, pop_counts: pl.DataFrame
) -> pl.DataFrame:
    """
    Format AE summary table with percentages for display.

    Parameters:
    -----------
    ae_summary : pl.DataFrame
        AE summary from create_ae_summary()
    pop_counts : pl.DataFrame
        Population counts for percentage calculations

    Returns:
    --------
    pl.DataFrame
        Formatted AE summary table ready for reporting
    """
    # Add population totals and calculate percentages
    formatted = ae_summary.join(
        pop_counts.select(["TRT01A", "N"]), on="TRT01A", how="left"
    ).with_columns(
        [
            # Fill missing counts with 0
            pl.col("n").fill_null(0),
            # Calculate percentage
            pl.when(pl.col("category") == "Participants in population")
            .then(None)  # No percentage for population row
            .otherwise((100.0 * pl.col("n") / pl.col("N")).round(1))
            .alias("pct"),
        ]
    )

    # Format display values
    ae_formatted = formatted.with_columns(
        [
            # Show counts as strings, including zeros
            pl.col("n").cast(str).alias("n_display"),
            # Format percentages with parentheses; blank out population row
            pl.when(pl.col("category") == "Participants in population")
            .then(pl.lit(""))
            .otherwise(pl.format("({})", pl.col("pct").fill_null(0).round(1).cast(str)))
            .alias("pct_display"),
        ]
    )

    return ae_formatted


def create_ae_by_soc_table(
    adae_safety: pl.DataFrame, pop_counts: pl.DataFrame, treatments: list[str]
) -> pl.DataFrame:
    """
    Create adverse events by System Organ Class table.

    Parameters:
    -----------
    adae_safety : pl.DataFrame
        Adverse events data for safety population
    pop_counts : pl.DataFrame
        Population counts by treatment group
    treatments : list[str]
        List of treatment names in display order

    Returns:
    --------
    pl.DataFrame
        Hierarchical table with SOC headers and specific AE terms
    """
    # Standardize AE term formatting
    ae_counts = (
        adae_safety.with_columns(
            [
                pl.col("AEDECOD").str.to_titlecase().alias("AEDECOD_STD"),
                pl.col("AEBODSYS").str.to_titlecase().alias("AEBODSYS_STD"),
            ]
        )
        .group_by(["TRT01A", "AEBODSYS_STD", "AEDECOD_STD"])
        .agg(n=pl.col("USUBJID").n_unique())
        .sort(["AEBODSYS_STD", "AEDECOD_STD", "TRT01A"])
    )

    # Initialize table with population counts
    table_data = [
        ["Participants in population"]
        + [str(pop_counts.filter(pl.col("TRT01A") == t)["N"][0]) for t in treatments],
        [""] * (len(treatments) + 1),  # Blank separator row
    ]

    # Build hierarchical structure: SOC -> Specific AE terms
    for soc in ae_counts["AEBODSYS_STD"].unique().sort():
        # Add SOC header row
        table_data.append([soc] + [""] * len(treatments))

        # Get all AE terms within this SOC
        soc_data = ae_counts.filter(pl.col("AEBODSYS_STD") == soc)

        # Add each specific AE term with counts
        for ae_term in soc_data["AEDECOD_STD"].unique().sort():
            row = [f"  {ae_term}"]  # Indent specific terms

            # Add counts for each treatment group
            for trt in treatments:
                count_data = soc_data.filter(
                    (pl.col("AEDECOD_STD") == ae_term) & (pl.col("TRT01A") == trt)
                )
                count = count_data["n"][0] if count_data.height > 0 else 0
                row.append(str(count))

            table_data.append(row)

    # Convert to Polars DataFrame
    return pl.DataFrame(
        table_data,
        schema=["System Organ Class / Preferred Term"] + treatments,
        orient="row",
    )
