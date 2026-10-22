#!/usr/bin/env python3
"""
Script to analyze and compare enrichment results from rloop-eclip-engine.
"""

import sys
from pathlib import Path

import polars as pl


def load_results(result_file: Path) -> dict:
    """Load enrichment results from Parquet file."""
    df = pl.read_parquet(result_file)
    return df.to_dicts()[0]


def print_comparison(results_list: list[tuple[str, dict]]) -> None:
    """Print formatted comparison table of multiple results."""
    print("\n" + "=" * 80)
    print("ENRICHMENT ANALYSIS COMPARISON")
    print("=" * 80)

    headers = ["Metric"] + [name for name, _ in results_list]

    metrics = [
        ("Dataset", lambda r, n: n),
        ("eCLIP Peaks", lambda r, n: f"{r['n_eclip_peaks']:,}"),
        ("R-loop Regions", lambda r, n: f"{r['n_drip_regions']:,}"),
        ("Observed Overlap (bp)", lambda r, n: f"{r['observed_overlap_bp']:,}"),
        ("Expected Overlap (bp)", lambda r, n: f"{r['expected_overlap_bp']:,.1f}"),
        ("Log2 Fold Enrichment", lambda r, n: f"{r['log2_fold_enrichment']:.3f}"),
        ("Fold Change", lambda r, n: f"{2**r['log2_fold_enrichment']:.2f}x"),
        ("Empirical p-value", lambda r, n: f"{r['empirical_pvalue']:.2e}"),
        ("Significance", lambda r, n: get_significance(r['empirical_pvalue'])),
        ("95% CI Lower", lambda r, n: f"{r['ci_95_lower']:,.1f}"),
        ("95% CI Upper", lambda r, n: f"{r['ci_95_upper']:,.1f}"),
        ("Permutations", lambda r, n: f"{r['n_iterations']:,}"),
        ("Runtime (sec)", lambda r, n: f"{r['runtime_seconds']:.1f}"),
    ]

    # Print table
    col_width = 30
    print(f"{'Metric':<{col_width}}", end="")
    for name, _ in results_list:
        print(f"{name:>20}", end="")
    print()
    print("-" * (col_width + 20 * len(results_list)))

    for metric_name, metric_fn in metrics:
        print(f"{metric_name:<{col_width}}", end="")
        for name, result in results_list:
            value = metric_fn(result, name)
            print(f"{value:>20}", end="")
        print()

    print("=" * 80 + "\n")


def get_significance(pvalue: float) -> str:
    """Return significance annotation based on p-value."""
    if pvalue < 0.001:
        return "***"
    elif pvalue < 0.01:
        return "**"
    elif pvalue < 0.05:
        return "*"
    else:
        return "n.s."


def print_interpretation(name: str, result: dict) -> None:
    """Print biological interpretation of results."""
    log2fc = result['log2_fold_enrichment']
    fc = 2 ** log2fc
    pval = result['empirical_pvalue']

    print(f"\nINTERPRETATION: {name}")
    print("-" * 60)

    if pval < 0.05 and log2fc > 0:
        print(f"✓ SIGNIFICANT ENRICHMENT detected (p = {pval:.2e})")
        print(f"  eCLIP peaks are {fc:.2f}-fold enriched at R-loop regions")
        print(f"  Observed {result['observed_overlap_bp']:,} bp overlap vs")
        print(f"  expected {result['expected_overlap_bp']:,.0f} bp by chance")
    elif pval < 0.05 and log2fc < 0:
        print(f"✓ SIGNIFICANT DEPLETION detected (p = {pval:.2e})")
        print(f"  eCLIP peaks are {1/fc:.2f}-fold depleted at R-loop regions")
    elif pval >= 0.05:
        print(f"✗ NO SIGNIFICANT ENRICHMENT (p = {pval:.2e})")
        print(f"  Overlap is consistent with random chance")

    print()


def main():
    if len(sys.argv) < 2:
        print("Usage: analyze_results.py <result1.parquet> [result2.parquet] ...")
        sys.exit(1)

    results_list = []
    for result_file in sys.argv[1:]:
        path = Path(result_file)
        if not path.exists():
            print(f"Error: File not found: {result_file}")
            continue

        result = load_results(path)
        name = path.stem
        results_list.append((name, result))

    if not results_list:
        print("Error: No valid result files provided")
        sys.exit(1)

    # Print comparison if multiple results
    if len(results_list) > 1:
        print_comparison(results_list)

    # Print detailed interpretation for each result
    for name, result in results_list:
        print_interpretation(name, result)


if __name__ == "__main__":
    main()
