#!/usr/bin/env python3
"""
Script to summarize BED file statistics for quality control.
"""

import sys
from pathlib import Path

import polars as pl


def summarize_bed(bed_file: Path) -> None:
    """Summarize BED file statistics."""
    print(f"\nSummary of: {bed_file.name}")
    print("=" * 70)

    # Auto-detect column count
    with open(bed_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                ncols = len(line.split('\t'))
                break

    # Generate column names
    col_names = ["chrom", "start", "end"] + [f"col{i}" for i in range(3, ncols)]

    # Read file
    df = pl.read_csv(
        bed_file,
        separator="\t",
        has_header=False,
        new_columns=col_names,
        comment_prefix="#",
    ).select(["chrom", "start", "end"])

    # Calculate interval widths
    df = df.with_columns(
        (pl.col("end") - pl.col("start")).alias("width")
    )

    # Basic statistics
    print(f"Total intervals: {df.height:,}")
    print(f"Total genomic coverage: {df['width'].sum():,} bp")

    # Chromosome distribution
    print("\nChromosome distribution:")
    chrom_counts = df.group_by("chrom").agg([
        pl.len().alias("count"),
        pl.col("width").sum().alias("total_bp")
    ]).sort("count", descending=True).head(10)

    print(f"  {'Chromosome':<15} {'Intervals':>12} {'Total bp':>15}")
    print("  " + "-" * 45)
    for row in chrom_counts.iter_rows(named=True):
        print(f"  {row['chrom']:<15} {row['count']:>12,} {row['total_bp']:>15,}")

    if df.height > chrom_counts.height:
        remaining = df.height - chrom_counts.height
        print(f"  ... and {remaining:,} more chromosomes/contigs")

    # Interval width statistics
    print("\nInterval width statistics:")
    width_stats = df.select([
        pl.col("width").min().alias("min"),
        pl.col("width").quantile(0.25).alias("q25"),
        pl.col("width").median().alias("median"),
        pl.col("width").quantile(0.75).alias("q75"),
        pl.col("width").max().alias("max"),
        pl.col("width").mean().alias("mean"),
    ])

    stats = width_stats.to_dicts()[0]
    print(f"  Min:     {stats['min']:>10,} bp")
    print(f"  Q1:      {stats['q25']:>10,.0f} bp")
    print(f"  Median:  {stats['median']:>10,.0f} bp")
    print(f"  Mean:    {stats['mean']:>10,.0f} bp")
    print(f"  Q3:      {stats['q75']:>10,.0f} bp")
    print(f"  Max:     {stats['max']:>10,} bp")

    # Check for potential issues
    print("\nQuality checks:")
    invalid = df.filter(pl.col("start") >= pl.col("end"))
    if invalid.height > 0:
        print(f"  ⚠  {invalid.height:,} invalid intervals (start >= end)")
    else:
        print(f"  ✓  All intervals valid (start < end)")

    # Check for unlocalized/random contigs
    unlocalized = df.filter(
        pl.col("chrom").str.contains("random") |
        pl.col("chrom").str.contains("Un_") |
        pl.col("chrom").str.contains("^GL") |
        pl.col("chrom").str.contains("^KI")
    )
    if unlocalized.height > 0:
        pct = 100 * unlocalized.height / df.height
        print(f"  ⚠  {unlocalized.height:,} intervals ({pct:.1f}%) on unlocalized contigs")
        print(f"     (will be filtered out during analysis)")
    else:
        print(f"  ✓  All intervals on standard chromosomes")

    print("=" * 70 + "\n")


def main():
    if len(sys.argv) < 2:
        print("Usage: summarize_bed.py <file1.bed> [file2.bed] ...")
        sys.exit(1)

    for bed_file in sys.argv[1:]:
        path = Path(bed_file)
        if not path.exists():
            print(f"Error: File not found: {bed_file}")
            continue

        summarize_bed(path)


if __name__ == "__main__":
    main()
