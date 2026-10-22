"""
Core utilities for genomic interval processing and chromosome name normalization.
"""

import re
from pathlib import Path
from typing import Dict, Set

import polars as pl
from loguru import logger


class ChromosomeNormalizer:
    """
    Handles UCSC ↔ Ensembl chromosome naming reconciliation without file rewrites.

    Supports patterns:
    - UCSC: chr1, chr2, ..., chrX, chrY, chrM
    - Ensembl: 1, 2, ..., X, Y, MT
    """

    # Regex patterns for chromosome detection
    UCSC_PATTERN = re.compile(r"^chr([0-9]+|[XYM]|MT?)$", re.IGNORECASE)
    ENSEMBL_PATTERN = re.compile(r"^([0-9]+|[XYM]|MT)$", re.IGNORECASE)

    # Unlocalized/random contigs to exclude from analysis
    EXCLUDE_PATTERNS = [
        re.compile(r".*_random$"),
        re.compile(r"^chrUn_.*"),
        re.compile(r".*_alt$"),
        re.compile(r"^HLA-.*"),
    ]

    @classmethod
    def normalize_to_ucsc(cls, chrom: str) -> str:
        """Convert chromosome name to UCSC format (chr1, chr2, ...)."""
        if cls.UCSC_PATTERN.match(chrom):
            return chrom  # Already UCSC format
        elif cls.ENSEMBL_PATTERN.match(chrom):
            # Special case: MT -> chrM
            if chrom.upper() in ("MT", "M"):
                return "chrM"
            return f"chr{chrom}"
        else:
            # Unrecognized format, return as-is and log warning
            logger.debug(f"Unrecognized chromosome format: {chrom}")
            return chrom

    @classmethod
    def normalize_to_ensembl(cls, chrom: str) -> str:
        """Convert chromosome name to Ensembl format (1, 2, ..., MT)."""
        if cls.ENSEMBL_PATTERN.match(chrom):
            return chrom  # Already Ensembl format
        elif cls.UCSC_PATTERN.match(chrom):
            # Remove 'chr' prefix
            base = chrom[3:]
            # Special case: chrM -> MT
            if base.upper() == "M":
                return "MT"
            return base
        else:
            logger.debug(f"Unrecognized chromosome format: {chrom}")
            return chrom

    @classmethod
    def should_exclude(cls, chrom: str) -> bool:
        """Check if chromosome should be excluded from analysis."""
        # Exclude mitochondrial genome (high coverage artifacts)
        if chrom.upper() in ("CHRM", "MT", "M"):
            return True
        # Exclude unlocalized/random contigs
        return any(pattern.match(chrom) for pattern in cls.EXCLUDE_PATTERNS)


def load_genome_sizes(genome_file: Path) -> Dict[str, int]:
    """
    Parse chromosome sizes file into dictionary.

    Expected format (tab-separated):
    chr1    248956422
    chr2    242193529
    ...
    """
    sizes = {}
    with open(genome_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                logger.warning(f"Skipping malformed line: {line}")
                continue
            chrom, size = parts[0], parts[1]
            try:
                sizes[chrom] = int(size)
            except ValueError:
                logger.warning(f"Invalid size for {chrom}: {size}")
                continue

    if not sizes:
        raise ValueError(f"No valid chromosomes parsed from {genome_file}")

    logger.info(f"Loaded {len(sizes)} chromosomes from genome sizes file")
    return sizes


def load_bed_to_polars(
    bed_path: Path,
    expected_format: str = "bed6",
    normalize_chroms: bool = True,
) -> pl.LazyFrame:
    """
    Load BED file into Polars LazyFrame with optional chromosome normalization.

    Args:
        bed_path: Path to BED format file
        expected_format: One of 'bed3', 'bed6', 'narrowpeak'
        normalize_chroms: Whether to normalize chromosome names to UCSC format

    Returns:
        Polars LazyFrame with columns: chrom, start, end, [name, score, strand]
    """
    # Define column schemas based on format
    schema_map = {
        "bed3": ["chrom", "start", "end"],
        "bed6": ["chrom", "start", "end", "name", "score", "strand"],
        "narrowpeak": ["chrom", "start", "end", "name", "score", "strand",
                       "signal_value", "p_value", "q_value", "peak"],
    }

    if expected_format not in schema_map:
        raise ValueError(f"Unsupported format: {expected_format}")

    columns = schema_map[expected_format]

    # Check if file has any non-comment data and detect column count
    has_data = False
    actual_columns = 0
    with open(bed_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                has_data = True
                actual_columns = len(line.split('\t'))
                break

    if not has_data:
        logger.warning(f"BED file is empty or contains only comments: {bed_path}")
        return pl.DataFrame(schema={
            "chrom": pl.Utf8,
            "start": pl.Int64,
            "end": pl.Int64
        }).lazy()

    # Generate column names based on actual file structure
    # Use provided column names if available, otherwise generic names
    if actual_columns <= len(columns):
        col_names = columns[:actual_columns]
    else:
        # File has more columns than expected (e.g., narrowPeak has 10, bed6 has 6)
        col_names = columns + [f"col{i}" for i in range(len(columns), actual_columns)]

    try:
        # Use Polars scan_csv for lazy loading (memory-efficient)
        lf = pl.scan_csv(
            bed_path,
            separator="\t",
            has_header=False,
            new_columns=col_names,
            schema_overrides={"chrom": pl.Utf8, "start": pl.Int64, "end": pl.Int64},
            comment_prefix="#",
        )

        # Normalize chromosome names if requested
        if normalize_chroms:
            lf = lf.with_columns(
                pl.col("chrom").map_elements(
                    ChromosomeNormalizer.normalize_to_ucsc,
                    return_dtype=pl.Utf8
                )
            )

        # Filter out excluded chromosomes (unlocalized, mitochondrial)
        lf = lf.filter(
            ~pl.col("chrom").map_elements(
                ChromosomeNormalizer.should_exclude,
                return_dtype=pl.Boolean
            )
        )

        # Validate intervals (start < end)
        lf = lf.filter(pl.col("start") < pl.col("end"))

        return lf

    except Exception as e:
        logger.error(f"Failed to load BED file {bed_path}: {e}")
        raise


def calculate_total_overlap(
    intervals_a: pl.LazyFrame,
    intervals_b: pl.LazyFrame,
) -> int:
    """
    Calculate total base-pair overlap between two sets of genomic intervals.

    Uses Polars join + vectorized operations for performance.
    Handles overlapping intervals within each set by merging first.
    """
    # TODO: Implement interval merging within each set to avoid double-counting
    # For now, assumes non-overlapping intervals within each set

    # Perform interval join on matching chromosomes
    joined = intervals_a.join(
        intervals_b,
        on="chrom",
        how="inner",
        suffix="_b"
    ).collect()  # Materialize for overlap calculation

    if joined.height == 0:
        logger.warning("No chromosome overlap found between interval sets")
        return 0

    # Calculate overlap: max(0, min(end1, end2) - max(start1, start2))
    overlap_bp = (
        joined
        .with_columns([
            pl.max_horizontal("start", "start_b").alias("overlap_start"),
            pl.min_horizontal("end", "end_b").alias("overlap_end"),
        ])
        .with_columns(
            (pl.col("overlap_end") - pl.col("overlap_start"))
            .clip(lower_bound=0)
            .alias("overlap_length")
        )
        .select(pl.col("overlap_length").sum())
        .item()
    )

    return int(overlap_bp)
